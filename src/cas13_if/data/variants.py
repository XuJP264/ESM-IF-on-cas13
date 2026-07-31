"""Cas13d variant mutation recovery, assay separation, and label auditing."""

from __future__ import annotations

import hashlib
import json
import re
from pathlib import Path
from typing import Any, cast

import pandas as pd
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import yaml

from cas13_if.provenance import atomic_write_text

SUBSTITUTION = re.compile(r"([A-Z])(\d+)([A-Z])")
DELETION = re.compile(r"del(\d+)-(\d+)")
INSERTION = re.compile(r"ins(\d+)_(\d+):([A-Z]+)")


def apply_mutation(wild_type: str, mutation: str) -> str:
    """Apply validated one-based substitutions/deletions/insertions."""
    sequence = wild_type.upper()
    parts = [part.strip() for part in re.split(r"[,+]", mutation) if part.strip()]
    substitutions: dict[int, str] = {}
    deletions: list[tuple[int, int]] = []
    insertion: tuple[int, int, str] | None = None
    for part in parts:
        substitution = SUBSTITUTION.fullmatch(part)
        deletion = DELETION.fullmatch(part)
        inserted = INSERTION.fullmatch(part)
        if substitution:
            before, raw_position, after = substitution.groups()
            position = int(raw_position)
            if position < 1 or position > len(sequence):
                raise ValueError(f"substitution outside WT sequence: {part}")
            if sequence[position - 1] != before:
                raise ValueError(
                    f"WT token mismatch for {part}: observed {sequence[position - 1]}"
                )
            substitutions[position - 1] = after
        elif deletion:
            start, end = map(int, deletion.groups())
            if start < 1 or end < start or end > len(sequence):
                raise ValueError(f"invalid deletion: {part}")
            deletions.append((start - 1, end))
        elif inserted:
            left, right, tokens = inserted.groups()
            left_index, right_index = int(left), int(right)
            if right_index != left_index + 1 or right_index > len(sequence) + 1:
                raise ValueError(f"invalid insertion junction: {part}")
            if insertion is not None:
                raise ValueError("only one insertion is supported per curated variant")
            insertion = left_index, right_index, tokens
        else:
            raise ValueError(f"unsupported mutation notation: {part}")
    if substitutions and (deletions or insertion):
        raise ValueError("mixed substitution/indel records must be separately curated")
    if substitutions:
        output = list(sequence)
        for index, token in substitutions.items():
            output[index] = token
        return "".join(output)
    if deletions:
        covered: set[int] = set()
        for start, end in deletions:
            interval = set(range(start, end))
            if covered.intersection(interval):
                raise ValueError("overlapping deletions are not allowed")
            covered.update(interval)
        return "".join(
            token for index, token in enumerate(sequence) if index not in covered
        )
    if insertion:
        left, _, tokens = insertion
        return sequence[:left] + tokens + sequence[left:]
    raise ValueError("mutation string is empty")


def activity_label(
    cis_activity: float | None,
    trans_activity: float | None,
    *,
    active_minimum: float,
    inactive_maximum: float,
    cis_retained_minimum: float,
    trans_reduced_maximum: float,
) -> str:
    if cis_activity is None:
        return "not_labeled"
    if (
        trans_activity is not None
        and cis_activity >= cis_retained_minimum
        and trans_activity <= trans_reduced_maximum
    ):
        return "cis-retained/trans-reduced"
    if cis_activity >= active_minimum:
        return "active"
    if cis_activity < inactive_maximum:
        return "inactive"
    return "partial"


def _yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("variant manifest root must be a mapping")
    return value


def build_variant_dataset(
    *, manifest_path: Path, scaffold_csv: Path, output_dir: Path
) -> dict[str, Any]:
    manifest = _yaml(manifest_path)
    scaffolds = pd.read_csv(scaffold_csv).set_index("scaffold_id")
    defaults = dict(manifest["defaults"])
    policy = manifest["label_policy"]
    source_by_id = {item["study_id"]: item for item in manifest["sources"]}
    rows: list[dict[str, Any]] = []
    for raw in manifest["variants"]:
        row = {**defaults, **raw}
        study = source_by_id[str(row["study_id"])]
        scaffold = str(row.get("scaffold", "CasRx"))
        wild_type = str(scaffolds.loc[scaffold, "full_natural_sequence"])
        mutant = apply_mutation(wild_type, str(row["mutation"]))
        cis = (
            float(row["cis_activity"]) if row.get("cis_activity") is not None else None
        )
        trans = (
            float(row["trans_activity"])
            if row.get("trans_activity") is not None
            else None
        )
        label = str(
            row.get("label")
            or activity_label(
                cis,
                trans,
                active_minimum=float(policy["active_minimum"]),
                inactive_maximum=float(policy["inactive_maximum"]),
                cis_retained_minimum=float(policy["cis_retained_minimum"]),
                trans_reduced_maximum=float(policy["trans_reduced_maximum"]),
            )
        )
        assay_type = str(row.get("assay_type", "in_vitro_cis_and_trans_cleavage"))
        key_text = "|".join(
            (
                str(row["study_id"]),
                str(row["variant_id"]),
                str(row["mutation"]),
                assay_type,
                str(row["comparability_group"]),
            )
        )
        rows.append(
            {
                "study_id": row["study_id"],
                "publication": study["publication"],
                "doi": study["doi"],
                "official_url": study["official_url"],
                "variant_id": row["variant_id"],
                "scaffold": scaffold,
                "WT_sequence": wild_type,
                "mutation": row["mutation"],
                "full_mutant_sequence": mutant,
                "residue_numbering_system": row["residue_numbering_system"],
                "numbering_unified": True,
                "assay_type": assay_type,
                "guide": row.get("guide"),
                "target": row.get("target"),
                "cis_activity": cis,
                "trans_activity": trans,
                "cell_knockdown": row.get("cell_knockdown"),
                "expression": row.get("expression"),
                "solubility": row.get("solubility"),
                "WT_normalized_activity": cis,
                "activity_scale": "WT=1" if cis is not None else None,
                "active_inactive_partial_label": label,
                "replicate_information": row.get("replicate_information"),
                "source_table_figure_supplement": row.get(
                    "source_table_figure_supplement"
                ),
                "extraction_method": row.get("extraction_method"),
                "numeric_is_approximate": "approximately"
                in str(row.get("extraction_method", "")),
                "confidence": row.get("confidence"),
                "comparability_group": row.get("comparability_group"),
                "notes": row.get("notes"),
                "dedup_key": hashlib.sha256(key_text.encode("utf-8")).hexdigest(),
                "source_evidence": "published_wet_lab_measurement_or_category",
                "is_mock": False,
                "processing_evidence_level": 0,
            }
        )
    frame = pd.DataFrame.from_records(rows)
    duplicate = frame["dedup_key"].duplicated(keep=False)
    if duplicate.any():
        raise ValueError(
            "duplicate variant records: "
            + ",".join(frame.loc[duplicate, "variant_id"].astype(str))
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pandas(frame, preserve_index=False)
    pq.write_table(table, output_dir / "cas13d_variant_activity.parquet")
    source_frame = pd.DataFrame.from_records(manifest["sources"])
    atomic_write_text(
        output_dir / "cas13d_variant_activity_sources.csv",
        source_frame.to_csv(index=False),
    )
    sensitivity_rows: list[dict[str, Any]] = []
    numeric = frame.loc[frame["cis_activity"].notna()]
    for active_minimum in (0.70, 0.80, 0.90):
        for inactive_maximum in (0.10, 0.20, 0.30):
            labels = [
                activity_label(
                    cast(float, row.cis_activity),
                    cast(float, row.trans_activity)
                    if pd.notna(row.trans_activity)
                    else None,
                    active_minimum=active_minimum,
                    inactive_maximum=inactive_maximum,
                    cis_retained_minimum=float(policy["cis_retained_minimum"]),
                    trans_reduced_maximum=float(policy["trans_reduced_maximum"]),
                )
                for row in numeric.itertuples()
            ]
            counts = pd.Series(labels).value_counts().to_dict()
            sensitivity_rows.append(
                {
                    "active_minimum": active_minimum,
                    "inactive_maximum": inactive_maximum,
                    **{str(key): int(value) for key, value in counts.items()},
                }
            )
    sensitivity = pd.DataFrame.from_records(sensitivity_rows).fillna(0)
    atomic_write_text(
        output_dir / "cas13d_variant_label_sensitivity.csv",
        sensitivity.to_csv(index=False),
    )
    missing = {
        column: int(frame[column].isna().sum())
        for column in (
            "cis_activity",
            "trans_activity",
            "cell_knockdown",
            "guide",
            "target",
            "expression",
            "solubility",
        )
    }
    summary = {
        "records": len(frame),
        "studies": int(frame["study_id"].nunique()),
        "unique_mutations": int(frame["mutation"].nunique()),
        "comparability_groups": int(frame["comparability_group"].nunique()),
        "numeric_cis_trans_records": int(
            (frame["cis_activity"].notna() & frame["trans_activity"].notna()).sum()
        ),
        "full_mutant_sequences_recovered": int(
            frame["full_mutant_sequence"].notna().sum()
        ),
        "missing_fields": missing,
        "labels": {
            str(key): int(value)
            for key, value in frame["active_inactive_partial_label"]
            .value_counts()
            .sort_index()
            .items()
        },
        "is_mock": False,
        "processing_evidence_level": 0,
    }
    atomic_write_text(
        output_dir / "cas13d_variant_activity_summary.json",
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    return summary
