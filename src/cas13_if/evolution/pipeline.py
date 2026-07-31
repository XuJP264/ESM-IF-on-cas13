"""Audited subtype-specific conservation table generation."""

from __future__ import annotations

import json
from dataclasses import asdict
from pathlib import Path
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from cas13_if.alignments.msa import read_aligned_fasta
from cas13_if.evolution.conservation import conservation_statistics
from cas13_if.provenance import atomic_write_text


def compute_subtype_conservation(
    *,
    msa_root: Path,
    output_dir: Path,
    identity_threshold: float,
    allowed_frequency: float,
) -> dict[str, Any]:
    if output_dir.exists():
        raise FileExistsError(
            f"refusing to overwrite conservation output: {output_dir}"
        )
    manifest_path = msa_root / "msa_manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"MSA manifest is missing: {manifest_path}")
    msa_manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    subtype_entries = msa_manifest.get("subtypes")
    if not isinstance(subtype_entries, dict):
        raise ValueError("MSA manifest subtypes must be a mapping")
    output_dir.mkdir(parents=True, exist_ok=False)
    subtype_summary: dict[str, Any] = {}
    for subtype, metadata in sorted(subtype_entries.items()):
        if not isinstance(metadata, dict) or metadata.get("status") != "success":
            subtype_summary[str(subtype)] = {
                "status": "not_run",
                "reason": (
                    metadata.get("reason")
                    if isinstance(metadata, dict)
                    else "invalid_msa_manifest_entry"
                ),
            }
            continue
        command = metadata.get("command")
        if not isinstance(command, list) or not command:
            raise ValueError(f"MSA command is missing for subtype {subtype}")
        input_path = Path(str(command[-1]))
        alignment_path = input_path.parent / "alignment.fasta"
        alignment = read_aligned_fasta(alignment_path)
        statistics = conservation_statistics(
            alignment,
            identity_threshold=identity_threshold,
            allowed_frequency=allowed_frequency,
        )
        rows = [
            {
                **asdict(item),
                "allowed_residues": list(item.allowed_residues),
                "subtype": str(subtype),
                "is_mock": False,
            }
            for item in statistics
        ]
        subtype_label = input_path.parent.name
        output_path = output_dir / f"{subtype_label}.parquet"
        pq.write_table(pa.Table.from_pylist(rows), output_path, compression="zstd")
        effective_count = statistics[0].effective_sequence_count if statistics else 0.0
        subtype_summary[str(subtype)] = {
            "status": "success",
            "sequences": alignment.n_sequences,
            "columns": alignment.n_columns,
            "effective_sequence_count": effective_count,
            "mean_gap_fraction": (
                sum(item.gap_fraction for item in statistics) / len(statistics)
                if statistics
                else None
            ),
            "mean_conservation": (
                sum(item.conservation for item in statistics) / len(statistics)
                if statistics
                else None
            ),
            "output": str(output_path),
        }
    manifest = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 0,
        "identity_threshold": identity_threshold,
        "allowed_frequency": allowed_frequency,
        "subtypes": subtype_summary,
        "warning": (
            "Position-level values are subtype-specific and must not be merged "
            "across incompatible Cas13 domain layouts."
        ),
    }
    atomic_write_text(
        output_dir / "conservation_manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return manifest
