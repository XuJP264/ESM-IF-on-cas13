"""Subtype-specific MAFFT orchestration with explicit sequence exclusions."""

from __future__ import annotations

import json
import shutil
import subprocess
from collections import Counter
from pathlib import Path
from typing import Any

import pyarrow.parquet as pq  # type: ignore[import-untyped]

from cas13_if.alignments.msa import read_aligned_fasta
from cas13_if.data.fasta import write_fasta
from cas13_if.provenance import atomic_write_text
from cas13_if.schemas import STANDARD_AA


def safe_subtype_label(subtype: str) -> str:
    label = "".join(
        character.lower() if character.isalnum() else "-"
        for character in subtype.strip()
    ).strip("-")
    if not label:
        raise ValueError("subtype label is empty")
    return label


def mafft_version(executable: str = "mafft") -> str:
    path = shutil.which(executable)
    if path is None:
        raise FileNotFoundError(f"MAFFT executable not found: {executable}")
    completed = subprocess.run(
        [path, "--version"], text=True, capture_output=True, check=False
    )
    output = (completed.stdout + "\n" + completed.stderr).strip()
    if completed.returncode != 0 or not output:
        raise RuntimeError(output or "mafft --version failed")
    return output.splitlines()[0]


def build_subtype_msas(
    *,
    exact_unique_path: Path,
    cluster_mapping_path: Path,
    output_dir: Path,
    executable: str,
    threads: int,
) -> dict[str, Any]:
    """Align one representative set per subtype and audit all exclusions."""
    if threads < 1:
        raise ValueError("MAFFT threads must be positive")
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite MSA output: {output_dir}")
    executable_path = shutil.which(executable)
    if executable_path is None:
        raise FileNotFoundError(f"MAFFT executable not found: {executable}")
    unique_rows = pq.read_table(exact_unique_path).to_pylist()
    cluster_rows = pq.read_table(
        cluster_mapping_path,
        columns=["sequence_sha256", "representative_sha256"],
    ).to_pylist()
    representative_by_member = {
        str(row["sequence_sha256"]): str(row["representative_sha256"])
        for row in cluster_rows
    }
    eligible_by_subtype_cluster: dict[str, dict[str, list[tuple[str, str]]]] = {}
    excluded: list[dict[str, Any]] = []
    for row in unique_rows:
        digest = str(row["sequence_sha256"])
        representative = representative_by_member.get(digest)
        if representative is None:
            excluded.append(
                {
                    "sequence_sha256": digest,
                    "reason": "missing_cluster_mapping",
                    "invalid_symbols": [],
                }
            )
            continue
        sequence = str(row["protein_sequence"]).upper()
        invalid = sorted(set(sequence).difference(STANDARD_AA))
        subtypes = row["subtypes"]
        if not isinstance(subtypes, list):
            excluded.append(
                {
                    "sequence_sha256": digest,
                    "reason": "subtypes_not_list",
                    "invalid_symbols": invalid,
                }
            )
            continue
        type_vi_subtypes = [
            str(subtype)
            for subtype in subtypes
            if str(subtype).upper().startswith("VI-")
        ]
        if not type_vi_subtypes:
            continue
        nonconflicting = int(row.get("nonconflicting_record_count", 1))
        complete = int(row.get("complete_record_count", 1))
        if nonconflicting < 1:
            excluded.append(
                {
                    "sequence_sha256": digest,
                    "subtypes": type_vi_subtypes,
                    "reason": "only_subtype_conflicting_records",
                    "invalid_symbols": invalid,
                }
            )
            continue
        if complete < 1:
            excluded.append(
                {
                    "sequence_sha256": digest,
                    "subtypes": type_vi_subtypes,
                    "reason": "no_explicitly_complete_record",
                    "invalid_symbols": invalid,
                }
            )
            continue
        if invalid:
            excluded.append(
                {
                    "sequence_sha256": digest,
                    "subtypes": type_vi_subtypes,
                    "reason": "noncanonical_amino_acid",
                    "invalid_symbols": invalid,
                }
            )
            continue
        for subtype in type_vi_subtypes:
            eligible_by_subtype_cluster.setdefault(subtype, {}).setdefault(
                representative, []
            ).append((digest, sequence))
    by_subtype = {
        subtype: [
            min(cluster_members, key=lambda item: item[0])
            for cluster_members in clusters.values()
        ]
        for subtype, clusters in eligible_by_subtype_cluster.items()
    }
    if not by_subtype:
        raise ValueError("no canonical Type VI representative sequences for MSA")

    output_dir.mkdir(parents=True, exist_ok=False)
    version = mafft_version(executable_path)
    summary: dict[str, Any] = {}
    for subtype, records in sorted(by_subtype.items()):
        label = safe_subtype_label(subtype)
        subtype_dir = output_dir / label
        subtype_dir.mkdir()
        input_fasta = subtype_dir / "representatives.fasta"
        write_fasta(sorted(records), input_fasta)
        if len(records) < 2:
            summary[subtype] = {
                "status": "not_run",
                "reason": "fewer_than_two_representatives",
                "input_sequences": len(records),
                "is_mock": False,
            }
            continue
        output_fasta = subtype_dir / "alignment.fasta"
        temporary = subtype_dir / ".alignment.fasta.part"
        command = [
            executable_path,
            "--auto",
            "--thread",
            str(threads),
            "--reorder",
            str(input_fasta),
        ]
        with temporary.open("w", encoding="utf-8") as stdout_handle:
            completed = subprocess.run(
                command,
                text=True,
                stdout=stdout_handle,
                stderr=subprocess.PIPE,
                check=False,
            )
        atomic_write_text(subtype_dir / "stderr.log", completed.stderr)
        if completed.returncode != 0:
            atomic_write_text(
                subtype_dir / "FAILED", f"exit_code={completed.returncode}\n"
            )
            raise RuntimeError(f"MAFFT failed for {subtype}; see {subtype_dir}")
        temporary.replace(output_fasta)
        alignment = read_aligned_fasta(output_fasta)
        summary[subtype] = {
            "status": "success",
            "input_sequences": len(records),
            "alignment_sequences": alignment.n_sequences,
            "alignment_columns": alignment.n_columns,
            "mafft_version": version,
            "command": command,
            "is_mock": False,
        }
    exclusion_counts = Counter(str(row["reason"]) for row in excluded)
    manifest = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 0,
        "selection": (
            "one eligible sequence per MMseqs2 70% cluster and Type VI subtype; "
            "requires at least one nonconflicting, explicitly complete record"
        ),
        "mafft_version": version,
        "subtypes": summary,
        "excluded_sequence_count": len(excluded),
        "exclusion_counts": dict(sorted(exclusion_counts.items())),
    }
    atomic_write_text(
        output_dir / "excluded_sequences.jsonl",
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in excluded),
    )
    atomic_write_text(
        output_dir / "msa_manifest.json",
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
    )
    return manifest
