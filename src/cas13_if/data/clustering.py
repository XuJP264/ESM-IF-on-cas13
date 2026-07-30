"""MMseqs2 clustering, cluster-level splits, and leakage gates."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any


@dataclass(frozen=True)
class MmseqsParameters:
    minimum_identity: float
    coverage: float = 0.8
    coverage_mode: int = 0
    cluster_mode: int = 2
    threads: int = 16

    def validate(self) -> None:
        if not 0 < self.minimum_identity <= 1:
            raise ValueError("minimum_identity must be in (0, 1]")
        if not 0 < self.coverage <= 1:
            raise ValueError("coverage must be in (0, 1]")
        if self.threads < 1:
            raise ValueError("threads must be positive")


def mmseqs_version(executable: str = "mmseqs") -> str:
    path = shutil.which(executable)
    if path is None:
        raise FileNotFoundError("MMseqs2 executable not found")
    result = subprocess.run(
        [path, "version"], text=True, capture_output=True, check=False
    )
    if result.returncode != 0:
        raise RuntimeError(result.stderr.strip() or "mmseqs version failed")
    return result.stdout.strip()


def run_mmseqs_clustering(
    fasta: Path,
    output_dir: Path,
    parameters: MmseqsParameters,
    *,
    executable: str = "mmseqs",
) -> dict[str, Any]:
    parameters.validate()
    path = shutil.which(executable)
    if path is None:
        raise FileNotFoundError("MMseqs2 executable not found")
    output_dir.mkdir(parents=True, exist_ok=False)
    prefix = output_dir / "clusters"
    temporary = output_dir / "tmp"
    command = [
        path,
        "easy-cluster",
        str(fasta),
        str(prefix),
        str(temporary),
        "--min-seq-id",
        str(parameters.minimum_identity),
        "-c",
        str(parameters.coverage),
        "--cov-mode",
        str(parameters.coverage_mode),
        "--cluster-mode",
        str(parameters.cluster_mode),
        "--threads",
        str(parameters.threads),
    ]
    result = subprocess.run(command, text=True, capture_output=True, check=False)
    (output_dir / "stdout.log").write_text(result.stdout, encoding="utf-8")
    (output_dir / "stderr.log").write_text(result.stderr, encoding="utf-8")
    metadata = {
        "tool": "MMseqs2",
        "version": mmseqs_version(path),
        "parameters": asdict(parameters),
        "command": command,
        "exit_code": result.returncode,
        "is_mock": False,
    }
    (output_dir / "metadata.json").write_text(
        json.dumps(metadata, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    if result.returncode != 0:
        raise RuntimeError(f"MMseqs2 clustering failed; see {output_dir}")
    cluster_file = Path(f"{prefix}_cluster.tsv")
    if not cluster_file.is_file():
        raise RuntimeError(f"MMseqs2 output missing: {cluster_file}")
    mapping = parse_cluster_tsv(cluster_file)
    summary = cluster_summary(mapping)
    (output_dir / "cluster_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return {"metadata": metadata, "mapping": mapping, "summary": summary}


def parse_cluster_tsv(path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 2 or not all(fields):
                raise ValueError(f"malformed MMseqs cluster row {line_number}")
            representative, member = fields
            if member in mapping and mapping[member] != representative:
                raise ValueError(f"member {member} appears in multiple clusters")
            mapping[member] = representative
    if not mapping:
        raise ValueError("empty MMseqs cluster mapping")
    return mapping


def cluster_summary(mapping: dict[str, str]) -> dict[str, Any]:
    clusters: dict[str, list[str]] = {}
    for member, representative in mapping.items():
        clusters.setdefault(representative, []).append(member)
    sizes = sorted(len(members) for members in clusters.values())
    return {
        "sequence_count": len(mapping),
        "cluster_count": len(clusters),
        "cluster_sizes": sizes,
        "minimum_cluster_size": min(sizes),
        "maximum_cluster_size": max(sizes),
    }


def assign_cluster_splits(
    mapping: dict[str, str],
    *,
    seed: int,
    train_fraction: float = 0.8,
    validation_fraction: float = 0.1,
) -> dict[str, str]:
    if not 0 < train_fraction < 1:
        raise ValueError("train_fraction must be in (0, 1)")
    if not 0 <= validation_fraction < 1 - train_fraction:
        raise ValueError("validation_fraction leaves no test set")
    representatives = sorted(set(mapping.values()))
    cluster_split: dict[str, str] = {}
    for representative in representatives:
        digest = hashlib.sha256(f"{seed}:{representative}".encode()).digest()
        fraction = int.from_bytes(digest[:8], "big") / 2**64
        if fraction < train_fraction:
            split = "train"
        elif fraction < train_fraction + validation_fraction:
            split = "validation"
        else:
            split = "test"
        cluster_split[representative] = split
    assignments = {
        member: cluster_split[representative]
        for member, representative in mapping.items()
    }
    assert_no_cluster_leakage(mapping, assignments)
    return assignments


def assert_no_cluster_leakage(
    mapping: dict[str, str], assignments: dict[str, str]
) -> None:
    missing = set(mapping).difference(assignments)
    if missing:
        raise ValueError(f"split assignments missing {len(missing)} sequence(s)")
    cluster_splits: dict[str, set[str]] = {}
    for member, representative in mapping.items():
        cluster_splits.setdefault(representative, set()).add(assignments[member])
    leaked = {
        representative: sorted(splits)
        for representative, splits in cluster_splits.items()
        if len(splits) != 1
    }
    if leaked:
        raise RuntimeError(
            "DATA LEAKAGE: clusters cross splits; affected="
            + json.dumps(leaked, sort_keys=True)
        )
