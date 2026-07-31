"""Real candidate QC and Atlas similarity search orchestration."""

from __future__ import annotations

import json
import statistics
import subprocess
from collections import Counter, defaultdict
from dataclasses import dataclass
from itertools import combinations
from pathlib import Path
from typing import Any

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

from cas13_if.data.fasta import write_fasta
from cas13_if.novelty.metrics import (
    composition_deviation,
    designed_position_identity,
    hydrophobicity_proxy,
    longest_homopolymer,
    low_complexity_windows,
    net_charge_proxy,
    sequence_identity,
    shannon_entropy,
    validate_sequence,
)
from cas13_if.provenance import atomic_write_text, sha256_file
from cas13_if.schemas import STANDARD_AA


@dataclass(frozen=True)
class NoveltyThresholds:
    maximum_parent_identity: float
    maximum_atlas_identity: float
    maximum_homopolymer_length: int
    maximum_low_complexity_windows: int
    minimum_designed_position_entropy: float
    low_complexity_window: int
    low_complexity_maximum_fraction: float


def load_benchmark_candidates(path: Path) -> list[dict[str, Any]]:
    """Read only compact fields from a trace-heavy benchmark JSONL."""
    rows: list[dict[str, Any]] = []
    seen: set[str] = set()
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip():
                continue
            payload = json.loads(line)
            candidate = payload.get("candidate")
            recovery = payload.get("recovery")
            if not isinstance(candidate, dict) or not isinstance(recovery, dict):
                raise ValueError(f"invalid candidate record at line {line_number}")
            candidate_id = str(candidate.get("candidate_id", ""))
            if not candidate_id or candidate_id in seen:
                raise ValueError(f"empty or duplicate candidate ID: {candidate_id!r}")
            seen.add(candidate_id)
            fixed_raw = candidate.get("fixed_positions", {})
            if not isinstance(fixed_raw, dict):
                raise ValueError(f"fixed_positions is not a mapping: {candidate_id}")
            fixed_positions = sorted(int(index) for index in fixed_raw)
            rows.append(
                {
                    "candidate_id": candidate_id,
                    "pdb_id": str(payload.get("pdb_id", "")),
                    "method": str(payload.get("method", "")),
                    "scaffold_id": str(candidate.get("scaffold_id", "")),
                    "sequence": str(candidate.get("sequence", "")),
                    "parent_sequence": str(candidate.get("parent_sequence", "")),
                    "temperature": float(candidate.get("temperature", 0.0)),
                    "seed": int(candidate.get("seed", 0)),
                    "fixed_positions": fixed_positions,
                    "fixed_position_violations": int(
                        recovery.get("fixed_position_violations", 0)
                    ),
                    "source_is_mock": bool(candidate.get("is_mock", False)),
                    "source_evidence_level": int(candidate.get("evidence_level", 0)),
                }
            )
    if not rows:
        raise ValueError(f"candidate JSONL is empty: {path}")
    if any(row["source_is_mock"] for row in rows):
        raise ValueError("real novelty analysis refuses mock candidates")
    return rows


def run_mmseqs_atlas_search(
    *,
    candidates: list[dict[str, Any]],
    atlas_fasta: Path,
    output_dir: Path,
    executable: Path,
    threads: int,
    sensitivity: float,
    minimum_query_coverage: float,
    maximum_evalue: float,
    maximum_sequences: int,
) -> tuple[dict[str, dict[str, Any]], list[str]]:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite novelty output: {output_dir}")
    if not executable.is_file():
        raise FileNotFoundError(f"MMseqs2 executable is missing: {executable}")
    if not atlas_fasta.is_file():
        raise FileNotFoundError(f"Atlas FASTA is missing: {atlas_fasta}")
    if threads < 1 or maximum_sequences < 1:
        raise ValueError("MMseqs2 threads and maximum sequences must be positive")
    if not 0 <= minimum_query_coverage <= 1:
        raise ValueError("minimum query coverage must be in [0, 1]")
    output_dir.mkdir(parents=True, exist_ok=False)
    query_fasta = output_dir / "candidates.fasta"
    write_fasta(
        ((str(row["candidate_id"]), str(row["sequence"])) for row in candidates),
        query_fasta,
    )
    alignment_path = output_dir / "atlas_alignments.tsv"
    temporary_dir = output_dir / "mmseqs_tmp"
    command = [
        str(executable),
        "easy-search",
        str(query_fasta),
        str(atlas_fasta),
        str(alignment_path),
        str(temporary_dir),
        "--format-output",
        "query,target,fident,alnlen,qcov,tcov,evalue,bits",
        "--threads",
        str(threads),
        "-s",
        str(sensitivity),
        "-e",
        str(maximum_evalue),
        "--cov-mode",
        "2",
        "-c",
        str(minimum_query_coverage),
        "--max-seqs",
        str(maximum_sequences),
        "--max-accept",
        str(maximum_sequences),
        "--sort-results",
        "1",
        "--remove-tmp-files",
        "1",
    ]
    completed = subprocess.run(command, text=True, capture_output=True, check=False)
    atomic_write_text(output_dir / "mmseqs_stdout.log", completed.stdout)
    atomic_write_text(output_dir / "mmseqs_stderr.log", completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError("MMseqs2 novelty search failed; see mmseqs_stderr.log")
    expected_ids = {str(row["candidate_id"]) for row in candidates}
    best: dict[str, dict[str, Any]] = {}
    if alignment_path.is_file():
        with alignment_path.open(encoding="utf-8") as handle:
            for line_number, line in enumerate(handle, start=1):
                fields = line.rstrip("\n").split("\t")
                if len(fields) != 8:
                    message = (
                        f"invalid MMseqs2 result at line {line_number}: "
                        f"{len(fields)} fields"
                    )
                    raise ValueError(message)
                query, target = fields[:2]
                if query not in expected_ids:
                    raise ValueError(f"unexpected MMseqs2 query ID: {query}")
                hit: dict[str, Any] = {
                    "target_sequence_sha256": target,
                    "identity": float(fields[2]),
                    "alignment_length": int(fields[3]),
                    "query_coverage": float(fields[4]),
                    "target_coverage": float(fields[5]),
                    "evalue": float(fields[6]),
                    "bits": float(fields[7]),
                }
                if not 0 <= hit["identity"] <= 1:
                    raise ValueError(
                        f"MMseqs2 fident is not a fraction in [0, 1]: {fields[2]}"
                    )
                previous = best.get(query)
                ranking = (
                    float(hit["identity"]),
                    float(hit["query_coverage"]),
                    float(hit["bits"]),
                )
                previous_ranking = (
                    (
                        float(previous["identity"]),
                        float(previous["query_coverage"]),
                        float(previous["bits"]),
                    )
                    if previous is not None
                    else (-1.0, -1.0, -1.0)
                )
                if ranking > previous_ranking:
                    best[query] = hit
    return best, command


def evaluate_candidate_novelty(
    candidates: list[dict[str, Any]],
    atlas_hits: dict[str, dict[str, Any]],
    thresholds: NoveltyThresholds,
) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    evaluated: list[dict[str, Any]] = []
    failure_counts: Counter[str] = Counter()
    sequences_by_scaffold: dict[str, list[str]] = defaultdict(list)
    for candidate in candidates:
        candidate_id = str(candidate["candidate_id"])
        sequence = str(candidate["sequence"])
        parent = str(candidate["parent_sequence"])
        invalid = sorted(set(sequence.upper()).difference(STANDARD_AA))
        if invalid:
            raise ValueError(f"candidate {candidate_id} has invalid symbols: {invalid}")
        sequence = validate_sequence(sequence)
        parent = validate_sequence(parent)
        fixed_positions = {int(index) for index in candidate["fixed_positions"]}
        designed_positions = set(range(len(sequence))).difference(fixed_positions)
        if not designed_positions:
            raise ValueError(f"candidate has no designed positions: {candidate_id}")
        designed_sequence = "".join(
            sequence[index] for index in sorted(designed_positions)
        )
        low_complexity = low_complexity_windows(
            sequence,
            window=thresholds.low_complexity_window,
            maximum_single_residue_fraction=(
                thresholds.low_complexity_maximum_fraction
            ),
        )
        hit = atlas_hits.get(candidate_id)
        atlas_identity = float(hit["identity"]) if hit is not None else None
        parent_identity = sequence_identity(sequence, parent)
        designed_identity = designed_position_identity(
            sequence, parent, designed_positions
        )
        designed_entropy = shannon_entropy(designed_sequence)
        homopolymer = longest_homopolymer(sequence)
        failures: list[str] = []
        if int(candidate["fixed_position_violations"]) != 0:
            failures.append("fixed_position_violation")
        if parent_identity > thresholds.maximum_parent_identity:
            failures.append("parent_identity_above_threshold")
        if atlas_identity is None:
            failures.append("no_atlas_hit_at_required_query_coverage")
        elif atlas_identity > thresholds.maximum_atlas_identity:
            failures.append("atlas_identity_above_threshold")
        if homopolymer > thresholds.maximum_homopolymer_length:
            failures.append("homopolymer_above_threshold")
        if len(low_complexity) > thresholds.maximum_low_complexity_windows:
            failures.append("low_complexity_windows_above_threshold")
        if designed_entropy < thresholds.minimum_designed_position_entropy:
            failures.append("designed_entropy_below_threshold")
        failure_counts.update(failures)
        passes = not failures
        evaluated.append(
            {
                **candidate,
                "sequence_length": len(sequence),
                "parent_identity": parent_identity,
                "designed_position_identity": designed_identity,
                "maximum_atlas_identity": atlas_identity,
                "maximum_atlas_hit": (
                    hit["target_sequence_sha256"] if hit is not None else None
                ),
                "maximum_atlas_query_coverage": (
                    hit["query_coverage"] if hit is not None else None
                ),
                "maximum_atlas_target_coverage": (
                    hit["target_coverage"] if hit is not None else None
                ),
                "shannon_entropy": shannon_entropy(sequence),
                "designed_position_entropy": designed_entropy,
                "longest_homopolymer": homopolymer,
                "low_complexity_window_count": len(low_complexity),
                "composition_deviation": composition_deviation(sequence, parent),
                "net_charge_proxy": net_charge_proxy(sequence),
                "hydrophobicity_proxy": hydrophobicity_proxy(sequence),
                "novelty_filter_failures": failures,
                "passes_level1_novelty": passes,
                "is_mock": False,
                "evidence_level": 1 if passes else 0,
            }
        )
        sequences_by_scaffold[str(candidate["scaffold_id"])].append(sequence)
    diversity: dict[str, Any] = {}
    for scaffold, sequences in sorted(sequences_by_scaffold.items()):
        identities = [
            sequence_identity(left, right) for left, right in combinations(sequences, 2)
        ]
        diversity[scaffold] = {
            "candidate_count": len(sequences),
            "pair_count": len(identities),
            "minimum_pairwise_identity": min(identities) if identities else None,
            "mean_pairwise_identity": (
                statistics.fmean(identities) if identities else None
            ),
            "maximum_pairwise_identity": max(identities) if identities else None,
        }
    passing = sum(bool(row["passes_level1_novelty"]) for row in evaluated)
    atlas_identities = [
        float(row["maximum_atlas_identity"])
        for row in evaluated
        if row["maximum_atlas_identity"] is not None
    ]
    summary = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level_max": 1 if passing else 0,
        "candidate_count": len(evaluated),
        "candidates_with_atlas_hit": len(atlas_identities),
        "passes_level1_novelty": passing,
        "filter_failure_counts": dict(sorted(failure_counts.items())),
        "maximum_observed_atlas_identity": (
            max(atlas_identities) if atlas_identities else None
        ),
        "minimum_observed_parent_identity": min(
            float(row["parent_identity"]) for row in evaluated
        ),
        "maximum_observed_parent_identity": max(
            float(row["parent_identity"]) for row in evaluated
        ),
        "diversity_by_scaffold": diversity,
        "thresholds": {
            "maximum_parent_identity": thresholds.maximum_parent_identity,
            "maximum_atlas_identity": thresholds.maximum_atlas_identity,
            "maximum_homopolymer_length": thresholds.maximum_homopolymer_length,
            "maximum_low_complexity_windows": (
                thresholds.maximum_low_complexity_windows
            ),
            "minimum_designed_position_entropy": (
                thresholds.minimum_designed_position_entropy
            ),
            "low_complexity_window": thresholds.low_complexity_window,
            "low_complexity_maximum_fraction": (
                thresholds.low_complexity_maximum_fraction
            ),
        },
        "claim_scope": (
            "Level 1 sequence-statistical novelty only for rows that pass every "
            "registered filter; no candidate is a validated or effective Cas13."
        ),
    }
    return evaluated, summary


def run_candidate_novelty_pipeline(
    *,
    candidate_jsonl: Path,
    atlas_fasta: Path,
    output_dir: Path,
    executable: Path,
    threads: int,
    sensitivity: float,
    minimum_query_coverage: float,
    maximum_evalue: float,
    maximum_sequences: int,
    thresholds: NoveltyThresholds,
) -> dict[str, Any]:
    candidates = load_benchmark_candidates(candidate_jsonl)
    hits, command = run_mmseqs_atlas_search(
        candidates=candidates,
        atlas_fasta=atlas_fasta,
        output_dir=output_dir,
        executable=executable,
        threads=threads,
        sensitivity=sensitivity,
        minimum_query_coverage=minimum_query_coverage,
        maximum_evalue=maximum_evalue,
        maximum_sequences=maximum_sequences,
    )
    rows, summary = evaluate_candidate_novelty(candidates, hits, thresholds)
    result_path = output_dir / "candidate_novelty.parquet"
    pq.write_table(pa.Table.from_pylist(rows), result_path, compression="zstd")
    summary.update(
        {
            "candidate_source": str(candidate_jsonl),
            "candidate_source_sha256": sha256_file(candidate_jsonl),
            "atlas_fasta": str(atlas_fasta),
            "atlas_fasta_sha256": sha256_file(atlas_fasta),
            "mmseqs_command": command,
            "mmseqs_result": str(output_dir / "atlas_alignments.tsv"),
        }
    )
    atomic_write_text(
        output_dir / "summary.json",
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    return summary
