#!/usr/bin/env python
"""Fail closed when the canonical matched-baseline handoff is inconsistent."""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

METHODS = {
    "matched_random_mutation",
    "msa_profile_sampling",
    "unconstrained_esm_if1",
    "catalytic_only_fixed_esm_if1",
    "conservation_constrained_esm_if1",
    "conservation_rna_contact_esm_if1",
    "proteinmpnn",
    "ligandmpnn",
    "esm_if1_ligandmpnn_consensus",
}


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--report-dir", type=Path, default=Path("reports/matched_baselines")
    )
    return parser.parse_args()


def _json(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"expected JSON object: {path}")
    return value


def _jsonl(path: Path) -> list[dict[str, Any]]:
    rows = [json.loads(line) for line in path.read_text(encoding="utf-8").splitlines()]
    if not all(isinstance(row, dict) for row in rows):
        raise ValueError(f"expected JSONL objects: {path}")
    return rows


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def verify(report_dir: Path) -> dict[str, Any]:
    required = {
        "methods_table.csv",
        "candidates.jsonl",
        "candidate_funnel.csv",
        "matched_statistics.csv",
        "per_region_metrics.csv",
        "failure_analysis.csv",
        "report.md",
        "report.html",
        "summary.json",
        "gpu_hpc_job_manifest.jsonl",
    }
    missing = sorted(name for name in required if not (report_dir / name).is_file())
    if missing:
        raise ValueError(f"missing canonical report files: {missing}")

    summary = _json(report_dir / "summary.json")
    candidates = _jsonl(report_dir / "candidates.jsonl")
    methods = _csv(report_dir / "methods_table.csv")
    statistics = _csv(report_dir / "matched_statistics.csv")
    regions = _csv(report_dir / "per_region_metrics.csv")
    failures = _csv(report_dir / "failure_analysis.csv")
    gpu_jobs = _jsonl(report_dir / "gpu_hpc_job_manifest.jsonl")

    if summary.get("is_mock") is not False or int(summary["mock_candidates"]) != 0:
        raise ValueError("formal summary contains mock output")
    if summary.get("identity_matching_passed") is not True:
        raise ValueError("identity matching did not pass")
    if int(summary["fixed_position_violations"]) != 0:
        raise ValueError("formal summary contains fixed-position violations")
    if len(candidates) != 18 or len(methods) != 9 or len(regions) != 18:
        raise ValueError("canonical 9-method x 2-seed matrix is incomplete")
    if {row["method"] for row in methods} != METHODS:
        raise ValueError("method table does not contain the registered methods")

    counts = Counter(str(row["method"]) for row in candidates)
    seeds: dict[str, set[int]] = defaultdict(set)
    fixed_hashes: set[str] = set()
    free_hashes: set[str] = set()
    for row in candidates:
        method = str(row["method"])
        seeds[method].add(int(row["seed_block"]))
        fixed_hashes.add(str(row["fixed_position_hash"]))
        free_hashes.add(str(row["free_position_hash"]))
        if bool(row["is_mock"]):
            raise ValueError(f"mock candidate: {row['candidate_id']}")
        if int(row["fixed_position_violations"]) != 0:
            raise ValueError(f"fixed-position violation: {row['candidate_id']}")
        parent_identity = float(row["parent_identity"])
        designed_identity = float(row["designed_position_identity"])
        if not 0.18 <= parent_identity <= 0.32:
            raise ValueError(f"parent identity outside registered interval: {method}")
        if not 0.18 <= designed_identity <= 0.32:
            raise ValueError(f"designed identity outside registered interval: {method}")
    if set(counts) != METHODS or set(counts.values()) != {2}:
        raise ValueError(f"unbalanced candidate counts: {dict(counts)}")
    if len({tuple(sorted(value)) for value in seeds.values()}) != 1:
        raise ValueError("methods do not share the same seed blocks")
    if len(fixed_hashes) != 1 or len(free_hashes) != 1:
        raise ValueError("methods do not share one fixed/free position mask")
    if fixed_hashes != {str(summary["fixed_position_hash"])}:
        raise ValueError("candidate and summary fixed-position hashes differ")
    if free_hashes != {str(summary["free_position_hash"])}:
        raise ValueError("candidate and summary free-position hashes differ")
    if not statistics or {row["independent_unit"] for row in statistics} != {"seed"}:
        raise ValueError("statistics are missing or use the wrong independent unit")

    report = (report_dir / "report.md").read_text(encoding="utf-8")
    for phrase in (
        "Level 1",
        "Level 2",
        "No mock candidate",
        "no sequence is described as an effective or validated Cas13",
    ):
        if phrase not in report:
            raise ValueError(f"report evidence boundary is missing: {phrase}")

    if len(gpu_jobs) != 1:
        raise ValueError("expected exactly one deferred GPU job")
    gpu = gpu_jobs[0]
    if gpu.get("status") != "not_run" or gpu.get("is_mock") is not False:
        raise ValueError("GPU extension must remain genuine but not_run")
    if gpu.get("expected_fixed_position_hash") != next(iter(fixed_hashes)):
        raise ValueError("GPU fixed-position hash differs")
    if gpu.get("expected_free_position_hash") != next(iter(free_hashes)):
        raise ValueError("GPU free-position hash differs")
    if not all(
        str(path).startswith("reports/matched_baselines_gpu/")
        for path in gpu["expected_outputs"]
    ):
        raise ValueError("GPU outputs are not isolated from the CPU report")

    return {
        "candidate_count": len(candidates),
        "method_count": len(methods),
        "seed_blocks": sorted(next(iter(seeds.values()))),
        "fixed_position_hash": next(iter(fixed_hashes)),
        "free_position_hash": next(iter(free_hashes)),
        "fixed_position_violations": 0,
        "mock_candidates": 0,
        "level1_pass_count": int(summary["level1_pass_count"]),
        "level2_scored_count": int(summary["level2_scored_count"]),
        "statistics_rows": len(statistics),
        "failure_rows": len(failures),
        "gpu_status": str(gpu["status"]),
    }


def main() -> int:
    result = verify(_arguments().report_dir.resolve())
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
