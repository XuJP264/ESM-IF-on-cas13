"""Deterministic identity matching and seed-level paired statistics."""

from __future__ import annotations

import hashlib
import itertools
import json
from collections import defaultdict
from typing import Any

import numpy as np

from cas13_if.novelty.metrics import designed_position_identity, sequence_identity
from cas13_if.statistics.resampling import benjamini_hochberg, bootstrap_mean


def position_set_hash(positions: set[int]) -> str:
    payload = json.dumps(sorted(positions), separators=(",", ":"))
    return hashlib.sha256(payload.encode()).hexdigest()


def add_identity_metrics(
    rows: list[dict[str, Any]],
    *,
    parent_sequence: str,
    designed_positions: set[int],
) -> list[dict[str, Any]]:
    output: list[dict[str, Any]] = []
    for row in rows:
        sequence = str(row["sequence"])
        output.append(
            {
                **row,
                "parent_identity": sequence_identity(sequence, parent_sequence),
                "designed_position_identity": designed_position_identity(
                    sequence, parent_sequence, designed_positions
                ),
            }
        )
    return output


def select_balanced_candidates(
    proposals: list[dict[str, Any]],
    *,
    methods: list[str],
    seed_blocks: list[int],
    minimum_parent_identity: float,
    maximum_parent_identity: float,
    minimum_designed_identity: float,
    maximum_designed_identity: float,
    target_identity: float,
) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    """Select exactly one row per method/seed from one common identity box."""
    if not (
        0 <= minimum_parent_identity <= maximum_parent_identity <= 1
        and 0 <= minimum_designed_identity <= maximum_designed_identity <= 1
        and 0 <= target_identity <= 1
    ):
        raise ValueError("invalid identity-matching bounds")
    grouped: dict[tuple[str, int], list[dict[str, Any]]] = defaultdict(list)
    failures: list[dict[str, Any]] = []
    for proposal in proposals:
        grouped[(str(proposal["method"]), int(proposal["seed_block"]))].append(proposal)
    selected: list[dict[str, Any]] = []
    for method in methods:
        for seed_block in seed_blocks:
            group = grouped.get((method, seed_block), [])
            eligible: list[dict[str, Any]] = []
            for row in group:
                parent_identity = float(row["parent_identity"])
                designed_identity = float(row["designed_position_identity"])
                within = (
                    minimum_parent_identity
                    <= parent_identity
                    <= maximum_parent_identity
                    and minimum_designed_identity
                    <= designed_identity
                    <= maximum_designed_identity
                )
                if within:
                    eligible.append(row)
                else:
                    failures.append(
                        {
                            "candidate_id": row["candidate_id"],
                            "method": method,
                            "seed_block": seed_block,
                            "stage": "identity_matching",
                            "reason": "outside_common_identity_interval",
                            "parent_identity": parent_identity,
                            "designed_position_identity": designed_identity,
                            "is_mock": bool(row.get("is_mock", False)),
                        }
                    )
            if not eligible:
                raise ValueError(
                    f"no identity-matched proposal for {method}/seed={seed_block}"
                )
            eligible.sort(
                key=lambda row: (
                    abs(float(row["parent_identity"]) - target_identity)
                    + abs(float(row["designed_position_identity"]) - target_identity),
                    str(row["candidate_id"]),
                )
            )
            chosen = eligible[0]
            selected.append(chosen)
            for row in eligible[1:]:
                failures.append(
                    {
                        "candidate_id": row["candidate_id"],
                        "method": method,
                        "seed_block": seed_block,
                        "stage": "identity_matching",
                        "reason": "eligible_not_selected_by_deterministic_tie_break",
                        "parent_identity": row["parent_identity"],
                        "designed_position_identity": row["designed_position_identity"],
                        "is_mock": bool(row.get("is_mock", False)),
                    }
                )
    return selected, failures


def exact_paired_permutation_p_value(differences: list[float]) -> float:
    """Two-sided exact sign-flip p-value for paired independent units."""
    if not differences:
        raise ValueError("paired permutation requires differences")
    observed = abs(float(np.mean(differences)))
    null = [
        abs(
            float(
                np.mean(
                    [
                        difference * sign
                        for difference, sign in zip(differences, signs, strict=True)
                    ]
                )
            )
        )
        for signs in itertools.product((-1.0, 1.0), repeat=len(differences))
    ]
    return sum(value >= observed - 1e-12 for value in null) / len(null)


def paired_seed_statistics(
    rows: list[dict[str, Any]],
    *,
    metrics: list[str],
    reference_method: str,
    bootstrap_replicates: int,
    confidence: float,
    seed: int,
) -> list[dict[str, Any]]:
    """Compare method means using seed blocks, never candidate pseudoreplicates."""
    by_method_seed = {(str(row["method"]), int(row["seed_block"])): row for row in rows}
    methods = sorted({str(row["method"]) for row in rows})
    seed_blocks = sorted({int(row["seed_block"]) for row in rows})
    results: list[dict[str, Any]] = []
    p_values: list[float] = []
    for method in methods:
        if method == reference_method:
            continue
        for metric in metrics:
            pairs = [
                (
                    by_method_seed.get((method, seed_block)),
                    by_method_seed.get((reference_method, seed_block)),
                )
                for seed_block in seed_blocks
            ]
            usable = [
                (left, right)
                for left, right in pairs
                if left is not None
                and right is not None
                and left.get(metric) is not None
                and right.get(metric) is not None
            ]
            if not usable:
                continue
            differences = [
                float(left[metric]) - float(right[metric]) for left, right in usable
            ]
            interval = bootstrap_mean(
                differences,
                replicates=bootstrap_replicates,
                confidence=confidence,
                seed=seed,
            )
            p_value = exact_paired_permutation_p_value(differences)
            p_values.append(p_value)
            results.append(
                {
                    "method": method,
                    "reference_method": reference_method,
                    "metric": metric,
                    "independent_unit": "seed",
                    "independent_unit_count": len(differences),
                    "paired_mean_difference": interval.estimate,
                    "bootstrap_ci_lower": interval.lower,
                    "bootstrap_ci_upper": interval.upper,
                    "bootstrap_replicates": bootstrap_replicates,
                    "confidence": confidence,
                    "unadjusted_permutation_p_value": p_value,
                    "inference_label": "low_power_descriptive",
                }
            )
    adjusted = benjamini_hochberg(p_values)
    for row, value in zip(results, adjusted, strict=True):
        row["benjamini_hochberg_p_value"] = value
    return results
