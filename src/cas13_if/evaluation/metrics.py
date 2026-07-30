"""Matched sequence and structural-region evaluation metrics."""

from __future__ import annotations

import math
from dataclasses import dataclass


@dataclass(frozen=True)
class RecoveryMetrics:
    overall: float
    designed_positions: float | None
    regions: dict[str, float | None]
    fixed_position_violations: int


def native_recovery(
    candidate: str,
    native: str,
    *,
    designed_positions: set[int] | None = None,
    fixed_positions: dict[int, str] | None = None,
    regions: dict[str, set[int]] | None = None,
) -> RecoveryMetrics:
    if len(candidate) != len(native) or not native:
        raise ValueError("candidate and nonempty native sequence lengths must match")
    length = len(native)
    _validate_positions(designed_positions or set(), length)
    for positions in (regions or {}).values():
        _validate_positions(positions, length)
    recovery = sum(left == right for left, right in zip(candidate, native, strict=True))
    designed = None
    if designed_positions:
        designed = sum(
            candidate[index] == native[index] for index in designed_positions
        ) / len(designed_positions)
    region_values = {
        name: (
            sum(candidate[index] == native[index] for index in positions)
            / len(positions)
            if positions
            else None
        )
        for name, positions in (regions or {}).items()
    }
    fixed = fixed_positions or {}
    violations = sum(
        index >= length or candidate[index] != token for index, token in fixed.items()
    )
    return RecoveryMetrics(
        overall=recovery / length,
        designed_positions=designed,
        regions=region_values,
        fixed_position_violations=violations,
    )


def perplexity(per_residue_log_probabilities: list[float]) -> float:
    if not per_residue_log_probabilities:
        raise ValueError("per-residue log probabilities cannot be empty")
    mean_log_probability = sum(per_residue_log_probabilities) / len(
        per_residue_log_probabilities
    )
    return math.exp(-mean_log_probability)


def _validate_positions(positions: set[int], length: int) -> None:
    if positions and (min(positions) < 0 or max(positions) >= length):
        raise ValueError("evaluation position outside sequence")
