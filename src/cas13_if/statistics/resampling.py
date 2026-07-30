"""Seeded paired and hierarchical resampling helpers."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass(frozen=True)
class ConfidenceInterval:
    estimate: float
    lower: float
    upper: float
    confidence: float
    replicates: int
    seed: int


def bootstrap_mean(
    values: list[float],
    *,
    replicates: int = 1000,
    confidence: float = 0.95,
    seed: int = 20260731,
) -> ConfidenceInterval:
    if not values:
        raise ValueError("bootstrap requires values at independent-unit level")
    if replicates < 100:
        raise ValueError("at least 100 bootstrap replicates are required")
    if not 0 < confidence < 1:
        raise ValueError("confidence must be in (0, 1)")
    observations = np.asarray(values, dtype=np.float64)
    rng = np.random.default_rng(seed)
    samples = rng.choice(
        observations, size=(replicates, len(observations)), replace=True
    )
    means = samples.mean(axis=1)
    tail = (1 - confidence) / 2
    return ConfidenceInterval(
        estimate=float(observations.mean()),
        lower=float(np.quantile(means, tail)),
        upper=float(np.quantile(means, 1 - tail)),
        confidence=confidence,
        replicates=replicates,
        seed=seed,
    )


def paired_effect(left: list[float], right: list[float]) -> float:
    if len(left) != len(right) or not left:
        raise ValueError("paired effects require equal nonempty independent units")
    return float(np.mean(np.asarray(left) - np.asarray(right)))


def benjamini_hochberg(p_values: list[float]) -> list[float]:
    if any(not 0 <= value <= 1 for value in p_values):
        raise ValueError("p-values must be in [0, 1]")
    count = len(p_values)
    if count == 0:
        return []
    order = np.argsort(p_values)
    adjusted = np.empty(count, dtype=np.float64)
    previous = 1.0
    for reverse_rank, index in enumerate(order[::-1], start=1):
        rank = count - reverse_rank + 1
        value = min(previous, p_values[int(index)] * count / rank, 1.0)
        adjusted[int(index)] = value
        previous = value
    return [float(value) for value in adjusted]
