"""Subtype-specific weighted conservation and entropy."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from cas13_if.alignments.msa import Alignment, sequence_weights
from cas13_if.schemas import STANDARD_AA

AA_ORDER = tuple(sorted(STANDARD_AA))


@dataclass(frozen=True)
class ColumnStatistics:
    column: int
    conservation: float
    entropy: float
    gap_fraction: float
    effective_sequence_count: float
    consensus: str | None
    allowed_residues: tuple[str, ...]
    weighted_frequencies: dict[str, float]
    unweighted_frequencies: dict[str, float]
    coverage: float


def conservation_statistics(
    alignment: Alignment,
    *,
    identity_threshold: float = 0.8,
    allowed_frequency: float = 0.05,
) -> list[ColumnStatistics]:
    if not 0 <= allowed_frequency <= 1:
        raise ValueError("allowed_frequency must be in [0, 1]")
    matrix = alignment.matrix()
    weights = sequence_weights(alignment, identity_threshold=identity_threshold)
    total_weight = float(weights.sum())
    statistics: list[ColumnStatistics] = []
    for column in range(alignment.n_columns):
        tokens = matrix[:, column]
        nongap = tokens != "-"
        nongap_weight = float(weights[nongap].sum())
        weighted: dict[str, float] = {}
        unweighted: dict[str, float] = {}
        for amino_acid in AA_ORDER:
            mask = tokens == amino_acid
            weighted[amino_acid] = (
                float(weights[mask].sum()) / nongap_weight if nongap_weight else 0.0
            )
            unweighted[amino_acid] = (
                float(mask.sum()) / int(nongap.sum()) if nongap.any() else 0.0
            )
        probabilities = np.array(list(weighted.values()), dtype=np.float64)
        positive = probabilities[probabilities > 0]
        entropy = float(-(positive * np.log(positive)).sum())
        consensus = (
            max(weighted, key=weighted.__getitem__) if nongap_weight > 0 else None
        )
        allowed = tuple(
            amino_acid
            for amino_acid, frequency in weighted.items()
            if frequency >= allowed_frequency and frequency > 0
        )
        statistics.append(
            ColumnStatistics(
                column=column,
                conservation=float(probabilities.max(initial=0.0)),
                entropy=entropy,
                gap_fraction=1.0 - float(nongap.mean()),
                effective_sequence_count=total_weight,
                consensus=consensus,
                allowed_residues=allowed,
                weighted_frequencies=weighted,
                unweighted_frequencies=unweighted,
                coverage=float(nongap.mean()),
            )
        )
    return statistics
