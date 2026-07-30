"""Sequence novelty, diversity, and inexpensive physicochemical proxies."""

from __future__ import annotations

import math
from collections import Counter
from itertools import combinations, pairwise

from cas13_if.schemas import STANDARD_AA

HYDROPHOBIC = frozenset("AVILMFWY")
POSITIVE = frozenset("KR")
NEGATIVE = frozenset("DE")


def validate_sequence(sequence: str) -> str:
    normalized = "".join(sequence.split()).upper()
    invalid = sorted(set(normalized).difference(STANDARD_AA))
    if not normalized or invalid:
        raise ValueError(f"invalid protein sequence; invalid={invalid}")
    return normalized


def sequence_identity(first: str, second: str) -> float:
    left = validate_sequence(first)
    right = validate_sequence(second)
    if len(left) != len(right):
        raise ValueError("identity requires equal-length unaligned sequences")
    return sum(a == b for a, b in zip(left, right, strict=True)) / len(left)


def designed_position_identity(
    candidate: str, parent: str, designed_positions: set[int]
) -> float:
    if not designed_positions:
        raise ValueError("designed_positions cannot be empty")
    if len(candidate) != len(parent):
        raise ValueError("candidate and parent lengths differ")
    if min(designed_positions) < 0 or max(designed_positions) >= len(parent):
        raise ValueError("designed position outside sequence")
    return sum(candidate[index] == parent[index] for index in designed_positions) / len(
        designed_positions
    )


def shannon_entropy(sequence: str) -> float:
    normalized = validate_sequence(sequence)
    counts = Counter(normalized)
    return -sum(
        (count / len(normalized)) * math.log(count / len(normalized))
        for count in counts.values()
    )


def longest_homopolymer(sequence: str) -> int:
    normalized = validate_sequence(sequence)
    longest = current = 1
    for previous, token in pairwise(normalized):
        current = current + 1 if token == previous else 1
        longest = max(longest, current)
    return longest


def low_complexity_windows(
    sequence: str,
    *,
    window: int = 12,
    maximum_single_residue_fraction: float = 0.5,
) -> list[tuple[int, int]]:
    normalized = validate_sequence(sequence)
    if window < 2 or window > len(normalized):
        raise ValueError("window must be between 2 and sequence length")
    flagged = []
    for start in range(len(normalized) - window + 1):
        segment = normalized[start : start + window]
        if max(Counter(segment).values()) / window > maximum_single_residue_fraction:
            flagged.append((start, start + window))
    return flagged


def composition_deviation(candidate: str, parent: str) -> float:
    first = Counter(validate_sequence(candidate))
    second = Counter(validate_sequence(parent))
    first_total = sum(first.values())
    second_total = sum(second.values())
    return 0.5 * sum(
        abs(first[amino_acid] / first_total - second[amino_acid] / second_total)
        for amino_acid in STANDARD_AA
    )


def net_charge_proxy(sequence: str) -> float:
    normalized = validate_sequence(sequence)
    return (
        sum(token in POSITIVE for token in normalized)
        - sum(token in NEGATIVE for token in normalized)
    ) / len(normalized)


def hydrophobicity_proxy(sequence: str) -> float:
    normalized = validate_sequence(sequence)
    return sum(token in HYDROPHOBIC for token in normalized) / len(normalized)


def pairwise_candidate_identity(sequences: list[str]) -> list[float]:
    return [
        sequence_identity(left, right) for left, right in combinations(sequences, 2)
    ]
