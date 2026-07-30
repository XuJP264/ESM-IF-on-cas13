"""Validated MSA parsing, sequence reweighting, and scaffold-column mapping."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from cas13_if.data.fasta import iter_fasta
from cas13_if.schemas import STANDARD_AA


class AlignmentError(ValueError):
    """Raised when an input is not a valid aligned FASTA."""


@dataclass(frozen=True)
class Alignment:
    identifiers: tuple[str, ...]
    sequences: tuple[str, ...]

    @property
    def n_sequences(self) -> int:
        return len(self.sequences)

    @property
    def n_columns(self) -> int:
        return len(self.sequences[0])

    def matrix(self) -> np.ndarray:
        return np.array([list(sequence) for sequence in self.sequences], dtype="U1")


def read_aligned_fasta(path: Path, *, alphabet: str = "protein") -> Alignment:
    records = list(iter_fasta(path))
    if len(records) < 2:
        raise AlignmentError("an MSA requires at least two sequences")
    lengths = {len(sequence) for _, sequence in records}
    if len(lengths) != 1:
        raise AlignmentError(
            "aligned FASTA rows must have identical lengths; unaligned FASTA rejected"
        )
    if lengths == {0}:
        raise AlignmentError("MSA rows cannot be empty")
    if alphabet == "protein":
        permitted = STANDARD_AA.union({"-", "."})
    elif alphabet == "rna":
        permitted = frozenset("ACGUTN-.")
    elif alphabet == "mixed":
        permitted = STANDARD_AA.union(set("U-."))
    else:
        raise AlignmentError(f"unknown MSA alphabet: {alphabet}")
    invalid = sorted(
        set().union(*(set(sequence) for _, sequence in records)).difference(permitted)
    )
    if invalid:
        raise AlignmentError(f"invalid {alphabet} MSA symbols: {invalid}")
    return Alignment(
        identifiers=tuple(identifier for identifier, _ in records),
        sequences=tuple(sequence.replace(".", "-") for _, sequence in records),
    )


def pairwise_identity(first: str, second: str) -> float:
    if len(first) != len(second):
        raise AlignmentError("identity requires aligned rows of equal length")
    comparable = [
        (left, right)
        for left, right in zip(first, second, strict=True)
        if left != "-" and right != "-"
    ]
    if not comparable:
        return 0.0
    return sum(left == right for left, right in comparable) / len(comparable)


def sequence_weights(
    alignment: Alignment,
    *,
    identity_threshold: float = 0.8,
) -> np.ndarray:
    """Return 1/neighborhood-size weights under aligned sequence identity."""
    if not 0 < identity_threshold <= 1:
        raise ValueError("identity threshold must be in (0, 1]")
    counts = np.ones(alignment.n_sequences, dtype=np.float64)
    for left in range(alignment.n_sequences):
        for right in range(left + 1, alignment.n_sequences):
            if (
                pairwise_identity(alignment.sequences[left], alignment.sequences[right])
                >= identity_threshold
            ):
                counts[left] += 1
                counts[right] += 1
    return 1.0 / counts


def map_ungapped_sequence_to_columns(
    aligned_reference: str, ungapped_sequence: str
) -> dict[int, int]:
    """Map zero-based biological indices to MSA columns, handling insertions/gaps."""
    normalized = ungapped_sequence.replace("-", "").upper()
    if aligned_reference.replace("-", "").upper() != normalized:
        raise AlignmentError("aligned reference does not match ungapped scaffold")
    mapping: dict[int, int] = {}
    biological_index = 0
    for column, token in enumerate(aligned_reference):
        if token != "-":
            mapping[biological_index] = column
            biological_index += 1
    return mapping
