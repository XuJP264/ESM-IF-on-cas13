"""Deterministic hard/soft/free constraint merging with reason retention."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

ConstraintClass = Literal["hard_fixed", "soft_constrained", "free"]


@dataclass(frozen=True)
class PositionEvidence:
    index: int
    source: str
    proposed_class: ConstraintClass
    reason: str
    manually_confirmed: bool = False


@dataclass(frozen=True)
class MergedConstraint:
    index: int
    final_class: ConstraintClass
    decision_reasons: tuple[str, ...]
    sources: tuple[str, ...]


def merge_position_evidence(
    *,
    length: int,
    evidence: list[PositionEvidence],
) -> list[MergedConstraint]:
    if length < 1:
        raise ValueError("length must be positive")
    grouped: dict[int, list[PositionEvidence]] = {index: [] for index in range(length)}
    for item in evidence:
        if item.index not in grouped:
            raise ValueError(f"constraint position {item.index} outside [0, {length})")
        if item.source == "hepn_regex" and item.proposed_class == "hard_fixed":
            if not item.manually_confirmed:
                raise ValueError(
                    "HEPN regex candidates require manual confirmation "
                    "before hard fixing"
                )
        grouped[item.index].append(item)
    rank: dict[ConstraintClass, int] = {
        "free": 0,
        "soft_constrained": 1,
        "hard_fixed": 2,
    }
    output: list[MergedConstraint] = []
    for index, items in grouped.items():
        final_class: ConstraintClass = "free"
        if items:
            final_class = max(
                (item.proposed_class for item in items),
                key=lambda item_class: rank[item_class],
            )
        output.append(
            MergedConstraint(
                index=index,
                final_class=final_class,
                decision_reasons=tuple(
                    sorted({item.reason for item in items})
                    or ["no_constraint_evidence"]
                ),
                sources=tuple(sorted({item.source for item in items})),
            )
        )
    return output
