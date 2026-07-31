"""Multi-state score aggregation and mapping-safe structural masks."""

from __future__ import annotations

import math
from collections.abc import Iterable, Mapping, Sequence
from dataclasses import dataclass
from itertools import combinations


@dataclass(frozen=True)
class MultiStateScore:
    single_state_score: dict[str, float]
    normalized_weights: dict[str, float]
    multi_state_mean_score: float
    multi_state_min_score: float
    multi_state_variance: float


@dataclass(frozen=True)
class StateResidue:
    full_index_0: int
    wild_type: str
    mapping_status: str


@dataclass(frozen=True)
class StateMask:
    state: str
    sequence: str
    hard_positions: frozenset[int]
    risk_positions: frozenset[int]


@dataclass(frozen=True)
class MultiStateMasks:
    state_intersection_hard_mask: frozenset[int]
    state_union_risk_mask: frozenset[int]
    state_variable_hinge_mask: frozenset[int]


def single_state_score(score: float) -> float:
    """Validate and return one state score."""
    value = float(score)
    if not math.isfinite(value):
        raise ValueError("state score must be finite")
    return value


def normalize_state_weights(
    states: Iterable[str], weights: Mapping[str, float] | None = None
) -> dict[str, float]:
    state_list = list(states)
    if not state_list or len(state_list) != len(set(state_list)):
        raise ValueError("states must be non-empty and unique")
    raw = (
        {state: 1.0 for state in state_list}
        if weights is None
        else {state: float(weights[state]) for state in state_list}
    )
    if weights is not None and set(weights) != set(state_list):
        raise ValueError("state weight keys must exactly match score states")
    if any(not math.isfinite(value) or value < 0 for value in raw.values()):
        raise ValueError("state weights must be finite and non-negative")
    total = sum(raw.values())
    if total <= 0:
        raise ValueError("at least one state weight must be positive")
    return {state: raw[state] / total for state in state_list}


def aggregate_multistate_scores(
    scores: Mapping[str, float],
    *,
    weights: Mapping[str, float] | None = None,
    required_states: set[str] | None = None,
) -> MultiStateScore:
    if not scores:
        raise ValueError("multi-state aggregation requires at least one state")
    if required_states is not None and set(scores) != required_states:
        missing = required_states.difference(scores)
        unexpected = set(scores).difference(required_states)
        raise ValueError(
            f"state set mismatch: missing={sorted(missing)} "
            f"unexpected={sorted(unexpected)}"
        )
    clean = {state: single_state_score(value) for state, value in scores.items()}
    normalized = normalize_state_weights(clean, weights)
    mean = sum(normalized[state] * clean[state] for state in clean)
    variance = sum(normalized[state] * (clean[state] - mean) ** 2 for state in clean)
    return MultiStateScore(
        single_state_score=clean,
        normalized_weights=normalized,
        multi_state_mean_score=mean,
        multi_state_min_score=min(clean.values()),
        multi_state_variance=variance,
    )


def _rank(values: Sequence[float]) -> list[float]:
    order = sorted(range(len(values)), key=lambda index: values[index])
    ranks = [0.0] * len(values)
    start = 0
    while start < len(order):
        end = start + 1
        while end < len(order) and values[order[end]] == values[order[start]]:
            end += 1
        average_rank = (start + end - 1) / 2.0
        for offset in range(start, end):
            ranks[order[offset]] = average_rank
        start = end
    return ranks


def _correlation(left: Sequence[float], right: Sequence[float]) -> float:
    left_mean = sum(left) / len(left)
    right_mean = sum(right) / len(right)
    numerator = sum(
        (a - left_mean) * (b - right_mean) for a, b in zip(left, right, strict=True)
    )
    left_scale = math.sqrt(sum((item - left_mean) ** 2 for item in left))
    right_scale = math.sqrt(sum((item - right_mean) ** 2 for item in right))
    if left_scale == 0 or right_scale == 0:
        return 1.0 if list(left) == list(right) else 0.0
    return numerator / (left_scale * right_scale)


def state_rank_consistency(
    scores_by_candidate: Mapping[str, Mapping[str, float]],
) -> float:
    """Mean pairwise Spearman rank correlation across structural states."""
    if len(scores_by_candidate) < 2:
        raise ValueError("rank consistency requires at least two candidates")
    candidates = sorted(scores_by_candidate)
    state_sets = [set(scores_by_candidate[item]) for item in candidates]
    if not state_sets or any(states != state_sets[0] for states in state_sets[1:]):
        raise ValueError("every candidate must have exactly the same states")
    states = sorted(state_sets[0])
    if len(states) < 2:
        raise ValueError("rank consistency requires at least two states")
    ranks = {
        state: _rank(
            [
                single_state_score(scores_by_candidate[item][state])
                for item in candidates
            ]
        )
        for state in states
    }
    correlations = [
        _correlation(ranks[left], ranks[right])
        for left, right in combinations(states, 2)
    ]
    return sum(correlations) / len(correlations)


def select_state_combination(
    available: Mapping[str, str], combination: str
) -> list[str]:
    """Resolve preregistered apo/binary/ternary state combinations."""
    categories = {
        "apo": sorted(key for key, value in available.items() if value == "apo"),
        "binary": sorted(key for key, value in available.items() if value == "binary"),
        "ternary": sorted(
            key for key, value in available.items() if value.startswith("ternary")
        ),
    }
    requested = {
        "apo only": ("apo",),
        "binary only": ("binary",),
        "ternary only": ("ternary",),
        "apo + binary": ("apo", "binary"),
        "binary + ternary": ("binary", "ternary"),
        "apo + binary + ternary": ("apo", "binary", "ternary"),
    }.get(combination)
    if requested is None:
        raise ValueError(f"unknown state combination: {combination}")
    missing = [category for category in requested if not categories[category]]
    if missing:
        raise ValueError(
            f"state combination {combination!r} is unavailable; missing={missing}"
        )
    return [state for category in requested for state in categories[category]]


def validate_residue_maps(
    maps: Mapping[str, Sequence[StateResidue]],
) -> tuple[str, ...]:
    if not maps:
        raise ValueError("residue maps cannot be empty")
    reference_state = sorted(maps)[0]
    reference = maps[reference_state]
    reference_indices = tuple(item.full_index_0 for item in reference)
    reference_sequence = tuple(item.wild_type for item in reference)
    if reference_indices != tuple(range(len(reference))):
        raise ValueError("reference residue map is not a complete zero-based sequence")
    for state, residues in maps.items():
        indices = tuple(item.full_index_0 for item in residues)
        sequence = tuple(item.wild_type for item in residues)
        if indices != reference_indices:
            raise ValueError(f"residue mapping indices disagree for state {state}")
        if sequence != reference_sequence:
            raise ValueError(f"wild-type tokens disagree for state {state}")
    return reference_sequence


def validate_fixed_tokens(
    sequences_by_state: Mapping[str, str], fixed_positions: Mapping[int, str]
) -> None:
    if not sequences_by_state:
        raise ValueError("fixed-token validation requires states")
    lengths = {len(sequence) for sequence in sequences_by_state.values()}
    if len(lengths) != 1:
        raise ValueError("state sequences have inconsistent lengths")
    length = next(iter(lengths))
    for position, expected in fixed_positions.items():
        if position < 0 or position >= length:
            raise ValueError(f"fixed position is outside state sequences: {position}")
        tokens = {sequence[position] for sequence in sequences_by_state.values()}
        if tokens != {expected.upper()}:
            raise ValueError(
                f"fixed token mismatch at {position}: observed={sorted(tokens)} "
                f"expected={expected.upper()}"
            )


def build_multistate_masks(states: Sequence[StateMask]) -> MultiStateMasks:
    if not states:
        raise ValueError("at least one state mask is required")
    if len({item.state for item in states}) != len(states):
        raise ValueError("state mask identifiers must be unique")
    sequences = {item.state: item.sequence for item in states}
    if len(set(sequences.values())) != 1:
        raise ValueError("state masks do not share an identical full parent sequence")
    length = len(states[0].sequence)
    for item in states:
        invalid = (
            set(item.hard_positions)
            .union(item.risk_positions)
            .difference(range(length))
        )
        if invalid:
            raise ValueError(f"state {item.state} has invalid mask positions {invalid}")
    hard = set(states[0].hard_positions)
    for item in states[1:]:
        hard.intersection_update(item.hard_positions)
    risk = set().union(*(item.risk_positions for item in states))
    variable = {
        position
        for position in range(length)
        if len(
            {
                (
                    position in item.hard_positions,
                    position in item.risk_positions,
                )
                for item in states
            }
        )
        > 1
    }
    return MultiStateMasks(
        state_intersection_hard_mask=frozenset(hard),
        state_union_risk_mask=frozenset(risk),
        state_variable_hinge_mask=frozenset(variable),
    )
