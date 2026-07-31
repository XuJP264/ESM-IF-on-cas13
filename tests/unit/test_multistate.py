import pytest

from cas13_if.scoring.multistate import (
    StateMask,
    StateResidue,
    aggregate_multistate_scores,
    build_multistate_masks,
    normalize_state_weights,
    select_state_combination,
    state_rank_consistency,
    validate_fixed_tokens,
    validate_residue_maps,
)


def test_identical_state_scores_and_normalized_weights() -> None:
    result = aggregate_multistate_scores(
        {"apo": -1.0, "binary": -1.0}, weights={"apo": 2.0, "binary": 2.0}
    )
    assert result.normalized_weights == {"apo": 0.5, "binary": 0.5}
    assert result.multi_state_mean_score == -1.0
    assert result.multi_state_min_score == -1.0
    assert result.multi_state_variance == 0.0


def test_missing_state_and_invalid_weights_fail() -> None:
    with pytest.raises(ValueError, match="missing"):
        aggregate_multistate_scores({"apo": -1.0}, required_states={"apo", "binary"})
    with pytest.raises(ValueError, match="positive"):
        normalize_state_weights(["apo"], {"apo": 0.0})


def test_state_rank_consistency_detects_reversal() -> None:
    scores = {
        "one": {"apo": 1.0, "binary": 2.0},
        "two": {"apo": 2.0, "binary": 1.0},
        "three": {"apo": 3.0, "binary": 0.0},
    }
    assert state_rank_consistency(scores) == pytest.approx(-1.0)


def test_state_combination_requires_available_category() -> None:
    available = {"6e9e": "binary", "6e9f": "ternary"}
    assert select_state_combination(available, "binary + ternary") == [
        "6e9e",
        "6e9f",
    ]
    with pytest.raises(ValueError, match="missing"):
        select_state_combination(available, "apo + binary")


def test_residue_mapping_and_fixed_token_inconsistency_stop() -> None:
    valid = [StateResidue(0, "A", "exact"), StateResidue(1, "C", "exact")]
    assert validate_residue_maps({"apo": valid, "binary": valid}) == ("A", "C")
    invalid = [StateResidue(0, "A", "exact"), StateResidue(1, "D", "exact")]
    with pytest.raises(ValueError, match="tokens disagree"):
        validate_residue_maps({"apo": valid, "binary": invalid})
    with pytest.raises(ValueError, match="fixed token mismatch"):
        validate_fixed_tokens({"apo": "AC", "binary": "AD"}, {1: "C"})


def test_intersection_union_and_variable_hinge_masks() -> None:
    masks = build_multistate_masks(
        [
            StateMask("apo", "ACDE", frozenset({0, 1}), frozenset({1, 2})),
            StateMask("binary", "ACDE", frozenset({0, 2}), frozenset({2, 3})),
        ]
    )
    assert masks.state_intersection_hard_mask == {0}
    assert masks.state_union_risk_mask == {1, 2, 3}
    assert masks.state_variable_hinge_mask == {1, 2, 3}


def test_state_masks_reject_parent_mismatch() -> None:
    with pytest.raises(ValueError, match="identical full parent"):
        build_multistate_masks(
            [
                StateMask("apo", "AC", frozenset(), frozenset()),
                StateMask("binary", "AD", frozenset(), frozenset()),
            ]
        )
