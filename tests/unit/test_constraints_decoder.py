import numpy as np
import pytest

from cas13_if.constraints.masks import PositionEvidence, merge_position_evidence
from cas13_if.generation.constrained_decoder import constrained_autoregressive_sample

ALPHABET = ("A", "C", "D")


def logits(prefix: str, index: int) -> np.ndarray:
    del prefix
    return np.array([0.0, float(index), -1.0])


def test_all_fixed_exactly_recovers_and_records_trace() -> None:
    expected = "ACDA"
    decoded = constrained_autoregressive_sample(
        length=len(expected),
        logits_function=logits,
        alphabet=ALPHABET,
        fixed_positions=dict(enumerate(expected)),
        seed=4,
    )
    assert decoded.sequence == expected
    assert decoded.fixed_position_violations == 0
    assert all(position.fixed for position in decoded.trace)
    assert all(len(position.logits) == len(ALPHABET) for position in decoded.trace)


def test_free_and_partial_sampling_are_seeded_and_causal() -> None:
    first = constrained_autoregressive_sample(
        length=8, logits_function=logits, alphabet=ALPHABET, seed=9
    )
    second = constrained_autoregressive_sample(
        length=8, logits_function=logits, alphabet=ALPHABET, seed=9
    )
    partial = constrained_autoregressive_sample(
        length=8,
        logits_function=logits,
        alphabet=ALPHABET,
        fixed_positions={2: "D", 7: "A"},
        allowed_residues={3: {"A", "C"}},
        seed=9,
    )
    assert first.sequence == second.sequence
    assert partial.sequence[2] == "D"
    assert partial.sequence[7] == "A"
    assert partial.fixed_position_violations == 0
    assert partial.semantics == "left_to_right_causal_hard_fixed"


def test_decoder_fails_immediately_on_invalid_constraints() -> None:
    with pytest.raises(ValueError, match="outside"):
        constrained_autoregressive_sample(
            length=2,
            logits_function=logits,
            alphabet=ALPHABET,
            fixed_positions={2: "A"},
        )
    with pytest.raises(ValueError, match="disallowed"):
        constrained_autoregressive_sample(
            length=2,
            logits_function=logits,
            alphabet=ALPHABET,
            fixed_positions={0: "A"},
            allowed_residues={0: {"C"}},
        )


def test_mask_merge_priority_and_manual_hepn_gate() -> None:
    merged = merge_position_evidence(
        length=3,
        evidence=[
            PositionEvidence(0, "conservation", "soft_constrained", "high C_i"),
            PositionEvidence(0, "manual", "hard_fixed", "catalytic"),
        ],
    )
    assert merged[0].final_class == "hard_fixed"
    assert merged[1].final_class == "free"
    with pytest.raises(ValueError, match="manual confirmation"):
        merge_position_evidence(
            length=2,
            evidence=[PositionEvidence(0, "hepn_regex", "hard_fixed", "RxxxxH motif")],
        )
