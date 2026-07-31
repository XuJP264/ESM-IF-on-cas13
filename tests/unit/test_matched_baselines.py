import pytest

from cas13_if.evaluation.matching import (
    add_identity_metrics,
    exact_paired_permutation_p_value,
    identity_matched_source_consensus,
    paired_seed_statistics,
    position_set_hash,
    proposal_seed,
    select_balanced_candidates,
)


def test_matching_is_deterministic_and_balanced_by_seed() -> None:
    proposals = []
    for method in ("a", "b"):
        for seed in (1, 2):
            for index, identity in enumerate((0.24, 0.27)):
                proposals.append(
                    {
                        "candidate_id": f"{method}-{seed}-{index}",
                        "method": method,
                        "seed_block": seed,
                        "parent_identity": identity,
                        "designed_position_identity": identity,
                        "is_mock": False,
                    }
                )
    selected, failures = select_balanced_candidates(
        proposals,
        methods=["a", "b"],
        seed_blocks=[1, 2],
        minimum_parent_identity=0.2,
        maximum_parent_identity=0.3,
        minimum_designed_identity=0.2,
        maximum_designed_identity=0.3,
        target_identity=0.25,
    )
    assert [row["candidate_id"] for row in selected] == [
        "a-1-0",
        "a-2-0",
        "b-1-0",
        "b-2-0",
    ]
    assert len(failures) == 4


def test_matching_fails_when_one_method_seed_has_no_common_identity() -> None:
    with pytest.raises(ValueError, match="no identity-matched"):
        select_balanced_candidates(
            [
                {
                    "candidate_id": "outside",
                    "method": "a",
                    "seed_block": 1,
                    "parent_identity": 0.9,
                    "designed_position_identity": 0.9,
                }
            ],
            methods=["a"],
            seed_blocks=[1],
            minimum_parent_identity=0.2,
            maximum_parent_identity=0.3,
            minimum_designed_identity=0.2,
            maximum_designed_identity=0.3,
            target_identity=0.25,
        )


def test_identity_hash_and_seed_level_statistics() -> None:
    assert position_set_hash({3, 1}) == position_set_hash({1, 3})
    rows = add_identity_metrics(
        [
            {
                "candidate_id": "x",
                "sequence": "ACDE",
                "method": "x",
                "seed_block": 1,
            }
        ],
        parent_sequence="ACDF",
        designed_positions={2, 3},
    )
    assert rows[0]["parent_identity"] == 0.75
    assert rows[0]["designed_position_identity"] == 0.5
    evaluated = [
        {"method": "ref", "seed_block": 1, "metric": 1.0},
        {"method": "ref", "seed_block": 2, "metric": 2.0},
        {"method": "test", "seed_block": 1, "metric": 2.0},
        {"method": "test", "seed_block": 2, "metric": 4.0},
    ]
    statistics = paired_seed_statistics(
        evaluated,
        metrics=["metric"],
        reference_method="ref",
        bootstrap_replicates=200,
        confidence=0.95,
        seed=7,
    )
    assert statistics[0]["independent_unit_count"] == 2
    assert statistics[0]["paired_mean_difference"] == 1.5
    assert exact_paired_permutation_p_value([1.0, 2.0]) == 0.5


def test_proposal_seeds_do_not_overlap_adjacent_seed_blocks() -> None:
    seeds = {
        proposal_seed(seed_block, proposal_index)
        for seed_block in (20260731, 20260732)
        for proposal_index in range(2)
    }
    assert seeds == {20260731, 20260732, 21260731, 21260732}
    with pytest.raises(ValueError, match="non-negative"):
        proposal_seed(1, -1)


def test_source_consensus_matches_identity_without_inventing_tokens() -> None:
    parent = "AAAAAAAAAA"
    first = "AACCCCCCCC"
    second = "CCADDDDDDD"
    sequence, metadata = identity_matched_source_consensus(
        parent_sequence=parent,
        first_sequence=first,
        second_sequence=second,
        first_confidences=[0.9] * len(parent),
        second_confidences=[0.8] * len(parent),
        target_identity=0.2,
    )
    assert sum(token == "A" for token in sequence) == 2
    assert metadata["achieved_parent_identity"] == 0.2
    assert all(
        token in {left, right}
        for token, left, right in zip(sequence, first, second, strict=True)
    )
    with pytest.raises(ValueError, match="identical non-zero lengths"):
        identity_matched_source_consensus(
            parent_sequence="AA",
            first_sequence="A",
            second_sequence="AA",
            first_confidences=[1.0],
            second_confidences=[1.0, 1.0],
            target_identity=0.5,
        )
