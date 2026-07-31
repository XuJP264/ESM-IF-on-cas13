import math

import pytest

from cas13_if.evaluation.metrics import native_recovery, perplexity
from cas13_if.novelty.metrics import (
    composition_deviation,
    designed_position_identity,
    hydrophobicity_proxy,
    longest_homopolymer,
    low_complexity_windows,
    net_charge_proxy,
    pairwise_candidate_identity,
    sequence_identity,
    shannon_entropy,
)
from cas13_if.novelty.pipeline import NoveltyThresholds, evaluate_candidate_novelty


def test_novelty_and_sequence_qc_metrics() -> None:
    assert sequence_identity("ACDE", "ACDF") == 0.75
    assert designed_position_identity("ACDE", "ACDF", {3}) == 0
    assert longest_homopolymer("AAAACD") == 4
    assert low_complexity_windows(
        "AAAAACDE", window=4, maximum_single_residue_fraction=0.5
    )
    assert shannon_entropy("ACDE") == math.log(4)
    assert composition_deviation("AAAA", "CCCC") == 1
    assert net_charge_proxy("KRDE") == 0
    assert hydrophobicity_proxy("AVDE") == 0.5
    assert pairwise_candidate_identity(["AC", "AD", "CD"]) == [0.5, 0.0, 0.5]


def test_recovery_by_region_and_perplexity() -> None:
    result = native_recovery(
        "ACDF",
        "ACDE",
        designed_positions={2, 3},
        fixed_positions={0: "A"},
        regions={"interface": {1, 3}, "empty": set()},
    )
    assert result.overall == 0.75
    assert result.designed_positions == 0.5
    assert result.regions["interface"] == 0.5
    assert result.regions["empty"] is None
    assert result.fixed_position_violations == 0
    assert perplexity([-math.log(2), -math.log(2)]) == pytest.approx(2)


def test_candidate_novelty_requires_atlas_hit_and_sequence_qc() -> None:
    candidates = [
        {
            "candidate_id": "candidate-a",
            "pdb_id": "fixture",
            "method": "fixture",
            "scaffold_id": "fixture-A",
            "sequence": "ACDEFGHIK",
            "parent_sequence": "ACDEFGHIL",
            "temperature": 1.0,
            "seed": 7,
            "fixed_positions": [0],
            "fixed_position_violations": 0,
            "source_is_mock": False,
            "source_evidence_level": 2,
        },
        {
            "candidate_id": "candidate-b",
            "pdb_id": "fixture",
            "method": "fixture",
            "scaffold_id": "fixture-A",
            "sequence": "AAAAAAAAA",
            "parent_sequence": "CCCCCCCCC",
            "temperature": 1.0,
            "seed": 8,
            "fixed_positions": [0],
            "fixed_position_violations": 0,
            "source_is_mock": False,
            "source_evidence_level": 2,
        },
    ]
    hits = {
        "candidate-a": {
            "target_sequence_sha256": "atlas-a",
            "identity": 0.5,
            "alignment_length": 9,
            "query_coverage": 1.0,
            "target_coverage": 1.0,
            "evalue": 1e-5,
            "bits": 50.0,
        }
    }
    rows, summary = evaluate_candidate_novelty(
        candidates,
        hits,
        NoveltyThresholds(
            maximum_parent_identity=0.95,
            maximum_atlas_identity=0.8,
            maximum_homopolymer_length=4,
            maximum_low_complexity_windows=0,
            minimum_designed_position_entropy=0.0,
            low_complexity_window=4,
            low_complexity_maximum_fraction=0.5,
        ),
    )
    assert rows[0]["passes_level1_novelty"] is True
    assert rows[0]["evidence_level"] == 1
    assert rows[1]["passes_level1_novelty"] is False
    assert (
        "no_atlas_hit_at_required_query_coverage" in rows[1]["novelty_filter_failures"]
    )
    assert "homopolymer_above_threshold" in rows[1]["novelty_filter_failures"]
    assert summary["passes_level1_novelty"] == 1
