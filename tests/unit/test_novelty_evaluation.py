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
