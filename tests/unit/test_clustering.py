from pathlib import Path

import pytest

from cas13_if.data.clustering import (
    MmseqsParameters,
    assert_no_cluster_leakage,
    assign_cluster_splits,
    cluster_summary,
    parse_cluster_tsv,
)


def test_cluster_mapping_summary_and_deterministic_split() -> None:
    mapping = parse_cluster_tsv(Path("tests/fixtures/clusters.tsv"))
    assert mapping["member1"] == "rep1"
    assert cluster_summary(mapping)["cluster_count"] == 2
    first = assign_cluster_splits(
        mapping, seed=7, train_fraction=0.5, validation_fraction=0.2
    )
    second = assign_cluster_splits(
        mapping, seed=7, train_fraction=0.5, validation_fraction=0.2
    )
    assert first == second
    assert first["rep1"] == first["member1"]


def test_leakage_gate_fails_closed() -> None:
    mapping = {"a": "a", "b": "a"}
    with pytest.raises(RuntimeError, match="DATA LEAKAGE"):
        assert_no_cluster_leakage(mapping, {"a": "train", "b": "test"})
    with pytest.raises(ValueError, match="missing"):
        assert_no_cluster_leakage(mapping, {"a": "train"})
    with pytest.raises(ValueError):
        MmseqsParameters(minimum_identity=0).validate()
