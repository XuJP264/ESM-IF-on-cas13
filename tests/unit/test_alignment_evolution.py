from pathlib import Path

import numpy as np
import pytest

from cas13_if.alignments.msa import (
    AlignmentError,
    map_ungapped_sequence_to_columns,
    read_aligned_fasta,
    sequence_weights,
)
from cas13_if.evolution.coevolution import (
    average_product_correction,
    bootstrap_top_pair_frequency,
    compute_mi_apc,
    permuted_cross_block_maxima,
)
from cas13_if.evolution.conservation import conservation_statistics


def test_alignment_validation_weights_and_mapping(tmp_path: Path) -> None:
    alignment = read_aligned_fasta(Path("tests/fixtures/protein_msa.fasta"))
    weights = sequence_weights(alignment)
    assert alignment.n_columns == 6
    assert weights.shape == (4,)
    assert np.all(weights > 0)
    assert map_ungapped_sequence_to_columns("AC-DE", "ACDE") == {
        0: 0,
        1: 1,
        2: 3,
        3: 4,
    }
    bad = tmp_path / "unaligned.fa"
    bad.write_text(">a\nACD\n>b\nACDE\n", encoding="utf-8")
    with pytest.raises(AlignmentError, match="identical lengths"):
        read_aligned_fasta(bad)


def test_conservation_has_gap_and_weighted_frequencies() -> None:
    alignment = read_aligned_fasta(Path("tests/fixtures/protein_msa.fasta"))
    statistics = conservation_statistics(alignment)
    assert statistics[0].conservation == pytest.approx(5 / 7)
    assert statistics[2].gap_fraction == 0.5
    assert statistics[0].consensus == "A"
    assert statistics[0].effective_sequence_count > 0


def test_mi_apc_bootstrap_and_permutation_are_seeded() -> None:
    alignment = read_aligned_fasta(
        Path("tests/fixtures/paired_msa.fasta"), alphabet="mixed"
    )
    result = compute_mi_apc(alignment)
    assert result.mutual_information.shape == (6, 6)
    assert np.allclose(np.diag(result.apc_corrected), 0)
    assert np.allclose(
        average_product_correction(result.mutual_information),
        result.apc_corrected,
    )
    first = bootstrap_top_pair_frequency(alignment, replicates=5, seed=7, top_n=2)
    second = bootstrap_top_pair_frequency(alignment, replicates=5, seed=7, top_n=2)
    assert np.array_equal(first, second)
    null = permuted_cross_block_maxima(alignment, split_column=3, replicates=4, seed=7)
    assert null.shape == (4,)
