"""Weighted MI/APC baselines, bootstrap stability, and permutation nulls."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from cas13_if.alignments.msa import Alignment, sequence_weights


@dataclass(frozen=True)
class CoevolutionResult:
    mutual_information: NDArray[np.float64]
    apc_corrected: NDArray[np.float64]
    weights: NDArray[np.float64]
    effective_sequence_count: float


def mutual_information_matrix(
    alignment: Alignment,
    *,
    weights: NDArray[np.float64] | None = None,
    pseudocount: float = 0.0,
) -> NDArray[np.float64]:
    if pseudocount < 0:
        raise ValueError("pseudocount must be non-negative")
    matrix = alignment.matrix()
    n_sequences, n_columns = matrix.shape
    if weights is None:
        weights = np.ones(n_sequences, dtype=np.float64)
    if weights.shape != (n_sequences,) or np.any(weights < 0):
        raise ValueError("weights must be a non-negative vector per MSA row")
    result = np.zeros((n_columns, n_columns), dtype=np.float64)
    for left in range(n_columns):
        for right in range(left + 1, n_columns):
            score = weighted_mutual_information(
                matrix[:, left],
                matrix[:, right],
                weights,
                pseudocount=pseudocount,
            )
            result[left, right] = result[right, left] = score
    return result


def weighted_mutual_information(
    left: NDArray[np.str_],
    right: NDArray[np.str_],
    weights: NDArray[np.float64],
    *,
    pseudocount: float = 0.0,
) -> float:
    valid = (left != "-") & (right != "-")
    if not valid.any():
        return 0.0
    left_valid = left[valid]
    right_valid = right[valid]
    weights_valid = weights[valid]
    left_symbols = sorted(set(left_valid.tolist()))
    right_symbols = sorted(set(right_valid.tolist()))
    total = float(weights_valid.sum()) + pseudocount * len(left_symbols) * len(
        right_symbols
    )
    if total <= 0:
        return 0.0
    joint = np.zeros((len(left_symbols), len(right_symbols)), dtype=np.float64)
    for left_index, left_symbol in enumerate(left_symbols):
        for right_index, right_symbol in enumerate(right_symbols):
            joint[left_index, right_index] = float(
                weights_valid[
                    (left_valid == left_symbol) & (right_valid == right_symbol)
                ].sum()
            )
    joint += pseudocount
    joint /= joint.sum()
    left_marginal = joint.sum(axis=1, keepdims=True)
    right_marginal = joint.sum(axis=0, keepdims=True)
    independent = left_marginal @ right_marginal
    positive = joint > 0
    contributions = joint[positive] * np.log(joint[positive] / independent[positive])
    return float(contributions.sum())


def average_product_correction(
    matrix: NDArray[np.float64],
) -> NDArray[np.float64]:
    if matrix.ndim != 2 or matrix.shape[0] != matrix.shape[1]:
        raise ValueError("APC requires a square matrix")
    corrected: NDArray[np.float64] = matrix.copy().astype(np.float64)
    size = matrix.shape[0]
    if size < 2:
        return corrected
    off_diagonal = ~np.eye(size, dtype=bool)
    global_mean = float(matrix[off_diagonal].mean())
    if global_mean == 0:
        np.fill_diagonal(corrected, 0.0)
        return corrected
    row_means = np.array(
        [float(np.delete(matrix[index], index).mean()) for index in range(size)]
    )
    corrected -= np.outer(row_means, row_means) / global_mean
    np.fill_diagonal(corrected, 0.0)
    return corrected


def compute_mi_apc(
    alignment: Alignment,
    *,
    identity_threshold: float = 0.8,
    pseudocount: float = 0.0,
) -> CoevolutionResult:
    weights = sequence_weights(alignment, identity_threshold=identity_threshold)
    mi = mutual_information_matrix(alignment, weights=weights, pseudocount=pseudocount)
    return CoevolutionResult(
        mutual_information=mi,
        apc_corrected=average_product_correction(mi),
        weights=weights,
        effective_sequence_count=float(weights.sum()),
    )


def bootstrap_top_pair_frequency(
    alignment: Alignment,
    *,
    replicates: int,
    seed: int,
    top_n: int,
) -> NDArray[np.float64]:
    if replicates < 1 or top_n < 1:
        raise ValueError("replicates and top_n must be positive")
    rng = np.random.default_rng(seed)
    counts = np.zeros((alignment.n_columns, alignment.n_columns), dtype=np.int64)
    for _ in range(replicates):
        indices = rng.integers(0, alignment.n_sequences, alignment.n_sequences)
        sampled = Alignment(
            identifiers=tuple(alignment.identifiers[index] for index in indices),
            sequences=tuple(alignment.sequences[index] for index in indices),
        )
        scores = compute_mi_apc(sampled).apc_corrected
        upper = np.triu_indices(alignment.n_columns, k=1)
        order = np.argsort(scores[upper])[::-1][:top_n]
        counts[upper[0][order], upper[1][order]] += 1
    counts += counts.T
    return counts / replicates


def permuted_cross_block_maxima(
    alignment: Alignment,
    *,
    split_column: int,
    replicates: int,
    seed: int,
) -> NDArray[np.float64]:
    if not 0 < split_column < alignment.n_columns:
        raise ValueError("split_column must separate nonempty protein/RNA blocks")
    rng = np.random.default_rng(seed)
    maxima = np.zeros(replicates, dtype=np.float64)
    protein = [sequence[:split_column] for sequence in alignment.sequences]
    rna = [sequence[split_column:] for sequence in alignment.sequences]
    for replicate in range(replicates):
        permutation = rng.permutation(alignment.n_sequences)
        permuted = Alignment(
            identifiers=alignment.identifiers,
            sequences=tuple(
                protein[index] + rna[int(permutation[index])]
                for index in range(alignment.n_sequences)
            ),
        )
        scores = compute_mi_apc(permuted).apc_corrected
        maxima[replicate] = float(scores[:split_column, split_column:].max())
    return maxima
