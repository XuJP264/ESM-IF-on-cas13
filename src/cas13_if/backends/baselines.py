"""Evolutionary-profile and matched-random sequence baselines."""

from __future__ import annotations

import math
from typing import Any

import numpy as np

from cas13_if.backends.base import InverseFoldingBackend
from cas13_if.schemas import (
    STANDARD_AA,
    BackendCapabilities,
    Candidate,
    EvidenceLevel,
    PositionTrace,
    SampleRequest,
    ScoreRequest,
    ScoreResult,
)

BASELINE_ALPHABET = tuple(sorted(STANDARD_AA))


class MsaProfileBackend(InverseFoldingBackend):
    """Sample independently from a validated aligned-position profile."""

    def __init__(
        self,
        frequencies: list[dict[str, float]],
        *,
        pseudocount: float = 1e-6,
    ) -> None:
        if pseudocount <= 0:
            raise ValueError("pseudocount must be positive")
        self.frequencies = frequencies
        self.pseudocount = pseudocount
        self._probabilities: np.ndarray | None = None

    def capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            scoring=True,
            sampling=True,
            hard_fixed=True,
            allowed_residue_filter=True,
            per_residue_probabilities=True,
        )

    def load(self) -> None:
        matrix = np.zeros(
            (len(self.frequencies), len(BASELINE_ALPHABET)), dtype=np.float64
        )
        for position, frequencies in enumerate(self.frequencies):
            unknown = set(frequencies).difference(STANDARD_AA)
            if unknown:
                raise ValueError(
                    f"profile position {position} has unknown residues: "
                    f"{sorted(unknown)}"
                )
            for token_index, token in enumerate(BASELINE_ALPHABET):
                value = float(frequencies.get(token, 0.0))
                if value < 0 or not math.isfinite(value):
                    raise ValueError(
                        f"invalid profile frequency at {position}/{token}: {value}"
                    )
                matrix[position, token_index] = value + self.pseudocount
            matrix[position] /= matrix[position].sum()
        self._probabilities = matrix

    def _require_loaded(self) -> np.ndarray:
        if self._probabilities is None:
            raise RuntimeError("MsaProfileBackend.load() must be called first")
        return self._probabilities

    def score(self, request: ScoreRequest) -> ScoreResult:
        probabilities = self._require_loaded()
        if len(request.sequence) != probabilities.shape[0]:
            raise ValueError("sequence length does not match MSA profile length")
        token_indices = {token: index for index, token in enumerate(BASELINE_ALPHABET)}
        per_residue = [
            math.log(float(probabilities[position, token_indices[token]]))
            for position, token in enumerate(request.sequence)
        ]
        mean = sum(per_residue) / len(per_residue)
        return ScoreResult(
            scaffold_id=request.scaffold_id,
            backend="msa_profile",
            sequence=request.sequence,
            conditional_log_likelihood=sum(per_residue),
            perplexity=math.exp(-mean),
            per_residue_log_probabilities=per_residue,
            is_mock=False,
            evidence_level=EvidenceLevel.IO_VALIDATED,
            metadata={
                "model": "independent_aligned_position_profile",
                "pseudocount": self.pseudocount,
                "not_an_inverse_folding_score": True,
            },
        )

    def sample(self, request: SampleRequest) -> list[Candidate]:
        probabilities = self._require_loaded()
        if len(request.parent_sequence) != probabilities.shape[0]:
            raise ValueError("parent length does not match MSA profile length")
        candidates: list[Candidate] = []
        for sample_index in range(request.count):
            seed = request.seed + sample_index
            rng = np.random.default_rng(seed)
            sequence: list[str] = []
            traces: list[PositionTrace] = []
            for position in range(len(request.parent_sequence)):
                position_probabilities = probabilities[position].copy()
                allowed = request.allowed_residues.get(position)
                if allowed is not None:
                    for token_index, token in enumerate(BASELINE_ALPHABET):
                        if token not in allowed:
                            position_probabilities[token_index] = 0.0
                    if position_probabilities.sum() <= 0:
                        raise ValueError(
                            f"profile has no probability mass allowed at {position}"
                        )
                    position_probabilities /= position_probabilities.sum()
                fixed = position in request.fixed_positions
                token = (
                    request.fixed_positions[position].upper()
                    if fixed
                    else str(rng.choice(BASELINE_ALPHABET, p=position_probabilities))
                )
                sequence.append(token)
                traces.append(
                    PositionTrace(
                        index=position,
                        logits=[
                            math.log(max(float(value), 1e-12))
                            for value in probabilities[position]
                        ],
                        probabilities=[
                            float(value) for value in position_probabilities
                        ],
                        selected_token=token,
                        fixed=fixed,
                        temperature=request.temperature,
                        seed=seed,
                    )
                )
            candidates.append(
                Candidate(
                    candidate_id=(
                        f"msa-profile-{request.scaffold_id}-{sample_index:04d}"
                    ),
                    scaffold_id=request.scaffold_id,
                    backend="msa_profile",
                    sequence="".join(sequence),
                    parent_sequence=request.parent_sequence,
                    seed=seed,
                    temperature=request.temperature,
                    is_mock=False,
                    evidence_level=EvidenceLevel.IO_VALIDATED,
                    fixed_positions=request.fixed_positions,
                    traces=traces,
                    metadata={
                        "profile_length": int(probabilities.shape[0]),
                        "independent_positions": True,
                        "evidence_level_note": (
                            "Level 1 requires downstream novelty evaluation"
                        ),
                    },
                )
            )
        return candidates

    def metadata(self) -> dict[str, Any]:
        return {
            "backend": "msa_profile",
            "loaded": self._probabilities is not None,
            "positions": len(self.frequencies),
            "pseudocount": self.pseudocount,
            "is_mock": False,
        }


class MatchedRandomMutationBackend(InverseFoldingBackend):
    """Mutate every free design position while preserving the registered mask."""

    def __init__(self) -> None:
        self._loaded = False

    def capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            scoring=False,
            sampling=True,
            hard_fixed=True,
            allowed_residue_filter=True,
            per_residue_probabilities=True,
        )

    def load(self) -> None:
        self._loaded = True

    def score(self, request: ScoreRequest) -> ScoreResult:
        del request
        raise NotImplementedError(
            "matched random mutation has no intrinsic sequence score"
        )

    def sample(self, request: SampleRequest) -> list[Candidate]:
        if not self._loaded:
            raise RuntimeError(
                "MatchedRandomMutationBackend.load() must be called first"
            )
        candidates: list[Candidate] = []
        for sample_index in range(request.count):
            seed = request.seed + sample_index
            rng = np.random.default_rng(seed)
            sequence: list[str] = []
            traces: list[PositionTrace] = []
            for position, parent_token in enumerate(request.parent_sequence):
                fixed = position in request.fixed_positions
                allowed = set(request.allowed_residues.get(position, STANDARD_AA))
                if fixed:
                    token = request.fixed_positions[position].upper()
                else:
                    mutation_choices = sorted(allowed.difference({parent_token}))
                    token = (
                        str(rng.choice(mutation_choices))
                        if mutation_choices
                        else parent_token
                    )
                probabilities = [
                    (
                        1.0 / len(allowed.difference({parent_token}))
                        if aa in allowed.difference({parent_token})
                        and allowed.difference({parent_token})
                        else 0.0
                    )
                    for aa in BASELINE_ALPHABET
                ]
                sequence.append(token)
                traces.append(
                    PositionTrace(
                        index=position,
                        logits=[
                            0.0 if probability > 0 else -1e9
                            for probability in probabilities
                        ],
                        probabilities=probabilities,
                        selected_token=token,
                        fixed=fixed,
                        temperature=request.temperature,
                        seed=seed,
                    )
                )
            candidates.append(
                Candidate(
                    candidate_id=(
                        f"matched-random-{request.scaffold_id}-{sample_index:04d}"
                    ),
                    scaffold_id=request.scaffold_id,
                    backend="matched_random_mutation",
                    sequence="".join(sequence),
                    parent_sequence=request.parent_sequence,
                    seed=seed,
                    temperature=request.temperature,
                    is_mock=False,
                    evidence_level=EvidenceLevel.IO_VALIDATED,
                    fixed_positions=request.fixed_positions,
                    traces=traces,
                    metadata={
                        "mutation_rule": "uniform_non_parent_at_every_free_position",
                        "evidence_level_note": (
                            "Level 1 requires downstream novelty evaluation"
                        ),
                    },
                )
            )
        return candidates

    def metadata(self) -> dict[str, Any]:
        return {
            "backend": "matched_random_mutation",
            "loaded": self._loaded,
            "is_mock": False,
            "intrinsic_score": False,
        }
