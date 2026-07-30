"""Deterministic mock backend for tests only."""

from __future__ import annotations

import math
import random
from typing import Any

from cas13_if.backends.base import InverseFoldingBackend
from cas13_if.schemas import (
    STANDARD_AA,
    BackendCapabilities,
    Candidate,
    EvidenceLevel,
    SampleRequest,
    ScoreRequest,
    ScoreResult,
)


class MockBackend(InverseFoldingBackend):
    """Test double that is permanently marked as mock and Level 0."""

    def __init__(self) -> None:
        self._loaded = False

    def capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            scoring=True,
            sampling=True,
            hard_fixed=True,
            allowed_residue_filter=True,
            per_residue_probabilities=True,
        )

    def load(self) -> None:
        self._loaded = True

    def _require_loaded(self) -> None:
        if not self._loaded:
            raise RuntimeError("MockBackend.load() must be called first")

    def score(self, request: ScoreRequest) -> ScoreResult:
        self._require_loaded()
        per_residue = [-math.log(20.0)] * len(request.sequence)
        return ScoreResult(
            scaffold_id=request.scaffold_id,
            backend="mock",
            sequence=request.sequence,
            conditional_log_likelihood=sum(per_residue),
            perplexity=20.0,
            per_residue_log_probabilities=per_residue,
            is_mock=True,
            evidence_level=EvidenceLevel.IO_VALIDATED,
            metadata={"purpose": "tests_only"},
        )

    def sample(self, request: SampleRequest) -> list[Candidate]:
        self._require_loaded()
        alphabet = sorted(STANDARD_AA)
        candidates: list[Candidate] = []
        for sample_index in range(request.count):
            rng = random.Random(request.seed + sample_index)
            tokens: list[str] = []
            for index in range(len(request.parent_sequence)):
                if index in request.fixed_positions:
                    token = request.fixed_positions[index].upper()
                else:
                    allowed = sorted(request.allowed_residues.get(index, STANDARD_AA))
                    token = rng.choice(allowed)
                tokens.append(token)
            candidates.append(
                Candidate(
                    candidate_id=f"mock-{request.scaffold_id}-{sample_index:04d}",
                    scaffold_id=request.scaffold_id,
                    backend="mock",
                    sequence="".join(tokens),
                    parent_sequence=request.parent_sequence,
                    seed=request.seed + sample_index,
                    temperature=request.temperature,
                    is_mock=True,
                    evidence_level=EvidenceLevel.IO_VALIDATED,
                    fixed_positions=request.fixed_positions,
                    metadata={"purpose": "tests_only", "alphabet": alphabet},
                )
            )
        return candidates

    def metadata(self) -> dict[str, Any]:
        return {
            "backend": "mock",
            "is_mock": True,
            "purpose": "unit_and_integration_tests_only",
            "loaded": self._loaded,
        }
