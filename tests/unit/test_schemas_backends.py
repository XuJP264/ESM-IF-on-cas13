import pytest
from pydantic import ValidationError

from cas13_if.backends.mock import MockBackend
from cas13_if.schemas import (
    Candidate,
    EvidenceLevel,
    SampleRequest,
    ScoreRequest,
)


def test_mock_backend_requires_load_and_is_deterministic() -> None:
    backend = MockBackend()
    request = SampleRequest(
        scaffold_id="x",
        structure_path="unused.pdb",
        parent_sequence="ACDE",
        count=2,
        seed=11,
        fixed_positions={0: "A"},
        allowed_residues={1: {"C"}},
    )
    with pytest.raises(RuntimeError, match="load"):
        backend.sample(request)
    backend.load()
    first = backend.sample(request)
    second = backend.sample(request)
    assert [candidate.sequence for candidate in first] == [
        candidate.sequence for candidate in second
    ]
    assert all(candidate.sequence.startswith("AC") for candidate in first)
    assert all(candidate.is_mock for candidate in first)
    score = backend.score(
        ScoreRequest(scaffold_id="x", structure_path="unused.pdb", sequence="ACDE")
    )
    assert score.perplexity == 20
    assert backend.metadata()["purpose"] == "unit_and_integration_tests_only"


def test_candidate_rejects_fixed_violation_and_mock_overclaim() -> None:
    with pytest.raises(ValidationError, match="fixed-position"):
        Candidate(
            candidate_id="bad",
            scaffold_id="x",
            backend="mock",
            sequence="CC",
            seed=1,
            temperature=1,
            is_mock=True,
            evidence_level=EvidenceLevel.IO_VALIDATED,
            fixed_positions={0: "A"},
        )
    with pytest.raises(ValidationError, match="Level 0"):
        Candidate(
            candidate_id="bad-level",
            scaffold_id="x",
            backend="mock",
            sequence="AC",
            seed=1,
            temperature=1,
            is_mock=True,
            evidence_level=EvidenceLevel.INVERSE_FOLDING_COMPATIBILITY,
        )


def test_sample_request_rejects_bad_positions() -> None:
    with pytest.raises(ValidationError, match="outside"):
        SampleRequest(
            scaffold_id="x",
            structure_path="unused",
            parent_sequence="AC",
            fixed_positions={2: "D"},
        )
