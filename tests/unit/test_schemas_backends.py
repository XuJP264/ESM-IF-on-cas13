from contextlib import nullcontext
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from pydantic import ValidationError

from cas13_if.backends.baselines import (
    MatchedRandomMutationBackend,
    MsaProfileBackend,
)
from cas13_if.backends.esm_if1 import (
    EsmIf1Backend,
    EsmIf1ConstrainedBackend,
    _sampling_request_digest,
    _sha256,
)
from cas13_if.backends.mock import MockBackend
from cas13_if.backends.mpnn import _proteinmpnn_slot_map, _restore_pdb_sequence
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
    with pytest.raises(ValidationError, match="disallowed"):
        SampleRequest(
            scaffold_id="x",
            structure_path="unused",
            parent_sequence="AC",
            fixed_positions={0: "A"},
            allowed_residues={0: {"C"}},
        )


def test_esm_backends_are_offline_and_declare_constraint_semantics(
    tmp_path: Path,
) -> None:
    standard = EsmIf1Backend(tmp_path / "missing.pt")
    constrained = EsmIf1ConstrainedBackend(tmp_path / "missing.pt")
    assert not standard.capabilities().hard_fixed
    assert constrained.capabilities().hard_fixed
    assert constrained.capabilities().protein_multichain
    assert not constrained.capabilities().rna_atomic_context
    with pytest.raises(FileNotFoundError, match="fetch_models"):
        standard.load()


def test_esm_sampling_request_digest_separates_conditions() -> None:
    base = SampleRequest(
        scaffold_id="6E9F-A",
        structure_path="/node-a/6e9f.cif",
        parent_sequence="ACDE",
        temperature=0.1,
        seed=7,
        fixed_positions={0: "A"},
        protein_chains=["A"],
    )
    equivalent_on_another_node = base.model_copy(
        update={"structure_path": "/node-b/6e9f.cif"}
    )
    another_temperature = base.model_copy(update={"temperature": 0.5})
    another_constraint = base.model_copy(update={"fixed_positions": {1: "C"}})
    digest = _sampling_request_digest("esm_if1_constrained", base)
    assert digest == _sampling_request_digest(
        "esm_if1_constrained", equivalent_on_another_node
    )
    assert digest != _sampling_request_digest(
        "esm_if1_constrained", another_temperature
    )
    assert digest != _sampling_request_digest("esm_if1_constrained", another_constraint)


def test_esm_backend_score_and_sample_contracts(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
) -> None:
    checkpoint = tmp_path / "checkpoint.pt"
    checkpoint.write_bytes(b"offline-checkpoint")
    backend = EsmIf1Backend(checkpoint, device="cpu")
    backend._model = object()
    backend._torch = SimpleNamespace(inference_mode=nullcontext)
    backend._device = "cpu"
    backend._checkpoint_sha256 = _sha256(checkpoint)
    coords = np.zeros((2, 3, 3), dtype=np.float32)
    monkeypatch.setattr(
        backend,
        "_load_conditioning_coords",
        lambda _path, _chains: (coords, "A", ["A"]),
    )
    monkeypatch.setattr(
        backend,
        "_get_sequence_loss",
        lambda _coords, _sequence: (
            np.asarray([1.0, 2.0]),
            np.asarray([False, False]),
        ),
    )

    result = backend.score(
        ScoreRequest(
            scaffold_id="fixture",
            structure_path="fixture.pdb",
            sequence="AC",
            protein_chains=["A"],
        )
    )
    assert result.conditional_log_likelihood == -3.0
    assert result.per_residue_log_probabilities == [-1.0, -2.0]
    assert result.metadata["rna_atomic_context"] is False

    monkeypatch.setattr(
        backend,
        "_decode",
        lambda **kwargs: ("AC", []),
    )
    request = SampleRequest(
        scaffold_id="fixture",
        structure_path="fixture.pdb",
        parent_sequence="AC",
        count=2,
        seed=5,
        protein_chains=["A"],
    )
    candidates = backend.sample(request)
    assert len(candidates) == 2
    assert len({candidate.candidate_id for candidate in candidates}) == 2
    assert [candidate.seed for candidate in candidates] == [5, 6]
    assert all(
        candidate.metadata["semantics"] == "left_to_right_causal"
        for candidate in candidates
    )

    constrained = request.model_copy(update={"fixed_positions": {0: "A"}})
    with pytest.raises(ValueError, match="does not accept residue constraints"):
        backend.sample(constrained)
    with pytest.raises(ValueError, match="sequence length"):
        backend.score(
            ScoreRequest(
                scaffold_id="fixture",
                structure_path="fixture.pdb",
                sequence="ACD",
                protein_chains=["A"],
            )
        )
    assert backend.metadata()["offline_local_checkpoint"] is True


def test_esm_coordinate_contracts(tmp_path: Path) -> None:
    backend = EsmIf1Backend(tmp_path / "checkpoint.pt", device="cpu")
    fixture = Path("tests/fixtures/minimal_complex.pdb")
    coordinates, target, conditioning = backend._load_conditioning_coords(
        str(fixture), ["A"]
    )
    assert target == "A"
    assert conditioning == ["A"]
    assert coordinates.shape == (2, 3, 3)
    assert np.isnan(coordinates[1, 2]).all()

    multichain = {"A": coordinates, "B": coordinates[:1]}
    concatenated = backend._concatenate(multichain, "A")
    assert concatenated.shape == (13, 3, 3)
    assert backend._target_length(multichain, "A") == 2
    assert backend._target_length(coordinates, "A") == 2

    with pytest.raises(ValueError, match="must explicitly name"):
        backend._load_conditioning_coords(str(fixture), [])
    with pytest.raises(FileNotFoundError, match="structure is missing"):
        backend._load_conditioning_coords(str(tmp_path / "missing.pdb"), ["A"])
    with pytest.raises(ValueError, match="has no residues"):
        backend._load_conditioning_coords(str(fixture), ["R"])


def test_profile_and_matched_random_baselines_are_seeded() -> None:
    profile = MsaProfileBackend(
        [
            {"A": 0.8, "C": 0.2},
            {"C": 1.0},
            {"D": 0.5, "E": 0.5},
        ]
    )
    profile.load()
    request = SampleRequest(
        scaffold_id="profile",
        structure_path="not_used",
        parent_sequence="ACD",
        count=2,
        seed=12,
        fixed_positions={1: "C"},
    )
    assert [item.sequence for item in profile.sample(request)] == [
        item.sequence for item in profile.sample(request)
    ]
    score = profile.score(
        ScoreRequest(
            scaffold_id="profile",
            structure_path="not_used",
            sequence="ACD",
        )
    )
    assert score.metadata["not_an_inverse_folding_score"]

    random_backend = MatchedRandomMutationBackend()
    random_backend.load()
    candidate = random_backend.sample(request)[0]
    assert candidate.sequence[1] == "C"
    assert candidate.sequence[0] != "A"
    assert candidate.sequence[2] != "D"
    with pytest.raises(NotImplementedError, match="no intrinsic"):
        random_backend.score(
            ScoreRequest(
                scaffold_id="random",
                structure_path="not_used",
                sequence="ACD",
            )
        )

    retained_backend = MatchedRandomMutationBackend(mutation_probability=0.0)
    retained_backend.load()
    retained = retained_backend.sample(request)[0]
    assert retained.sequence == "ACD"


def test_mpnn_adapter_restores_fixed_identity_and_insertion_mapping(
    tmp_path: Path,
) -> None:
    restored_path = tmp_path / "restored.pdb"
    sequence, residue_keys = _restore_pdb_sequence(
        source=Path("tests/fixtures/minimal_complex.pdb"),
        destination=restored_path,
        chain="A",
        fixed_positions={0: "R"},
    )
    assert sequence == "RG"
    assert residue_keys[1].insertion_code == "A"
    slots = _proteinmpnn_slot_map(residue_keys)
    assert slots[residue_keys[0]] == 0
    assert slots[residue_keys[1]] == 1
