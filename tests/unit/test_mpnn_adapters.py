import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest

from cas13_if.backends.mpnn import LigandMpnnBackend, ProteinMpnnBackend
from cas13_if.schemas import SampleRequest, ScoreRequest


def _assets(tmp_path: Path) -> dict[str, Path]:
    protein_upstream = tmp_path / "ProteinMPNN"
    ligand_upstream = tmp_path / "LigandMPNN"
    protein_upstream.mkdir()
    ligand_upstream.mkdir()
    (protein_upstream / "protein_mpnn_run.py").write_text("# fixture\n")
    (ligand_upstream / "run.py").write_text("# fixture\n")
    files = {
        "protein_upstream": protein_upstream,
        "ligand_upstream": ligand_upstream,
        "protein": tmp_path / "v_48_020.pt",
        "ligand": tmp_path / "ligand.pt",
        "soluble": tmp_path / "soluble.pt",
        "python": tmp_path / "python",
    }
    for key in ("protein", "ligand", "soluble", "python"):
        files[key].write_bytes(b"fixture")
    return files


def _request() -> SampleRequest:
    return SampleRequest(
        scaffold_id="fixture",
        structure_path="tests/fixtures/minimal_complex.pdb",
        parent_sequence="RG",
        count=1,
        temperature=1.0,
        seed=7,
        fixed_positions={0: "R"},
        protein_chains=["A"],
    )


def test_proteinmpnn_adapter_runs_fixed_fixture(monkeypatch, tmp_path: Path) -> None:
    assets = _assets(tmp_path)

    def fake_run(command, **_kwargs):
        output = Path(command[command.index("--out_folder") + 1])
        pdb = Path(command[command.index("--pdb_path") + 1])
        (output / "seqs").mkdir(parents=True)
        (output / "probs").mkdir(parents=True)
        (output / f"seqs/{pdb.stem}.fa").write_text(
            ">native\nRG\n>sample\nRG\n", encoding="utf-8"
        )
        probabilities = np.ones((1, 2, 21), dtype=np.float32) / 21
        np.savez(output / f"probs/{pdb.stem}.npz", probs=probabilities)
        return SimpleNamespace(returncode=0, stdout="fixture", stderr="")

    monkeypatch.setattr("cas13_if.backends.mpnn.subprocess.run", fake_run)
    backend = ProteinMpnnBackend(
        upstream=assets["protein_upstream"],
        checkpoint=assets["protein"],
        python_executable=assets["python"],
        device="cpu",
    )
    assert backend.capabilities().hard_fixed
    with pytest.raises(RuntimeError, match="load"):
        backend.sample(_request())
    backend.load()
    candidate = backend.sample(_request())[0]
    assert candidate.sequence == "RG"
    assert candidate.fixed_positions == {0: "R"}
    assert candidate.metadata["fixed_upstream_slots_1"] == [1]
    assert len(candidate.metadata["selected_token_probabilities"]) == 2
    assert backend.metadata()["checkpoint_sha256"]
    with pytest.raises(NotImplementedError, match="not ESM"):
        backend.score(ScoreRequest(scaffold_id="x", structure_path="x", sequence="RG"))
    with pytest.raises(ValueError, match="filters"):
        backend.sample(_request().model_copy(update={"allowed_residues": {1: {"G"}}}))


def test_ligandmpnn_adapter_runs_fixed_fixture(monkeypatch, tmp_path: Path) -> None:
    assets = _assets(tmp_path)

    def fake_run(command, **_kwargs):
        if "-c" in command:
            return SimpleNamespace(
                returncode=0,
                stdout=json.dumps({"shape": [2, 21], "selected": [0.8, 0.7]}),
                stderr="",
            )
        output = Path(command[command.index("--out_folder") + 1])
        pdb = Path(command[command.index("--pdb_path") + 1])
        (output / "seqs").mkdir(parents=True)
        (output / "stats").mkdir(parents=True)
        (output / f"seqs/{pdb.stem}.fa").write_text(
            ">native\nRG\n>sample\nRG\n", encoding="utf-8"
        )
        (output / f"stats/{pdb.stem}.pt").write_bytes(b"fixture")
        return SimpleNamespace(returncode=0, stdout="fixture", stderr="")

    monkeypatch.setattr("cas13_if.backends.mpnn.subprocess.run", fake_run)
    backend = LigandMpnnBackend(
        upstream=assets["ligand_upstream"],
        ligand_checkpoint=assets["ligand"],
        protein_checkpoint=assets["protein"],
        soluble_checkpoint=assets["soluble"],
        python_executable=assets["python"],
        rna_context_chains=["R"],
        device="cpu",
    )
    assert backend.capabilities().rna_atomic_context
    backend.load()
    candidate = backend.sample(_request())[0]
    assert candidate.sequence == "RG"
    assert candidate.metadata["rna_context_chains"] == ["R"]
    assert candidate.metadata["fixed_pdb_residue_identifiers"] == ["A1"]
    assert backend.metadata()["checkpoint_sha256"]
    with pytest.raises(NotImplementedError, match="not ESM"):
        backend.score(ScoreRequest(scaffold_id="x", structure_path="x", sequence="RG"))
