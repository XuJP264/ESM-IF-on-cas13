"""Offline adapters for pinned ProteinMPNN and LigandMPNN command runners."""

from __future__ import annotations

import hashlib
import json
import os
import subprocess
import tempfile
from pathlib import Path
from typing import Any

import numpy as np

from cas13_if.backends.base import InverseFoldingBackend
from cas13_if.data.fasta import iter_fasta
from cas13_if.provenance import atomic_write_text, sha256_file
from cas13_if.schemas import (
    BackendCapabilities,
    Candidate,
    EvidenceLevel,
    SampleRequest,
    ScoreRequest,
    ScoreResult,
)
from cas13_if.structures.parser import parse_structure, protein_chain_sequence

AA_TO_THREE = {
    "A": "ALA",
    "C": "CYS",
    "D": "ASP",
    "E": "GLU",
    "F": "PHE",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "K": "LYS",
    "L": "LEU",
    "M": "MET",
    "N": "ASN",
    "P": "PRO",
    "Q": "GLN",
    "R": "ARG",
    "S": "SER",
    "T": "THR",
    "V": "VAL",
    "W": "TRP",
    "Y": "TYR",
}


def _candidate_digest(backend: str, request: SampleRequest, seed: int) -> str:
    payload = {
        "backend": backend,
        "scaffold_id": request.scaffold_id,
        "structure_path": str(Path(request.structure_path).name),
        "parent_sha256": hashlib.sha256(
            request.parent_sequence.encode("ascii")
        ).hexdigest(),
        "fixed_positions": sorted(request.fixed_positions.items()),
        "protein_chains": request.protein_chains,
        "temperature": request.temperature,
        "seed": seed,
    }
    return hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:12]


def _subprocess_environment(device: str) -> dict[str, str]:
    environment = os.environ.copy()
    environment["PYTHONNOUSERSITE"] = "1"
    if device == "cpu":
        environment["CUDA_VISIBLE_DEVICES"] = ""
    return environment


def _restore_pdb_sequence(
    *,
    source: Path,
    destination: Path,
    chain: str,
    fixed_positions: dict[int, str],
) -> tuple[str, list[Any]]:
    atoms = parse_structure(source)
    deposited, residue_keys = protein_chain_sequence(atoms, chain)
    if any(index >= len(residue_keys) for index in fixed_positions):
        raise ValueError("fixed position exceeds strict coordinate sequence")
    replacement_by_residue = {
        (
            residue_keys[index].residue_number,
            residue_keys[index].insertion_code,
        ): AA_TO_THREE[token.upper()]
        for index, token in fixed_positions.items()
    }
    lines: list[str] = []
    for line in source.read_text(encoding="utf-8").splitlines(keepends=True):
        if line.startswith(("ATOM  ", "HETATM")) and line[21:22].strip() == chain:
            try:
                residue_number = int(line[22:26])
            except ValueError:
                lines.append(line)
                continue
            insertion_code = line[26:27].strip()
            replacement = replacement_by_residue.get((residue_number, insertion_code))
            if replacement is not None:
                line = f"{line[:17]}{replacement:>3}{line[20:]}"
        lines.append(line)
    atomic_write_text(destination, "".join(lines))
    restored_atoms = parse_structure(destination)
    restored, restored_keys = protein_chain_sequence(restored_atoms, chain)
    expected = list(deposited)
    for index, token in fixed_positions.items():
        expected[index] = token.upper()
    original_ids = [
        (key.chain_id, key.residue_number, key.insertion_code) for key in residue_keys
    ]
    restored_ids = [
        (key.chain_id, key.residue_number, key.insertion_code) for key in restored_keys
    ]
    if restored != "".join(expected) or restored_ids != original_ids:
        raise RuntimeError("temporary biological-residue restoration changed mapping")
    return restored, residue_keys


def _fasta_sequences(path: Path) -> list[str]:
    sequences = [
        sequence.replace("/", "").replace(":", "") for _, sequence in iter_fasta(path)
    ]
    if len(sequences) < 2:
        raise RuntimeError(f"upstream output lacks native and sample rows: {path}")
    return sequences


def _proteinmpnn_slot_map(residue_keys: list[Any]) -> dict[Any, int]:
    by_number: dict[int, list[Any]] = {}
    for key in residue_keys:
        by_number.setdefault(int(key.residue_number), []).append(key)
    first = min(by_number)
    last = max(by_number)
    mapping: dict[Any, int] = {}
    slot = 0
    for residue_number in range(first, last + 1):
        keys = sorted(
            by_number.get(residue_number, []), key=lambda item: item.insertion_code
        )
        if not keys:
            slot += 1
            continue
        for key in keys:
            mapping[key] = slot
            slot += 1
    return mapping


class ProteinMpnnBackend(InverseFoldingBackend):
    """Protein-backbone-only ProteinMPNN adapter with hard-fixed restoration."""

    backend_name = "proteinmpnn"

    def __init__(
        self,
        *,
        upstream: Path,
        checkpoint: Path,
        python_executable: Path,
        device: str = "cpu",
    ) -> None:
        self.upstream = upstream.resolve()
        self.checkpoint = checkpoint.resolve()
        self.python_executable = python_executable.resolve()
        self.device = device
        self._loaded = False

    def capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            scoring=False,
            sampling=True,
            hard_fixed=True,
            per_residue_probabilities=True,
        )

    def load(self) -> None:
        for path in (
            self.upstream / "protein_mpnn_run.py",
            self.checkpoint,
            self.python_executable,
        ):
            if not path.is_file():
                raise FileNotFoundError(f"required ProteinMPNN asset missing: {path}")
        self._loaded = True

    def score(self, request: ScoreRequest) -> ScoreResult:
        del request
        raise NotImplementedError(
            "ProteinMPNN intrinsic scores are not ESM conditional log-likelihoods"
        )

    def sample(self, request: SampleRequest) -> list[Candidate]:
        if not self._loaded:
            raise RuntimeError("ProteinMpnnBackend.load() must be called first")
        if request.allowed_residues:
            raise ValueError(
                "ProteinMPNN adapter does not silently approximate filters"
            )
        if len(request.protein_chains) != 1:
            raise ValueError(
                "ProteinMPNN comparison requires one explicit design chain"
            )
        source = Path(request.structure_path)
        chain = request.protein_chains[0]
        candidates: list[Candidate] = []
        for sample_index in range(request.count):
            seed = request.seed + sample_index
            with tempfile.TemporaryDirectory(prefix="cas13-if-proteinmpnn-") as tmp:
                root = Path(tmp)
                restored_path = root / f"{source.stem}_restored.pdb"
                restored, residue_keys = _restore_pdb_sequence(
                    source=source,
                    destination=restored_path,
                    chain=chain,
                    fixed_positions=request.fixed_positions,
                )
                if restored != request.parent_sequence:
                    raise ValueError(
                        "ProteinMPNN biological-restored PDB sequence differs "
                        "from parent"
                    )
                slot_map = _proteinmpnn_slot_map(residue_keys)
                fixed_slots = [
                    slot_map[residue_keys[index]] + 1
                    for index in sorted(request.fixed_positions)
                ]
                fixed_path = root / "fixed_positions.jsonl"
                fixed_payload = {restored_path.stem: {chain: fixed_slots}}
                atomic_write_text(
                    fixed_path, json.dumps(fixed_payload, sort_keys=True) + "\n"
                )
                output = root / "output"
                command = [
                    str(self.python_executable),
                    str(self.upstream / "protein_mpnn_run.py"),
                    "--pdb_path",
                    str(restored_path),
                    "--pdb_path_chains",
                    chain,
                    "--out_folder",
                    str(output),
                    "--path_to_model_weights",
                    str(self.checkpoint.parent),
                    "--model_name",
                    self.checkpoint.stem,
                    "--fixed_positions_jsonl",
                    str(fixed_path),
                    "--num_seq_per_target",
                    "1",
                    "--batch_size",
                    "1",
                    "--sampling_temp",
                    str(request.temperature),
                    "--seed",
                    str(seed),
                    "--save_probs",
                    "1",
                    "--suppress_print",
                    "1",
                ]
                completed = subprocess.run(
                    command,
                    cwd=self.upstream,
                    env=_subprocess_environment(self.device),
                    text=True,
                    capture_output=True,
                    check=False,
                )
                if completed.returncode != 0:
                    raise RuntimeError(
                        "ProteinMPNN failed: "
                        f"exit={completed.returncode} stderr={completed.stderr[-2000:]}"
                    )
                fasta = output / f"seqs/{restored_path.stem}.fa"
                upstream_native, upstream_sample = _fasta_sequences(fasta)[-2:]
                if upstream_native.replace("X", "") != request.parent_sequence:
                    raise RuntimeError(
                        "ProteinMPNN native tensor does not match parent"
                    )
                sequence = upstream_sample.replace("X", "")
                if len(sequence) != len(request.parent_sequence):
                    raise RuntimeError("ProteinMPNN resolved sample length mismatch")
                probability_path = next((output / "probs").glob("*.npz"))
                with np.load(probability_path, allow_pickle=True) as probabilities:
                    matrix = np.asarray(probabilities["probs"])[0]
                    resolved_slots = [
                        index
                        for index, token in enumerate(upstream_native)
                        if token != "X"
                    ]
                    alphabet = "ACDEFGHIKLMNPQRSTVWYX"
                    selected_probabilities = [
                        float(matrix[slot, alphabet.index(token)])
                        for slot, token in zip(resolved_slots, sequence, strict=True)
                    ]
                digest = _candidate_digest(self.backend_name, request, seed)
                candidates.append(
                    Candidate(
                        candidate_id=(
                            f"{self.backend_name}-{request.scaffold_id}-{digest}-0000"
                        ),
                        scaffold_id=request.scaffold_id,
                        backend=self.backend_name,
                        sequence=sequence,
                        parent_sequence=request.parent_sequence,
                        seed=seed,
                        temperature=request.temperature,
                        is_mock=False,
                        evidence_level=EvidenceLevel.INVERSE_FOLDING_COMPATIBILITY,
                        fixed_positions=request.fixed_positions,
                        metadata={
                            "implementation": "dauparas/ProteinMPNN",
                            "checkpoint_sha256": sha256_file(self.checkpoint),
                            "device": self.device,
                            "protein_backbone_only": True,
                            "rna_atomic_context": False,
                            "biological_identity_restored_in_temporary_pdb": True,
                            "upstream_tensor_length": len(upstream_native),
                            "unresolved_slots": upstream_native.count("X"),
                            "fixed_upstream_slots_1": fixed_slots,
                            "selected_token_probabilities": selected_probabilities,
                            "upstream_stdout_tail": completed.stdout[-1000:],
                        },
                    )
                )
        return candidates

    def metadata(self) -> dict[str, Any]:
        return {
            "backend": self.backend_name,
            "loaded": self._loaded,
            "checkpoint": str(self.checkpoint),
            "checkpoint_sha256": (
                sha256_file(self.checkpoint) if self.checkpoint.is_file() else None
            ),
            "device": self.device,
            "is_mock": False,
        }


class LigandMpnnBackend(InverseFoldingBackend):
    """LigandMPNN adapter retaining declared RNA chains as atomic context."""

    backend_name = "ligandmpnn"

    def __init__(
        self,
        *,
        upstream: Path,
        ligand_checkpoint: Path,
        protein_checkpoint: Path,
        soluble_checkpoint: Path,
        python_executable: Path,
        rna_context_chains: list[str],
        device: str = "cpu",
    ) -> None:
        self.upstream = upstream.resolve()
        self.ligand_checkpoint = ligand_checkpoint.resolve()
        self.protein_checkpoint = protein_checkpoint.resolve()
        self.soluble_checkpoint = soluble_checkpoint.resolve()
        self.python_executable = python_executable.resolve()
        self.rna_context_chains = list(rna_context_chains)
        self.device = device
        self._loaded = False

    def capabilities(self) -> BackendCapabilities:
        return BackendCapabilities(
            scoring=False,
            sampling=True,
            rna_atomic_context=True,
            hard_fixed=True,
            per_residue_probabilities=True,
        )

    def load(self) -> None:
        for path in (
            self.upstream / "run.py",
            self.ligand_checkpoint,
            self.protein_checkpoint,
            self.soluble_checkpoint,
            self.python_executable,
        ):
            if not path.is_file():
                raise FileNotFoundError(f"required LigandMPNN asset missing: {path}")
        self._loaded = True

    def score(self, request: ScoreRequest) -> ScoreResult:
        del request
        raise NotImplementedError(
            "LigandMPNN intrinsic scores are not ESM conditional log-likelihoods"
        )

    def sample(self, request: SampleRequest) -> list[Candidate]:
        if not self._loaded:
            raise RuntimeError("LigandMpnnBackend.load() must be called first")
        if request.allowed_residues:
            raise ValueError("LigandMPNN adapter does not silently approximate filters")
        if len(request.protein_chains) != 1:
            raise ValueError("LigandMPNN comparison requires one design protein chain")
        source = Path(request.structure_path)
        chain = request.protein_chains[0]
        candidates: list[Candidate] = []
        for sample_index in range(request.count):
            seed = request.seed + sample_index
            with tempfile.TemporaryDirectory(prefix="cas13-if-ligandmpnn-") as tmp:
                root = Path(tmp)
                restored_path = root / f"{source.stem}_restored.pdb"
                restored, residue_keys = _restore_pdb_sequence(
                    source=source,
                    destination=restored_path,
                    chain=chain,
                    fixed_positions=request.fixed_positions,
                )
                if restored != request.parent_sequence:
                    raise ValueError(
                        "LigandMPNN biological-restored PDB sequence differs "
                        "from parent"
                    )
                fixed_identifiers = [
                    f"{chain}{residue_keys[index].residue_number}"
                    f"{residue_keys[index].insertion_code}"
                    for index in sorted(request.fixed_positions)
                ]
                output = root / "output"
                command = [
                    str(self.python_executable),
                    str(self.upstream / "run.py"),
                    "--model_type",
                    "ligand_mpnn",
                    "--checkpoint_ligand_mpnn",
                    str(self.ligand_checkpoint),
                    "--checkpoint_protein_mpnn",
                    str(self.protein_checkpoint),
                    "--checkpoint_soluble_mpnn",
                    str(self.soluble_checkpoint),
                    "--pdb_path",
                    str(restored_path),
                    "--out_folder",
                    str(output),
                    "--parse_these_chains_only",
                    ",".join([chain, *self.rna_context_chains]),
                    "--chains_to_design",
                    chain,
                    "--fixed_residues",
                    " ".join(fixed_identifiers),
                    "--batch_size",
                    "1",
                    "--number_of_batches",
                    "1",
                    "--temperature",
                    str(request.temperature),
                    "--seed",
                    str(seed),
                    "--ligand_mpnn_use_atom_context",
                    "1",
                    "--save_stats",
                    "1",
                    "--verbose",
                    "0",
                ]
                completed = subprocess.run(
                    command,
                    cwd=self.upstream,
                    env=_subprocess_environment(self.device),
                    text=True,
                    capture_output=True,
                    check=False,
                )
                if completed.returncode != 0:
                    raise RuntimeError(
                        "LigandMPNN failed: "
                        f"exit={completed.returncode} stderr={completed.stderr[-2000:]}"
                    )
                fasta = output / f"seqs/{restored_path.stem}.fa"
                sequence = _fasta_sequences(fasta)[-1]
                if len(sequence) != len(request.parent_sequence):
                    raise RuntimeError("LigandMPNN sample length mismatch")
                stats_path = next((output / "stats").glob("*.pt"))
                torch_load_script = (
                    "import json,torch;d=torch.load(r'"
                    + str(stats_path)
                    + "',map_location='cpu',weights_only=True);"
                    "p=d['sampling_probs'][0];s=d['generated_sequences'][0];"
                    "print(json.dumps({'shape':list(p.shape),"
                    "'selected':[float(p[i,int(s[i])]) for i in range(len(s))]}))"
                )
                stats_completed = subprocess.run(
                    [str(self.python_executable), "-c", torch_load_script],
                    env=_subprocess_environment(self.device),
                    text=True,
                    capture_output=True,
                    check=False,
                )
                if stats_completed.returncode != 0:
                    raise RuntimeError(
                        "failed to inspect LigandMPNN stats: "
                        f"{stats_completed.stderr[-2000:]}"
                    )
                stats = json.loads(stats_completed.stdout)
                selected_probabilities = [float(value) for value in stats["selected"]]
                if len(selected_probabilities) != len(sequence):
                    raise RuntimeError("LigandMPNN per-residue probability mismatch")
                digest = _candidate_digest(self.backend_name, request, seed)
                candidates.append(
                    Candidate(
                        candidate_id=(
                            f"{self.backend_name}-{request.scaffold_id}-{digest}-0000"
                        ),
                        scaffold_id=request.scaffold_id,
                        backend=self.backend_name,
                        sequence=sequence,
                        parent_sequence=request.parent_sequence,
                        seed=seed,
                        temperature=request.temperature,
                        is_mock=False,
                        evidence_level=EvidenceLevel.INVERSE_FOLDING_COMPATIBILITY,
                        fixed_positions=request.fixed_positions,
                        metadata={
                            "implementation": "dauparas/LigandMPNN",
                            "checkpoint_sha256": sha256_file(self.ligand_checkpoint),
                            "device": self.device,
                            "rna_atomic_context": True,
                            "rna_context_chains": self.rna_context_chains,
                            "biological_identity_restored_in_temporary_pdb": True,
                            "fixed_pdb_residue_identifiers": fixed_identifiers,
                            "selected_token_probabilities": selected_probabilities,
                            "stats_probability_shape": stats["shape"],
                            "upstream_stdout_tail": completed.stdout[-1000:],
                        },
                    )
                )
        return candidates

    def metadata(self) -> dict[str, Any]:
        return {
            "backend": self.backend_name,
            "loaded": self._loaded,
            "checkpoint": str(self.ligand_checkpoint),
            "checkpoint_sha256": (
                sha256_file(self.ligand_checkpoint)
                if self.ligand_checkpoint.is_file()
                else None
            ),
            "device": self.device,
            "rna_context_chains": self.rna_context_chains,
            "is_mock": False,
        }
