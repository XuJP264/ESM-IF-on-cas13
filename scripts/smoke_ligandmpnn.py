#!/usr/bin/env python
"""Run pinned LigandMPNN with Cas13 protein design and RNA atom context."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path
from typing import Any

import torch

from cas13_if.provenance import sha256_file
from cas13_if.structures.parser import parse_structure, protein_chain_sequence


def _last_sequence(path: Path) -> str:
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    ]
    if len(lines) < 2:
        raise RuntimeError(f"LigandMPNN did not produce native + sampled FASTA: {path}")
    return lines[-1].replace(":", "")


def _validate_rna_context(upstream: Path, structure: Path) -> dict[str, Any]:
    sys.path.insert(0, str(upstream))
    data_utils = __import__("data_utils")
    protein_dict, _, other_atoms, _, _ = data_utils.parse_PDB(
        str(structure),
        device="cpu",
        chains=["A", "B", "C"],
        parse_all_atoms=False,
    )
    if other_atoms is None:
        raise RuntimeError("LigandMPNN parser removed all RNA context atoms")
    chain_ids = [str(value) for value in other_atoms.getChids()]
    residue_names = [str(value) for value in other_atoms.getResnames()]
    atom_count = int(other_atoms.numAtoms())
    context_chains = sorted(set(chain_ids))
    if atom_count <= 0 or not {"B", "C"}.issubset(context_chains):
        raise RuntimeError(
            f"RNA chains B/C are absent from atomic context: {context_chains}"
        )
    if set(str(value) for value in protein_dict["chain_letters"]) != {"A"}:
        raise RuntimeError("design protein parsing did not isolate chain A")
    return {
        "atom_count": atom_count,
        "chains": context_chains,
        "residue_names": sorted(set(residue_names)),
    }


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    upstream = repo / "third_party/LigandMPNN"
    checkpoint = repo / "models/ligandmpnn/ligandmpnn_v_32_010_25.pt"
    protein_checkpoint = repo / "models/ligandmpnn/proteinmpnn_v_48_020.pt"
    structure = repo / "data/experimental_structures/6e9f.pdb"
    for path in (
        upstream / "run.py",
        checkpoint,
        protein_checkpoint,
        structure,
    ):
        if not path.is_file():
            raise FileNotFoundError(f"required local asset is missing: {path}")

    context = _validate_rna_context(upstream, structure)
    atoms = parse_structure(structure)
    native, residue_keys = protein_chain_sequence(atoms, "A")
    first_residue = residue_keys[0]
    fixed_identifier = f"A{first_residue.residue_number}{first_residue.insertion_code}"
    with tempfile.TemporaryDirectory(
        prefix="ligandmpnn-smoke-", dir=repo / "results"
    ) as temporary:
        output = Path(temporary)
        command = [
            sys.executable,
            str(upstream / "run.py"),
            "--model_type",
            "ligand_mpnn",
            "--checkpoint_ligand_mpnn",
            str(checkpoint),
            "--checkpoint_protein_mpnn",
            str(protein_checkpoint),
            "--pdb_path",
            str(structure),
            "--out_folder",
            str(output),
            "--parse_these_chains_only",
            "A,B,C",
            "--chains_to_design",
            "A",
            "--fixed_residues",
            fixed_identifier,
            "--batch_size",
            "1",
            "--number_of_batches",
            "1",
            "--temperature",
            "0.1",
            "--seed",
            "20260731",
            "--ligand_mpnn_use_atom_context",
            "1",
            "--save_stats",
            "1",
            "--verbose",
            "0",
        ]
        started = time.perf_counter()
        completed = subprocess.run(
            command,
            cwd=upstream,
            text=True,
            capture_output=True,
            check=False,
        )
        elapsed = time.perf_counter() - started
        if completed.returncode != 0:
            raise RuntimeError(
                "LigandMPNN failed with exit code "
                f"{completed.returncode}\nSTDOUT:\n{completed.stdout}\n"
                f"STDERR:\n{completed.stderr}"
            )
        fasta = output / "seqs/6e9f.fa"
        sampled = _last_sequence(fasta)
        stats_files = sorted((output / "stats").glob("*.pt"))
        backbones = sorted((output / "backbones").glob("*.pdb"))
        if len(sampled) != len(native):
            raise RuntimeError(
                f"sample length {len(sampled)} differs from chain A "
                f"length {len(native)}"
            )
        if sampled[0] != native[0]:
            raise RuntimeError(
                f"fixed residue {fixed_identifier} was not preserved: "
                f"{sampled[0]} != {native[0]}"
            )
        if not stats_files or not backbones:
            raise RuntimeError("LigandMPNN did not save stats and backbone outputs")

    result = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 2,
        "implementation_commit": "26ec57ac976ade5379920dbd43c7f97a91cf82de",
        "checkpoint": {
            "path": str(checkpoint.relative_to(repo)),
            "sha256": sha256_file(checkpoint),
            "size_bytes": checkpoint.stat().st_size,
        },
        "runtime": {
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
            "elapsed_seconds": elapsed,
        },
        "cas13_6e9f_a": {
            "native_length": len(native),
            "sample_length": len(sampled),
            "sample": sampled,
            "design_chain": "A",
            "context_chains": ["B", "C"],
            "rna_atomic_context_validated": True,
            "parsed_nonprotein_context": context,
            "fixed_residue": fixed_identifier,
            "fixed_residue_preserved": True,
            "per_residue_statistics_saved": True,
        },
    }
    output_path = repo / "artifacts/system/ligandmpnn_real_smoke.json"
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
