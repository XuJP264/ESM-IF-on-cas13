#!/usr/bin/env python
"""Run pinned LigandMPNN with Cas13 protein design and RNA atom context."""

from __future__ import annotations

import argparse
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


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--pdb-id", default="6e9f")
    parser.add_argument("--protein-chain", default="A")
    parser.add_argument("--context-chains", default="B,C")
    parser.add_argument(
        "--model-type",
        choices=("ligand_mpnn", "protein_mpnn", "soluble_mpnn"),
        default="ligand_mpnn",
    )
    parser.add_argument("--seed", type=int, default=20260731)
    parser.add_argument("--output-path", type=Path)
    return parser.parse_args()


def _last_sequence(path: Path) -> str:
    lines = [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    ]
    if len(lines) < 2:
        raise RuntimeError(f"LigandMPNN did not produce native + sampled FASTA: {path}")
    return lines[-1].replace(":", "")


def _validate_rna_context(
    upstream: Path,
    structure: Path,
    *,
    protein_chain: str,
    context_chains: list[str],
) -> dict[str, Any]:
    expected_context_chains = set(context_chains)
    sys.path.insert(0, str(upstream))
    data_utils = __import__("data_utils")
    protein_dict, _, other_atoms, _, _ = data_utils.parse_PDB(
        str(structure),
        device="cpu",
        chains=[protein_chain, *context_chains],
        parse_all_atoms=False,
    )
    if other_atoms is None:
        raise RuntimeError("LigandMPNN parser removed all RNA context atoms")
    chain_ids = [str(value) for value in other_atoms.getChids()]
    residue_names = [str(value) for value in other_atoms.getResnames()]
    atom_count = int(other_atoms.numAtoms())
    parsed_context_chains = sorted(set(chain_ids))
    if atom_count <= 0 or not expected_context_chains.issubset(parsed_context_chains):
        raise RuntimeError(
            "RNA chains are absent from atomic context: "
            f"expected={sorted(expected_context_chains)} "
            f"parsed={parsed_context_chains}"
        )
    if set(str(value) for value in protein_dict["chain_letters"]) != {protein_chain}:
        raise RuntimeError(
            f"design protein parsing did not isolate chain {protein_chain}"
        )
    return {
        "atom_count": atom_count,
        "chains": parsed_context_chains,
        "residue_names": sorted(set(residue_names)),
    }


def main() -> int:
    arguments = _arguments()
    repo = Path(__file__).resolve().parents[1]
    upstream = repo / "third_party/LigandMPNN"
    checkpoints = {
        "ligand_mpnn": repo / "models/ligandmpnn/ligandmpnn_v_32_010_25.pt",
        "protein_mpnn": repo / "models/ligandmpnn/proteinmpnn_v_48_020.pt",
        "soluble_mpnn": repo / "models/ligandmpnn/solublempnn_v_48_020.pt",
    }
    model_type = str(arguments.model_type)
    checkpoint = checkpoints[model_type]
    protein_checkpoint = repo / "models/ligandmpnn/proteinmpnn_v_48_020.pt"
    soluble_checkpoint = repo / "models/ligandmpnn/solublempnn_v_48_020.pt"
    ligand_checkpoint = repo / "models/ligandmpnn/ligandmpnn_v_32_010_25.pt"
    pdb_id = str(arguments.pdb_id).lower()
    protein_chain = str(arguments.protein_chain)
    requested_context_chains = [
        value for value in str(arguments.context_chains).split(",") if value
    ]
    context_chains = requested_context_chains if model_type == "ligand_mpnn" else []
    structure = repo / f"data/experimental_structures/{pdb_id}.pdb"
    for path in (
        upstream / "run.py",
        checkpoint,
        protein_checkpoint,
        soluble_checkpoint,
        ligand_checkpoint,
        structure,
    ):
        if not path.is_file():
            raise FileNotFoundError(f"required local asset is missing: {path}")

    context = (
        _validate_rna_context(
            upstream,
            structure,
            protein_chain=protein_chain,
            context_chains=context_chains,
        )
        if model_type == "ligand_mpnn"
        else {"atom_count": 0, "chains": [], "residue_names": []}
    )
    atoms = parse_structure(structure)
    native, residue_keys = protein_chain_sequence(atoms, protein_chain)
    first_residue = residue_keys[0]
    fixed_identifier = (
        f"{protein_chain}{first_residue.residue_number}{first_residue.insertion_code}"
    )
    with tempfile.TemporaryDirectory(
        prefix="ligandmpnn-smoke-", dir=repo / "results"
    ) as temporary:
        output = Path(temporary)
        command = [
            sys.executable,
            str(upstream / "run.py"),
            "--model_type",
            model_type,
            "--checkpoint_ligand_mpnn",
            str(ligand_checkpoint),
            "--checkpoint_protein_mpnn",
            str(protein_checkpoint),
            "--checkpoint_soluble_mpnn",
            str(soluble_checkpoint),
            "--pdb_path",
            str(structure),
            "--out_folder",
            str(output),
            "--parse_these_chains_only",
            ",".join([protein_chain, *context_chains]),
            "--chains_to_design",
            protein_chain,
            "--fixed_residues",
            fixed_identifier,
            "--batch_size",
            "1",
            "--number_of_batches",
            "1",
            "--temperature",
            "0.1",
            "--seed",
            str(arguments.seed),
            "--ligand_mpnn_use_atom_context",
            "1" if model_type == "ligand_mpnn" else "0",
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
        fasta = output / f"seqs/{structure.stem}.fa"
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
        "model_type": model_type,
        "runtime": {
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
            "elapsed_seconds": elapsed,
        },
        f"cas13_{pdb_id}_{protein_chain.lower()}": {
            "native_length": len(native),
            "sample_length": len(sampled),
            "sample": sampled,
            "design_chain": protein_chain,
            "context_chains": context_chains,
            "rna_atomic_context_validated": model_type == "ligand_mpnn",
            "parsed_nonprotein_context": context,
            "fixed_residue": fixed_identifier,
            "fixed_residue_preserved": True,
            "per_residue_statistics_saved": True,
        },
    }
    output_path = arguments.output_path or (
        repo
        / (
            "artifacts/system/ligandmpnn_real_smoke.json"
            if pdb_id == "6e9f" and model_type == "ligand_mpnn"
            else (
                "artifacts/system/"
                f"ligandmpnn_{model_type}_{pdb_id}_{protein_chain.lower()}_real_smoke.json"
            )
        )
    )
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
