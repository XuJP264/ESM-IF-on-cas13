#!/usr/bin/env python
"""Run the pinned ProteinMPNN implementation on experimental Cas13 backbone."""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import torch

from cas13_if.provenance import sha256_file
from cas13_if.structures.parser import parse_structure, protein_chain_sequence


def _sequences(path: Path) -> list[str]:
    sequences = [
        line.strip().replace("/", "")
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith(">")
    ]
    if len(sequences) < 2:
        raise RuntimeError(
            f"ProteinMPNN did not produce native + sampled FASTA: {path}"
        )
    return sequences


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    upstream = repo / "third_party/ProteinMPNN"
    checkpoint = repo / "models/proteinmpnn/v_48_020.pt"
    structure = repo / "data/experimental_structures/6e9f.pdb"
    for path in (upstream / "protein_mpnn_run.py", checkpoint, structure):
        if not path.is_file():
            raise FileNotFoundError(f"required local asset is missing: {path}")

    atoms = parse_structure(structure)
    native, _ = protein_chain_sequence(atoms, "A")
    with tempfile.TemporaryDirectory(
        prefix="proteinmpnn-smoke-", dir=repo / "results"
    ) as temporary:
        output = Path(temporary)
        command = [
            sys.executable,
            str(upstream / "protein_mpnn_run.py"),
            "--pdb_path",
            str(structure),
            "--pdb_path_chains",
            "A",
            "--out_folder",
            str(output),
            "--path_to_model_weights",
            str(checkpoint.parent),
            "--model_name",
            "v_48_020",
            "--num_seq_per_target",
            "1",
            "--batch_size",
            "1",
            "--sampling_temp",
            "0.1",
            "--seed",
            "20260731",
            "--save_probs",
            "1",
            "--suppress_print",
            "1",
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
                "ProteinMPNN failed with exit code "
                f"{completed.returncode}\nSTDOUT:\n{completed.stdout}\n"
                f"STDERR:\n{completed.stderr}"
            )
        fasta = output / "seqs/6e9f.fa"
        upstream_native, sampled = _sequences(fasta)[0], _sequences(fasta)[-1]
        probability_files = sorted((output / "probs").glob("*.npz"))
        if len(sampled) != len(upstream_native):
            raise RuntimeError(
                f"sample length {len(sampled)} differs from upstream tensor "
                f"length {len(upstream_native)}"
            )
        if upstream_native.replace("X", "") != native:
            raise RuntimeError(
                "ProteinMPNN resolved sequence does not match the strict "
                "coordinate-residue parser"
            )
        if not probability_files:
            raise RuntimeError("ProteinMPNN --save_probs produced no NPZ file")

    result = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 2,
        "implementation_commit": "8907e6671bfbfc92303b5f79c4b5e6ce47cdef57",
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
            "resolved_coordinate_length": len(native),
            "upstream_tensor_length": len(upstream_native),
            "unresolved_internal_slots": upstream_native.count("X"),
            "sample_length": len(sampled),
            "sample": sampled,
            "protein_backbone_only": True,
            "rna_atomic_context": False,
            "per_residue_probabilities_saved": True,
            "missing_coordinate_positions_are_masked_by_upstream": True,
        },
    }
    output_path = repo / "artifacts/system/proteinmpnn_real_smoke.json"
    output_path.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
