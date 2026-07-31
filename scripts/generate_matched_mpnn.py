#!/usr/bin/env python
"""Generate real ProteinMPNN and RNA-context LigandMPNN proposals."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import yaml

from cas13_if.backends.mpnn import LigandMpnnBackend, ProteinMpnnBackend
from cas13_if.provenance import atomic_write_text
from cas13_if.schemas import SampleRequest
from cas13_if.structures.parser import parse_structure, protein_chain_sequence


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def _load(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"configuration is not a mapping: {path}")
    return value


def _repo_path(repo: Path, value: Any) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else repo / path


def main() -> int:
    arguments = _arguments()
    repo = Path(__file__).resolve().parents[1]
    config = _load(arguments.config)
    inputs = config["inputs"]
    models = config["models"]
    sampling = config["sampling"]
    execution = config["execution"]
    environments = config["environments"]
    structure = _repo_path(repo, inputs["structure_pdb"])
    chain = str(inputs["chain_id"])
    deposited, residue_keys = protein_chain_sequence(parse_structure(structure), chain)
    functional = _load(_repo_path(repo, inputs["functional_manifest"]))
    entries = functional["structures"]["6E9F"]["residues"]
    by_pdb = {
        (key.residue_number, key.insertion_code): index
        for index, key in enumerate(residue_keys)
    }
    fixed = {
        by_pdb[
            (int(entry["pdb_residue_number"]), str(entry.get("insertion_code") or ""))
        ]: str(entry["biological_amino_acid"]).upper()
        for entry in entries
    }
    parent_tokens = list(deposited)
    for index, token in fixed.items():
        parent_tokens[index] = token
    parent = "".join(parent_tokens)
    python = _repo_path(repo, environments["ligandmpnn_python"])
    device = str(execution.get("device", "cpu"))
    backends = {
        "proteinmpnn": ProteinMpnnBackend(
            upstream=_repo_path(repo, models["proteinmpnn_upstream"]),
            checkpoint=_repo_path(repo, models["proteinmpnn_checkpoint"]),
            python_executable=python,
            device=device,
        ),
        "ligandmpnn": LigandMpnnBackend(
            upstream=_repo_path(repo, models["ligandmpnn_upstream"]),
            ligand_checkpoint=_repo_path(repo, models["ligandmpnn_checkpoint"]),
            protein_checkpoint=_repo_path(
                repo, models["ligandmpnn_protein_checkpoint"]
            ),
            soluble_checkpoint=_repo_path(
                repo, models["ligandmpnn_soluble_checkpoint"]
            ),
            python_executable=python,
            rna_context_chains=[
                str(value)
                for value in [
                    *inputs["crrna_chains"],
                    *inputs["target_rna_chains"],
                ]
            ],
            device=device,
        ),
    }
    rows: list[dict[str, Any]] = []
    for method, backend in backends.items():
        backend.load()
        temperature = float(sampling["temperatures"][method])
        for seed_block_value in sampling["seed_blocks"]:
            seed_block = int(seed_block_value)
            for proposal_index in range(int(sampling["proposals_per_seed"])):
                actual_seed = seed_block + proposal_index
                candidate = backend.sample(
                    SampleRequest(
                        scaffold_id="6E9F-A",
                        structure_path=str(structure),
                        parent_sequence=parent,
                        count=1,
                        temperature=temperature,
                        seed=actual_seed,
                        fixed_positions=fixed,
                        protein_chains=[chain],
                    )
                )[0]
                payload = candidate.model_dump(mode="json")
                payload["candidate_id"] = (
                    f"{method}-seed{seed_block}-proposal{proposal_index}-"
                    f"{payload['candidate_id']}"
                )
                payload["metadata"] = {
                    **payload["metadata"],
                    "method": method,
                    "seed_block": seed_block,
                    "proposal_index": proposal_index,
                    "actual_model_seed": actual_seed,
                }
                rows.append(
                    {
                        "method": method,
                        "seed_block": seed_block,
                        "proposal_index": proposal_index,
                        "candidate": payload,
                    }
                )
    atomic_write_text(
        arguments.output,
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
    )
    summary = {
        "is_mock": False,
        "proposal_count": len(rows),
        "fixed_positions": fixed,
        "backends": {name: backend.metadata() for name, backend in backends.items()},
    }
    atomic_write_text(
        arguments.output.with_suffix(".summary.json"),
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
