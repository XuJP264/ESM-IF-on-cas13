#!/usr/bin/env python
"""Run one genuine hard-fixed Cas13 ESM-IF1 sample on CPU/GPU."""

from __future__ import annotations

import json
import time
from pathlib import Path

import torch
import yaml

from cas13_if.backends.esm_if1 import EsmIf1ConstrainedBackend
from cas13_if.provenance import sha256_file
from cas13_if.schemas import SampleRequest
from cas13_if.structures.parser import parse_structure, protein_chain_sequence


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    checkpoint = repo / "models/esm_if1/esm_if1_gvp4_t16_142M_UR50.pt"
    structure = repo / "data/experimental_structures/6e9f.cif"
    functional_path = repo / "data/manifests/cas13_functional_residues.yaml"
    functional = yaml.safe_load(functional_path.read_text(encoding="utf-8"))
    entries = functional["structures"]["6E9F"]["residues"]
    atoms = parse_structure(structure)
    deposited, keys = protein_chain_sequence(atoms, "A")
    number_to_index = {key.residue_number: index for index, key in enumerate(keys)}
    fixed = {
        number_to_index[int(entry["pdb_residue_number"])]: str(
            entry["biological_amino_acid"]
        )
        for entry in entries
    }
    biological_tokens = list(deposited)
    for index, token in fixed.items():
        biological_tokens[index] = token
    biological = "".join(biological_tokens)
    device = "cuda" if torch.cuda.is_available() else "cpu"
    backend = EsmIf1ConstrainedBackend(checkpoint, device=device)
    backend.load()
    started = time.perf_counter()
    candidate = backend.sample(
        SampleRequest(
            scaffold_id="6E9F-A-biological-restored",
            structure_path=str(structure),
            parent_sequence=biological,
            count=1,
            temperature=0.1,
            seed=20260731,
            fixed_positions=fixed,
            protein_chains=["A"],
        )
    )[0]
    elapsed = time.perf_counter() - started
    violations = sum(
        candidate.sequence[index] != token for index, token in fixed.items()
    )
    if violations:
        raise RuntimeError(f"hard-fixed sampling produced {violations} violation(s)")
    result = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 2,
        "checkpoint_sha256": sha256_file(checkpoint),
        "structure": "6E9F",
        "chain": "A",
        "sequence_length": len(candidate.sequence),
        "device": device,
        "temperature": candidate.temperature,
        "seed": candidate.seed,
        "fixed_positions": fixed,
        "fixed_position_violations": violations,
        "elapsed_seconds": elapsed,
        "candidate": candidate.model_dump(mode="json"),
        "claim_scope": "Level 2 computational candidate; not a validated Cas13",
    }
    output = repo / "artifacts/system/esm_if1_constrained_6e9f_smoke.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(
        json.dumps(
            {key: value for key, value in result.items() if key != "candidate"},
            indent=2,
            sort_keys=True,
        )
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
