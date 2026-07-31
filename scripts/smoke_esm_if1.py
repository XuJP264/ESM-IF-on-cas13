#!/usr/bin/env python
"""Run genuine ESM-IF1 toy and Cas13 inference from local assets."""

from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import torch

from cas13_if.backends.esm_if1 import EsmIf1Backend, EsmIf1ConstrainedBackend
from cas13_if.provenance import sha256_file
from cas13_if.schemas import SampleRequest, ScoreRequest
from cas13_if.structures.parser import parse_structure, protein_chain_sequence


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    checkpoint = repo / "models/esm_if1/esm_if1_gvp4_t16_142M_UR50.pt"
    fixture = repo / "tests/fixtures/minimal_complex.pdb"
    cas13 = repo / "data/experimental_structures/6e9f.cif"
    cas13a = repo / "data/experimental_structures/5xwp.cif"
    for path in (checkpoint, fixture, cas13, cas13a):
        if not path.is_file():
            raise FileNotFoundError(f"required local asset is missing: {path}")

    device = "cuda" if torch.cuda.is_available() else "cpu"
    started = time.perf_counter()
    backend = EsmIf1Backend(checkpoint, device=device)
    backend.load()
    load_seconds = time.perf_counter() - started

    fixture_atoms = parse_structure(fixture)
    fixture_sequence, _ = protein_chain_sequence(fixture_atoms, "A")
    toy_score = backend.score(
        ScoreRequest(
            scaffold_id="minimal-complex",
            structure_path=str(fixture),
            sequence=fixture_sequence,
            protein_chains=["A"],
        )
    )
    toy_sample = backend.sample(
        SampleRequest(
            scaffold_id="minimal-complex",
            structure_path=str(fixture),
            parent_sequence=fixture_sequence,
            count=1,
            temperature=0.1,
            seed=20260731,
            protein_chains=["A"],
        )
    )[0]

    constrained = EsmIf1ConstrainedBackend(checkpoint, device=device)
    # Reuse the already loaded model to avoid duplicating 142M parameters.
    constrained.__dict__.update(backend.__dict__)
    constrained_candidate = constrained.sample(
        SampleRequest(
            scaffold_id="minimal-complex",
            structure_path=str(fixture),
            parent_sequence=fixture_sequence,
            count=1,
            temperature=0.5,
            seed=20260731,
            fixed_positions={
                index: token for index, token in enumerate(fixture_sequence)
            },
            protein_chains=["A"],
        )
    )[0]

    cas13_atoms = parse_structure(cas13)
    cas13_sequence, _ = protein_chain_sequence(cas13_atoms, "A")
    cas13_started = time.perf_counter()
    cas13_score = backend.score(
        ScoreRequest(
            scaffold_id="6E9F-A",
            structure_path=str(cas13),
            sequence=cas13_sequence,
            protein_chains=["A"],
        )
    )
    cas13_seconds = time.perf_counter() - cas13_started
    cas13a_atoms = parse_structure(cas13a)
    cas13a_sequence, _ = protein_chain_sequence(cas13a_atoms, "A")
    cas13a_started = time.perf_counter()
    cas13a_score = backend.score(
        ScoreRequest(
            scaffold_id="5XWP-A",
            structure_path=str(cas13a),
            sequence=cas13a_sequence,
            protein_chains=["A"],
        )
    )
    cas13a_seconds = time.perf_counter() - cas13a_started

    result = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 2,
        "checkpoint": {
            "path": str(checkpoint.relative_to(repo)),
            "sha256": sha256_file(checkpoint),
            "size_bytes": checkpoint.stat().st_size,
        },
        "runtime": {
            "device": device,
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
            "python": sys.version,
            "platform": platform.platform(),
            "load_seconds": load_seconds,
        },
        "toy": {
            "sequence": fixture_sequence,
            "conditional_log_likelihood": toy_score.conditional_log_likelihood,
            "perplexity": toy_score.perplexity,
            "sample": toy_sample.sequence,
            "fixed_all_recovery": constrained_candidate.sequence == fixture_sequence,
            "fixed_position_violations": sum(
                constrained_candidate.sequence[index] != token
                for index, token in enumerate(fixture_sequence)
            ),
        },
        "cas13_6e9f_a": {
            "length": len(cas13_sequence),
            "conditional_log_likelihood": cas13_score.conditional_log_likelihood,
            "mean_conditional_log_likelihood": cas13_score.metadata[
                "mean_conditional_log_likelihood"
            ],
            "perplexity": cas13_score.perplexity,
            "inference_seconds": cas13_seconds,
            "protein_chains": ["A"],
            "rna_chains_excluded_from_esm_if1": ["B", "C"],
        },
        "cas13_5xwp_a": {
            "length": len(cas13a_sequence),
            "conditional_log_likelihood": cas13a_score.conditional_log_likelihood,
            "mean_conditional_log_likelihood": cas13a_score.metadata[
                "mean_conditional_log_likelihood"
            ],
            "perplexity": cas13a_score.perplexity,
            "inference_seconds": cas13a_seconds,
            "protein_chains": ["A"],
            "rna_chains_excluded_from_esm_if1": ["C", "D"],
        },
    }
    output = repo / "artifacts/system/esm_if1_real_smoke.json"
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    print(json.dumps(result, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
