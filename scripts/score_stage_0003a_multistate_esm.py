#!/usr/bin/env python
"""Score prepared Stage-0003A state projections with genuine offline ESM-IF1."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cas13_if.backends.esm_if1 import EsmIf1Backend
from cas13_if.provenance import atomic_write_text
from cas13_if.schemas import ScoreRequest


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--jobs", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--checkpoint", type=Path, required=True)
    parser.add_argument("--device", default="cpu")
    args = parser.parse_args()
    jobs = [
        json.loads(line)
        for line in args.jobs.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    backend = EsmIf1Backend(args.checkpoint, device=args.device)
    backend.load()
    output: list[dict[str, object]] = []
    for job in jobs:
        score = backend.score(
            ScoreRequest(
                scaffold_id=str(job["scaffold_id"]),
                structure_path=str(job["structure_path"]),
                sequence=str(job["sequence"]),
                protein_chains=[str(job["protein_chain"])],
                seed=20260801,
            )
        )
        output.append(
            {
                "job_id": job["job_id"],
                "pdb_id": job["pdb_id"],
                "scaffold_id": job["scaffold_id"],
                "state": job["state"],
                "resolved_positions": len(str(job["sequence"])),
                "conditional_log_likelihood": score.conditional_log_likelihood,
                "mean_log_likelihood_per_resolved_residue": score.metadata[
                    "mean_conditional_log_likelihood"
                ],
                "perplexity": score.perplexity,
                "checkpoint_sha256": score.metadata["checkpoint_sha256"],
                "device": score.metadata["device"],
                "is_mock": False,
                "evidence_level": 2,
            }
        )
    atomic_write_text(
        args.output,
        "".join(json.dumps(item, sort_keys=True) + "\n" for item in output),
    )
    print(json.dumps({"scored": len(output), "device": backend.metadata()["device"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
