#!/usr/bin/env python
"""Score all selected matched candidates with one offline ESM-IF1 model load."""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

import yaml

from cas13_if.backends.esm_if1 import EsmIf1Backend
from cas13_if.provenance import atomic_write_text
from cas13_if.schemas import ScoreRequest


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--candidates", type=Path, required=True)
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
    checkpoint = _repo_path(repo, config["models"]["esm_if1_checkpoint"])
    structure = _repo_path(repo, config["inputs"]["structure_pdb"])
    chain = str(config["inputs"]["chain_id"])
    backend = EsmIf1Backend(
        checkpoint, device=str(config["execution"].get("device", "cpu"))
    )
    backend.load()
    rows: list[dict[str, Any]] = []
    for line in arguments.candidates.read_text(encoding="utf-8").splitlines():
        if not line.strip():
            continue
        candidate = json.loads(line)
        score = backend.score(
            ScoreRequest(
                scaffold_id="6E9F-A",
                structure_path=str(structure),
                sequence=str(candidate["sequence"]),
                protein_chains=[chain],
                seed=int(candidate["actual_model_seed"]),
            )
        )
        rows.append(
            {
                "candidate_id": candidate["candidate_id"],
                "score": score.model_dump(mode="json"),
            }
        )
    atomic_write_text(
        arguments.output,
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
    )
    print(json.dumps({"scored": len(rows), "device": backend.metadata()["device"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
