#!/usr/bin/env python
"""CLI wrapper for the Stage-0003A real structure atlas."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from cas13_if.structures.atlas import build_structure_atlas


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--config", type=Path, default=Path("configs/stage_0003a_structures.yaml")
    )
    parser.add_argument("--output", type=Path, default=Path("reports/stage_0003a"))
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    summary = build_structure_atlas(
        repo=repo,
        config_path=(
            repo / args.config if not args.config.is_absolute() else args.config
        ),
        output_dir=(
            repo / args.output if not args.output.is_absolute() else args.output
        ),
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
