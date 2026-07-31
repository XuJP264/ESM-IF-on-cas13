#!/usr/bin/env python
"""Build the curated, assay-stratified Cas13d variant dataset."""

from __future__ import annotations

import json
from pathlib import Path

from cas13_if.data.variants import build_variant_dataset


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    summary = build_variant_dataset(
        manifest_path=repo / "data/manifests/cas13d_variant_activity_curated.yaml",
        scaffold_csv=repo / "reports/stage_0003a/scaffolds.csv",
        output_dir=repo / "data/processed",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
