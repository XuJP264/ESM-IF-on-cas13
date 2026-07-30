"""Mock prediction result writer for fixture E2E only."""

from __future__ import annotations

import json
from pathlib import Path

from cas13_if.refold.providers import PredictionJob


def write_mock_predictions(jobs: list[PredictionJob], output_root: Path) -> None:
    for job in jobs:
        if not job.is_mock:
            raise ValueError("mock writer refuses non-mock jobs")
        job_dir = output_root / job.candidate_id
        job_dir.mkdir(parents=True, exist_ok=False)
        (job_dir / "prediction.cif").write_text("data_mock\n#\n", encoding="utf-8")
        (job_dir / "pae.json").write_text(
            json.dumps({"predicted_aligned_error": [[0.0]], "is_mock": True}) + "\n",
            encoding="utf-8",
        )
        result = {
            "candidate_id": job.candidate_id,
            "provider": job.provider,
            "mean_plddt": 50.0,
            "structure_path": "prediction.cif",
            "pae_path": "pae.json",
            "is_mock": True,
        }
        (job_dir / "result.json").write_text(
            json.dumps(result, indent=2, sort_keys=True) + "\n",
            encoding="utf-8",
        )
