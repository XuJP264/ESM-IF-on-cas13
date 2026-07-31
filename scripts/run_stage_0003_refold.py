#!/usr/bin/env python
"""H100-only dispatcher for an explicitly configured prediction-site adapter."""

from __future__ import annotations

import argparse
import json
import os
import subprocess
from pathlib import Path
from typing import Any

import yaml

from cas13_if.provenance import atomic_write_text


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    return parser.parse_args()


def _load(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("refold configuration root must be a mapping")
    return value


def _gpu_memory_mib() -> int:
    completed = subprocess.run(
        [
            "nvidia-smi",
            "--query-gpu=memory.total",
            "--format=csv,noheader,nounits",
        ],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError("nvidia-smi failed; refusing structure prediction")
    values = [
        int(float(line)) for line in completed.stdout.splitlines() if line.strip()
    ]
    if not values:
        raise RuntimeError("nvidia-smi returned no visible GPUs")
    return max(values)


def main() -> int:
    arguments = _arguments()
    repo = Path(__file__).resolve().parents[1]
    config_path = arguments.config
    if not config_path.is_absolute():
        config_path = repo / config_path
    config = _load(config_path)
    minimum = int(config["execution"]["minimum_gpu_memory_mib"])
    measured = _gpu_memory_mib()
    if measured < minimum:
        raise RuntimeError(
            f"Stage 0003 refold requires >= {minimum} MiB visible GPU memory; "
            f"measured {measured} MiB. Local RTX 4060 execution is forbidden."
        )
    task_root = repo / str(config["execution"]["task_root"])
    checksums = task_root / "manifests/SHA256SUMS"
    manifest = task_root / "manifests/all_jobs.jsonl"
    if not checksums.is_file() or not manifest.is_file():
        raise FileNotFoundError("Stage 0003 manifests are missing")
    checked = subprocess.run(
        ["sha256sum", "--check", str(checksums)],
        cwd=task_root,
        text=True,
        capture_output=True,
        check=False,
    )
    if checked.returncode != 0:
        raise RuntimeError(
            f"Stage 0003 manifest hash failure: {checked.stdout[-2000:]}"
        )
    variable = str(config["execution"]["site_adapter_environment_variable"])
    adapter_value = os.environ.get(variable)
    if not adapter_value:
        raise RuntimeError(
            f"{variable} is not configured. Install the chosen predictor backends "
            "on the H100 node and set this to an executable site adapter accepting "
            "--manifest, --task-root, and --config. No prediction was started."
        )
    adapter = Path(adapter_value).expanduser().resolve()
    if not adapter.is_file() or not os.access(adapter, os.X_OK):
        raise RuntimeError(f"configured site adapter is not executable: {adapter}")
    command = [
        str(adapter),
        "--manifest",
        str(manifest),
        "--task-root",
        str(task_root),
        "--config",
        str(config_path),
    ]
    dispatch = {
        "command": command,
        "gpu_memory_mib": measured,
        "manifest": str(manifest),
        "manifest_line_count": sum(
            bool(line.strip())
            for line in manifest.read_text(encoding="utf-8").splitlines()
        ),
        "is_mock": False,
    }
    atomic_write_text(
        task_root / "manifests/last_dispatch.json",
        json.dumps(dispatch, indent=2, sort_keys=True) + "\n",
    )
    completed = subprocess.run(command, cwd=repo, check=False)
    return int(completed.returncode)


if __name__ == "__main__":
    raise SystemExit(main())
