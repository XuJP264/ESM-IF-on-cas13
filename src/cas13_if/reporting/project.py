"""Aggregate real and mock artifacts into an evidence-scoped project report."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any

from cas13_if.provenance import atomic_write_text


def _read_json(path: Path) -> Any:
    return json.loads(path.read_text(encoding="utf-8"))


def _latest(paths: list[Path]) -> Path | None:
    return max(paths, key=lambda path: path.stat().st_mtime) if paths else None


def build_project_report(
    *,
    repo_root: Path,
    output_dir: Path,
) -> dict[str, Any]:
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite report output: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)
    runs_root = repo_root / "results/runs"
    completed_runs = (
        sorted(
            run_dir
            for run_dir in runs_root.iterdir()
            if run_dir.is_dir()
            and ((run_dir / "SUCCESS").is_file() or (run_dir / "FAILED").is_file())
        )
        if runs_root.is_dir()
        else []
    )
    run_inventory: list[dict[str, Any]] = []
    for run_dir in completed_runs:
        git = _read_json(run_dir / "git.json")
        metrics = _read_json(run_dir / "metrics.json")
        failures = [
            json.loads(line)
            for line in (run_dir / "failures.jsonl")
            .read_text(encoding="utf-8")
            .splitlines()
            if line.strip()
        ]
        run_inventory.append(
            {
                "run_id": run_dir.name,
                "status": "SUCCESS" if (run_dir / "SUCCESS").is_file() else "FAILED",
                "is_mock": bool(metrics.get("is_mock")),
                "git_commit": git.get("commit"),
                "failure_count": len(failures),
            }
        )
    _write_json(output_dir / "run_inventory.json", run_inventory)

    artifacts: dict[str, Any] = {}
    declared = {
        "machine_audit": repo_root / "artifacts/system/hardware.json",
        "structure_funnel": (
            repo_root / "data/manifests/experimental_structure_funnel.json"
        ),
        "esm_if1_smoke": repo_root / "artifacts/system/esm_if1_real_smoke.json",
        "proteinmpnn_smoke": (
            repo_root / "artifacts/system/proteinmpnn_real_smoke.json"
        ),
        "ligandmpnn_smoke": (repo_root / "artifacts/system/ligandmpnn_real_smoke.json"),
        "atlas_funnel": (repo_root / "data/processed/atlas/v1.0/data_funnel.json"),
        "clustering": (
            repo_root / "data/processed/atlas/v1.0/clusters/clustering_summary.json"
        ),
        "msa": repo_root / "data/processed/atlas/v1.0/msa/msa_manifest.json",
        "conservation": (
            repo_root
            / "data/processed/atlas/v1.0/conservation/conservation_manifest.json"
        ),
        "matched_multimodel_benchmark": (
            repo_root / "reports/tables/matched_multimodel_benchmark.json"
        ),
        "candidate_novelty": repo_root / "reports/tables/candidate_novelty.json",
        "real_refold": repo_root / "reports/tables/real_refold_summary.json",
    }
    for name, path in declared.items():
        artifacts[name] = (
            {
                "status": "available",
                "path": str(path.relative_to(repo_root)),
                "data": _read_json(path),
            }
            if path.is_file()
            else {"status": "not_run", "path": str(path.relative_to(repo_root))}
        )
    atlas_data = artifacts.get("atlas_funnel", {}).get("data", {})
    if (
        isinstance(atlas_data, dict)
        and atlas_data.get("high_confidence_pairs") == 0
        and int(atlas_data.get("ambiguous_pairs", 0)) > 0
    ):
        pairing_blocker = {
            "status": "blocked",
            "is_mock": False,
            "reason_code": "atlas_repeat_orientation_unavailable",
            "high_confidence_pairs": 0,
            "ambiguous_pairs": int(atlas_data["ambiguous_pairs"]),
            "reason": (
                "Atlas v1.0 does not provide a recoverable direct-repeat "
                "orientation for these records; ambiguous strands are excluded."
            ),
        }
        for artifact_name in ("paired_msa_real", "mi_apc_real", "formal_dca_real"):
            artifacts[artifact_name] = dict(pairing_blocker)
    else:
        for artifact_name in ("paired_msa_real", "mi_apc_real", "formal_dca_real"):
            artifacts[artifact_name] = {"status": "not_run"}
    benchmark_path = _latest(
        [
            run_dir / "benchmark/summary.json"
            for run_dir in completed_runs
            if (run_dir / "SUCCESS").is_file()
            and (run_dir / "benchmark/summary.json").is_file()
        ]
    )
    artifacts["experimental_benchmark"] = (
        {
            "status": "available",
            "path": str(benchmark_path.relative_to(repo_root)),
            "data": _read_json(benchmark_path),
        }
        if benchmark_path is not None
        else {"status": "not_run"}
    )
    _write_json(output_dir / "artifact_inventory.json", artifacts)

    available_real = sorted(
        name
        for name, value in artifacts.items()
        if value.get("status") == "available"
        and not bool(value.get("data", {}).get("is_mock", False))
    )
    not_run = sorted(
        name for name, value in artifacts.items() if value.get("status") == "not_run"
    )
    blocked = sorted(
        name for name, value in artifacts.items() if value.get("status") == "blocked"
    )
    summary = {
        "schema_version": "1.0",
        "is_mock": False,
        "maximum_evidence_level": 2
        if "experimental_benchmark" in available_real
        else 0,
        "available_real_artifacts": available_real,
        "not_run_artifacts": not_run,
        "blocked_artifacts": blocked,
        "completed_run_count": len(run_inventory),
        "successful_run_count": sum(
            row["status"] == "SUCCESS" for row in run_inventory
        ),
        "failed_run_count": sum(row["status"] == "FAILED" for row in run_inventory),
        "claim_scope": (
            "Computational evidence Levels 0-2 only. No result in this report "
            "is wet-lab validated or an effective Cas13."
        ),
    }
    _write_json(output_dir / "summary.json", summary)

    available_lines = "\n".join(f"- `{name}`" for name in available_real) or "- None"
    blocked_lines = "\n".join(f"- `{name}`" for name in blocked) or "- None"
    not_run_lines = "\n".join(f"- `{name}`" for name in not_run) or "- None"
    inventory_json = json.dumps(run_inventory, indent=2, sort_keys=True)
    markdown = f"""# ESM-IF-on-Cas13 project report

**Result class:** REAL AGGREGATION (individual mock runs remain labeled)

**Maximum evidence level:** {summary["maximum_evidence_level"]}

{summary["claim_scope"]}

## Available real artifacts

{available_lines}

## Explicitly not run

{not_run_lines}

## Blocked by declared source data

{blocked_lines}

## Run audit

- Completed: {summary["completed_run_count"]}
- Successful: {summary["successful_run_count"]}
- Failed: {summary["failed_run_count"]}

```json
{inventory_json}
```

The machine-readable scientific values and their `is_mock` fields are retained
in `artifact_inventory.json`; this page does not invent missing performance
numbers.
"""
    atomic_write_text(output_dir / "report.md", markdown)
    html_document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8">
<title>ESM-IF-on-Cas13 project report</title>
<style>body{{font-family:system-ui;max-width:1000px;margin:2rem auto;padding:0 1rem}}
pre{{background:#f4f4f4;padding:1rem;overflow:auto}}.warning{{color:#9b1c1c}}</style>
</head><body><h1>ESM-IF-on-Cas13 project report</h1>
<p><strong>Maximum evidence level:</strong> {summary["maximum_evidence_level"]}</p>
<p class="warning">{html.escape(summary["claim_scope"])}</p>
<h2>Available real artifacts</h2><pre>{html.escape(chr(10).join(available_real))}</pre>
<h2>Explicitly not run</h2><pre>{html.escape(chr(10).join(not_run))}</pre>
<h2>Blocked by declared source data</h2><pre>{html.escape(chr(10).join(blocked))}</pre>
<h2>Run audit</h2><pre>{html.escape(inventory_json)}</pre>
</body></html>
"""
    atomic_write_text(output_dir / "report.html", html_document)
    return summary


def _write_json(path: Path, value: Any) -> None:
    atomic_write_text(path, json.dumps(value, indent=2, sort_keys=True) + "\n")
