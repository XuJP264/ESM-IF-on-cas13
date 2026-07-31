from json import dumps, loads
from pathlib import Path

import pytest

from cas13_if.reporting.project import build_project_report
from cas13_if.reporting.reports import render_run_report
from cas13_if.statistics.resampling import (
    benjamini_hochberg,
    bootstrap_mean,
    paired_effect,
)


def test_resampling_and_adjustment_are_deterministic() -> None:
    first = bootstrap_mean([1, 2, 3], replicates=100, seed=5)
    second = bootstrap_mean([1, 2, 3], replicates=100, seed=5)
    assert first == second
    assert first.estimate == 2
    assert paired_effect([2, 4], [1, 1]) == 2
    adjusted = benjamini_hochberg([0.01, 0.04, 0.2])
    assert adjusted == pytest.approx([0.03, 0.06, 0.2])


def test_mock_report_is_visibly_labeled(tmp_path: Path) -> None:
    markdown = tmp_path / "report.md"
    html = tmp_path / "report.html"
    render_run_report(
        title="Fixture",
        evidence_level=0,
        is_mock=True,
        metrics={"fixed_violations": 0},
        failures=[],
        markdown_path=markdown,
        html_path=html,
    )
    assert "MOCK — TESTS ONLY" in markdown.read_text(encoding="utf-8")
    assert "supports no scientific claim" in html.read_text(encoding="utf-8")


def test_project_report_marks_missing_artifacts_not_run(tmp_path: Path) -> None:
    output = tmp_path / "output"
    summary = build_project_report(repo_root=tmp_path, output_dir=output)
    inventory = loads((output / "artifact_inventory.json").read_text(encoding="utf-8"))
    assert summary["maximum_evidence_level"] == 0
    assert inventory["experimental_benchmark"]["status"] == "not_run"
    assert "does not invent" in (output / "report.md").read_text(encoding="utf-8")


def test_project_report_ignores_benchmark_summary_from_failed_run(
    tmp_path: Path,
) -> None:
    runs = tmp_path / "results/runs"
    for run_id, status, value in (
        ("successful", "SUCCESS", 1),
        ("failed-later", "FAILED", 999),
    ):
        run = runs / run_id
        (run / "benchmark").mkdir(parents=True)
        (run / status).write_text("\n", encoding="utf-8")
        (run / "git.json").write_text(dumps({"commit": "abc"}), encoding="utf-8")
        (run / "metrics.json").write_text(dumps({"is_mock": False}), encoding="utf-8")
        (run / "failures.jsonl").write_text("", encoding="utf-8")
        (run / "benchmark/summary.json").write_text(
            dumps({"is_mock": False, "value": value}), encoding="utf-8"
        )
    output = tmp_path / "report"
    build_project_report(repo_root=tmp_path, output_dir=output)
    inventory = loads((output / "artifact_inventory.json").read_text(encoding="utf-8"))
    benchmark = inventory["experimental_benchmark"]
    assert benchmark["data"]["value"] == 1
    assert "successful" in benchmark["path"]
