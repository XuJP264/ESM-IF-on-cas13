from pathlib import Path

import pytest

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
