"""Minimal evidence-aware Markdown and HTML report rendering."""

from __future__ import annotations

import html
import json
from pathlib import Path
from typing import Any


def render_run_report(
    *,
    title: str,
    evidence_level: int,
    is_mock: bool,
    metrics: dict[str, Any],
    failures: list[dict[str, Any]],
    markdown_path: Path,
    html_path: Path,
) -> None:
    if not 0 <= evidence_level <= 4:
        raise ValueError("evidence level must be between 0 and 4")
    label = "MOCK — TESTS ONLY" if is_mock else "REAL"
    caveat = (
        "This report contains mock outputs and supports no scientific claim."
        if is_mock
        else (
            "This computational report does not establish Level 4 or an "
            "effective/validated Cas13."
        )
    )
    metrics_json = json.dumps(metrics, indent=2, sort_keys=True)
    failures_json = json.dumps(failures, indent=2, sort_keys=True)
    markdown = (
        f"# {title}\n\n"
        f"**Result class:** {label}\n\n"
        f"**Maximum evidence level:** {evidence_level}\n\n"
        f"{caveat}\n\n"
        "## Metrics\n\n"
        f"```json\n{metrics_json}\n```\n\n"
        "## Failures\n\n"
        f"```json\n{failures_json}\n```\n"
    )
    markdown_path.parent.mkdir(parents=True, exist_ok=True)
    markdown_path.write_text(markdown, encoding="utf-8")
    html_document = f"""<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>{html.escape(title)}</title>
<style>body{{font-family:system-ui;max-width:960px;margin:2rem auto;padding:0 1rem}}
.mock{{color:#a00;font-weight:700}}pre{{background:#f5f5f5;padding:1rem;overflow:auto}}</style>
</head><body>
<h1>{html.escape(title)}</h1>
<p class="{"mock" if is_mock else "real"}">Result class: {html.escape(label)}</p>
<p>Maximum evidence level: {evidence_level}</p>
<p>{html.escape(caveat)}</p>
<h2>Metrics</h2><pre>{html.escape(metrics_json)}</pre>
<h2>Failures</h2><pre>{html.escape(failures_json)}</pre>
</body></html>
"""
    html_path.parent.mkdir(parents=True, exist_ok=True)
    html_path.write_text(html_document, encoding="utf-8")
