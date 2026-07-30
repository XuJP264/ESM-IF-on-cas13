import json
from pathlib import Path

from typer.testing import CliRunner

from cas13_if import cli


def _minimal_repo(root: Path) -> Path:
    for relative in (
        "AGENTS.md",
        ".agent/PLANS.md",
        "pyproject.toml",
        "workflow/Snakefile",
        "docs/STATUS.md",
    ):
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("fixture\n", encoding="utf-8")
    hardware = root / "artifacts/system/hardware.json"
    hardware.parent.mkdir(parents=True)
    hardware.write_text(
        json.dumps({"is_mock": True, "fixture": True}) + "\n", encoding="utf-8"
    )
    config = root / "config.yaml"
    config.write_text(
        "experiment:\n  name: cli-fixture\n  seed: 20260731\n",
        encoding="utf-8",
    )
    return config


def test_cli_help_fetch_and_fixture_preflight(
    tmp_path: Path, monkeypatch: object
) -> None:
    monkeypatch.setattr(cli, "_repo_root", lambda: tmp_path)  # type: ignore[attr-defined]
    runner = CliRunner()
    help_result = runner.invoke(cli.app, ["--help"])
    assert help_result.exit_code == 0
    assert "preflight" in help_result.output
    fetch_result = runner.invoke(cli.app, ["fetch", "models"])
    assert fetch_result.exit_code == 0
    assert "fetch_models.sh" in fetch_result.output
    invalid_fetch = runner.invoke(cli.app, ["fetch", "unknown"])
    assert invalid_fetch.exit_code == 2
    config = _minimal_repo(tmp_path)
    preflight = runner.invoke(
        cli.app, ["preflight", "--config", str(config), "--fixture"]
    )
    assert preflight.exit_code == 0, preflight.output
    assert '"fixed_position_violations": 0' in preflight.output
    run_dirs = list((tmp_path / "results/runs").iterdir())
    assert len(run_dirs) == 1
    assert (run_dirs[0] / "SUCCESS").is_file()
