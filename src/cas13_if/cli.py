"""Command-line interface for auditable Cas13 inverse-folding workflows."""

from __future__ import annotations

import json
import shutil
import sys
from pathlib import Path
from typing import Annotated, Any

import typer
from rich.console import Console

from cas13_if import __version__
from cas13_if.backends.mock import MockBackend
from cas13_if.config import ConfigDict, ConfigError, load_config
from cas13_if.provenance import RunExistsError, RunRecorder
from cas13_if.schemas import SampleRequest

app = typer.Typer(
    name="cas13-if",
    help="Auditable structural, evolutionary, and RNA-aware Cas13 inverse folding.",
    no_args_is_help=True,
)
console = Console()
ConfigOption = Annotated[
    Path,
    typer.Option(
        "--config",
        "-c",
        exists=True,
        file_okay=True,
        dir_okay=False,
        readable=True,
        help="YAML configuration to resolve and record.",
    ),
]


def _repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def _experiment(config: ConfigDict, fallback: str) -> str:
    experiment = config.get("experiment", {})
    if isinstance(experiment, dict):
        value = experiment.get("name")
        if isinstance(value, str) and value:
            return value
    return fallback


def _start_run(
    command_name: str,
    config_path: Path,
    *,
    is_mock: bool,
) -> tuple[ConfigDict, RunRecorder]:
    config = load_config(config_path)
    config.setdefault("execution", {})
    execution = config["execution"]
    if not isinstance(execution, dict):
        raise ConfigError("execution configuration must be a mapping")
    execution["command"] = command_name
    execution["is_mock"] = is_mock
    recorder = RunRecorder(
        root=_repo_root() / "results/runs",
        experiment=_experiment(config, command_name),
        resolved_config=config,
        command=sys.argv,
        repo_root=_repo_root(),
        is_mock=is_mock,
    )
    console.print(f"run_dir={recorder.run_dir}")
    return config, recorder


def _fail_not_run(command_name: str, config: ConfigOption) -> None:
    try:
        _, recorder = _start_run(command_name, config, is_mock=False)
        message = (
            f"{command_name} prerequisites are not available or its production "
            "stage has not been completed; see failures.jsonl"
        )
        recorder.record_failure("precondition", message)
        recorder.finish(success=False)
        console.print(f"[red]not_run:[/red] {message}")
    except (ConfigError, RunExistsError) as exc:
        console.print(f"[red]error:[/red] {exc}")
    raise typer.Exit(code=2)


@app.callback()
def main(
    version: Annotated[
        bool,
        typer.Option("--version", help="Show package version and exit."),
    ] = False,
) -> None:
    if version:
        console.print(__version__)
        raise typer.Exit()


@app.command()
def preflight(
    config: ConfigOption,
    fixture: Annotated[
        bool, typer.Option(help="Run deterministic fixture backend I/O validation.")
    ] = False,
) -> None:
    """Validate configuration, repository layout, tools, and fixture backend I/O."""
    try:
        resolved, recorder = _start_run("preflight", config, is_mock=fixture)
        required_files = [
            "AGENTS.md",
            ".agent/PLANS.md",
            "pyproject.toml",
            "workflow/Snakefile",
            "docs/STATUS.md",
        ]
        missing = [
            relative
            for relative in required_files
            if not (_repo_root() / relative).is_file()
        ]
        tool_status = {
            name: shutil.which(name) is not None
            for name in ("git", "conda", "docker", "apptainer", "mmseqs")
        }
        metrics: dict[str, Any] = {
            "required_files": len(required_files),
            "missing_files": missing,
            "tools": tool_status,
            "config_keys": sorted(resolved),
        }
        if fixture:
            backend = MockBackend()
            backend.load()
            candidates = backend.sample(
                SampleRequest(
                    scaffold_id="fixture",
                    structure_path="tests/fixtures/minimal_complex.pdb",
                    parent_sequence="ACDEFGHIK",
                    count=2,
                    temperature=1.0,
                    seed=20260731,
                    fixed_positions={0: "A", 8: "K"},
                )
            )
            metrics["fixture_candidates"] = len(candidates)
            metrics["fixed_position_violations"] = sum(
                candidate.sequence[0] != "A" or candidate.sequence[8] != "K"
                for candidate in candidates
            )
            metrics["is_mock"] = True
        success = not missing and metrics.get("fixed_position_violations", 0) == 0
        if not success:
            recorder.record_failure("preflight", json.dumps(metrics, sort_keys=True))
        recorder.finish(success=success, metrics=metrics)
        console.print_json(json.dumps(metrics))
        if not success:
            raise typer.Exit(code=1)
    except (ConfigError, RunExistsError, ValueError) as exc:
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


@app.command()
def fetch(
    asset: Annotated[
        str,
        typer.Argument(help="references, third-party, models, atlas, or structures"),
    ],
) -> None:
    """Show the explicit offline-safe fetch script for an asset family."""
    scripts = {
        "references": "scripts/fetch_references.sh",
        "third-party": "scripts/fetch_third_party.sh",
        "models": "scripts/fetch_models.sh",
        "atlas": "scripts/fetch_atlas.sh",
        "structures": "scripts/fetch_experimental_structures.sh",
    }
    if asset not in scripts:
        console.print(f"[red]unknown asset:[/red] {asset}")
        raise typer.Exit(code=2)
    console.print(f"bash {scripts[asset]}")


@app.command("inspect-atlas")
def inspect_atlas(config: ConfigOption) -> None:
    """Inspect an Atlas source without materializing it in memory."""
    _fail_not_run("inspect-atlas", config)


@app.command("build-dataset")
def build_dataset(config: ConfigOption) -> None:
    """Build normalized Atlas and high-confidence Cas13 tables."""
    _fail_not_run("build-dataset", config)


@app.command()
def cluster(config: ConfigOption) -> None:
    """Cluster Cas13 sequences at registered MMseqs2 thresholds."""
    _fail_not_run("cluster", config)


@app.command("build-msa")
def build_msa(config: ConfigOption) -> None:
    """Build and validate subtype-specific alignments."""
    _fail_not_run("build-msa", config)


@app.command()
def conservation(config: ConfigOption) -> None:
    """Compute weighted subtype-specific conservation and entropy."""
    _fail_not_run("conservation", config)


@app.command("build-paired-msa")
def build_paired_msa(config: ConfigOption) -> None:
    """Build high-confidence paired Cas13/direct-repeat alignments."""
    _fail_not_run("build-paired-msa", config)


@app.command()
def coevolution(
    config: ConfigOption,
    fixture: Annotated[
        bool, typer.Option(help="Use deterministic fixture paired alignment.")
    ] = False,
) -> None:
    """Compute MI/APC or ingest a declared formal direct-coupling result."""
    if not fixture:
        _fail_not_run("coevolution", config)
    _fail_not_run("coevolution-fixture", config)


@app.command("inspect-structure")
def inspect_structure(config: ConfigOption) -> None:
    """Run strict protein/RNA structure QC and mapping."""
    _fail_not_run("inspect-structure", config)


@app.command("annotate-contacts")
def annotate_contacts(config: ConfigOption) -> None:
    """Annotate direct and second-shell RNA contacts."""
    _fail_not_run("annotate-contacts", config)


@app.command("build-mask")
def build_mask(config: ConfigOption) -> None:
    """Merge functional evidence into hard/soft/free positions."""
    _fail_not_run("build-mask", config)


@app.command()
def score(config: ConfigOption) -> None:
    """Score declared sequences with an offline local backend."""
    _fail_not_run("score", config)


@app.command()
def sample(config: ConfigOption) -> None:
    """Sample candidates with a declared offline local backend."""
    _fail_not_run("sample", config)


@app.command("sequence-qc")
def sequence_qc(config: ConfigOption) -> None:
    """Compute sequence validity, complexity, and composition metrics."""
    _fail_not_run("sequence-qc", config)


@app.command()
def novelty(config: ConfigOption) -> None:
    """Search candidates against the full declared Atlas sequence resource."""
    _fail_not_run("novelty", config)


@app.command()
def benchmark(config: ConfigOption) -> None:
    """Run the registered matched experimental-structure benchmark."""
    _fail_not_run("benchmark", config)


@app.command("export-refold")
def export_refold(config: ConfigOption) -> None:
    """Export deterministic structure-prediction manifests and shards."""
    _fail_not_run("export-refold", config)


@app.command("ingest-refold")
def ingest_refold(config: ConfigOption) -> None:
    """Ingest and audit declared structure-prediction outputs."""
    _fail_not_run("ingest-refold", config)


@app.command()
def report(config: ConfigOption) -> None:
    """Generate Markdown/HTML reports with evidence labels."""
    _fail_not_run("report", config)


@app.command()
def bundle(config: ConfigOption) -> None:
    """Create a manifest-based GPU/HPC transfer bundle."""
    _fail_not_run("bundle", config)


@app.command()
def verify(config: ConfigOption) -> None:
    """Verify hashes, manifests, output schemas, and leakage gates."""
    _fail_not_run("verify", config)


if __name__ == "__main__":
    app()
