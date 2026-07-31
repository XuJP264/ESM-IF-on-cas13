"""Command-line interface for auditable Cas13 inverse-folding workflows."""

from __future__ import annotations

import json
import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Annotated, Any

import numpy as np
import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import typer
from rich.console import Console

from cas13_if import __version__
from cas13_if.alignments.msa import read_aligned_fasta
from cas13_if.alignments.pipeline import build_subtype_msas
from cas13_if.alignments.scaffold_mapping import build_scaffold_mapping
from cas13_if.backends.mock import MockBackend
from cas13_if.config import ConfigDict, ConfigError, load_config
from cas13_if.data.atlas import iter_json_array, process_atlas
from cas13_if.data.clustering import (
    MmseqsParameters,
    assert_no_cluster_leakage,
    assign_cluster_splits,
    run_mmseqs_clustering,
)
from cas13_if.data.fasta import write_fasta
from cas13_if.evolution.coevolution import (
    bootstrap_top_pair_frequency,
    compute_mi_apc,
    permuted_cross_block_maxima,
)
from cas13_if.evolution.pipeline import compute_subtype_conservation
from cas13_if.novelty.pipeline import NoveltyThresholds, run_candidate_novelty_pipeline
from cas13_if.provenance import (
    RunExistsError,
    RunRecorder,
    atomic_write_text,
    sha256_file,
)
from cas13_if.reporting.project import build_project_report
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


def _mapping(config: ConfigDict, key: str) -> dict[str, Any]:
    value = config.get(key)
    if not isinstance(value, dict):
        raise ConfigError(f"{key} configuration must be a mapping")
    return value


def _path(value: Any, *, key: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ConfigError(f"{key} must be a non-empty path")
    path = Path(value)
    return path if path.is_absolute() else _repo_root() / path


def _file_entry(path: Path) -> dict[str, Any]:
    return {
        "path": str(path.relative_to(_repo_root())),
        "size_bytes": path.stat().st_size,
        "sha256": sha256_file(path),
    }


def _record_inputs(recorder: RunRecorder, files: list[Path]) -> None:
    atomic_write_text(
        recorder.run_dir / "input_manifest.json",
        json.dumps(
            {
                "files": [_file_entry(path) for path in files],
                "is_mock": recorder.is_mock,
            },
            indent=2,
            sort_keys=True,
        )
        + "\n",
    )


def _tree_files(path: Path) -> list[Path]:
    return sorted(item for item in path.rglob("*") if item.is_file())


def _value_counts(values: Any) -> dict[str, int]:
    counts: dict[str, int] = {}
    for value in values:
        key = str(value)
        counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


def _subtype_counts(rows: list[dict[str, Any]]) -> dict[str, int]:
    counts: dict[str, int] = {}
    for row in rows:
        subtypes = row.get("subtypes")
        if not isinstance(subtypes, list):
            continue
        for subtype in subtypes:
            key = str(subtype)
            counts[key] = counts.get(key, 0) + 1
    return dict(sorted(counts.items()))


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
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("inspect-atlas", config, is_mock=False)
        atlas = _mapping(resolved, "atlas")
        source = _path(atlas.get("input_file"), key="atlas.input_file")
        if not source.is_file():
            raise FileNotFoundError(f"Atlas input is missing: {source}")
        _record_inputs(recorder, [source])
        operons = 0
        type_vi = 0
        cas13_annotations = 0
        for raw in iter_json_array(source):
            operons += 1
            if not isinstance(raw, dict):
                continue
            subtype = str((raw.get("summary") or {}).get("subtype") or "")
            type_vi += int(subtype.upper().startswith("VI-"))
            cas13_annotations += sum(
                "cas13"
                in " ".join(
                    str(entry.get(key, ""))
                    for key in ("gene_name", "hmm_name", "annotation")
                ).lower()
                or "c2c2"
                in " ".join(
                    str(entry.get(key, ""))
                    for key in ("gene_name", "hmm_name", "annotation")
                ).lower()
                for entry in (raw.get("cas") or [])
                if isinstance(entry, dict)
            )
        metrics = {
            "is_mock": False,
            "evidence_level": 0,
            "operons": operons,
            "type_vi_operons": type_vi,
            "cas13_annotations": cas13_annotations,
        }
        recorder.finish(success=True, metrics=metrics)
        console.print_json(json.dumps(metrics))
    except (ConfigError, RunExistsError, FileNotFoundError, ValueError) as exc:
        if recorder is not None:
            recorder.record_failure("inspect-atlas", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


@app.command("build-dataset")
def build_dataset(config: ConfigOption) -> None:
    """Build normalized Atlas and high-confidence Cas13 tables."""
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("build-dataset", config, is_mock=False)
        atlas = _mapping(resolved, "atlas")
        source = _path(atlas.get("input_file"), key="atlas.input_file")
        output_dir = _path(atlas.get("output_dir"), key="atlas.output_dir")
        batch_size = int(atlas.get("batch_size", 10_000))
        _record_inputs(recorder, [source])
        funnel = process_atlas(source, output_dir, batch_size=batch_size)
        outputs = [
            _file_entry(path)
            for path in sorted(output_dir.iterdir())
            if path.is_file() and path.name != "cas13_exact_dedup.sqlite"
        ]
        recorder.finish(success=True, metrics=funnel, outputs=outputs)
        console.print_json(json.dumps(funnel))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("build-dataset", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


@app.command()
def cluster(config: ConfigOption) -> None:
    """Cluster Cas13 sequences at registered MMseqs2 thresholds."""
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("cluster", config, is_mock=False)
        clustering = _mapping(resolved, "clustering")
        source = _path(clustering.get("input_file"), key="clustering.input_file")
        output_root = _path(clustering.get("output_dir"), key="clustering.output_dir")
        executable = _path(clustering.get("executable"), key="clustering.executable")
        if output_root.exists():
            raise FileExistsError(
                f"refusing to overwrite clustering output: {output_root}"
            )
        if not executable.is_file():
            raise FileNotFoundError(f"MMseqs2 executable is missing: {executable}")
        _record_inputs(recorder, [source])
        table = pq.read_table(source)
        rows = table.select(
            ["sequence_sha256", "protein_sequence", "protein_length", "subtypes"]
        ).to_pylist()
        if not rows:
            raise ValueError("exact-unique Cas13 table is empty")
        output_root.mkdir(parents=True, exist_ok=False)
        fasta = output_root / "cas13_exact_unique.fasta"
        write_fasta(
            (
                (str(row["sequence_sha256"]), str(row["protein_sequence"]))
                for row in rows
            ),
            fasta,
        )
        identities = clustering.get("identities")
        if not isinstance(identities, list) or not identities:
            raise ConfigError("clustering.identities must be a non-empty list")
        seed = int(_mapping(resolved, "experiment").get("seed", 20260731))
        strict_identity = float(
            _mapping(resolved, "splits").get("strict_identity", 0.4)
        )
        summaries: dict[str, Any] = {}
        output_files: list[Path] = [fasta]
        for identity_value in identities:
            identity = float(identity_value)
            label = f"identity_{round(identity * 100):03d}"
            parameters = MmseqsParameters(
                minimum_identity=identity,
                coverage=float(clustering.get("coverage", 0.8)),
                coverage_mode=int(clustering.get("coverage_mode", 0)),
                cluster_mode=int(clustering.get("cluster_mode", 2)),
                threads=int(clustering.get("threads", 16)),
            )
            result = run_mmseqs_clustering(
                fasta,
                output_root / label,
                parameters,
                executable=str(executable),
            )
            mapping = result["mapping"]
            expected = {str(row["sequence_sha256"]) for row in rows}
            missing = expected.difference(mapping)
            unexpected = set(mapping).difference(expected)
            if missing or unexpected:
                raise RuntimeError(
                    "MMseqs2 mapping mismatch: "
                    f"missing={len(missing)}, unexpected={len(unexpected)}"
                )
            split_assignments = (
                assign_cluster_splits(mapping, seed=seed)
                if identity == strict_identity
                else {}
            )
            if split_assignments:
                assert_no_cluster_leakage(mapping, split_assignments)
            mapping_rows = [
                {
                    "sequence_sha256": str(row["sequence_sha256"]),
                    "representative_sha256": mapping[str(row["sequence_sha256"])],
                    "split": split_assignments.get(str(row["sequence_sha256"])),
                    "protein_length": int(row["protein_length"]),
                    "subtypes": row["subtypes"],
                }
                for row in rows
            ]
            mapping_path = output_root / label / "cluster_mapping.parquet"
            pq.write_table(pa.Table.from_pylist(mapping_rows), mapping_path)
            output_files.extend(
                [
                    output_root / label / "metadata.json",
                    output_root / label / "cluster_summary.json",
                    mapping_path,
                ]
            )
            summaries[label] = {
                **result["summary"],
                "subtype_counts": _subtype_counts(mapping_rows),
                "split_counts": _value_counts(split_assignments.values()),
                "leakage_gate": "passed" if split_assignments else "not_applicable",
            }
        summary_path = output_root / "clustering_summary.json"
        atomic_write_text(
            summary_path, json.dumps(summaries, indent=2, sort_keys=True) + "\n"
        )
        output_files.append(summary_path)
        recorder.finish(
            success=True,
            metrics={
                "is_mock": False,
                "evidence_level": 0,
                "thresholds": summaries,
            },
            outputs=[_file_entry(path) for path in output_files],
        )
        console.print_json(json.dumps(summaries))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("cluster", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


@app.command("build-msa")
def build_msa(config: ConfigOption) -> None:
    """Build and validate subtype-specific alignments."""
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("build-msa", config, is_mock=False)
        msa = _mapping(resolved, "msa")
        source = _path(msa.get("input_file"), key="msa.input_file")
        cluster_mapping = _path(msa.get("cluster_mapping"), key="msa.cluster_mapping")
        output_dir = _path(msa.get("output_dir"), key="msa.output_dir")
        executable = _path(msa.get("executable"), key="msa.executable")
        _record_inputs(recorder, [source, cluster_mapping])
        manifest = build_subtype_msas(
            exact_unique_path=source,
            cluster_mapping_path=cluster_mapping,
            output_dir=output_dir,
            executable=str(executable),
            threads=int(msa.get("threads", 16)),
            minimum_protein_length=int(msa.get("minimum_protein_length", 1)),
            maximum_protein_length=(
                int(msa["maximum_protein_length"])
                if msa.get("maximum_protein_length") is not None
                else None
            ),
        )
        recorder.finish(
            success=True,
            metrics=manifest,
            outputs=[_file_entry(path) for path in _tree_files(output_dir)],
        )
        console.print_json(json.dumps(manifest))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("build-msa", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


@app.command()
def conservation(config: ConfigOption) -> None:
    """Compute weighted subtype-specific conservation and entropy."""
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("conservation", config, is_mock=False)
        settings = _mapping(resolved, "conservation")
        msa_dir = _path(settings.get("msa_dir"), key="conservation.msa_dir")
        output_dir = _path(settings.get("output_dir"), key="conservation.output_dir")
        _record_inputs(recorder, [msa_dir / "msa_manifest.json"])
        manifest = compute_subtype_conservation(
            msa_root=msa_dir,
            output_dir=output_dir,
            identity_threshold=float(settings.get("sequence_identity_threshold", 0.8)),
            allowed_frequency=float(settings.get("allowed_residue_frequency", 0.05)),
            constraint_minimum_column_coverage=float(
                settings.get("constraint_minimum_column_coverage", 0.8)
            ),
        )
        recorder.finish(
            success=True,
            metrics=manifest,
            outputs=[_file_entry(path) for path in _tree_files(output_dir)],
        )
        console.print_json(json.dumps(manifest))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("conservation", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


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
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("coevolution-fixture", config, is_mock=True)
        fixture_path = _repo_root() / "data/fixtures/paired_msa.fasta"
        _record_inputs(recorder, [fixture_path])
        alignment = read_aligned_fasta(fixture_path, alphabet="mixed")
        split_column = 3
        seed = int(_mapping(resolved, "experiment").get("seed", 20260731))
        result = compute_mi_apc(alignment)
        bootstrap = bootstrap_top_pair_frequency(
            alignment,
            replicates=20,
            seed=seed,
            top_n=2,
        )
        permutation = permuted_cross_block_maxima(
            alignment,
            split_column=split_column,
            replicates=20,
            seed=seed,
        )
        cross_block = result.apc_corrected[:split_column, split_column:]
        ranked = np.argsort(cross_block, axis=None)[::-1]
        top_pairs = [
            {
                "protein_column_0": int(np.unravel_index(index, cross_block.shape)[0]),
                "rna_column_0": int(
                    split_column + np.unravel_index(index, cross_block.shape)[1]
                ),
                "apc_mi": float(cross_block.flat[index]),
            }
            for index in ranked[:3]
        ]
        output_dir = recorder.run_dir / "coevolution"
        output_dir.mkdir()
        matrices_path = output_dir / "mi_apc_bootstrap.npz"
        np.savez_compressed(
            matrices_path,
            mutual_information=result.mutual_information,
            apc_corrected=result.apc_corrected,
            bootstrap_frequency=bootstrap,
            permutation_cross_block_maxima=permutation,
        )
        summary = {
            "schema_version": "1.0",
            "is_mock": True,
            "evidence_level": 0,
            "alignment_sequences": alignment.n_sequences,
            "alignment_columns": alignment.n_columns,
            "split_column": split_column,
            "effective_sequence_count": result.effective_sequence_count,
            "bootstrap_replicates": 20,
            "permutation_replicates": 20,
            "top_cross_block_pairs": top_pairs,
            "formal_dca_status": "not_run",
            "claim_scope": (
                "Deterministic fixture smoke only; MI/APC is not DCA and "
                "supports no scientific Cas13/direct-repeat claim."
            ),
        }
        summary_path = output_dir / "summary.json"
        atomic_write_text(
            summary_path,
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
        )
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[_file_entry(path) for path in _tree_files(output_dir)],
        )
        console.print_json(json.dumps(summary))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("coevolution-fixture", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


@app.command("inspect-structure")
def inspect_structure(config: ConfigOption) -> None:
    """Run strict protein/RNA structure QC and mapping."""
    _fail_not_run("inspect-structure", config)


@app.command("map-scaffold")
def map_scaffold(config: ConfigOption) -> None:
    """Map a resolved structure chain through full scaffold and subtype MSA."""
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("map-scaffold", config, is_mock=False)
        settings = _mapping(resolved, "mapping")
        structure = _path(settings.get("structure"), key="mapping.structure")
        entity = _path(settings.get("protein_entity"), key="mapping.protein_entity")
        msa = _path(settings.get("msa"), key="mapping.msa")
        conservation_path = _path(
            settings.get("conservation"), key="mapping.conservation"
        )
        output_dir = _path(settings.get("output_dir"), key="mapping.output_dir")
        executable = _path(
            settings.get("mafft_executable"), key="mapping.mafft_executable"
        )
        inputs = [structure, entity, msa, conservation_path]
        _record_inputs(recorder, inputs)
        summary = build_scaffold_mapping(
            structure_path=structure,
            entity_path=entity,
            msa_path=msa,
            conservation_path=conservation_path,
            output_dir=output_dir,
            chain_id=str(settings.get("chain_id", "A")),
            subtype=str(settings.get("subtype", "VI-D")),
            mafft_executable=str(executable),
            threads=int(settings.get("threads", 1)),
            minimum_conservation_coverage=float(
                settings.get("minimum_conservation_coverage", 0.8)
            ),
        )
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[_file_entry(path) for path in _tree_files(output_dir)],
        )
        console.print_json(json.dumps(summary))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("map-scaffold", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


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
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("novelty", config, is_mock=False)
        inputs = _mapping(resolved, "inputs")
        search = _mapping(resolved, "search")
        filters = _mapping(resolved, "filters")
        report = _mapping(resolved, "report")
        candidate_jsonl = _path(
            inputs.get("candidate_jsonl"), key="inputs.candidate_jsonl"
        )
        atlas_fasta = _path(inputs.get("atlas_fasta"), key="inputs.atlas_fasta")
        executable = _path(search.get("executable"), key="search.executable")
        canonical_summary = _path(
            report.get("canonical_summary"), key="report.canonical_summary"
        )
        if canonical_summary.exists():
            raise FileExistsError(
                f"refusing to overwrite canonical novelty summary: {canonical_summary}"
            )
        _record_inputs(recorder, [candidate_jsonl, atlas_fasta, executable])
        output_dir = recorder.run_dir / "novelty"
        summary = run_candidate_novelty_pipeline(
            candidate_jsonl=candidate_jsonl,
            atlas_fasta=atlas_fasta,
            output_dir=output_dir,
            executable=executable,
            threads=int(search.get("threads", 16)),
            sensitivity=float(search.get("sensitivity", 7.5)),
            minimum_query_coverage=float(search.get("minimum_query_coverage", 0.8)),
            maximum_evalue=float(search.get("maximum_evalue", 1000.0)),
            maximum_sequences=int(search.get("maximum_sequences", 5000)),
            thresholds=NoveltyThresholds(
                maximum_parent_identity=float(filters["max_parent_identity"]),
                maximum_atlas_identity=float(filters["max_atlas_identity"]),
                maximum_homopolymer_length=int(filters["max_homopolymer_length"]),
                maximum_low_complexity_windows=int(
                    filters["max_low_complexity_windows"]
                ),
                minimum_designed_position_entropy=float(
                    filters["min_designed_position_entropy"]
                ),
                low_complexity_window=int(
                    _mapping(filters, "low_complexity")["window"]
                ),
                low_complexity_maximum_fraction=float(
                    _mapping(filters, "low_complexity")["max_single_residue_fraction"]
                ),
            ),
        )
        atomic_write_text(
            canonical_summary,
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
        )
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[_file_entry(path) for path in _tree_files(output_dir)]
            + [_file_entry(canonical_summary)],
        )
        console.print_json(json.dumps(summary))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
        KeyError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("novelty", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


@app.command()
def benchmark(config: ConfigOption) -> None:
    """Run the registered matched experimental-structure benchmark."""
    recorder: RunRecorder | None = None
    try:
        resolved, recorder = _start_run("benchmark", config, is_mock=False)
        model = _mapping(resolved, "model")
        structure = _mapping(resolved, "structures")
        environment = _path(model.get("environment"), key="model.environment")
        python = environment / "bin/python"
        script = _repo_root() / "scripts/benchmark_esm_if1.py"
        checkpoint = _path(model.get("checkpoint"), key="model.checkpoint")
        structure_manifest = _path(structure.get("manifest"), key="structures.manifest")
        functional_manifest = _path(
            structure.get("functional_manifest"),
            key="structures.functional_manifest",
        )
        inputs = [config, checkpoint, structure_manifest, functional_manifest]
        for pdb_id in structure.get("pdb_ids", []):
            inputs.append(
                _repo_root()
                / "data/experimental_structures"
                / f"{str(pdb_id).lower()}.cif"
            )
        for path in [python, script, *inputs]:
            if not path.is_file():
                raise FileNotFoundError(f"benchmark prerequisite is missing: {path}")
        _record_inputs(recorder, inputs)
        output_dir = recorder.run_dir / "benchmark"
        command = [
            str(python),
            str(script),
            "--config",
            str(config),
            "--output-dir",
            str(output_dir),
        ]
        completed = subprocess.run(
            command,
            cwd=_repo_root(),
            text=True,
            capture_output=True,
            check=False,
            env={**os.environ, "PYTHONNOUSERSITE": "1"},
        )
        atomic_write_text(recorder.run_dir / "stdout.log", completed.stdout)
        atomic_write_text(recorder.run_dir / "stderr.log", completed.stderr)
        if completed.returncode != 0:
            raise RuntimeError(
                "experimental benchmark failed; see stdout.log and stderr.log"
            )
        summary_path = output_dir / "summary.json"
        summary = json.loads(summary_path.read_text(encoding="utf-8"))
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[_file_entry(path) for path in _tree_files(output_dir)],
        )
        console.print_json(json.dumps(summary))
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("benchmark", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


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
    recorder: RunRecorder | None = None
    try:
        _, recorder = _start_run("report", config, is_mock=False)
        output_dir = recorder.run_dir / "report"
        summary = build_project_report(
            repo_root=_repo_root(),
            output_dir=output_dir,
        )
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[_file_entry(path) for path in _tree_files(output_dir)],
        )
        pointer = _repo_root() / "reports/latest.txt"
        atomic_write_text(
            pointer,
            str(output_dir.relative_to(_repo_root())) + "\n",
        )
        console.print_json(json.dumps(summary))
        console.print(f"report={output_dir / 'report.html'}")
    except (
        ConfigError,
        RunExistsError,
        FileNotFoundError,
        FileExistsError,
        RuntimeError,
        ValueError,
        OSError,
    ) as exc:
        if recorder is not None:
            recorder.record_failure("report", str(exc))
            recorder.finish(success=False)
        console.print(f"[red]error:[/red] {exc}")
        raise typer.Exit(code=2) from exc


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
