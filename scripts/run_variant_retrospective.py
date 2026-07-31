#!/usr/bin/env python
"""Run a real, descriptive CasRx point-variant retrospective benchmark."""

from __future__ import annotations

import json
import os
import re
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]
from scipy.stats import spearmanr

from cas13_if.config import ConfigDict, load_config
from cas13_if.data.fasta import iter_fasta, write_fasta
from cas13_if.provenance import RunRecorder, atomic_write_text, sha256_file
from cas13_if.scoring.multistate import (
    aggregate_multistate_scores,
    state_rank_consistency,
)

POINT_MUTATION = re.compile(r"([A-Z])(\d+)([A-Z])")


def _repo_path(repo: Path, value: Any) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else repo / path


def _jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    atomic_write_text(
        path, "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows)
    )


def _point_positions(mutation: str) -> list[int] | None:
    parts = mutation.split(",")
    matches = [POINT_MUTATION.fullmatch(part.strip()) for part in parts]
    if any(item is None for item in matches):
        return None
    return [int(item.group(2)) for item in matches if item is not None]


def _project(mapping: pd.DataFrame, sequence: str) -> str:
    coordinate = mapping.loc[mapping["coordinate_index_0"].notna()].copy()
    coordinate["coordinate_index_0"] = coordinate["coordinate_index_0"].astype(int)
    coordinate = coordinate.sort_values("coordinate_index_0")
    if coordinate["coordinate_index_0"].tolist() != list(range(len(coordinate))):
        raise ValueError("coordinate projection is not contiguous")
    return "".join(
        sequence[int(index)] for index in coordinate["full_scaffold_index_0"]
    )


def _isolated_environment() -> dict[str, str]:
    environment = os.environ.copy()
    environment["PYTHONNOUSERSITE"] = "1"
    return environment


def _run(command: list[str], *, cwd: Path, stdout: Path, stderr: Path) -> None:
    completed = subprocess.run(
        command,
        cwd=cwd,
        env=_isolated_environment(),
        text=True,
        capture_output=True,
        check=False,
    )
    atomic_write_text(stdout, completed.stdout)
    atomic_write_text(stderr, completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(
            f"subprocess failed ({completed.returncode}): {' '.join(command)}; "
            f"stderr={completed.stderr[-2000:]}"
        )


def _proteinmpnn_scores(
    *,
    repo: Path,
    config: ConfigDict,
    records: list[dict[str, Any]],
    work: Path,
) -> dict[str, dict[str, float]]:
    model = _repo_path(repo, config["models"]["proteinmpnn_checkpoint"])
    upstream = _repo_path(repo, config["models"]["proteinmpnn_upstream"])
    python = _repo_path(repo, config["models"]["mpnn_python"])
    structure = _repo_path(repo, config["inputs"]["primary_structure"])
    output: dict[str, dict[str, float]] = {}
    for record in records:
        identifier = str(record["score_id"])
        root = work / f"proteinmpnn-{identifier}"
        root.mkdir()
        fasta = root / "target.fasta"
        write_fasta([(identifier, str(record["full_sequence"]))], fasta)
        out = root / "output"
        command = [
            str(python),
            str(upstream / "protein_mpnn_run.py"),
            "--pdb_path",
            str(structure),
            "--pdb_path_chains",
            str(config["inputs"]["protein_chain"]),
            "--out_folder",
            str(out),
            "--path_to_model_weights",
            str(model.parent),
            "--model_name",
            model.stem,
            "--score_only",
            "1",
            "--path_to_fasta",
            str(fasta),
            "--num_seq_per_target",
            str(config["execution"]["proteinmpnn_decoding_orders"]),
            "--batch_size",
            "1",
            "--seed",
            str(config["execution"]["seed"]),
            "--suppress_print",
            "1",
        ]
        _run(
            command,
            cwd=upstream,
            stdout=root / "stdout.log",
            stderr=root / "stderr.log",
        )
        score_file = next((out / "score_only").glob("*_fasta_1.npz"))
        with np.load(score_file) as values:
            scores = np.asarray(values["score"], dtype=float)
        output[identifier] = {
            "proteinmpnn_mean_nll": float(scores.mean()),
            "proteinmpnn_nll_std": float(scores.std()),
            "proteinmpnn_decoding_orders": int(scores.size),
        }
    return output


def _stats_score(
    python: Path, stats_path: Path, expected_sequence: str
) -> dict[str, Any]:
    script = (
        "import json,torch;d=torch.load(r'"
        + str(stats_path)
        + "',map_location='cpu',weights_only=True);"
        "lp=d['log_probs'][0];s=d['generated_sequences'][0];"
        "m=d['mask']*d['chain_mask'];sel=lp[range(len(s)),s];"
        "print(json.dumps({'sequence':''.join('ACDEFGHIKLMNPQRSTVWYX'[int(x)] "
        "for x in s),'mean_log_prob':float((sel*m).sum()/m.sum()),"
        "'positions':int(m.sum())}))"
    )
    completed = subprocess.run(
        [str(python), "-c", script],
        env=_isolated_environment(),
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        raise RuntimeError(f"LigandMPNN stats read failed: {completed.stderr}")
    result = json.loads(completed.stdout)
    if result["sequence"] != expected_sequence:
        raise RuntimeError("forced LigandMPNN sequence differs from requested variant")
    return result


def _ligandmpnn_scores(
    *,
    repo: Path,
    config: ConfigDict,
    records: list[dict[str, Any]],
    mapping: pd.DataFrame,
    work: Path,
) -> dict[str, dict[str, Any]]:
    upstream = _repo_path(repo, config["models"]["ligandmpnn_upstream"])
    python = _repo_path(repo, config["models"]["mpnn_python"])
    structure = _repo_path(repo, config["inputs"]["primary_structure"])
    coordinate = mapping.loc[mapping["coordinate_index_0"].notna()].copy()
    coordinate["coordinate_index_0"] = coordinate["coordinate_index_0"].astype(int)
    coordinate = coordinate.sort_values("coordinate_index_0")
    output: dict[str, dict[str, Any]] = {}
    bias_value = float(config["execution"]["ligandmpnn_forcing_bias"])
    for record in records:
        identifier = str(record["score_id"])
        root = work / f"ligandmpnn-{identifier}"
        root.mkdir()
        projected = str(record["projected_9M31"])
        bias = {
            f"A{int(row.pdb_residue_number)}"
            + (
                str(row.pdb_insertion_code) if pd.notna(row.pdb_insertion_code) else ""
            ): {projected[int(row.coordinate_index_0)]: bias_value}
            for row in coordinate.itertuples()
        }
        bias_path = root / "forcing_bias.json"
        atomic_write_text(bias_path, json.dumps(bias, sort_keys=True) + "\n")
        out = root / "output"
        command = [
            str(python),
            str(upstream / "run.py"),
            "--model_type",
            "ligand_mpnn",
            "--checkpoint_ligand_mpnn",
            str(_repo_path(repo, config["models"]["ligandmpnn_checkpoint"])),
            "--checkpoint_protein_mpnn",
            str(_repo_path(repo, config["models"]["ligandmpnn_protein_checkpoint"])),
            "--checkpoint_soluble_mpnn",
            str(_repo_path(repo, config["models"]["ligandmpnn_soluble_checkpoint"])),
            "--pdb_path",
            str(structure),
            "--out_folder",
            str(out),
            "--parse_these_chains_only",
            ",".join(
                [
                    str(config["inputs"]["protein_chain"]),
                    *[str(item) for item in config["inputs"]["rna_context_chains"]],
                ]
            ),
            "--chains_to_design",
            str(config["inputs"]["protein_chain"]),
            "--bias_AA_per_residue",
            str(bias_path),
            "--batch_size",
            "1",
            "--number_of_batches",
            "1",
            "--temperature",
            "1.0",
            "--seed",
            str(config["execution"]["seed"]),
            "--ligand_mpnn_use_atom_context",
            "1",
            "--save_stats",
            "1",
            "--verbose",
            "0",
        ]
        _run(
            command,
            cwd=upstream,
            stdout=root / "stdout.log",
            stderr=root / "stderr.log",
        )
        fasta = next((out / "seqs").glob("*.fa"))
        fasta_records = list(iter_fasta(fasta))
        if not fasta_records or fasta_records[-1][1] != projected:
            raise RuntimeError("LigandMPNN output did not preserve forced target")
        headers = [
            line[1:]
            for line in fasta.read_text(encoding="utf-8").splitlines()
            if line.startswith(">")
        ]
        # The upstream FASTA records context configuration on the input record
        # and sampling statistics on generated records.  Require its explicit
        # context attestation without assuming which record is last.
        if not any("use_ligand_context=True" in header for header in headers):
            raise RuntimeError("LigandMPNN did not report atomic context")
        stats_path = next((out / "stats").glob("*.pt"))
        stats = _stats_score(python, stats_path, projected)
        output[identifier] = {
            "ligandmpnn_mean_log_probability": float(stats["mean_log_prob"]),
            "ligandmpnn_scored_positions": int(stats["positions"]),
            "ligandmpnn_rna_atomic_context": True,
            "ligandmpnn_forcing_semantics": (
                "target-only sampling bias; reported score uses raw pre-bias log_probs"
            ),
        }
    return output


def _bootstrap_mean(
    values: np.ndarray, *, replicates: int, seed: int
) -> tuple[float, float]:
    if values.size == 0:
        return float("nan"), float("nan")
    generator = np.random.default_rng(seed)
    estimates = np.asarray(
        [
            generator.choice(values, size=values.size, replace=True).mean()
            for _ in range(replicates)
        ]
    )
    return float(np.quantile(estimates, 0.025)), float(np.quantile(estimates, 0.975))


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    config_path = repo / "configs/stage_0003a_retrospective.yaml"
    config = load_config(config_path)
    recorder = RunRecorder(
        root=repo / "results/runs",
        experiment=str(config["experiment"]["name"]),
        resolved_config=config,
        command=[sys.executable, str(Path(__file__).resolve())],
        repo_root=repo,
        is_mock=False,
    )
    work = recorder.run_dir / "work"
    work.mkdir()
    report_root = repo / "reports/stage_0003a"
    try:
        variants = pq.read_table(
            _repo_path(repo, config["inputs"]["variants"])
        ).to_pandas()
        variants = variants.loc[
            variants["comparability_group"].eq("chen2025_in_vitro_cis_trans_figure6")
        ].copy()
        variants["point_positions"] = variants["mutation"].map(_point_positions)
        eligible = variants.loc[variants["point_positions"].notna()].copy()
        mappings = {
            pdb_id: pd.read_csv(
                report_root / f"residue_mapping/{pdb_id.lower()}/mapping.csv"
            )
            for pdb_id in config["inputs"]["scaffold_states"]
        }
        wild_type = str(eligible.iloc[0]["WT_sequence"])
        records: list[dict[str, Any]] = [
            {
                "score_id": "WT",
                "variant_id": "WT",
                "full_sequence": wild_type,
                "mutation": "WT",
            }
        ]
        records.extend(
            {
                "score_id": str(row.variant_id).replace("-", "_"),
                "variant_id": str(row.variant_id),
                "full_sequence": str(row.full_mutant_sequence),
                "mutation": str(row.mutation),
            }
            for row in eligible.itertuples()
        )
        for record in records:
            for pdb_id, mapping in mappings.items():
                record[f"projected_{pdb_id}"] = _project(
                    mapping, str(record["full_sequence"])
                )
        esm_jobs = [
            {
                "job_id": f"{record['score_id']}-{pdb_id}",
                "pdb_id": pdb_id,
                "scaffold_id": "CasRx",
                "state": "binary" if pdb_id == "9M31" else "ternary",
                "structure_path": str(
                    repo / f"data/experimental_structures/{pdb_id.lower()}.cif"
                ),
                "protein_chain": "A",
                "sequence": record[f"projected_{pdb_id}"],
            }
            for record in records
            for pdb_id in config["inputs"]["scaffold_states"]
        ]
        esm_jobs_path = work / "esm_jobs.jsonl"
        esm_scores_path = work / "esm_scores.jsonl"
        _jsonl(esm_jobs_path, esm_jobs)
        _run(
            [
                str(_repo_path(repo, config["models"]["esm_if1_python"])),
                str(repo / "scripts/score_stage_0003a_multistate_esm.py"),
                "--jobs",
                str(esm_jobs_path),
                "--output",
                str(esm_scores_path),
                "--checkpoint",
                str(_repo_path(repo, config["models"]["esm_if1_checkpoint"])),
                "--device",
                str(config["execution"]["device"]),
            ],
            cwd=repo,
            stdout=work / "esm.stdout.log",
            stderr=work / "esm.stderr.log",
        )
        esm_scores = {
            str(row["job_id"]): row
            for row in (
                json.loads(line)
                for line in esm_scores_path.read_text(encoding="utf-8").splitlines()
                if line.strip()
            )
        }
        protein_scores = _proteinmpnn_scores(
            repo=repo, config=config, records=records, work=work
        )
        ligand_scores = _ligandmpnn_scores(
            repo=repo,
            config=config,
            records=records,
            mapping=mappings["9M31"],
            work=work,
        )
        wt_esm = {
            pdb_id: float(
                esm_scores[f"WT-{pdb_id}"]["mean_log_likelihood_per_resolved_residue"]
            )
            for pdb_id in config["inputs"]["scaffold_states"]
        }
        wt_protein = protein_scores["WT"]["proteinmpnn_mean_nll"]
        wt_ligand = ligand_scores["WT"]["ligandmpnn_mean_log_probability"]
        evaluated: list[dict[str, Any]] = []
        rank_scores: dict[str, dict[str, float]] = {}
        mapping = mappings["9M31"].set_index("biological_index_1")
        eligible_by_id = eligible.set_index("variant_id")
        for record in records[1:]:
            variant_id = str(record["variant_id"])
            activity = eligible_by_id.loc[variant_id]
            state_values = {
                pdb_id: float(
                    esm_scores[f"{record['score_id']}-{pdb_id}"][
                        "mean_log_likelihood_per_resolved_residue"
                    ]
                )
                for pdb_id in config["inputs"]["scaffold_states"]
            }
            rank_scores[variant_id] = state_values
            aggregate = aggregate_multistate_scores(state_values)
            positions = [int(item) for item in activity["point_positions"]]
            annotated = mapping.loc[positions]
            protein = protein_scores[str(record["score_id"])]
            ligand = ligand_scores[str(record["score_id"])]
            evaluated.append(
                {
                    "variant_id": variant_id,
                    "mutation": activity["mutation"],
                    "label": activity["active_inactive_partial_label"],
                    "cis_activity": activity["cis_activity"],
                    "trans_activity": activity["trans_activity"],
                    "numeric_is_approximate": activity["numeric_is_approximate"],
                    "mutated_positions": ";".join(map(str, positions)),
                    "mean_mutated_position_conservation": float(
                        annotated["conservation"].mean()
                    ),
                    "RNA_contact_mutations": int(
                        annotated["RNA_contact"]
                        .map(lambda value: str(value).lower() == "true")
                        .sum()
                    ),
                    "RNA_second_shell_mutations": int(
                        annotated["RNA_second_shell"]
                        .map(lambda value: str(value).lower() == "true")
                        .sum()
                    ),
                    "protein_core_mutations": int(
                        annotated["protein_core"]
                        .map(lambda value: str(value).lower() == "true")
                        .sum()
                    ),
                    "esm_binary_mean_log_probability": state_values["9M31"],
                    "esm_ternary_mean_log_probability": state_values["9M8Q"],
                    "esm_binary_delta_vs_wt": state_values["9M31"] - wt_esm["9M31"],
                    "esm_ternary_delta_vs_wt": state_values["9M8Q"] - wt_esm["9M8Q"],
                    "multi_state_mean_score": aggregate.multi_state_mean_score,
                    "multi_state_min_score": aggregate.multi_state_min_score,
                    "multi_state_variance": aggregate.multi_state_variance,
                    **protein,
                    "proteinmpnn_delta_nll_vs_wt": protein["proteinmpnn_mean_nll"]
                    - wt_protein,
                    **ligand,
                    "ligandmpnn_delta_log_probability_vs_wt": ligand[
                        "ligandmpnn_mean_log_probability"
                    ]
                    - wt_ligand,
                    "is_mock": False,
                    "evidence_level": 2,
                }
            )
        frame = pd.DataFrame.from_records(evaluated)
        rank_consistency = state_rank_consistency(rank_scores)
        correlations: list[dict[str, Any]] = []
        numeric = frame.loc[frame["cis_activity"].notna()].copy()
        metrics = [
            "esm_binary_delta_vs_wt",
            "esm_ternary_delta_vs_wt",
            "multi_state_min_score",
            "proteinmpnn_delta_nll_vs_wt",
            "ligandmpnn_delta_log_probability_vs_wt",
            "mean_mutated_position_conservation",
        ]
        for endpoint in ("cis_activity", "trans_activity"):
            for metric in metrics:
                coefficient = float(
                    spearmanr(numeric[metric], numeric[endpoint]).statistic
                )
                correlations.append(
                    {
                        "endpoint": endpoint,
                        "metric": metric,
                        "spearman_rho": coefficient,
                        "n": len(numeric),
                        "p_value": None,
                        "inference": "descriptive_only",
                    }
                )
        label_rows: list[dict[str, Any]] = []
        replicates = int(config["statistics"]["bootstrap_replicates"])
        for label, group in numeric.groupby("label"):
            values = group["multi_state_min_score"].to_numpy(dtype=float)
            low, high = _bootstrap_mean(
                values,
                replicates=replicates,
                seed=int(config["execution"]["seed"]),
            )
            label_rows.append(
                {
                    "label": label,
                    "n": len(group),
                    "multi_state_min_mean": float(values.mean()),
                    "bootstrap_ci_low": low,
                    "bootstrap_ci_high": high,
                    "unit": "variant within one assay group",
                }
            )
        output = report_root / "variant_retrospective.csv"
        correlation_path = report_root / "variant_retrospective_correlations.csv"
        labels_path = report_root / "variant_retrospective_label_summary.csv"
        atomic_write_text(output, frame.to_csv(index=False))
        atomic_write_text(
            correlation_path, pd.DataFrame(correlations).to_csv(index=False)
        )
        atomic_write_text(labels_path, pd.DataFrame(label_rows).to_csv(index=False))
        excluded = variants.loc[
            variants["point_positions"].isna(),
            [
                "variant_id",
                "mutation",
                "active_inactive_partial_label",
            ],
        ].copy()
        excluded["reason"] = (
            "indel changes sequence length; same-backbone score invalid"
        )
        excluded_path = report_root / "variant_retrospective_exclusions.csv"
        atomic_write_text(excluded_path, excluded.to_csv(index=False))
        summary = {
            "source_records_in_comparability_group": len(variants),
            "point_variants_scored": len(frame),
            "indels_excluded": len(excluded),
            "real_esm_state_scores": len(esm_scores),
            "real_proteinmpnn_scores": len(protein_scores),
            "real_ligandmpnn_scores": len(ligand_scores),
            "state_rank_consistency": rank_consistency,
            "correlations_are_descriptive": True,
            "significance_tests": "not_run_insufficient_class_counts",
            "activity_values_are_approximate_graph_readings": True,
            "fixed_or_forced_sequence_violations": 0,
            "rna_atomic_context_verified": True,
            "is_mock": False,
            "evidence_level": 2,
        }
        summary_path = report_root / "variant_retrospective_summary.json"
        atomic_write_text(
            summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        outputs = [output, correlation_path, labels_path, excluded_path, summary_path]
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[
                {"path": str(path.relative_to(repo)), "sha256": sha256_file(path)}
                for path in outputs
            ],
        )
        print(json.dumps({**summary, "run_dir": str(recorder.run_dir)}, indent=2))
        return 0
    except Exception as error:
        recorder.record_failure("variant_retrospective", str(error))
        recorder.finish(success=False)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
