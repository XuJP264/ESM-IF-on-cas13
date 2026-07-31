#!/usr/bin/env python
"""Build real multi-state masks, run native ESM scoring, and aggregate states."""

from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from cas13_if.config import load_config
from cas13_if.provenance import RunRecorder, atomic_write_text, sha256_file
from cas13_if.scoring.multistate import (
    StateMask,
    StateResidue,
    aggregate_multistate_scores,
    build_multistate_masks,
    select_state_combination,
    validate_fixed_tokens,
    validate_residue_maps,
)


def _yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"YAML root is not a mapping: {path}")
    return value


def _bool(series: pd.Series) -> pd.Series:
    if series.dtype == bool:
        return series
    return series.map(lambda value: str(value).lower() == "true")


def _jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    atomic_write_text(
        path, "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows)
    )


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    config_path = repo / "configs/stage_0003a_multistate.yaml"
    config = load_config(config_path)
    report_root = repo / str(config["inputs"]["report_root"])
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
    try:
        structure_config = _yaml(repo / str(config["inputs"]["structure_config"]))
        scaffold_frame = pd.read_csv(report_root / "scaffolds.csv").set_index(
            "scaffold_id"
        )
        jobs: list[dict[str, Any]] = []
        masks_summary: list[dict[str, Any]] = []
        mask_dir = report_root / "multistate_masks"
        mask_dir.mkdir(parents=True, exist_ok=True)
        for scaffold in structure_config["scaffolds"]:
            scaffold_id = str(scaffold["scaffold_id"])
            parent = str(scaffold_frame.loc[scaffold_id, "full_natural_sequence"])
            state_masks: list[StateMask] = []
            residue_maps: dict[str, list[StateResidue]] = {}
            for state in scaffold["states"]:
                pdb_id = str(state["pdb_id"]).upper()
                mapping = pd.read_csv(
                    report_root / f"residue_mapping/{pdb_id.lower()}/mapping.csv"
                )
                mapping = mapping.sort_values("full_scaffold_index_0")
                residues = [
                    StateResidue(
                        full_index_0=int(row.full_scaffold_index_0),
                        wild_type=str(row.full_scaffold_amino_acid),
                        mapping_status=str(row.mapping_status),
                    )
                    for row in mapping.itertuples()
                ]
                residue_maps[pdb_id] = residues
                high = mapping["mapping_confidence"].eq("high")
                hepn = mapping["HEPN_annotation"].notna()
                core = _bool(mapping["protein_core"].fillna(False))
                direct = _bool(mapping["RNA_contact"])
                second = _bool(mapping["RNA_second_shell"])
                breaks = _bool(mapping["chain_break_adjacent"])
                gate_failure = ~_bool(mapping["mapping_gate_passed"])
                hard = set(
                    mapping.loc[hepn | (high & direct), "full_scaffold_index_0"].astype(
                        int
                    )
                )
                risk = set(
                    mapping.loc[
                        gate_failure | core | direct | second | breaks,
                        "full_scaffold_index_0",
                    ].astype(int)
                )
                state_masks.append(
                    StateMask(
                        state=pdb_id,
                        sequence=parent,
                        hard_positions=frozenset(hard),
                        risk_positions=frozenset(risk),
                    )
                )
                coordinate = mapping.loc[
                    mapping["coordinate_index_0"].notna()
                ].sort_values("coordinate_index_0")
                coordinate_indices = (
                    coordinate["coordinate_index_0"].astype(int).tolist()
                )
                if coordinate_indices != list(range(len(coordinate_indices))):
                    raise ValueError(f"coordinate mapping is incomplete for {pdb_id}")
                projected = "".join(coordinate["full_scaffold_amino_acid"].astype(str))
                jobs.append(
                    {
                        "job_id": f"native-{scaffold_id}-{pdb_id}",
                        "pdb_id": pdb_id,
                        "scaffold_id": scaffold_id,
                        "state": str(state["state"]),
                        "structure_path": str(
                            repo / f"data/experimental_structures/{pdb_id.lower()}.cif"
                        ),
                        "protein_chain": str(state["protein_chain"]),
                        "sequence": projected,
                        "full_parent_sha256": __import__("hashlib")
                        .sha256(parent.encode("ascii"))
                        .hexdigest(),
                        "is_mock": False,
                    }
                )
            validate_residue_maps(residue_maps)
            masks = build_multistate_masks(state_masks)
            fixed = {
                index: parent[index] for index in masks.state_intersection_hard_mask
            }
            validate_fixed_tokens(
                {item.state: item.sequence for item in state_masks}, fixed
            )
            payload = {
                "scaffold_id": scaffold_id,
                "states": [item.state for item in state_masks],
                "full_parent_length": len(parent),
                "state_intersection_hard_mask": sorted(
                    masks.state_intersection_hard_mask
                ),
                "state_union_risk_mask": sorted(masks.state_union_risk_mask),
                "state_variable_hinge_mask": sorted(masks.state_variable_hinge_mask),
                "fixed_tokens": fixed,
                "is_mock": False,
                "evidence_level": 0,
            }
            atomic_write_text(
                mask_dir / f"{scaffold_id.lower()}.json",
                json.dumps(payload, indent=2, sort_keys=True) + "\n",
            )
            masks_summary.append(
                {
                    "scaffold_id": scaffold_id,
                    "state_count": len(state_masks),
                    "intersection_hard_positions": len(
                        masks.state_intersection_hard_mask
                    ),
                    "union_risk_positions": len(masks.state_union_risk_mask),
                    "variable_hinge_positions": len(masks.state_variable_hinge_mask),
                    "fixed_position_violations": 0,
                    "is_mock": False,
                }
            )
        jobs_path = work / "esm_jobs.jsonl"
        scores_path = work / "esm_scores.jsonl"
        _jsonl(jobs_path, jobs)
        command = [
            str(repo / str(config["models"]["esm_if1_python"])),
            str(repo / "scripts/score_stage_0003a_multistate_esm.py"),
            "--jobs",
            str(jobs_path),
            "--output",
            str(scores_path),
            "--checkpoint",
            str(repo / str(config["models"]["esm_if1_checkpoint"])),
            "--device",
            str(config["execution"]["device"]),
        ]
        environment = os.environ.copy()
        environment["PYTHONNOUSERSITE"] = "1"
        completed = subprocess.run(
            command,
            cwd=repo,
            env=environment,
            text=True,
            capture_output=True,
            check=False,
        )
        atomic_write_text(recorder.run_dir / "stdout.log", completed.stdout)
        atomic_write_text(recorder.run_dir / "stderr.log", completed.stderr)
        if completed.returncode != 0:
            raise RuntimeError(
                f"ESM multi-state worker failed ({completed.returncode}): "
                f"{completed.stderr[-2000:]}"
            )
        scores = [
            json.loads(line)
            for line in scores_path.read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        if any(bool(row["is_mock"]) for row in scores):
            raise ValueError("mock score entered the real multi-state result")
        combinations: list[dict[str, Any]] = []
        failures: list[dict[str, Any]] = []
        for scaffold_id, grouped in pd.DataFrame(scores).groupby("scaffold_id"):
            available = {
                str(row.pdb_id): str(row.state) for row in grouped.itertuples()
            }
            score_by_pdb = {
                str(row.pdb_id): float(row.mean_log_likelihood_per_resolved_residue)
                for row in grouped.itertuples()
            }
            for combination in config["state_combinations"]:
                try:
                    selected = select_state_combination(available, str(combination))
                except ValueError as error:
                    failures.append(
                        {
                            "scaffold_id": scaffold_id,
                            "combination": combination,
                            "status": "not_applicable",
                            "reason": str(error),
                            "is_mock": False,
                        }
                    )
                    continue
                aggregate = aggregate_multistate_scores(
                    {state: score_by_pdb[state] for state in selected}
                )
                combinations.append(
                    {
                        "scaffold_id": scaffold_id,
                        "combination": combination,
                        "pdb_ids": ";".join(selected),
                        "single_state_score": json.dumps(
                            aggregate.single_state_score, sort_keys=True
                        ),
                        "normalized_weights": json.dumps(
                            aggregate.normalized_weights, sort_keys=True
                        ),
                        "multi_state_mean_score": aggregate.multi_state_mean_score,
                        "multi_state_min_score": aggregate.multi_state_min_score,
                        "multi_state_variance": aggregate.multi_state_variance,
                        "score_normalization": config["experiment"][
                            "normalized_score_policy"
                        ],
                        "state_rank_consistency": None,
                        "rank_consistency_reason": "one native sequence only",
                        "is_mock": False,
                        "evidence_level": 2,
                    }
                )
        score_frame = pd.DataFrame.from_records(scores)
        combination_frame = pd.DataFrame.from_records(combinations)
        failure_frame = pd.DataFrame.from_records(failures)
        mask_frame = pd.DataFrame.from_records(masks_summary)
        atomic_write_text(
            report_root / "multistate_native_scores.csv",
            score_frame.to_csv(index=False),
        )
        atomic_write_text(
            report_root / "multistate_aggregates.csv",
            combination_frame.to_csv(index=False),
        )
        atomic_write_text(
            report_root / "multistate_failures.csv", failure_frame.to_csv(index=False)
        )
        atomic_write_text(
            report_root / "multistate_masks.csv", mask_frame.to_csv(index=False)
        )
        summary = {
            "real_state_scores": len(scores),
            "real_multistate_aggregates": len(combinations),
            "not_applicable_combinations": len(failures),
            "scaffolds": len(masks_summary),
            "fixed_position_violations": int(
                mask_frame["fixed_position_violations"].sum()
            ),
            "device": str(config["execution"]["device"]),
            "score_policy": config["experiment"]["normalized_score_policy"],
            "is_mock": False,
            "evidence_level": 2,
        }
        atomic_write_text(
            report_root / "multistate_summary.json",
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
        )
        outputs = [
            report_root / name
            for name in (
                "multistate_native_scores.csv",
                "multistate_aggregates.csv",
                "multistate_failures.csv",
                "multistate_masks.csv",
                "multistate_summary.json",
            )
        ]
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
        recorder.record_failure("stage_0003a_multistate", str(error))
        recorder.finish(success=False)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
