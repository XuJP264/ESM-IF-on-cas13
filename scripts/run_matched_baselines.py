#!/usr/bin/env python
"""Run the preregistered small real VI-D matched-baseline comparison."""

from __future__ import annotations

import argparse
import html
import json
import os
import platform
import subprocess
import sys
from collections import Counter
from pathlib import Path
from statistics import fmean
from typing import Any

import matplotlib.pyplot as plt
import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]
import yaml

from cas13_if.backends.baselines import (
    BASELINE_ALPHABET,
    MatchedRandomMutationBackend,
    MsaProfileBackend,
)
from cas13_if.config import ConfigDict, load_config
from cas13_if.evaluation.matching import (
    add_identity_metrics,
    paired_seed_statistics,
    position_set_hash,
    select_balanced_candidates,
)
from cas13_if.evaluation.metrics import native_recovery
from cas13_if.evaluation.regions import build_structure_regions
from cas13_if.novelty.metrics import pairwise_candidate_identity
from cas13_if.novelty.pipeline import (
    NoveltyThresholds,
    evaluate_candidate_novelty,
    run_mmseqs_atlas_search,
)
from cas13_if.provenance import (
    RunRecorder,
    atomic_write_text,
    git_metadata,
    sha256_file,
)
from cas13_if.schemas import Candidate, EvidenceLevel, SampleRequest
from cas13_if.structures.parser import parse_structure, protein_chain_sequence

METHODS = [
    "matched_random_mutation",
    "msa_profile_sampling",
    "unconstrained_esm_if1",
    "catalytic_only_fixed_esm_if1",
    "conservation_constrained_esm_if1",
    "conservation_rna_contact_esm_if1",
    "proteinmpnn",
    "ligandmpnn",
    "esm_if1_ligandmpnn_consensus",
]


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    return parser.parse_args()


def _deep_merge(base: ConfigDict, override: ConfigDict) -> ConfigDict:
    merged: ConfigDict = json.loads(json.dumps(base))
    for key, value in override.items():
        if key == "defaults_from":
            continue
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _resolved_config(repo: Path, path: Path) -> ConfigDict:
    config = load_config(path)
    defaults = config.get("defaults_from")
    if defaults is None:
        return config
    base_path = Path(str(defaults))
    if not base_path.is_absolute():
        base_path = repo / base_path
    return _deep_merge(load_config(base_path), config)


def _repo_path(repo: Path, value: Any) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else repo / path


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"YAML root is not a mapping: {path}")
    return value


def _parent_and_fixed(
    structure: Path, chain: str, functional_manifest: Path
) -> tuple[str, dict[int, str], list[Any]]:
    deposited, residue_keys = protein_chain_sequence(parse_structure(structure), chain)
    entries = _load_yaml(functional_manifest)["structures"]["6E9F"]["residues"]
    by_pdb = {
        (key.residue_number, key.insertion_code): index
        for index, key in enumerate(residue_keys)
    }
    fixed = {
        by_pdb[
            (int(entry["pdb_residue_number"]), str(entry.get("insertion_code") or ""))
        ]: str(entry["biological_amino_acid"]).upper()
        for entry in entries
    }
    parent = list(deposited)
    for index, token in fixed.items():
        parent[index] = token
    return "".join(parent), fixed, residue_keys


def _run_subprocess(
    command: list[str], *, repo: Path, log_prefix: Path, device: str
) -> None:
    environment = os.environ.copy()
    environment["PYTHONNOUSERSITE"] = "1"
    if device == "cpu":
        environment["CUDA_VISIBLE_DEVICES"] = ""
    completed = subprocess.run(
        command,
        cwd=repo,
        env=environment,
        text=True,
        capture_output=True,
        check=False,
    )
    atomic_write_text(log_prefix.with_suffix(".stdout.log"), completed.stdout)
    atomic_write_text(log_prefix.with_suffix(".stderr.log"), completed.stderr)
    if completed.returncode != 0:
        raise RuntimeError(
            f"subprocess failed ({completed.returncode}): {' '.join(command)}; "
            f"see {log_prefix}.stderr.log"
        )


def _read_proposals(path: Path) -> list[dict[str, Any]]:
    rows: list[dict[str, Any]] = []
    for line_number, line in enumerate(
        path.read_text(encoding="utf-8").splitlines(), start=1
    ):
        if not line.strip():
            continue
        raw = json.loads(line)
        candidate = raw["candidate"]
        if bool(candidate["is_mock"]):
            raise ValueError(f"mock proposal rejected at {path}:{line_number}")
        fixed = {
            int(index): str(token)
            for index, token in candidate["fixed_positions"].items()
        }
        rows.append(
            {
                **candidate,
                "fixed_positions": fixed,
                "method": str(raw["method"]),
                "seed_block": int(raw["seed_block"]),
                "proposal_index": int(raw["proposal_index"]),
                "actual_model_seed": int(candidate["metadata"]["actual_model_seed"]),
            }
        )
    return rows


def _local_baseline_proposals(
    *,
    config: ConfigDict,
    parent: str,
    fixed: dict[int, str],
    mapping_path: Path,
    conservation_path: Path,
) -> list[dict[str, Any]]:
    mapping = pd.read_csv(mapping_path)
    mapping = mapping.loc[mapping["coordinate_index_0"].notna()].copy()
    mapping["coordinate_index_0"] = mapping["coordinate_index_0"].astype(int)
    mapping = mapping.sort_values("coordinate_index_0")
    if mapping["coordinate_index_0"].tolist() != list(range(len(parent))):
        raise ValueError("mapping does not cover every coordinate position")
    conservation = {
        int(row["column"]): row for row in pq.read_table(conservation_path).to_pylist()
    }
    frequencies = [
        conservation[int(row.msa_column_0)]["weighted_frequencies"]
        for row in mapping.itertuples()
    ]
    sampling = config["sampling"]
    backends = {
        "matched_random_mutation": MatchedRandomMutationBackend(
            mutation_probability=float(sampling["random_mutation_probability"])
        ),
        "msa_profile_sampling": MsaProfileBackend(frequencies),
    }
    rows: list[dict[str, Any]] = []
    for method, backend in backends.items():
        backend.load()
        temperature = (
            1.0
            if method == "matched_random_mutation"
            else float(sampling["temperatures"]["msa_profile_sampling"])
        )
        for seed_block_value in sampling["seed_blocks"]:
            seed_block = int(seed_block_value)
            for proposal_index in range(int(sampling["proposals_per_seed"])):
                actual_seed = seed_block + proposal_index
                candidate = backend.sample(
                    SampleRequest(
                        scaffold_id="6E9F-A",
                        structure_path="profile_or_random_no_structure_input",
                        parent_sequence=parent,
                        count=1,
                        temperature=temperature,
                        seed=actual_seed,
                        fixed_positions=fixed,
                        protein_chains=["A"],
                    )
                )[0]
                selected_confidences = [
                    trace.probabilities[BASELINE_ALPHABET.index(trace.selected_token)]
                    for trace in candidate.traces
                ]
                payload = candidate.model_dump(mode="json")
                payload["candidate_id"] = (
                    f"{method}-seed{seed_block}-proposal{proposal_index}-"
                    f"{payload['candidate_id']}"
                )
                payload["backend"] = method
                payload["traces"] = []
                payload["metadata"] = {
                    **payload["metadata"],
                    "method": method,
                    "seed_block": seed_block,
                    "proposal_index": proposal_index,
                    "actual_model_seed": actual_seed,
                    "selected_token_probabilities": selected_confidences,
                }
                rows.append(
                    {
                        **payload,
                        "fixed_positions": {
                            int(index): str(token)
                            for index, token in payload["fixed_positions"].items()
                        },
                        "method": method,
                        "seed_block": seed_block,
                        "proposal_index": proposal_index,
                        "actual_model_seed": actual_seed,
                    }
                )
    return rows


def _consensus_proposals(proposals: list[dict[str, Any]]) -> list[dict[str, Any]]:
    by_key = {
        (str(row["method"]), int(row["seed_block"]), int(row["proposal_index"])): row
        for row in proposals
    }
    rows: list[dict[str, Any]] = []
    keys = sorted(
        (seed, proposal)
        for method, seed, proposal in by_key
        if method == "unconstrained_esm_if1"
    )
    for seed_block, proposal_index in keys:
        esm = by_key[("unconstrained_esm_if1", seed_block, proposal_index)]
        ligand = by_key[("ligandmpnn", seed_block, proposal_index)]
        esm_confidence = esm["metadata"]["selected_token_probabilities"]
        ligand_confidence = ligand["metadata"]["selected_token_probabilities"]
        tokens: list[str] = []
        agreements = 0
        selected_confidences: list[float] = []
        for index, (esm_token, ligand_token) in enumerate(
            zip(esm["sequence"], ligand["sequence"], strict=True)
        ):
            if esm_token == ligand_token:
                token = esm_token
                agreements += 1
                confidence = max(esm_confidence[index], ligand_confidence[index])
            elif esm_confidence[index] >= ligand_confidence[index]:
                token = esm_token
                confidence = esm_confidence[index]
            else:
                token = ligand_token
                confidence = ligand_confidence[index]
            tokens.append(token)
            selected_confidences.append(float(confidence))
        sequence = "".join(tokens)
        fixed = {
            int(index): str(token) for index, token in esm["fixed_positions"].items()
        }
        candidate = Candidate(
            candidate_id=(
                "esm-if1-ligandmpnn-consensus-"
                f"seed{seed_block}-proposal{proposal_index}"
            ),
            scaffold_id="6E9F-A",
            backend="esm_if1_ligandmpnn_consensus",
            sequence=sequence,
            parent_sequence=str(esm["parent_sequence"]),
            seed=int(esm["actual_model_seed"]),
            temperature=float(esm["temperature"]),
            is_mock=False,
            evidence_level=EvidenceLevel.INVERSE_FOLDING_COMPATIBILITY,
            fixed_positions=fixed,
            metadata={
                "method": "esm_if1_ligandmpnn_consensus",
                "seed_block": seed_block,
                "proposal_index": proposal_index,
                "actual_model_seed": int(esm["actual_model_seed"]),
                "source_candidate_ids": [esm["candidate_id"], ligand["candidate_id"]],
                "source_model_exact_agreement": agreements / len(sequence),
                "disagreement_rule": "higher_selected_token_probability",
                "selected_token_probabilities": selected_confidences,
            },
        ).model_dump(mode="json")
        rows.append(
            {
                **candidate,
                "fixed_positions": fixed,
                "method": "esm_if1_ligandmpnn_consensus",
                "seed_block": seed_block,
                "proposal_index": proposal_index,
                "actual_model_seed": int(esm["actual_model_seed"]),
            }
        )
    return rows


def _write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    atomic_write_text(
        path, "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows)
    )


def _method_table(
    candidate_rows: list[dict[str, Any]],
    proposals: list[dict[str, Any]],
    *,
    fixed_hash: str,
    free_hash: str,
    fixed_count: int,
    free_count: int,
) -> pd.DataFrame:
    proposal_counts = Counter(str(row["method"]) for row in proposals)
    rows: list[dict[str, Any]] = []
    for method in METHODS:
        selected = [row for row in candidate_rows if row["method"] == method]
        atlas = [
            float(row["maximum_atlas_identity"])
            for row in selected
            if row["maximum_atlas_identity"] is not None
        ]
        rows.append(
            {
                "method": method,
                "is_mock": False,
                "actual_proposals": proposal_counts[method],
                "actual_selected_candidates": len(selected),
                "seed_count": len({int(row["seed_block"]) for row in selected}),
                "fixed_position_count": fixed_count,
                "free_position_count": free_count,
                "fixed_position_hash": fixed_hash,
                "free_position_hash": free_hash,
                "mean_conditional_log_likelihood": fmean(
                    float(row["conditional_log_likelihood"]) for row in selected
                ),
                "mean_perplexity": fmean(float(row["perplexity"]) for row in selected),
                "mean_parent_identity": fmean(
                    float(row["parent_identity"]) for row in selected
                ),
                "mean_designed_position_identity": fmean(
                    float(row["designed_position_identity"]) for row in selected
                ),
                "mean_maximum_atlas_identity_among_coverage_hits": (
                    fmean(atlas) if atlas else None
                ),
                "atlas_coverage_hit_count": len(atlas),
                "fixed_position_violations": sum(
                    int(row["fixed_position_violations"]) for row in selected
                ),
                "level1_pass_count": sum(
                    bool(row["passes_level1_novelty"]) for row in selected
                ),
                "mean_buried_core_recovery": fmean(
                    float(row["buried_core_recovery"]) for row in selected
                ),
                "mean_rna_interface_recovery": fmean(
                    float(row["rna_interface_recovery"]) for row in selected
                ),
                "mean_rna_second_shell_recovery": fmean(
                    float(row["rna_second_shell_recovery"]) for row in selected
                ),
                "mean_hepn_region_recovery": fmean(
                    float(row["hepn_region_recovery"]) for row in selected
                ),
                "mean_model_agreement": fmean(
                    float(row["model_agreement"]) for row in selected
                ),
                "candidate_diversity_one_minus_mean_pairwise_identity": fmean(
                    float(row["candidate_diversity"]) for row in selected
                ),
                "evidence_boundary": (
                    "Level 1 only if novelty gates pass; Level 2 from ESM score; "
                    "no functional validity claim"
                ),
            }
        )
    return pd.DataFrame(rows)


def _render_figures(
    output_dir: Path, candidates: pd.DataFrame, funnel: pd.DataFrame
) -> None:
    figure_dir = output_dir / "figures"
    figure_dir.mkdir(parents=True, exist_ok=False)
    method_order = METHODS
    labels = [method.replace("_", "\n") for method in method_order]
    for metric, filename, ylabel in (
        ("parent_identity", "matched_identity.png", "Parent identity"),
        (
            "conditional_log_likelihood",
            "conditional_log_likelihood.png",
            "ESM-IF1 conditional log-likelihood",
        ),
    ):
        figure, axis = plt.subplots(figsize=(13, 5))
        for index, method in enumerate(method_order):
            values = candidates.loc[candidates["method"].eq(method), metric]
            axis.scatter([index] * len(values), values, s=35)
        axis.set_xticks(range(len(labels)), labels, fontsize=7)
        axis.set_ylabel(ylabel)
        axis.grid(axis="y", alpha=0.25)
        figure.tight_layout()
        figure.savefig(figure_dir / filename, dpi=180)
        plt.close(figure)
    region_means = candidates.groupby("method", sort=False)[
        [
            "buried_core_recovery",
            "rna_interface_recovery",
            "rna_second_shell_recovery",
            "hepn_region_recovery",
        ]
    ].mean()
    region_means = region_means.reindex(method_order)
    axis = region_means.plot(kind="bar", figsize=(13, 5))
    axis.set_ylabel("Recovery")
    axis.set_xticklabels(labels, rotation=0, fontsize=7)
    axis.legend(fontsize=8)
    axis.figure.tight_layout()
    axis.figure.savefig(figure_dir / "per_region_recovery.png", dpi=180)
    plt.close(axis.figure)
    selected = funnel.loc[funnel["stage"].isin(["proposed", "selected", "level1_pass"])]
    pivot = selected.pivot(index="method", columns="stage", values="count").reindex(
        method_order
    )
    axis = pivot.plot(kind="bar", figsize=(13, 5))
    axis.set_ylabel("Candidates")
    axis.set_xticklabels(labels, rotation=0, fontsize=7)
    axis.figure.tight_layout()
    axis.figure.savefig(figure_dir / "candidate_funnel.png", dpi=180)
    plt.close(axis.figure)


def _report_text(
    methods: pd.DataFrame,
    candidates: pd.DataFrame,
    statistics: pd.DataFrame,
    funnel: pd.DataFrame,
    summary: dict[str, Any],
) -> str:
    identity = candidates.groupby("method")[
        ["parent_identity", "designed_position_identity"]
    ].agg(["min", "mean", "max"])
    return f"""# Matched VI-D baseline report

This is a **real, preregistered small CPU comparison** on the 6E9F/EsCas13d
coordinate scaffold. No mock candidate is included. The independent unit is
the seed block (n={summary["seed_count"]}), not each generated candidate; with
two seeds, confidence intervals and p-values are low-power descriptive checks,
not evidence of a method winner.

Level 1 is assigned only to individual rows that pass every sequence novelty
and QC gate. Level 2 denotes genuine ESM-IF1 compatibility scoring against the
experimental backbone. Neither level establishes functional activity. There
is no Level 3 refold evidence in this report and no Level 4 wet-lab evidence;
no sequence is described as an effective or validated Cas13.

All methods use the identical hard-fixed hash `{summary["fixed_position_hash"]}`
and free-position hash `{summary["free_position_hash"]}`. The four biological
HEPN R/H residues are restored and fixed in every method. Conservation and RNA
information alter allowed-token proposal support only inside the common free
set. `catalytic_only_fixed_esm_if1` is a technical replicate of the common-mask
ESM condition, as preregistered, rather than a distinct biological treatment.

## Method-level complete metric table

{methods.to_markdown(index=False)}

## Matched identity distribution

{identity.to_markdown()}

## Candidate funnel

{funnel.to_markdown(index=False)}

## Paired seed-level statistics

{statistics.to_markdown(index=False)}

Missing Atlas hits at 80% query coverage are failed closed and are not treated
as proof of extreme novelty. Every proposal excluded by identity matching and
every selected-candidate QC failure is retained in `failure_analysis.csv`.
"""


def _file_entries(repo: Path, paths: list[Path]) -> list[dict[str, Any]]:
    return [
        {
            "path": str(path.relative_to(repo)),
            "size_bytes": path.stat().st_size,
            "sha256": sha256_file(path),
        }
        for path in paths
    ]


def main() -> int:
    arguments = _arguments()
    repo = Path(__file__).resolve().parents[1]
    config_path = arguments.config.resolve()
    config = _resolved_config(repo, config_path)
    execution = config["execution"]
    git = git_metadata(repo)
    if bool(execution.get("require_clean_git", True)) and bool(git["dirty"]):
        raise RuntimeError("formal matched baseline requires a clean git worktree")
    recorder = RunRecorder(
        root=repo / "results/runs",
        experiment=str(config["experiment"]["name"]),
        resolved_config=config,
        command=sys.argv,
        repo_root=repo,
        is_mock=False,
    )
    report_dir = _repo_path(repo, config["outputs"]["canonical_report_dir"])
    work_dir = recorder.run_dir / "matched_baselines"
    work_dir.mkdir(parents=True, exist_ok=False)
    try:
        if report_dir.exists():
            raise FileExistsError(f"refusing to overwrite report: {report_dir}")
        inputs = config["inputs"]
        models = config["models"]
        environments = config["environments"]
        sampling = config["sampling"]
        structure_pdb = _repo_path(repo, inputs["structure_pdb"])
        structure_cif = _repo_path(repo, inputs["structure_cif"])
        mapping_path = _repo_path(repo, inputs["mapping_csv"])
        conservation_path = _repo_path(repo, inputs["conservation_parquet"])
        functional_path = _repo_path(repo, inputs["functional_manifest"])
        atlas_fasta = _repo_path(repo, inputs["atlas_fasta"])
        model_paths = [
            _repo_path(repo, models[key])
            for key in (
                "esm_if1_checkpoint",
                "proteinmpnn_checkpoint",
                "ligandmpnn_checkpoint",
                "ligandmpnn_protein_checkpoint",
                "ligandmpnn_soluble_checkpoint",
            )
        ]
        input_paths = [
            config_path,
            structure_pdb,
            structure_cif,
            mapping_path,
            conservation_path,
            functional_path,
            atlas_fasta,
            *model_paths,
        ]
        atomic_write_text(
            recorder.run_dir / "input_manifest.json",
            json.dumps(
                {"files": _file_entries(repo, input_paths), "is_mock": False},
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
        resolved_path = work_dir / "resolved_execution_config.yaml"
        atomic_write_text(resolved_path, yaml.safe_dump(config, sort_keys=True))
        parent, fixed, _ = _parent_and_fixed(
            structure_pdb, str(inputs["chain_id"]), functional_path
        )
        fixed_set = set(fixed)
        free_set = set(range(len(parent))).difference(fixed_set)
        fixed_hash = position_set_hash(fixed_set)
        free_hash = position_set_hash(free_set)
        constraints = config["constraints"]
        regions, region_annotations = build_structure_regions(
            structure_path=structure_cif,
            protein_chain=str(inputs["chain_id"]),
            crrna_chains={str(value) for value in inputs["crrna_chains"]},
            target_rna_chains={str(value) for value in inputs["target_rna_chains"]},
            hepn_positions=fixed_set,
            direct_cutoff=float(constraints["direct_rna_cutoff_angstrom"]),
            second_shell_cutoff=float(constraints["second_shell_cutoff_angstrom"]),
            buried_rsa_threshold=float(constraints["buried_rsa_threshold"]),
        )
        _write_jsonl(work_dir / "region_annotations.jsonl", region_annotations)
        device = str(execution.get("device", "cpu"))
        esm_proposals_path = work_dir / "esm_proposals.jsonl"
        mpnn_proposals_path = work_dir / "mpnn_proposals.jsonl"
        _run_subprocess(
            [
                str(_repo_path(repo, environments["esm_if1_python"])),
                str(repo / "scripts/generate_matched_esm.py"),
                "--config",
                str(resolved_path),
                "--output",
                str(esm_proposals_path),
            ],
            repo=repo,
            log_prefix=work_dir / "esm_generation",
            device=device,
        )
        _run_subprocess(
            [
                str(_repo_path(repo, environments["ligandmpnn_python"])),
                str(repo / "scripts/generate_matched_mpnn.py"),
                "--config",
                str(resolved_path),
                "--output",
                str(mpnn_proposals_path),
            ],
            repo=repo,
            log_prefix=work_dir / "mpnn_generation",
            device=device,
        )
        proposals = [
            *_read_proposals(esm_proposals_path),
            *_read_proposals(mpnn_proposals_path),
            *_local_baseline_proposals(
                config=config,
                parent=parent,
                fixed=fixed,
                mapping_path=mapping_path,
                conservation_path=conservation_path,
            ),
        ]
        proposals.extend(_consensus_proposals(proposals))
        expected_proposals = (
            len(METHODS)
            * len(sampling["seed_blocks"])
            * int(sampling["proposals_per_seed"])
        )
        if len(proposals) != expected_proposals:
            raise ValueError(
                f"proposal count {len(proposals)} != expected {expected_proposals}"
            )
        for proposal in proposals:
            if bool(proposal["is_mock"]):
                raise ValueError("formal proposal set contains mock output")
            proposal_fixed = {
                int(index): str(token)
                for index, token in proposal["fixed_positions"].items()
            }
            if proposal_fixed != fixed:
                raise ValueError(f"fixed mask differs for {proposal['candidate_id']}")
            if any(
                proposal["sequence"][index] != token for index, token in fixed.items()
            ):
                raise ValueError(
                    f"fixed-position violation: {proposal['candidate_id']}"
                )
            proposal["fixed_position_hash"] = fixed_hash
            proposal["free_position_hash"] = free_hash
        proposals = add_identity_metrics(
            proposals, parent_sequence=parent, designed_positions=free_set
        )
        matching = config["matching"]
        selected, matching_failures = select_balanced_candidates(
            proposals,
            methods=METHODS,
            seed_blocks=[int(value) for value in sampling["seed_blocks"]],
            minimum_parent_identity=float(matching["minimum_parent_identity"]),
            maximum_parent_identity=float(matching["maximum_parent_identity"]),
            minimum_designed_identity=float(matching["minimum_designed_identity"]),
            maximum_designed_identity=float(matching["maximum_designed_identity"]),
            target_identity=float(matching["target_identity"]),
        )
        selected_path = work_dir / "selected_for_esm_score.jsonl"
        compact_selected = [
            {
                "candidate_id": row["candidate_id"],
                "sequence": row["sequence"],
                "actual_model_seed": row["actual_model_seed"],
            }
            for row in selected
        ]
        _write_jsonl(selected_path, compact_selected)
        score_path = work_dir / "esm_scores.jsonl"
        _run_subprocess(
            [
                str(_repo_path(repo, environments["esm_if1_python"])),
                str(repo / "scripts/score_matched_esm.py"),
                "--config",
                str(resolved_path),
                "--candidates",
                str(selected_path),
                "--output",
                str(score_path),
            ],
            repo=repo,
            log_prefix=work_dir / "esm_scoring",
            device=device,
        )
        scores = {
            row["candidate_id"]: row["score"]
            for row in (
                json.loads(line)
                for line in score_path.read_text(encoding="utf-8").splitlines()
                if line.strip()
            )
        }
        novelty_input = [
            {
                "candidate_id": row["candidate_id"],
                "pdb_id": "6E9F",
                "method": row["method"],
                "scaffold_id": "6E9F-A",
                "sequence": row["sequence"],
                "parent_sequence": parent,
                "temperature": row["temperature"],
                "seed": row["actual_model_seed"],
                "fixed_positions": sorted(fixed),
                "fixed_position_violations": 0,
                "source_is_mock": False,
                "source_evidence_level": row["evidence_level"],
            }
            for row in selected
        ]
        novelty_config = config["novelty"]
        atlas_hits, mmseqs_command = run_mmseqs_atlas_search(
            candidates=novelty_input,
            atlas_fasta=atlas_fasta,
            output_dir=work_dir / "atlas_novelty",
            executable=_repo_path(repo, novelty_config["mmseqs_executable"]),
            threads=int(novelty_config["threads"]),
            sensitivity=float(novelty_config["sensitivity"]),
            minimum_query_coverage=float(novelty_config["minimum_query_coverage"]),
            maximum_evalue=float(novelty_config["maximum_evalue"]),
            maximum_sequences=int(novelty_config["maximum_sequences"]),
        )
        novelty_rows, novelty_summary = evaluate_candidate_novelty(
            novelty_input,
            atlas_hits,
            NoveltyThresholds(
                maximum_parent_identity=float(matching["maximum_parent_identity"]),
                maximum_atlas_identity=float(novelty_config["maximum_atlas_identity"]),
                maximum_homopolymer_length=int(
                    novelty_config["maximum_homopolymer_length"]
                ),
                maximum_low_complexity_windows=int(
                    novelty_config["maximum_low_complexity_windows"]
                ),
                minimum_designed_position_entropy=float(
                    novelty_config["minimum_designed_position_entropy"]
                ),
                low_complexity_window=int(novelty_config["low_complexity_window"]),
                low_complexity_maximum_fraction=float(
                    novelty_config["low_complexity_maximum_fraction"]
                ),
            ),
        )
        novelty_by_id = {row["candidate_id"]: row for row in novelty_rows}
        evaluated: list[dict[str, Any]] = []
        for row in selected:
            recovery = native_recovery(
                row["sequence"],
                parent,
                designed_positions=free_set,
                fixed_positions=fixed,
                regions=regions,
            )
            score = scores[row["candidate_id"]]
            novelty = novelty_by_id[row["candidate_id"]]
            evaluated.append(
                {
                    **row,
                    "conditional_log_likelihood": score["conditional_log_likelihood"],
                    "perplexity": score["perplexity"],
                    "per_residue_log_probabilities": score[
                        "per_residue_log_probabilities"
                    ],
                    "fixed_position_violations": recovery.fixed_position_violations,
                    "maximum_atlas_identity": novelty["maximum_atlas_identity"],
                    "maximum_atlas_hit": novelty["maximum_atlas_hit"],
                    "atlas_query_coverage": novelty["maximum_atlas_query_coverage"],
                    "low_complexity_failure": int(
                        "low_complexity_windows_above_threshold"
                        in novelty["novelty_filter_failures"]
                    ),
                    "homopolymer_failure": int(
                        "homopolymer_above_threshold"
                        in novelty["novelty_filter_failures"]
                    ),
                    "novelty_filter_failures": novelty["novelty_filter_failures"],
                    "passes_level1_novelty": novelty["passes_level1_novelty"],
                    "buried_core_recovery": recovery.regions["buried_core"],
                    "rna_interface_recovery": recovery.regions["rna_interface"],
                    "rna_second_shell_recovery": recovery.regions["rna_second_shell"],
                    "hepn_region_recovery": recovery.regions["hepn_region"],
                    "crrna_interface_recovery": recovery.regions["crrna_interface"],
                    "target_rna_interface_recovery": recovery.regions[
                        "target_rna_interface"
                    ],
                }
            )
        consensus_by_seed = {
            int(row["seed_block"]): str(row["sequence"])
            for row in evaluated
            if row["method"] == "esm_if1_ligandmpnn_consensus"
        }
        diversity_by_method = {
            method: (
                1.0 - fmean(identities)
                if (
                    identities := pairwise_candidate_identity(
                        [
                            str(row["sequence"])
                            for row in evaluated
                            if row["method"] == method
                        ]
                    )
                )
                else 0.0
            )
            for method in METHODS
        }
        for row in evaluated:
            row["model_agreement"] = sum(
                left == right
                for left, right in zip(
                    row["sequence"],
                    consensus_by_seed[int(row["seed_block"])],
                    strict=True,
                )
            ) / len(parent)
            row["candidate_diversity"] = diversity_by_method[row["method"]]
        if any(int(row["fixed_position_violations"]) != 0 for row in evaluated):
            raise ValueError("selected candidates contain fixed-position violations")
        if {row["fixed_position_hash"] for row in evaluated} != {fixed_hash}:
            raise ValueError("selected methods do not share one fixed-position hash")
        if {row["free_position_hash"] for row in evaluated} != {free_hash}:
            raise ValueError("selected methods do not share one free-position hash")

        report_dir.mkdir(parents=True, exist_ok=False)
        candidate_rows: list[dict[str, Any]] = []
        for row in evaluated:
            metadata = dict(row["metadata"])
            metadata.pop("selected_token_probabilities", None)
            canonical = {
                key: value
                for key, value in row.items()
                if key not in {"traces", "per_residue_log_probabilities", "metadata"}
            }
            canonical["metadata"] = metadata
            canonical["is_mock"] = False
            canonical["evidence_level_generation"] = row["evidence_level"]
            canonical["evidence_level_compatibility"] = 2
            canonical["claim_label"] = "computational_candidate"
            candidate_rows.append(canonical)
        _write_jsonl(report_dir / "candidates.jsonl", candidate_rows)
        candidate_frame = pd.DataFrame(candidate_rows)
        per_region = candidate_frame[
            [
                "candidate_id",
                "method",
                "seed_block",
                "buried_core_recovery",
                "rna_interface_recovery",
                "crrna_interface_recovery",
                "target_rna_interface_recovery",
                "rna_second_shell_recovery",
                "hepn_region_recovery",
            ]
        ]
        atomic_write_text(
            report_dir / "per_region_metrics.csv", per_region.to_csv(index=False)
        )
        methods = _method_table(
            candidate_rows,
            proposals,
            fixed_hash=fixed_hash,
            free_hash=free_hash,
            fixed_count=len(fixed_set),
            free_count=len(free_set),
        )
        atomic_write_text(report_dir / "methods_table.csv", methods.to_csv(index=False))
        proposal_counts = Counter(str(row["method"]) for row in proposals)
        eligible_counts = Counter(
            str(row["method"])
            for row in proposals
            if float(matching["minimum_parent_identity"])
            <= float(row["parent_identity"])
            <= float(matching["maximum_parent_identity"])
            and float(matching["minimum_designed_identity"])
            <= float(row["designed_position_identity"])
            <= float(matching["maximum_designed_identity"])
        )
        funnel_rows: list[dict[str, Any]] = []
        for method in METHODS:
            method_rows = [row for row in candidate_rows if row["method"] == method]
            stages = {
                "proposed": proposal_counts[method],
                "identity_interval_eligible": eligible_counts[method],
                "selected": len(method_rows),
                "fixed_zero_violation": sum(
                    int(row["fixed_position_violations"]) == 0 for row in method_rows
                ),
                "low_complexity_pass": sum(
                    int(row["low_complexity_failure"]) == 0 for row in method_rows
                ),
                "homopolymer_pass": sum(
                    int(row["homopolymer_failure"]) == 0 for row in method_rows
                ),
                "atlas_coverage_hit": sum(
                    row["maximum_atlas_identity"] is not None for row in method_rows
                ),
                "level1_pass": sum(
                    bool(row["passes_level1_novelty"]) for row in method_rows
                ),
                "level2_esm_scored": len(method_rows),
            }
            funnel_rows.extend(
                {"method": method, "stage": stage, "count": count}
                for stage, count in stages.items()
            )
        funnel = pd.DataFrame(funnel_rows)
        atomic_write_text(
            report_dir / "candidate_funnel.csv", funnel.to_csv(index=False)
        )
        failure_rows = list(matching_failures)
        for row in candidate_rows:
            for reason in row["novelty_filter_failures"]:
                failure_rows.append(
                    {
                        "candidate_id": row["candidate_id"],
                        "method": row["method"],
                        "seed_block": row["seed_block"],
                        "stage": "sequence_novelty_qc",
                        "reason": reason,
                        "parent_identity": row["parent_identity"],
                        "designed_position_identity": row["designed_position_identity"],
                        "is_mock": False,
                    }
                )
        failure_frame = pd.DataFrame(failure_rows)
        atomic_write_text(
            report_dir / "failure_analysis.csv",
            failure_frame.to_csv(index=False),
        )
        statistics_config = config["statistics"]
        statistic_metrics = [
            "conditional_log_likelihood",
            "perplexity",
            "parent_identity",
            "designed_position_identity",
            "maximum_atlas_identity",
            "fixed_position_violations",
            "low_complexity_failure",
            "homopolymer_failure",
            "buried_core_recovery",
            "rna_interface_recovery",
            "rna_second_shell_recovery",
            "hepn_region_recovery",
            "model_agreement",
            "candidate_diversity",
        ]
        statistics_rows = paired_seed_statistics(
            candidate_rows,
            metrics=statistic_metrics,
            reference_method="unconstrained_esm_if1",
            bootstrap_replicates=int(statistics_config["bootstrap_replicates"]),
            confidence=float(statistics_config["confidence"]),
            seed=int(config["experiment"]["seed"]),
        )
        statistics_frame = pd.DataFrame(statistics_rows)
        atomic_write_text(
            report_dir / "matched_statistics.csv",
            statistics_frame.to_csv(index=False),
        )
        summary = {
            "schema_version": "1.0",
            "is_mock": False,
            "evidence_level_max": 2,
            "candidate_count": len(candidate_rows),
            "methods": len(METHODS),
            "seed_count": len(sampling["seed_blocks"]),
            "candidate_count_per_method": {
                method: sum(row["method"] == method for row in candidate_rows)
                for method in METHODS
            },
            "fixed_position_count": len(fixed_set),
            "free_position_count": len(free_set),
            "fixed_position_hash": fixed_hash,
            "free_position_hash": free_hash,
            "fixed_position_violations": sum(
                int(row["fixed_position_violations"]) for row in candidate_rows
            ),
            "identity_matching_passed": True,
            "mock_candidates": sum(bool(row["is_mock"]) for row in candidate_rows),
            "level1_pass_count": sum(
                bool(row["passes_level1_novelty"]) for row in candidate_rows
            ),
            "level2_scored_count": len(candidate_rows),
            "novelty_parameters": novelty_summary["thresholds"],
            "mmseqs_command": mmseqs_command,
            "runtime": {
                "device_requested": device,
                "python": sys.version,
                "platform": platform.platform(),
            },
            "claim_scope": (
                "Level 1 novelty only for passing rows and Level 2 ESM-IF1 "
                "compatibility; no functional validity claim"
            ),
        }
        atomic_write_text(
            report_dir / "summary.json",
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
        )
        gpu_jobs = [
            {
                "job_id": "vi-d-matched-baselines-gpu-extension",
                "task": "matched-baselines",
                "config": "configs/matched_baselines_gpu.yaml",
                "seed_blocks": list(range(20260731, 20260741)),
                "command": (
                    "bash scripts/launch_gpu_tmux.sh "
                    "configs/matched_baselines_gpu.yaml matched-baselines"
                ),
                "status": "not_run",
                "is_mock": False,
                "reason": "large GPU extension intentionally deferred",
                "expected_outputs": [
                    "reports/matched_baselines/methods_table.csv",
                    "reports/matched_baselines/candidates.jsonl",
                    "reports/matched_baselines/matched_statistics.csv",
                ],
            }
        ]
        _write_jsonl(report_dir / "gpu_hpc_job_manifest.jsonl", gpu_jobs)
        _render_figures(report_dir, candidate_frame, funnel)
        report = _report_text(
            methods, candidate_frame, statistics_frame, funnel, summary
        )
        atomic_write_text(report_dir / "report.md", report)
        escaped_summary = html.escape(json.dumps(summary, indent=2, sort_keys=True))
        report_html = f"""<!doctype html><html><head><meta charset="utf-8">
<title>Matched VI-D baseline report</title><style>body{{font-family:system-ui;
margin:2rem;max-width:1500px}}table{{border-collapse:collapse;font-size:11px}}
th,td{{border:1px solid #ccd3db;padding:4px}}th{{background:#eef2f6}}</style>
</head><body><h1>Matched VI-D baseline report</h1>
<p><strong>Real small CPU matrix; evidence maximum Level 2.</strong> Level 1
is sequence novelty only; Level 2 is inverse-folding compatibility. No
functional validity claim.</p><h2>Summary</h2><pre>{escaped_summary}</pre>
<h2>Methods and metrics</h2>{methods.to_html(index=False, escape=True)}
<h2>Candidate funnel</h2>{funnel.to_html(index=False, escape=True)}
<h2>Paired seed statistics</h2>{statistics_frame.to_html(index=False, escape=True)}
</body></html>\n"""
        atomic_write_text(report_dir / "report.html", report_html)
        outputs = sorted(path for path in report_dir.rglob("*") if path.is_file())
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=_file_entries(repo, outputs),
        )
        print(json.dumps({"run_dir": str(recorder.run_dir), **summary}, indent=2))
        return 0
    except Exception as exc:
        recorder.record_failure("matched-baselines", str(exc))
        recorder.finish(success=False)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
