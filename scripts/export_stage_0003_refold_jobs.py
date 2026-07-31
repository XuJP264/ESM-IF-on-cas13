#!/usr/bin/env python
"""Export deterministic Stage-0003 Level-3 jobs without running predictors."""

from __future__ import annotations

import csv
import json
import sys
from io import StringIO
from pathlib import Path
from typing import Any

import yaml

from cas13_if.provenance import RunRecorder, atomic_write_text, sha256_file
from cas13_if.refold.stage3 import Stage3RefoldJob, export_stage3_jobs, make_stage3_job
from cas13_if.structures.parser import (
    group_residues,
    parse_structure,
    residue_polymer_type,
)


def _load(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"YAML root is not a mapping: {path}")
    return value


def _csv(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _jsonl(path: Path) -> list[dict[str, Any]]:
    return [
        json.loads(line)
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]


def _rna_sequence(structure: Path, chain: str) -> str:
    residues = group_residues(parse_structure(structure))
    tokens: list[str] = []
    for key in residues:
        if key.chain_id != chain or residue_polymer_type(key.residue_name) != "rna":
            continue
        token = key.residue_name.strip().upper()
        aliases = {"RA": "A", "RC": "C", "RG": "G", "RU": "U"}
        token = aliases.get(token, token)
        if token not in {"A", "C", "G", "U"}:
            raise ValueError(
                f"noncanonical RNA residue {key.residue_name!r} in {structure}/{chain}"
            )
        tokens.append(token)
    if not tokens:
        raise ValueError(f"no RNA sequence in {structure}/{chain}")
    return "".join(tokens)


def _project_pilot_to_full(
    row: dict[str, Any], mapping: list[dict[str, str]], full_parent: str
) -> str:
    coordinate_rows = sorted(
        (entry for entry in mapping if entry["coordinate_index_0"].strip()),
        key=lambda entry: int(float(entry["coordinate_index_0"])),
    )
    sequence = str(row["sequence"])
    if len(sequence) != len(coordinate_rows):
        raise ValueError("Stage-0002 pilot length does not match 6E9F mapping")
    output = list(full_parent)
    for token, entry in zip(sequence, coordinate_rows, strict=True):
        output[int(entry["full_scaffold_index_0"])] = token
    return "".join(output)


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    config_path = repo / "configs/stage_0003_refold.yaml"
    config = _load(config_path)
    task_root = repo / str(config["execution"]["task_root"])
    recorder = RunRecorder(
        root=repo / "results/runs",
        experiment="stage-0003-refold-job-export",
        resolved_config=config,
        command=[sys.executable, str(Path(__file__).resolve())],
        repo_root=repo,
        is_mock=False,
    )
    try:
        states = {
            row["pdb_id"]: row
            for row in _csv(repo / str(config["inputs"]["state_table"]))
        }
        representatives = {
            "EsCas13d": "6E9F",
            "UrCas13d": "6IV9",
            "DjCas13d": "9M33",
            "CasRx": "9M8Q",
        }
        contexts: dict[str, dict[str, Any]] = {}
        for scaffold_id, pdb_id in representatives.items():
            state = states[pdb_id]
            structure = repo / f"data/experimental_structures/{pdb_id.lower()}.cif"
            crrna_chain = str(state["crrna_chains"]).split(";")[0]
            target_value = str(state["target_rna_chains"])
            target_chain = (
                target_value.split(";")[0]
                if target_value and target_value.lower() != "nan"
                else None
            )
            contexts[scaffold_id] = {
                "pdb_id": pdb_id,
                "crrna": _rna_sequence(structure, crrna_chain),
                "target_rna": (
                    _rna_sequence(structure, target_chain) if target_chain else None
                ),
                "full_parent": state["full_natural_sequence"],
            }
        inventory: list[dict[str, Any]] = []
        mapping = _csv(repo / "reports/stage_0003a/residue_mapping/6e9f/mapping.csv")
        pilot_rows = _jsonl(repo / str(config["inputs"]["pilot_candidates"]))
        if len(pilot_rows) != 18:
            raise ValueError(
                f"expected 18 Stage-0002 pilot candidates, got {len(pilot_rows)}"
            )
        for row in pilot_rows:
            if bool(row["is_mock"]):
                raise ValueError("mock Stage-0002 candidate entered export")
            inventory.append(
                {
                    "candidate_id": f"pilot-{row['candidate_id']}",
                    "parent_scaffold": "EsCas13d",
                    "sequence": _project_pilot_to_full(
                        row, mapping, str(contexts["EsCas13d"]["full_parent"])
                    ),
                    "source": "stage_0002_pilot",
                    "source_candidate_id": row["candidate_id"],
                    "source_is_mock": False,
                }
            )
        local_rows = _jsonl(
            repo / str(config["inputs"]["local_multiscaffold_candidates"])
        )
        for row in local_rows:
            if bool(row["is_mock"]):
                raise ValueError("mock local multi-scaffold candidate entered export")
            inventory.append(
                {
                    "candidate_id": f"multiscaffold-{row['candidate_id']}",
                    "parent_scaffold": row["scaffold_id"],
                    "sequence": row["sequence"],
                    "source": "stage_0003a_local_real_smoke",
                    "source_candidate_id": row["candidate_id"],
                    "source_is_mock": False,
                }
            )
        for scaffold_id, context in contexts.items():
            inventory.append(
                {
                    "candidate_id": f"WT-{scaffold_id}",
                    "parent_scaffold": scaffold_id,
                    "sequence": context["full_parent"],
                    "source": "natural_parent_control",
                    "source_candidate_id": f"WT-{scaffold_id}",
                    "source_is_mock": False,
                }
            )
        candidate_ids = [str(row["candidate_id"]) for row in inventory]
        if len(candidate_ids) != len(set(candidate_ids)):
            raise ValueError("candidate inventory contains duplicate identifiers")
        jobs: list[Stage3RefoldJob] = []
        shards = int(config["execution"]["shard_count"])
        for candidate in inventory:
            context = contexts[str(candidate["parent_scaffold"])]
            for seed_value in config["execution"]["seeds"]:
                seed = int(seed_value)
                for state_name in ("monomer", "binary", "ternary"):
                    if state_name == "ternary" and context["target_rna"] is None:
                        continue
                    for backend in config["backends"][state_name]:
                        jobs.append(
                            make_stage3_job(
                                candidate_id=str(candidate["candidate_id"]),
                                parent_scaffold=str(candidate["parent_scaffold"]),
                                sequence=str(candidate["sequence"]),
                                crrna=(
                                    str(context["crrna"])
                                    if state_name != "monomer"
                                    else None
                                ),
                                target_rna=(
                                    str(context["target_rna"])
                                    if state_name == "ternary"
                                    else None
                                ),
                                state=state_name,  # type: ignore[arg-type]
                                seed=seed,
                                msa_policy=str(config["policies"]["msa"]),
                                template_policy=str(config["policies"]["template"]),
                                recycles=int(config["execution"]["recycles"]),
                                backend=str(backend),  # type: ignore[arg-type]
                                shards=shards,
                            )
                        )
        summary = export_stage3_jobs(jobs, task_root, shards=shards)
        inventory_handle = StringIO()
        writer = csv.DictWriter(
            inventory_handle,
            fieldnames=[
                "candidate_id",
                "parent_scaffold",
                "source",
                "source_candidate_id",
                "sequence_length",
                "sequence_sha256",
                "source_is_mock",
            ],
        )
        writer.writeheader()
        for row in inventory:
            sequence = str(row["sequence"])
            import hashlib

            writer.writerow(
                {
                    **{key: row[key] for key in writer.fieldnames if key in row},
                    "sequence_length": len(sequence),
                    "sequence_sha256": hashlib.sha256(
                        sequence.encode("ascii")
                    ).hexdigest(),
                }
            )
        inventory_path = task_root / "manifests/candidate_inventory.csv"
        atomic_write_text(inventory_path, inventory_handle.getvalue())
        summary.update(
            {
                "pilot_candidate_count": len(pilot_rows),
                "local_real_multiscaffold_candidate_count": len(local_rows),
                "wt_control_count": len(contexts),
                "template_policy": config["policies"]["template"],
                "large_prediction_execution": "not_run",
            }
        )
        summary_path = task_root / "manifests/summary.json"
        atomic_write_text(
            summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        files = sorted(path for path in task_root.rglob("*") if path.is_file())
        checksum_path = task_root / "manifests/SHA256SUMS"
        atomic_write_text(
            checksum_path,
            "".join(
                f"{sha256_file(path)}  {path.relative_to(task_root)}\n"
                for path in files
                if path != checksum_path
            ),
        )
        report_summary = repo / "reports/stage_0003a/gpu_jobs_summary.json"
        atomic_write_text(
            report_summary, json.dumps(summary, indent=2, sort_keys=True) + "\n"
        )
        recorder.finish(
            success=True,
            metrics=summary,
            outputs=[
                {
                    "path": str(report_summary.relative_to(repo)),
                    "sha256": sha256_file(report_summary),
                },
                {
                    "path": str(summary_path.relative_to(repo)),
                    "sha256": sha256_file(summary_path),
                },
                {
                    "path": str(checksum_path.relative_to(repo)),
                    "sha256": sha256_file(checksum_path),
                },
            ],
        )
        print(json.dumps({**summary, "run_dir": str(recorder.run_dir)}, indent=2))
        return 0
    except Exception as error:
        recorder.record_failure("export_stage_0003_refold_jobs", str(error))
        recorder.finish(success=False)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
