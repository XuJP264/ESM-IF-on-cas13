#!/usr/bin/env python
"""Run the bounded real Stage-0003A multi-scaffold generation smoke."""

from __future__ import annotations

import csv
import json
import math
import sys
from pathlib import Path
from typing import Any

import yaml

from cas13_if.backends.esm_if1 import TRACE_ALPHABET, EsmIf1ConstrainedBackend
from cas13_if.backends.mpnn import LigandMpnnBackend, ProteinMpnnBackend
from cas13_if.evaluation.matching import identity_matched_source_consensus
from cas13_if.novelty.metrics import (
    designed_position_identity,
    longest_homopolymer,
    low_complexity_windows,
    sequence_identity,
)
from cas13_if.provenance import RunRecorder, atomic_write_text, sha256_file
from cas13_if.schemas import SampleRequest, ScoreRequest


def _load(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"configuration root is not a mapping: {path}")
    return value


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _truth(value: Any) -> bool:
    return str(value).strip().lower() == "true"


def _integer(value: str) -> int:
    parsed = float(value)
    if not math.isfinite(parsed) or not parsed.is_integer():
        raise ValueError(f"expected integer-valued mapping field, received {value!r}")
    return int(parsed)


def _probabilities(candidate: Any) -> list[float]:
    if candidate.traces:
        return [
            float(trace.probabilities[TRACE_ALPHABET.index(trace.selected_token)])
            for trace in candidate.traces
        ]
    values = candidate.metadata.get("selected_token_probabilities")
    if not isinstance(values, list) or len(values) != len(candidate.sequence):
        raise RuntimeError("candidate lacks aligned selected-token probabilities")
    return [float(value) for value in values]


def _region_recovery(
    sequence: str,
    parent: str,
    mapping: list[dict[str, str]],
    field: str,
) -> float | None:
    positions = {
        int(row["full_scaffold_index_0"])
        for row in mapping
        if (
            (
                bool(str(row.get(field, "")).strip())
                and str(row.get(field, "")).strip().lower() != "nan"
                if field == "HEPN_annotation"
                else _truth(row.get(field, ""))
            )
            and row["coordinate_index_0"].strip()
        )
    }
    if not positions:
        return None
    return sum(sequence[index] == parent[index] for index in positions) / len(positions)


def _candidate_record(
    *,
    candidate: Any,
    method: str,
    scaffold_id: str,
    pdb_id: str,
    full_parent: str,
    coordinate_to_full: list[int],
    full_fixed: dict[int, str],
    mapping: list[dict[str, str]],
    score: Any,
    qc: dict[str, Any],
    source_metadata: dict[str, Any] | None = None,
) -> dict[str, Any]:
    full_tokens = list(full_parent)
    for coordinate_index, full_index in enumerate(coordinate_to_full):
        full_tokens[full_index] = candidate.sequence[coordinate_index]
    full_sequence = "".join(full_tokens)
    violations = sum(
        full_sequence[index] != token for index, token in full_fixed.items()
    )
    designable = set(coordinate_to_full).difference(full_fixed)
    low_complexity = low_complexity_windows(
        full_sequence,
        window=int(qc["low_complexity_window"]),
        maximum_single_residue_fraction=float(qc["low_complexity_maximum_fraction"]),
    )
    parent_low_complexity = low_complexity_windows(
        full_parent,
        window=int(qc["low_complexity_window"]),
        maximum_single_residue_fraction=float(qc["low_complexity_maximum_fraction"]),
    )
    newly_introduced_low_complexity = sorted(
        set(low_complexity).difference(parent_low_complexity)
    )
    return {
        "candidate_id": f"{method}-{scaffold_id}-{pdb_id}-seed{candidate.seed}",
        "scaffold_id": scaffold_id,
        "pdb_id": pdb_id,
        "method": method,
        "seed": int(candidate.seed),
        "sequence": full_sequence,
        "coordinate_sequence": candidate.sequence,
        "parent_sequence": full_parent,
        "coordinate_parent_sequence": candidate.parent_sequence,
        "coordinate_to_full_index_0": coordinate_to_full,
        "fixed_positions": full_fixed,
        "designable_positions": sorted(designable),
        "parent_identity": sequence_identity(full_sequence, full_parent),
        "designed_position_identity": designed_position_identity(
            full_sequence, full_parent, designable
        ),
        "fixed_position_violations": violations,
        "conditional_log_likelihood": score.conditional_log_likelihood,
        "mean_conditional_log_likelihood": score.metadata[
            "mean_conditional_log_likelihood"
        ],
        "perplexity": score.perplexity,
        "core_recovery": _region_recovery(
            full_sequence, full_parent, mapping, "protein_core"
        ),
        "rna_interface_recovery": _region_recovery(
            full_sequence, full_parent, mapping, "RNA_contact"
        ),
        "rna_second_shell_recovery": _region_recovery(
            full_sequence, full_parent, mapping, "RNA_second_shell"
        ),
        "hepn_recovery": _region_recovery(
            full_sequence, full_parent, mapping, "HEPN_annotation"
        ),
        "longest_homopolymer": longest_homopolymer(full_sequence),
        "homopolymer_failure": longest_homopolymer(full_sequence)
        > int(qc["maximum_homopolymer_length"]),
        "low_complexity_window_count": len(low_complexity),
        "parent_low_complexity_window_count": len(parent_low_complexity),
        "new_low_complexity_window_count": len(newly_introduced_low_complexity),
        "low_complexity_failure": bool(newly_introduced_low_complexity),
        "rna_atomic_context": bool(candidate.metadata.get("rna_atomic_context", False)),
        "source_model_metadata": source_metadata or candidate.metadata,
        "is_mock": False,
        "evidence_level": 2,
        "claim_scope": (
            "real inverse-folding compatibility smoke; not functional evidence"
        ),
    }


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    config_path = repo / "configs/stage_0003a_local_smoke.yaml"
    config = _load(config_path)
    output = repo / str(config["outputs"]["directory"])
    if output.exists():
        raise FileExistsError(f"refusing to overwrite canonical smoke output: {output}")
    recorder = RunRecorder(
        root=repo / "results/runs",
        experiment=str(config["experiment"]["name"]),
        resolved_config=config,
        command=[sys.executable, str(Path(__file__).resolve())],
        repo_root=repo,
        is_mock=False,
    )
    try:
        report_root = repo / str(config["inputs"]["report_root"])
        state_rows = {
            str(row["pdb_id"]): row for row in _rows(report_root / "states.csv")
        }
        esm = EsmIf1ConstrainedBackend(
            repo / str(config["models"]["esm_if1_checkpoint"]),
            device=str(config["execution"]["device"]),
        )
        esm.load()
        records: list[dict[str, Any]] = []
        failures: list[dict[str, Any]] = []
        method_counts: dict[str, int] = {}
        for representative in config["inputs"]["representatives"]:
            scaffold_id = str(representative["scaffold_id"])
            pdb_id = str(representative["pdb_id"]).upper()
            state = state_rows[pdb_id]
            full_parent = str(state["full_natural_sequence"])
            mapping = _rows(
                report_root / f"residue_mapping/{pdb_id.lower()}/mapping.csv"
            )
            coordinate_rows = sorted(
                (row for row in mapping if row["coordinate_index_0"].strip()),
                key=lambda row: _integer(row["coordinate_index_0"]),
            )
            coordinate_indices = [
                _integer(row["coordinate_index_0"]) for row in coordinate_rows
            ]
            if coordinate_indices != list(range(len(coordinate_rows))):
                raise ValueError(f"non-contiguous coordinate mapping for {pdb_id}")
            coordinate_to_full = [
                int(row["full_scaffold_index_0"]) for row in coordinate_rows
            ]
            coordinate_parent = "".join(
                str(row["full_scaffold_amino_acid"]) for row in coordinate_rows
            )
            mask_path = report_root / f"multistate_masks/{scaffold_id.lower()}.json"
            mask = json.loads(mask_path.read_text(encoding="utf-8"))
            full_fixed = {
                int(index): str(token) for index, token in mask["fixed_tokens"].items()
            }
            full_to_coordinate = {
                full_index: coordinate_index
                for coordinate_index, full_index in enumerate(coordinate_to_full)
            }
            unresolved_fixed = sorted(set(full_fixed).difference(full_to_coordinate))
            if unresolved_fixed:
                raise ValueError(
                    f"hard-fixed positions lack coordinates for {pdb_id}: "
                    f"{unresolved_fixed[:10]}"
                )
            coordinate_fixed = {
                full_to_coordinate[index]: token for index, token in full_fixed.items()
            }
            eligible = [
                row
                for row in coordinate_rows
                if _truth(row["mapping_gate_passed"])
                and float(row["msa_coverage"])
                >= float(config["sampling"]["minimum_msa_coverage"])
                and float(row["conservation"])
                >= float(config["sampling"]["minimum_conservation"])
                and _integer(row["coordinate_index_0"]) not in coordinate_fixed
            ]
            eligible.sort(
                key=lambda row: (
                    -float(row["conservation"]),
                    -float(row["msa_coverage"]),
                    _integer(row["coordinate_index_0"]),
                )
            )
            conservation_allowed = {
                _integer(row["coordinate_index_0"]): set(
                    str(row["allowed_residues"]).split(";")
                )
                for row in eligible[
                    : int(config["sampling"]["maximum_conservation_positions"])
                ]
            }
            rna_eligible = [
                row
                for row in coordinate_rows
                if _truth(row["mapping_gate_passed"])
                and _truth(row["RNA_second_shell"])
                and not _truth(row["RNA_contact"])
                and _integer(row["coordinate_index_0"]) not in coordinate_fixed
            ]
            rna_eligible.sort(
                key=lambda row: (
                    float(row["minimum_rna_distance"]),
                    _integer(row["coordinate_index_0"]),
                )
            )
            rna_allowed = {
                _integer(row["coordinate_index_0"]): {
                    str(row["full_scaffold_amino_acid"]),
                    str(row["msa_consensus"]),
                }
                for row in rna_eligible[
                    : int(config["sampling"]["maximum_rna_second_shell_positions"])
                ]
            }
            allowed_by_method = {
                "common_safety_mask_esm_if1": {},
                "conservation_esm_if1": conservation_allowed,
                "conservation_rna_esm_if1": {
                    **conservation_allowed,
                    **rna_allowed,
                },
            }
            rna_chains = [
                token
                for field in ("crrna_chains", "target_rna_chains")
                for token in str(state[field]).split(";")
                if token and token.lower() != "nan"
            ]
            if not rna_chains or int(float(state["rna_atom_count"])) <= 0:
                raise ValueError(f"representative state lacks RNA context: {pdb_id}")
            protein = ProteinMpnnBackend(
                upstream=repo / str(config["models"]["proteinmpnn_upstream"]),
                checkpoint=repo / str(config["models"]["proteinmpnn_checkpoint"]),
                python_executable=repo
                / str(config["environments"]["ligandmpnn_python"]),
                device=str(config["execution"]["device"]),
            )
            ligand = LigandMpnnBackend(
                upstream=repo / str(config["models"]["ligandmpnn_upstream"]),
                ligand_checkpoint=repo / str(config["models"]["ligandmpnn_checkpoint"]),
                protein_checkpoint=repo
                / str(config["models"]["ligandmpnn_protein_checkpoint"]),
                soluble_checkpoint=repo
                / str(config["models"]["ligandmpnn_soluble_checkpoint"]),
                python_executable=repo
                / str(config["environments"]["ligandmpnn_python"]),
                rna_context_chains=rna_chains,
                device=str(config["execution"]["device"]),
            )
            protein.load()
            ligand.load()
            structure_cif = repo / f"data/experimental_structures/{pdb_id.lower()}.cif"
            structure_pdb = repo / f"data/experimental_structures/{pdb_id.lower()}.pdb"
            for seed_value in config["execution"]["seeds"]:
                seed = int(seed_value)
                source_candidates: dict[str, Any] = {}
                for method, allowed in allowed_by_method.items():
                    candidate = esm.sample(
                        SampleRequest(
                            scaffold_id=f"{scaffold_id}-{pdb_id}",
                            structure_path=str(structure_cif),
                            parent_sequence=coordinate_parent,
                            count=1,
                            temperature=float(
                                config["sampling"]["temperature"]["esm_if1"]
                            ),
                            seed=seed,
                            fixed_positions=coordinate_fixed,
                            allowed_residues=allowed,
                            protein_chains=[str(state["protein_chain"])],
                        )
                    )[0]
                    source_candidates[method] = candidate
                for method, backend in (
                    ("proteinmpnn", protein),
                    ("ligandmpnn", ligand),
                ):
                    candidate = backend.sample(
                        SampleRequest(
                            scaffold_id=f"{scaffold_id}-{pdb_id}",
                            structure_path=str(structure_pdb),
                            parent_sequence=coordinate_parent,
                            count=1,
                            temperature=float(
                                config["sampling"]["temperature"][method]
                            ),
                            seed=seed,
                            fixed_positions=coordinate_fixed,
                            protein_chains=[str(state["protein_chain"])],
                        )
                    )[0]
                    source_candidates[method] = candidate
                if not bool(
                    source_candidates["ligandmpnn"].metadata.get(
                        "rna_atomic_context", False
                    )
                ):
                    raise RuntimeError(f"LigandMPNN RNA context gate failed: {pdb_id}")
                first = source_candidates["conservation_rna_esm_if1"]
                second = source_candidates["ligandmpnn"]
                consensus_sequence, consensus_metadata = (
                    identity_matched_source_consensus(
                        parent_sequence=coordinate_parent,
                        first_sequence=first.sequence,
                        second_sequence=second.sequence,
                        first_confidences=_probabilities(first),
                        second_confidences=_probabilities(second),
                        target_identity=float(
                            config["sampling"]["consensus_target_coordinate_identity"]
                        ),
                    )
                )
                consensus = first.model_copy(
                    update={
                        "candidate_id": (
                            f"esm_ligand_consensus-{scaffold_id}-{pdb_id}-seed{seed}"
                        ),
                        "backend": "esm_ligand_consensus",
                        "sequence": consensus_sequence,
                        "traces": [],
                        "metadata": {
                            **consensus_metadata,
                            "rna_atomic_context": True,
                            "source_candidate_ids": [
                                first.candidate_id,
                                second.candidate_id,
                            ],
                            "semantics": "derived_from_two_real_model_proposals",
                        },
                    }
                )
                source_candidates["esm_ligand_consensus"] = consensus
                for method, candidate in source_candidates.items():
                    if any(
                        candidate.sequence[index] != token
                        for index, token in coordinate_fixed.items()
                    ):
                        raise RuntimeError(
                            f"fixed-position violation in {method}/{pdb_id}/seed{seed}"
                        )
                    score = esm.score(
                        ScoreRequest(
                            scaffold_id=f"{scaffold_id}-{pdb_id}",
                            structure_path=str(structure_cif),
                            sequence=candidate.sequence,
                            protein_chains=[str(state["protein_chain"])],
                            seed=seed,
                        )
                    )
                    record = _candidate_record(
                        candidate=candidate,
                        method=method,
                        scaffold_id=scaffold_id,
                        pdb_id=pdb_id,
                        full_parent=full_parent,
                        coordinate_to_full=coordinate_to_full,
                        full_fixed=full_fixed,
                        mapping=mapping,
                        score=score,
                        qc=config["qc"],
                        source_metadata={
                            **candidate.metadata,
                            "coordinate_fixed_position_count": len(coordinate_fixed),
                            "full_fixed_position_count": len(full_fixed),
                            "conservation_allowed_position_count": len(
                                conservation_allowed
                            ),
                            "rna_allowed_position_count": len(rna_allowed),
                            "rna_context_chains": rna_chains,
                            "rna_atom_count": int(float(state["rna_atom_count"])),
                        },
                    )
                    records.append(record)
                    method_counts[method] = method_counts.get(method, 0) + 1
        if any(row["is_mock"] for row in records):
            raise RuntimeError("mock candidate entered real smoke output")
        if any(int(row["fixed_position_violations"]) for row in records):
            raise RuntimeError("fixed-position violations entered real smoke output")
        expected = (
            len(config["inputs"]["representatives"])
            * len(config["execution"]["seeds"])
            * 6
        )
        if len(records) != expected:
            raise RuntimeError(f"candidate count {len(records)} != expected {expected}")
        output.mkdir(parents=True, exist_ok=False)
        atomic_write_text(
            output / "candidates.jsonl",
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in records),
        )
        fields = [
            "candidate_id",
            "scaffold_id",
            "pdb_id",
            "method",
            "seed",
            "parent_identity",
            "designed_position_identity",
            "conditional_log_likelihood",
            "mean_conditional_log_likelihood",
            "perplexity",
            "fixed_position_violations",
            "core_recovery",
            "rna_interface_recovery",
            "rna_second_shell_recovery",
            "hepn_recovery",
            "rna_atomic_context",
            "homopolymer_failure",
            "low_complexity_failure",
            "is_mock",
            "evidence_level",
        ]
        from io import StringIO

        handle = StringIO()
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(records)
        atomic_write_text(output / "metrics.csv", handle.getvalue())
        summary = {
            "real_candidate_count": len(records),
            "scaffold_count": len({row["scaffold_id"] for row in records}),
            "method_counts": method_counts,
            "seed_count": len(config["execution"]["seeds"]),
            "fixed_position_violations": 0,
            "ligandmpnn_candidates_with_rna_context": sum(
                row["method"] == "ligandmpnn" and row["rna_atomic_context"]
                for row in records
            ),
            "mock_candidate_count": 0,
            "execution_failure_count": len(failures),
            "sequence_qc_failure_count": sum(
                bool(row["homopolymer_failure"] or row["low_complexity_failure"])
                for row in records
            ),
            "evidence_level_max": 2,
            "claim_scope": "bounded real execution smoke, not final statistics",
        }
        atomic_write_text(
            output / "failures.jsonl",
            "".join(json.dumps(row, sort_keys=True) + "\n" for row in failures),
        )
        atomic_write_text(
            output / "summary.json",
            json.dumps(summary, indent=2, sort_keys=True) + "\n",
        )
        outputs = [
            output / "candidates.jsonl",
            output / "metrics.csv",
            output / "failures.jsonl",
            output / "summary.json",
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
        recorder.record_failure("stage_0003a_local_smoke", str(error))
        recorder.finish(success=False)
        raise


if __name__ == "__main__":
    raise SystemExit(main())
