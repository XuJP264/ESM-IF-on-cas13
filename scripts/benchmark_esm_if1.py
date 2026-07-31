#!/usr/bin/env python
"""Run the genuine experimental-structure ESM-IF1 benchmark."""

from __future__ import annotations

import argparse
import json
import platform
import sys
import time
from dataclasses import asdict
from pathlib import Path
from typing import Any

import torch
import yaml

from cas13_if.backends.esm_if1 import EsmIf1Backend, EsmIf1ConstrainedBackend
from cas13_if.config import load_config
from cas13_if.evaluation.metrics import native_recovery
from cas13_if.provenance import atomic_write_text, sha256_file
from cas13_if.reporting.reports import render_run_report
from cas13_if.schemas import SampleRequest, ScoreRequest
from cas13_if.structures.contacts import annotate_rna_contacts
from cas13_if.structures.parser import (
    Atom,
    ResidueKey,
    parse_structure,
    protein_chain_sequence,
    residue_polymer_type,
)
from cas13_if.structures.sasa import relative_solvent_accessibility


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    return parser.parse_args()


def _load_yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"YAML root must be a mapping: {path}")
    return value


def _repo_path(repo: Path, value: Any, *, key: str) -> Path:
    if not isinstance(value, str) or not value:
        raise ValueError(f"{key} must be a non-empty path")
    path = Path(value)
    return path if path.is_absolute() else repo / path


def _selected_atoms(
    atoms: list[Atom], *, protein_chain: str, rna_chains: set[str]
) -> list[Atom]:
    return [
        atom
        for atom in atoms
        if (
            atom.residue.chain_id == protein_chain
            and residue_polymer_type(atom.residue.residue_name) == "protein"
        )
        or (
            atom.residue.chain_id in rna_chains
            and residue_polymer_type(atom.residue.residue_name) == "rna"
        )
    ]


def _functional_positions(
    functional_entry: dict[str, Any],
    residue_keys: list[ResidueKey],
    deposited_sequence: str,
) -> tuple[dict[int, str], list[dict[str, Any]]]:
    number_to_index = {
        (key.residue_number, key.insertion_code): index
        for index, key in enumerate(residue_keys)
    }
    fixed: dict[int, str] = {}
    rows: list[dict[str, Any]] = []
    residue_entries = functional_entry.get("residues")
    if not isinstance(residue_entries, list):
        raise ValueError("functional residue manifest lacks residue list")
    for raw in residue_entries:
        if not isinstance(raw, dict):
            raise ValueError("functional residue entry must be a mapping")
        number = int(raw["pdb_residue_number"])
        insertion_code = str(raw.get("insertion_code") or "")
        index = number_to_index.get((number, insertion_code))
        if index is None:
            raise ValueError(
                f"functional residue {number}{insertion_code} is unresolved"
            )
        deposited = deposited_sequence[index]
        declared_deposited = str(raw["deposited_amino_acid"]).upper()
        if deposited != declared_deposited:
            raise ValueError(
                f"deposited residue mismatch at {number}: "
                f"{deposited} != {declared_deposited}"
            )
        biological = str(raw["biological_amino_acid"]).upper()
        fixed[index] = biological
        rows.append(
            {
                **raw,
                "index_0": index,
                "coordinate_residue_name": residue_keys[index].residue_name,
                "mapping_validated": True,
            }
        )
    return fixed, rows


def _regions(
    *,
    atoms: list[Atom],
    structure_path: Path,
    protein_chain: str,
    cr_rna_chains: set[str],
    target_rna_chains: set[str],
    residue_keys: list[ResidueKey],
    hepn_positions: set[int],
    direct_cutoff: float,
    second_shell_cutoff: float,
    buried_rsa_threshold: float,
) -> tuple[dict[str, set[int]], list[dict[str, Any]]]:
    relevant_atoms = _selected_atoms(
        atoms,
        protein_chain=protein_chain,
        rna_chains=cr_rna_chains.union(target_rna_chains),
    )
    annotations = annotate_rna_contacts(
        relevant_atoms,
        direct_cutoff=direct_cutoff,
        second_shell_cutoff=second_shell_cutoff,
    )
    index_by_key = {key: index for index, key in enumerate(residue_keys)}
    rsa = relative_solvent_accessibility(
        structure_path,
        chain_id=protein_chain,
    )
    regions: dict[str, set[int]] = {
        "buried": set(),
        "surface": set(),
        "rna_interface": set(),
        "crrna_interface": set(),
        "target_rna_interface": set(),
        "rna_second_shell": set(),
        "non_rna_interface": set(range(len(residue_keys))),
        "hepn_catalytic_positions": set(hepn_positions),
    }
    rows: list[dict[str, Any]] = []
    for annotation in annotations:
        index = index_by_key.get(annotation.protein_residue)
        if index is None:
            continue
        relative_sasa = rsa.get(annotation.protein_residue)
        if relative_sasa is None:
            raise ValueError(f"SASA mapping missing for {annotation.protein_residue}")
        if relative_sasa < buried_rsa_threshold:
            regions["buried"].add(index)
        else:
            regions["surface"].add(index)
        if annotation.direct_rna_contact:
            regions["rna_interface"].add(index)
            regions["non_rna_interface"].discard(index)
        if set(annotation.contacted_rna_chains).intersection(cr_rna_chains):
            regions["crrna_interface"].add(index)
        if set(annotation.contacted_rna_chains).intersection(target_rna_chains):
            regions["target_rna_interface"].add(index)
        if annotation.second_shell:
            regions["rna_second_shell"].add(index)
        rows.append(
            {
                "index_0": index,
                "biological_index_1": index + 1,
                "pdb_residue_number": annotation.protein_residue.residue_number,
                "insertion_code": annotation.protein_residue.insertion_code,
                "residue_name": annotation.protein_residue.residue_name,
                "relative_sasa": relative_sasa,
                "burial": (
                    "buried" if relative_sasa < buried_rsa_threshold else "surface"
                ),
                "minimum_rna_distance": annotation.minimum_rna_distance,
                "direct_rna_contact": annotation.direct_rna_contact,
                "second_shell": annotation.second_shell,
                "contacted_rna_chains": annotation.contacted_rna_chains,
                "hepn_catalytic_position": index in hepn_positions,
            }
        )
    if len(rows) != len(residue_keys):
        raise ValueError(
            f"region annotation count {len(rows)} != sequence length "
            f"{len(residue_keys)}"
        )
    return regions, rows


def _region_score_means(
    per_residue_log_probabilities: list[float],
    regions: dict[str, set[int]],
) -> dict[str, float | None]:
    return {
        name: (
            sum(per_residue_log_probabilities[index] for index in positions)
            / len(positions)
            if positions
            else None
        )
        for name, positions in regions.items()
    }


def _jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    atomic_write_text(
        path,
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
    )


def main() -> int:
    arguments = _arguments()
    repo = Path(__file__).resolve().parents[1]
    config = load_config(arguments.config)
    output_dir = arguments.output_dir.resolve()
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite benchmark: {output_dir}")
    output_dir.mkdir(parents=True, exist_ok=False)

    model_config = config.get("model")
    structure_config = config.get("structures")
    sampling_config = config.get("sampling")
    if not all(
        isinstance(value, dict)
        for value in (model_config, structure_config, sampling_config)
    ):
        raise ValueError("model, structures, and sampling must be mappings")
    assert isinstance(model_config, dict)
    assert isinstance(structure_config, dict)
    assert isinstance(sampling_config, dict)
    checkpoint = _repo_path(
        repo, model_config.get("checkpoint"), key="model.checkpoint"
    )
    structure_manifest_path = _repo_path(
        repo,
        structure_config.get("manifest"),
        key="structures.manifest",
    )
    functional_manifest_path = _repo_path(
        repo,
        structure_config.get("functional_manifest"),
        key="structures.functional_manifest",
    )
    structure_manifest = _load_yaml(structure_manifest_path)
    functional_manifest = _load_yaml(functional_manifest_path)
    structure_entries = {
        str(entry["pdb_id"]).upper(): entry
        for entry in structure_manifest.get("structures", [])
        if isinstance(entry, dict)
    }
    functional_entries = functional_manifest.get("structures")
    if not isinstance(functional_entries, dict):
        raise ValueError("functional manifest structures must be a mapping")

    device = str(model_config.get("device", "auto"))
    backend = EsmIf1Backend(checkpoint, device=device)
    load_started = time.perf_counter()
    backend.load()
    load_seconds = time.perf_counter() - load_started
    constrained = EsmIf1ConstrainedBackend(checkpoint, device=device)
    constrained.__dict__.update(backend.__dict__)

    pdb_ids = structure_config.get("pdb_ids")
    temperatures = sampling_config.get("temperatures")
    methods = sampling_config.get("methods", ["unconstrained", "catalytic_only_fixed"])
    if (
        not isinstance(pdb_ids, list)
        or not isinstance(temperatures, list)
        or not isinstance(methods, list)
    ):
        raise ValueError("pdb_ids, temperatures, and methods must be lists")
    sample_count = int(sampling_config.get("samples_per_condition", 1))
    seed = int(config.get("experiment", {}).get("seed", 20260731))
    direct_cutoff = float(structure_config.get("rna_contact_cutoff_angstrom", 5))
    shell_cutoff = float(structure_config.get("second_shell_cutoff_angstrom", 8))
    buried_threshold = float(structure_config.get("buried_rsa_threshold", 0.2))

    score_rows: list[dict[str, Any]] = []
    position_rows: list[dict[str, Any]] = []
    structure_summaries: dict[str, Any] = {}
    candidate_partial = output_dir / "candidates.jsonl.part"
    candidate_handle = candidate_partial.open("x", encoding="utf-8")
    benchmark_started = time.perf_counter()
    for pdb_id_value in pdb_ids:
        pdb_id = str(pdb_id_value).upper()
        entry = structure_entries.get(pdb_id)
        functional = functional_entries.get(pdb_id)
        if not isinstance(entry, dict) or not isinstance(functional, dict):
            raise ValueError(f"manifest entry missing for {pdb_id}")
        chain = str(entry["selected_design_chain"])
        structure_path = _repo_path(
            repo, entry["files"]["mmcif"]["path"], key=f"{pdb_id}.mmcif"
        )
        atoms = parse_structure(structure_path)
        deposited_sequence, residue_keys = protein_chain_sequence(atoms, chain)
        fixed_catalytic, functional_rows = _functional_positions(
            functional,
            residue_keys,
            deposited_sequence,
        )
        biological_sequence_tokens = list(deposited_sequence)
        for index, token in fixed_catalytic.items():
            biological_sequence_tokens[index] = token
        biological_sequence = "".join(biological_sequence_tokens)
        cr_rna_chains = {str(value) for value in entry["crrna_chains"]}
        target_rna_chains = {str(value) for value in entry["target_rna_chains"]}
        regions, annotations = _regions(
            atoms=atoms,
            structure_path=structure_path,
            protein_chain=chain,
            cr_rna_chains=cr_rna_chains,
            target_rna_chains=target_rna_chains,
            residue_keys=residue_keys,
            hepn_positions=set(fixed_catalytic),
            direct_cutoff=direct_cutoff,
            second_shell_cutoff=shell_cutoff,
            buried_rsa_threshold=buried_threshold,
        )
        _jsonl(
            output_dir / f"{pdb_id.lower()}_position_annotations.jsonl",
            annotations,
        )
        atomic_write_text(
            output_dir / f"{pdb_id.lower()}_functional_mapping.json",
            json.dumps(functional_rows, indent=2, sort_keys=True) + "\n",
        )

        for sequence_label, sequence in (
            ("deposited_construct", deposited_sequence),
            ("biological_catalytic_restored", biological_sequence),
        ):
            score = backend.score(
                ScoreRequest(
                    scaffold_id=f"{pdb_id}-{chain}",
                    structure_path=str(structure_path),
                    sequence=sequence,
                    protein_chains=[chain],
                )
            )
            score_rows.append(
                {
                    "pdb_id": pdb_id,
                    "sequence_label": sequence_label,
                    **score.model_dump(mode="json"),
                    "region_mean_log_probabilities": _region_score_means(
                        score.per_residue_log_probabilities, regions
                    ),
                }
            )
            for index, log_probability in enumerate(
                score.per_residue_log_probabilities
            ):
                position_rows.append(
                    {
                        "pdb_id": pdb_id,
                        "sequence_label": sequence_label,
                        "index_0": index,
                        "pdb_residue_number": residue_keys[index].residue_number,
                        "insertion_code": residue_keys[index].insertion_code,
                        "amino_acid": sequence[index],
                        "log_probability": log_probability,
                        "probability": float(torch.exp(torch.tensor(log_probability))),
                        "is_hepn": index in regions["hepn_catalytic_positions"],
                        "is_rna_interface": index in regions["rna_interface"],
                    }
                )

        method_summaries: dict[str, Any] = {}
        for method in methods:
            method_name = str(method)
            if method_name == "unconstrained":
                method_backend = backend
                fixed_positions: dict[int, str] = {}
            elif method_name == "catalytic_only_fixed":
                method_backend = constrained
                fixed_positions = dict(fixed_catalytic)
            elif method_name == "catalytic_rna_fixed":
                method_backend = constrained
                fixed_positions = {
                    index: biological_sequence[index]
                    for index in regions["rna_interface"].union(fixed_catalytic)
                }
            else:
                raise ValueError(f"unknown benchmark method: {method_name}")
            method_recoveries: list[float] = []
            method_violations = 0
            for temperature_value in temperatures:
                temperature = float(temperature_value)
                candidates = method_backend.sample(
                    SampleRequest(
                        scaffold_id=f"{pdb_id}-{chain}",
                        structure_path=str(structure_path),
                        parent_sequence=biological_sequence,
                        count=sample_count,
                        temperature=temperature,
                        seed=seed,
                        fixed_positions=fixed_positions,
                        protein_chains=[chain],
                    )
                )
                for candidate in candidates:
                    recovery = native_recovery(
                        candidate.sequence,
                        biological_sequence,
                        designed_positions=set(
                            range(len(biological_sequence))
                        ).difference(fixed_positions),
                        fixed_positions=fixed_positions,
                        regions=regions,
                    )
                    method_recoveries.append(recovery.overall)
                    method_violations += recovery.fixed_position_violations
                    candidate_score = backend.score(
                        ScoreRequest(
                            scaffold_id=f"{pdb_id}-{chain}",
                            structure_path=str(structure_path),
                            sequence=candidate.sequence,
                            protein_chains=[chain],
                        )
                    )
                    candidate_payload = candidate.model_dump(mode="json")
                    candidate_payload["candidate_id"] = (
                        f"{pdb_id.lower()}-{method_name}-t{temperature:g}-"
                        f"{candidate_payload['candidate_id']}"
                    )
                    candidate_handle.write(
                        json.dumps(
                            {
                                "pdb_id": pdb_id,
                                "method": method_name,
                                "candidate": candidate_payload,
                                "recovery": asdict(recovery),
                                "score": candidate_score.model_dump(mode="json"),
                            },
                            sort_keys=True,
                        )
                        + "\n"
                    )
                    candidate_handle.flush()
            method_summaries[method_name] = {
                "candidate_count": len(method_recoveries),
                "mean_overall_recovery": (
                    sum(method_recoveries) / len(method_recoveries)
                    if method_recoveries
                    else None
                ),
                "fixed_position_violations": method_violations,
                "fixed_position_count": len(fixed_positions),
            }
        structure_summaries[pdb_id] = {
            "subtype": entry["cas13_subtype"],
            "chain": chain,
            "coordinate_length": len(deposited_sequence),
            "construct_differs_from_biological_at": sorted(
                index
                for index in fixed_catalytic
                if deposited_sequence[index] != biological_sequence[index]
            ),
            "region_sizes": {
                name: len(positions) for name, positions in regions.items()
            },
            "methods": method_summaries,
            "rna_chains_used_for_annotation_only": sorted(
                cr_rna_chains.union(target_rna_chains)
            ),
            "rna_passed_to_esm_if1": False,
        }
    elapsed = time.perf_counter() - benchmark_started
    candidate_handle.close()
    candidate_partial.replace(output_dir / "candidates.jsonl")
    _jsonl(output_dir / "scores.jsonl", score_rows)
    _jsonl(output_dir / "per_position_scores.jsonl", position_rows)
    summary = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 2,
        "claim_scope": (
            "Level 2 inverse-folding compatibility only; no candidate is a "
            "validated or effective Cas13"
        ),
        "checkpoint": {
            "path": str(checkpoint.relative_to(repo)),
            "sha256": sha256_file(checkpoint),
            "size_bytes": checkpoint.stat().st_size,
        },
        "runtime": {
            "device": backend.metadata()["device"],
            "torch": torch.__version__,
            "cuda_runtime": torch.version.cuda,
            "cuda_available": torch.cuda.is_available(),
            "gpu": torch.cuda.get_device_name(0) if torch.cuda.is_available() else None,
            "python": sys.version,
            "platform": platform.platform(),
            "model_load_seconds": load_seconds,
            "benchmark_seconds": elapsed,
        },
        "sampling": {
            "temperatures": [float(value) for value in temperatures],
            "samples_per_condition": sample_count,
            "seed": seed,
            "methods": methods,
            "causal_constraint_semantics": True,
            "future_fixed_tokens_visible": False,
        },
        "structures": structure_summaries,
    }
    atomic_write_text(
        output_dir / "summary.json",
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    render_run_report(
        title="Experimental Cas13 ESM-IF1 benchmark",
        evidence_level=2,
        is_mock=False,
        metrics=summary,
        failures=[],
        markdown_path=output_dir / "benchmark.md",
        html_path=output_dir / "benchmark.html",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
