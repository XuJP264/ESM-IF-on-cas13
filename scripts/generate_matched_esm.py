#!/usr/bin/env python
"""Generate the four real ESM-IF1 proposal conditions for the VI-D matrix."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import yaml

from cas13_if.backends.esm_if1 import TRACE_ALPHABET, EsmIf1ConstrainedBackend
from cas13_if.provenance import atomic_write_text
from cas13_if.schemas import SampleRequest
from cas13_if.structures.contacts import annotate_rna_contacts
from cas13_if.structures.parser import (
    parse_structure,
    protein_chain_sequence,
    residue_polymer_type,
)


def _arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--config", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    return parser.parse_args()


def _load(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"configuration is not a mapping: {path}")
    return value


def _repo_path(repo: Path, value: Any) -> Path:
    path = Path(str(value))
    return path if path.is_absolute() else repo / path


def _mapping_rows(path: Path) -> list[dict[str, str]]:
    """Read the mapping CSV without adding analysis-only dependencies.

    This runner executes inside the intentionally minimal legacy ESM-IF1
    environment, so tabular selection stays in the standard library.
    """

    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def _parent_and_fixed(
    *, repo: Path, structure: Path, chain: str, functional_manifest: Path
) -> tuple[str, dict[int, str], list[Any]]:
    atoms = parse_structure(structure)
    deposited, residue_keys = protein_chain_sequence(atoms, chain)
    manifest = _load(functional_manifest)
    entries = manifest["structures"]["6E9F"]["residues"]
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
    del repo
    return "".join(parent), fixed, residue_keys


def _rna_contact_positions(
    *, structure: Path, chain: str, rna_chains: set[str], residue_keys: list[Any]
) -> list[tuple[int, float]]:
    atoms = parse_structure(structure)
    selected = [
        atom
        for atom in atoms
        if (
            atom.residue.chain_id == chain
            and residue_polymer_type(atom.residue.residue_name) == "protein"
        )
        or (
            atom.residue.chain_id in rna_chains
            and residue_polymer_type(atom.residue.residue_name) == "rna"
        )
    ]
    by_key = {key: index for index, key in enumerate(residue_keys)}
    return sorted(
        (
            (by_key[item.protein_residue], float(item.minimum_rna_distance))
            for item in annotate_rna_contacts(selected)
            if item.direct_rna_contact and item.minimum_rna_distance is not None
        ),
        key=lambda value: (value[1], value[0]),
    )


def main() -> int:
    arguments = _arguments()
    repo = Path(__file__).resolve().parents[1]
    config = _load(arguments.config)
    inputs = config["inputs"]
    models = config["models"]
    sampling = config["sampling"]
    constraints = config["constraints"]
    execution = config["execution"]
    structure = _repo_path(repo, inputs["structure_pdb"])
    chain = str(inputs["chain_id"])
    parent, fixed, residue_keys = _parent_and_fixed(
        repo=repo,
        structure=structure,
        chain=chain,
        functional_manifest=_repo_path(repo, inputs["functional_manifest"]),
    )
    resolved = [
        row
        for row in _mapping_rows(_repo_path(repo, inputs["mapping_csv"]))
        if row["coordinate_index_0"].strip()
    ]
    eligible = [
        row
        for row in resolved
        if row["mapping_confidence"] == str(constraints["minimum_mapping_confidence"])
        and float(row["msa_coverage"]) >= float(constraints["minimum_msa_coverage"])
        and float(row["conservation"]) >= float(constraints["minimum_conservation"])
        and int(row["coordinate_index_0"]) not in fixed
    ]
    eligible.sort(
        key=lambda row: (
            -float(row["conservation"]),
            -float(row["msa_coverage"]),
            int(row["coordinate_index_0"]),
        )
    )
    conservation_rows = eligible[
        : int(constraints["maximum_conservation_biased_positions"])
    ]
    conservation_allowed = {
        int(row["coordinate_index_0"]): set(row["allowed_residues"].split(";"))
        for row in conservation_rows
    }
    rna_ranked = _rna_contact_positions(
        structure=structure,
        chain=chain,
        rna_chains={
            str(value)
            for value in [*inputs["crrna_chains"], *inputs["target_rna_chains"]]
        },
        residue_keys=residue_keys,
    )
    maximum_rna = int(constraints["maximum_rna_contact_biased_positions"])
    mapping_by_coordinate = {int(row["coordinate_index_0"]): row for row in resolved}
    rna_allowed: dict[int, set[str]] = {}
    for index, _ in rna_ranked:
        if index in fixed or len(rna_allowed) >= maximum_rna:
            continue
        row = mapping_by_coordinate[index]
        tokens = {parent[index], row["msa_consensus"]}
        rna_allowed[index] = {token for token in tokens if len(token) == 1}

    checkpoint = _repo_path(repo, models["esm_if1_checkpoint"])
    backend = EsmIf1ConstrainedBackend(
        checkpoint, device=str(execution.get("device", "cpu"))
    )
    backend.load()
    methods = {
        "unconstrained_esm_if1": {},
        "catalytic_only_fixed_esm_if1": {},
        "conservation_constrained_esm_if1": conservation_allowed,
        "conservation_rna_contact_esm_if1": {
            **conservation_allowed,
            **rna_allowed,
        },
    }
    rows: list[dict[str, Any]] = []
    temperature = float(sampling["temperatures"]["esm_if1"])
    for method, allowed in methods.items():
        for seed_block_value in sampling["seed_blocks"]:
            seed_block = int(seed_block_value)
            for proposal_index in range(int(sampling["proposals_per_seed"])):
                actual_seed = seed_block + proposal_index
                candidate = backend.sample(
                    SampleRequest(
                        scaffold_id="6E9F-A",
                        structure_path=str(structure),
                        parent_sequence=parent,
                        count=1,
                        temperature=temperature,
                        seed=actual_seed,
                        fixed_positions=fixed,
                        allowed_residues=allowed,
                        protein_chains=[chain],
                    )
                )[0]
                selected_confidence = [
                    trace.probabilities[TRACE_ALPHABET.index(trace.selected_token)]
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
                    "allowed_position_count": len(allowed),
                    "selected_token_probabilities": selected_confidence,
                }
                rows.append(
                    {
                        "method": method,
                        "seed_block": seed_block,
                        "proposal_index": proposal_index,
                        "candidate": payload,
                    }
                )
    atomic_write_text(
        arguments.output,
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
    )
    summary = {
        "is_mock": False,
        "proposal_count": len(rows),
        "methods": {method: len(allowed) for method, allowed in methods.items()},
        "conservation_biased_positions": sorted(conservation_allowed),
        "rna_contact_biased_positions": sorted(rna_allowed),
        "fixed_positions": fixed,
        "device": backend.metadata()["device"],
    }
    atomic_write_text(
        arguments.output.with_suffix(".summary.json"),
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(summary, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
