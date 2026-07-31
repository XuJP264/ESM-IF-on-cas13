"""Build the Stage-0003A experimental Cas13d scaffold/state atlas."""

from __future__ import annotations

import html
import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from cas13_if.alignments.scaffold_mapping import (
    coordinate_index_mapping,
    global_index_mapping,
)
from cas13_if.provenance import atomic_write_text, sha256_file
from cas13_if.structures.parser import (
    group_residues,
    parse_structure,
    protein_chain_sequence,
    residue_polymer_type,
    structure_qc,
)


@dataclass(frozen=True)
class ScaffoldSource:
    scaffold_id: str
    sequence: str
    confidence: str
    strategy: str
    accession: str | None
    note: str


def _yaml(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"configuration root is not a mapping: {path}")
    return value


def _entity(path: Path) -> dict[str, Any]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError(f"entity JSON root is not a mapping: {path}")
    if value.get("entity_poly", {}).get("rcsb_entity_polymer_type") != "Protein":
        raise ValueError(f"declared entity is not protein: {path}")
    return value


def _entity_sequence(entity: dict[str, Any]) -> str:
    return "".join(entity["entity_poly"]["pdbx_seq_one_letter_code_can"].split())


def _fasta_sequence(path: Path) -> str:
    lines = path.read_text(encoding="utf-8").splitlines()
    sequence = "".join(line.strip() for line in lines if not line.startswith(">"))
    if not sequence:
        raise ValueError(f"FASTA contains no sequence: {path}")
    return sequence.upper()


def _restore(sequence: str, restorations: dict[Any, Any]) -> str:
    restored = list(sequence)
    for raw_position, raw_token in restorations.items():
        position = int(raw_position)
        if position < 1 or position > len(restored):
            raise ValueError(f"restoration position {position} is outside sequence")
        restored[position - 1] = str(raw_token).upper()
    return "".join(restored)


def _natural_source(repo: Path, declaration: dict[str, Any]) -> ScaffoldSource:
    first_state = declaration["states"][0]
    entity_path = repo / (
        "data/experimental_structures/"
        f"{str(first_state['pdb_id']).lower()}.entity_{first_state['protein_entity']}.json"
    )
    deposited = _entity_sequence(_entity(entity_path))
    strategy = str(declaration["natural_sequence_strategy"])
    if strategy == "uniprot_fasta":
        sequence = _fasta_sequence(repo / str(declaration["natural_sequence_source"]))
    elif strategy == "rcsb_reference_restore":
        reference_length = int(declaration["natural_reference_length"])
        sequence = _restore(
            deposited[:reference_length], first_state.get("catalytic_restorations", {})
        )
    elif strategy == "deposited_primary_paper_restore":
        mutation = (
            _entity(entity_path).get("rcsb_polymer_entity", {}).get("pdbx_mutation")
        )
        if mutation not in {None, "", "?"}:
            raise ValueError(
                f"{declaration['scaffold_id']} declares no mutation but "
                f"RCSB says {mutation}"
            )
        sequence = _restore(deposited, declaration.get("natural_restorations", {}))
    else:
        raise ValueError(f"unknown natural sequence strategy: {strategy}")
    return ScaffoldSource(
        scaffold_id=str(declaration["scaffold_id"]),
        sequence=sequence,
        confidence=str(declaration["natural_sequence_confidence"]),
        strategy=strategy,
        accession=(
            str(declaration["natural_sequence_accession"])
            if declaration.get("natural_sequence_accession")
            else None
        ),
        note=str(declaration.get("natural_sequence_note", "")),
    )


def aligned_identity(left: str, right: str) -> float:
    """Return identity over aligned residue pairs from a deterministic global map."""
    mapping = global_index_mapping(left, right)
    pairs = [
        (index, query_index)
        for index, query_index in enumerate(mapping.reference_to_query)
        if query_index is not None
    ]
    if not pairs:
        raise ValueError("sequences have no aligned residue pairs")
    return sum(left[i] == right[j] for i, j in pairs) / len(pairs)


def _alternate_locations(path: Path) -> int:
    count = 0
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(("ATOM  ", "HETATM")) and line[16:17].strip() not in {
            "",
            "A",
        }:
            count += 1
    return count


def _citation(entry: dict[str, Any]) -> dict[str, Any]:
    primary: dict[str, Any] = next(
        (
            item
            for item in entry.get("citation", [])
            if item.get("rcsb_is_primary") == "Y" or item.get("id") == "primary"
        ),
        {},
    )
    return {
        "title": primary.get("title"),
        "doi": primary.get("pdbx_database_id_DOI"),
        "year": primary.get("year"),
        "venue": primary.get("journal_abbrev"),
    }


def _state_record(
    repo: Path,
    scaffold: dict[str, Any],
    source: ScaffoldSource,
    state: dict[str, Any],
    sources: dict[str, Any],
) -> dict[str, Any]:
    pdb_id = str(state["pdb_id"]).upper()
    lower = pdb_id.lower()
    cif_path = repo / f"data/experimental_structures/{lower}.cif"
    pdb_path = repo / f"data/experimental_structures/{lower}.pdb"
    entry_path = repo / f"data/experimental_structures/{lower}.entry.json"
    entity_path = repo / (
        f"data/experimental_structures/{lower}.entity_{state['protein_entity']}.json"
    )
    for path in (cif_path, pdb_path, entry_path, entity_path):
        if not path.is_file():
            raise FileNotFoundError(f"missing structure source: {path}")
    entry = json.loads(entry_path.read_text(encoding="utf-8"))
    entity = _entity(entity_path)
    deposited_full = _entity_sequence(entity)
    atoms = parse_structure(cif_path)
    qc = structure_qc(atoms)
    chain = str(state["protein_chain"])
    coordinate_sequence, residue_keys = protein_chain_sequence(atoms, chain)
    natural_to_coordinate, coordinate_mapping_strategy = coordinate_index_mapping(
        source.sequence, coordinate_sequence, residue_keys
    )
    mapped = [
        (index, coordinate)
        for index, coordinate in enumerate(natural_to_coordinate.reference_to_query)
        if coordinate is not None
    ]
    substitutions = [
        {
            "natural_index_1": index + 1,
            "natural": source.sequence[index],
            "coordinate": coordinate_sequence[coordinate],
        }
        for index, coordinate in mapped
        if source.sequence[index] != coordinate_sequence[coordinate]
    ]
    unresolved = [
        index + 1
        for index, coordinate in enumerate(natural_to_coordinate.reference_to_query)
        if coordinate is None
    ]
    selected_rna = set(state.get("crrna_chains", [])).union(
        state.get("target_rna_chains", [])
    )
    rna_atom_count = sum(
        atom.residue.chain_id in selected_rna
        and residue_polymer_type(atom.residue.residue_name) == "rna"
        for atom in atoms
    )
    residues = group_residues(atoms)
    missing_backbone = [
        {
            "chain": key.chain_id,
            "residue_number": key.residue_number,
            "insertion_code": key.insertion_code,
            "missing": sorted({"N", "CA", "C"}.difference({a.name for a in values})),
        }
        for key, values in residues.items()
        if key.chain_id == chain
        and residue_polymer_type(key.residue_name) == "protein"
        and not {"N", "CA", "C"}.issubset({a.name for a in values})
    ]
    resolution = entry.get("rcsb_entry_info", {}).get("resolution_combined") or []
    declared_mutation = entity.get("rcsb_polymer_entity", {}).get("pdbx_mutation")
    return {
        "pdb_id": pdb_id,
        "scaffold_id": source.scaffold_id,
        "subtype": "VI-D",
        "state": str(state["state"]),
        "protein_chain": chain,
        "crrna_chains": ";".join(state.get("crrna_chains", [])),
        "target_rna_chains": ";".join(state.get("target_rna_chains", [])),
        "resolution_angstrom": resolution[0] if resolution else None,
        "experimental_method": entry.get("exptl", [{}])[0].get("method"),
        "publication": json.dumps(_citation(entry), sort_keys=True),
        "full_natural_sequence": source.sequence,
        "full_natural_length": len(source.sequence),
        "natural_sequence_confidence": source.confidence,
        "deposited_full_sequence": deposited_full,
        "deposited_full_length": len(deposited_full),
        "declared_construct_mutation": declared_mutation,
        "paper_declared_restorations": json.dumps(
            state.get("catalytic_restorations", {}), sort_keys=True
        ),
        "coordinate_sequence": coordinate_sequence,
        "coordinate_length": len(coordinate_sequence),
        "coordinate_mapping_strategy": coordinate_mapping_strategy,
        "coordinate_residue_first": residue_keys[0].residue_number,
        "coordinate_residue_last": residue_keys[-1].residue_number,
        "missing_natural_positions": json.dumps(unresolved),
        "natural_coordinate_substitutions": json.dumps(substitutions, sort_keys=True),
        "missing_backbone": json.dumps(missing_backbone, sort_keys=True),
        "chain_break_count": len(qc.chain_breaks),
        "alternate_location_atoms_discarded": _alternate_locations(pdb_path),
        "modified_residue_count": len(qc.modified_residues),
        "rna_atom_count": rna_atom_count,
        "coordinates_are_finite": qc.coordinates_are_finite,
        "experimental_evidence": "experimental_structure",
        "inclusion": True,
        "inclusion_reason": (
            "preregistered Cas13d scaffold-state with official coordinates"
        ),
        "download_url": str(sources["rcsb_files"]).format(
            pdb_id=pdb_id, extension="cif"
        ),
        "accessed_at": "2026-08-01",
        "license": str(sources["archive_license"]),
        "cif_sha256": sha256_file(cif_path),
        "pdb_sha256": sha256_file(pdb_path),
        "entry_sha256": sha256_file(entry_path),
        "entity_sha256": sha256_file(entity_path),
        "is_mock": False,
        "evidence_level": 0,
    }


def build_structure_atlas(
    *, repo: Path, config_path: Path, output_dir: Path
) -> dict[str, Any]:
    """Build real scaffold/state CSVs and an auditable HTML atlas."""
    config = _yaml(config_path)
    scaffolds = config.get("scaffolds")
    if not isinstance(scaffolds, list) or len(scaffolds) < 3:
        raise ValueError("at least three scaffold declarations are required")
    output_dir.mkdir(parents=True, exist_ok=True)
    sources = [_natural_source(repo, item) for item in scaffolds]
    source_by_id = {item.scaffold_id: item for item in sources}
    identities = {
        left.scaffold_id: {
            right.scaffold_id: aligned_identity(left.sequence, right.sequence)
            for right in sources
        }
        for left in sources
    }
    scaffold_rows = [
        {
            "scaffold_id": source.scaffold_id,
            "aliases": ";".join(
                next(
                    item.get("aliases", [])
                    for item in scaffolds
                    if item["scaffold_id"] == source.scaffold_id
                )
            ),
            "organism": next(
                item["organism"]
                for item in scaffolds
                if item["scaffold_id"] == source.scaffold_id
            ),
            "subtype": "VI-D",
            "full_natural_sequence": source.sequence,
            "full_natural_length": len(source.sequence),
            "natural_sequence_strategy": source.strategy,
            "natural_sequence_confidence": source.confidence,
            "natural_sequence_accession": source.accession,
            "natural_sequence_note": source.note,
            "state_count": len(
                next(
                    item["states"]
                    for item in scaffolds
                    if item["scaffold_id"] == source.scaffold_id
                )
            ),
            "parent_identities": json.dumps(
                identities[source.scaffold_id], sort_keys=True
            ),
            "nearest_other_parent_identity": max(
                value
                for key, value in identities[source.scaffold_id].items()
                if key != source.scaffold_id
            ),
            "is_mock": False,
            "evidence_level": 0,
        }
        for source in sources
    ]
    state_rows: list[dict[str, Any]] = []
    for scaffold in scaffolds:
        source = source_by_id[str(scaffold["scaffold_id"])]
        for state in scaffold["states"]:
            state_rows.append(
                _state_record(repo, scaffold, source, state, config["sources"])
            )
    if len(state_rows) < 6:
        raise ValueError("at least six real scaffold-state units are required")
    scaffold_frame = pd.DataFrame.from_records(scaffold_rows)
    state_frame = pd.DataFrame.from_records(state_rows)
    atomic_write_text(output_dir / "scaffolds.csv", scaffold_frame.to_csv(index=False))
    atomic_write_text(output_dir / "states.csv", state_frame.to_csv(index=False))
    summary = {
        "schema_version": "1.0",
        "generated_at": "2026-08-01",
        "scaffolds": len(scaffold_rows),
        "state_units": len(state_rows),
        "states": state_frame["state"].value_counts().sort_index().to_dict(),
        "all_coordinates_finite": bool(state_frame["coordinates_are_finite"].all()),
        "all_rna_states_have_atoms": bool(
            (state_frame.loc[state_frame["state"] != "apo", "rna_atom_count"] > 0).all()
        ),
        "config_sha256": sha256_file(config_path),
        "is_mock": False,
        "evidence_level": 0,
    }
    atomic_write_text(
        output_dir / "structure_atlas_summary.json",
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
    )
    html_text = (
        """<!doctype html><html><head><meta charset="utf-8">
<title>Stage 0003A Cas13d structure atlas</title><style>
body{font-family:system-ui;margin:2rem;color:#17202a}table{border-collapse:collapse;
font-size:12px}th,td{border:1px solid #ccd3db;padding:4px 6px;max-width:28rem;
overflow-wrap:anywhere}th{background:#eef2f6;position:sticky;top:0}</style></head><body>
<h1>Stage 0003A experimental Cas13d structure atlas</h1>
<p><strong>Evidence Level 0.</strong> Real RCSB coordinates and provenance; no
functional or wet-lab claim. Natural-parent and deposited construct sequences
are separate fields.</p><h2>Summary</h2><pre>"""
        + html.escape(json.dumps(summary, indent=2, sort_keys=True))
        + "</pre><h2>Scaffolds</h2>"
        + scaffold_frame.drop(columns=["full_natural_sequence"]).to_html(
            index=False, escape=True
        )
        + "<h2>States</h2>"
        + state_frame.drop(
            columns=[
                "full_natural_sequence",
                "deposited_full_sequence",
                "coordinate_sequence",
            ]
        ).to_html(index=False, escape=True)
        + "</body></html>\n"
    )
    atomic_write_text(output_dir / "structure_atlas.html", html_text)
    return summary
