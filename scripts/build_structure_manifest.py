#!/usr/bin/env python
"""Build a strict, hash-addressed manifest for downloaded Cas13 structures."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import yaml
from Bio.Align import PairwiseAligner

from cas13_if.provenance import sha256_file
from cas13_if.structures.parser import (
    parse_structure,
    protein_chain_sequence,
    structure_qc,
)

STRUCTURES: dict[str, dict[str, Any]] = {
    "6E9F": {
        "subtype": "VI-D",
        "state": "crRNA_target_ternary",
        "selected_design_chain": "A",
        "crrna_chains": ["B"],
        "target_rna_chains": ["C"],
        "role": "primary_benchmark",
    },
    "5XWP": {
        "subtype": "VI-A",
        "state": "crRNA_target_ternary",
        "selected_design_chain": "A",
        "crrna_chains": ["C"],
        "target_rna_chains": ["D"],
        "symmetry_related_copies": {"protein": ["B"], "crrna": ["E"], "target": ["F"]},
        "role": "primary_benchmark",
    },
    "6E9E": {
        "subtype": "VI-D",
        "state": "crRNA_binary",
        "selected_design_chain": "A",
        "crrna_chains": ["B"],
        "target_rna_chains": [],
        "role": "matched_multistate_context",
    },
    "5XWY": {
        "subtype": "VI-A",
        "state": "crRNA_binary",
        "selected_design_chain": "A",
        "crrna_chains": ["B"],
        "target_rna_chains": [],
        "role": "matched_multistate_context",
    },
}

RELATED_QUERY = [
    {
        "pdb_id": "6IV8",
        "subtype": "VI-D",
        "state": "pre-crRNA_binary",
        "query_status": "identified_not_in_bootstrap_benchmark",
    },
    {
        "pdb_id": "6IV9",
        "subtype": "VI-D",
        "state": "crRNA_binary",
        "query_status": "identified_not_in_bootstrap_benchmark",
    },
    {
        "pdb_id": "6AAY",
        "subtype": "VI-B",
        "state": "crRNA_binary",
        "query_status": "identified_not_in_bootstrap_benchmark",
    },
    {
        "pdb_id": "7OS0",
        "subtype": "VI-A",
        "state": "crRNA_binary",
        "query_status": "identified_not_in_bootstrap_benchmark",
    },
    {
        "pdb_id": "7VTI",
        "subtype": "Cas13bt3",
        "state": "crRNA_binary",
        "query_status": "identified_not_in_bootstrap_benchmark",
    },
    {
        "pdb_id": "7VTN",
        "subtype": "Cas13bt3",
        "state": "crRNA_target_ternary",
        "query_status": "identified_not_in_bootstrap_benchmark",
    },
]


def _protein_entity(entity_files: list[Path]) -> dict[str, Any]:
    for path in entity_files:
        entity = json.loads(path.read_text(encoding="utf-8"))
        if entity["entity_poly"]["rcsb_entity_polymer_type"] == "Protein":
            return entity
    raise ValueError("entry has no protein polymer entity")


def _sequence_mapping(canonical: str, observed: str) -> dict[str, Any]:
    aligner = PairwiseAligner()
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(canonical, observed)[0]
    aligned_pairs = 0
    identical = 0
    for canonical_block, observed_block in zip(
        alignment.aligned[0], alignment.aligned[1], strict=True
    ):
        block_length = min(
            int(canonical_block[1] - canonical_block[0]),
            int(observed_block[1] - observed_block[0]),
        )
        aligned_pairs += block_length
        for offset in range(block_length):
            if (
                canonical[int(canonical_block[0]) + offset]
                == observed[int(observed_block[0]) + offset]
            ):
                identical += 1
    return {
        "seqres_length": len(canonical),
        "coordinate_length": len(observed),
        "unmodeled_or_unmapped_count": max(0, len(canonical) - len(observed)),
        "aligned_coordinate_count": aligned_pairs,
        "aligned_identity": identical / aligned_pairs if aligned_pairs else 0.0,
        "mapping_confidence": (
            "high" if aligned_pairs and identical / aligned_pairs >= 0.99 else "review"
        ),
    }


def _alternate_location_count(pdb_path: Path) -> int:
    count = 0
    with pdb_path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(("ATOM  ", "HETATM")):
                alternate = line[16:17].strip()
                if alternate not in {"", "A"}:
                    count += 1
    return count


def _record(repo: Path, pdb_id: str, declaration: dict[str, Any]) -> dict[str, Any]:
    lower = pdb_id.lower()
    root = repo / "data/experimental_structures"
    cif_path = root / f"{lower}.cif"
    pdb_path = root / f"{lower}.pdb"
    entry_path = root / f"{lower}.entry.json"
    entity_files = sorted(root.glob(f"{lower}.entity_*.json"))
    for path in (cif_path, pdb_path, entry_path):
        if not path.is_file():
            raise FileNotFoundError(f"missing structure asset: {path}")
    entry = json.loads(entry_path.read_text(encoding="utf-8"))
    protein_entity = _protein_entity(entity_files)
    atoms = parse_structure(cif_path)
    qc = structure_qc(atoms)
    design_chain = str(declaration["selected_design_chain"])
    observed_sequence, residue_keys = protein_chain_sequence(atoms, design_chain)
    canonical_sequence = (
        protein_entity["entity_poly"]["pdbx_seq_one_letter_code_can"]
        .replace("\n", "")
        .replace(" ", "")
    )
    entry_info = entry["rcsb_entry_info"]
    citation = next(
        (
            item
            for item in entry.get("citation", [])
            if item.get("rcsb_is_primary") == "Y" or item.get("id") == "primary"
        ),
        {},
    )
    resolution_values = entry_info.get("resolution_combined") or []
    breaks = [
        {
            "left": {
                "chain": left.chain_id,
                "residue_number": left.residue_number,
                "insertion_code": left.insertion_code,
            },
            "right": {
                "chain": right.chain_id,
                "residue_number": right.residue_number,
                "insertion_code": right.insertion_code,
            },
            "ca_distance_angstrom": distance,
        }
        for left, right, distance in qc.chain_breaks
    ]
    missing_backbone = [
        {
            "chain": key.chain_id,
            "residue_number": key.residue_number,
            "insertion_code": key.insertion_code,
            "missing_atoms": list(missing),
        }
        for key, missing in qc.missing_backbone
    ]
    record = {
        "pdb_id": pdb_id,
        "title": entry["struct"]["title"],
        "publication": {
            "title": citation.get("title"),
            "doi": citation.get("pdbx_database_id_DOI"),
            "year": citation.get("year"),
            "venue": citation.get("journal_abbrev"),
        },
        "experimental_method": entry["exptl"][0]["method"],
        "resolution_angstrom": resolution_values[0] if resolution_values else None,
        "protein_chains": list(qc.protein_chains),
        "rna_chains": list(qc.rna_chains),
        "target_rna_chains": declaration["target_rna_chains"],
        "crrna_chains": declaration["crrna_chains"],
        "cas13_subtype": declaration["subtype"],
        "selected_design_chain": design_chain,
        "selected_coordinate_sequence": observed_sequence,
        "selected_residue_range": {
            "first": {
                "number": residue_keys[0].residue_number,
                "insertion_code": residue_keys[0].insertion_code,
            },
            "last": {
                "number": residue_keys[-1].residue_number,
                "insertion_code": residue_keys[-1].insertion_code,
            },
        },
        "seqres_mapping": _sequence_mapping(canonical_sequence, observed_sequence),
        "uniprot_ids": protein_entity["rcsb_polymer_entity_container_identifiers"].get(
            "uniprot_ids", []
        ),
        "missing_backbone_atoms": missing_backbone,
        "chain_breaks": breaks,
        "modified_residues": [
            {
                "chain": key.chain_id,
                "residue_number": key.residue_number,
                "insertion_code": key.insertion_code,
                "residue_name": key.residue_name,
            }
            for key in qc.modified_residues
        ],
        "alternate_location_atoms_discarded": _alternate_location_count(pdb_path),
        "coordinates_are_finite": qc.coordinates_are_finite,
        "biological_assembly": 1,
        "assembly_count": entry_info.get("assembly_count"),
        "state": declaration["state"],
        "role": declaration["role"],
        "inclusion": True,
        "inclusion_reason": (
            "primary preregistered experimental complex"
            if declaration["role"] == "primary_benchmark"
            else "same-study matched structural state"
        ),
        "files": {
            "mmcif": {
                "path": str(cif_path.relative_to(repo)),
                "size_bytes": cif_path.stat().st_size,
                "sha256": sha256_file(cif_path),
            },
            "pdb": {
                "path": str(pdb_path.relative_to(repo)),
                "size_bytes": pdb_path.stat().st_size,
                "sha256": sha256_file(pdb_path),
            },
            "rcsb_entry_metadata": {
                "path": str(entry_path.relative_to(repo)),
                "size_bytes": entry_path.stat().st_size,
                "sha256": sha256_file(entry_path),
            },
        },
    }
    if "symmetry_related_copies" in declaration:
        record["symmetry_related_copies"] = declaration["symmetry_related_copies"]
    return record


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    records = [
        _record(repo, pdb_id, declaration) for pdb_id, declaration in STRUCTURES.items()
    ]
    manifest = {
        "schema_version": "1.0",
        "generated_at": "2026-07-31",
        "source": "RCSB PDB official download and Data API",
        "is_mock": False,
        "evidence_level": 0,
        "structures": records,
        "related_structure_query": RELATED_QUERY,
    }
    manifest_path = repo / "data/manifests/experimental_structures.yaml"
    manifest_path.write_text(
        yaml.safe_dump(manifest, sort_keys=False, width=100),
        encoding="utf-8",
    )
    funnel = {
        "queried_related_entries": len(records) + len(RELATED_QUERY),
        "downloaded_entries": len(records),
        "primary_benchmark_entries": sum(
            record["role"] == "primary_benchmark" for record in records
        ),
        "matched_state_entries": sum(
            record["role"] == "matched_multistate_context" for record in records
        ),
        "included_after_qc": sum(
            record["inclusion"]
            and record["coordinates_are_finite"]
            and not record["missing_backbone_atoms"]
            for record in records
        ),
        "is_mock": False,
        "evidence_level": 0,
    }
    funnel_path = repo / "data/manifests/experimental_structure_funnel.json"
    funnel_path.write_text(
        json.dumps(funnel, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    rows = []
    for record in records:
        mapping = record["seqres_mapping"]
        rows.append(
            "| {pdb_id} | {cas13_subtype} | {state} | {experimental_method} | "
            "{resolution} | {protein_chains} | {rna_chains} | {design_chain} | "
            "{coord_len}/{seqres_len} | {identity:.4f} | {breaks} | {status} |".format(
                pdb_id=record["pdb_id"],
                cas13_subtype=record["cas13_subtype"],
                state=record["state"],
                experimental_method=record["experimental_method"],
                resolution=record["resolution_angstrom"],
                protein_chains=",".join(record["protein_chains"]),
                rna_chains=",".join(record["rna_chains"]),
                design_chain=record["selected_design_chain"],
                coord_len=mapping["coordinate_length"],
                seqres_len=mapping["seqres_length"],
                identity=mapping["aligned_identity"],
                breaks=len(record["chain_breaks"]),
                status="included" if record["inclusion"] else "excluded",
            )
        )
    table_header = (
        "| PDB | Subtype | State | Method | Resolution (Å) | Protein chains | "
        "RNA chains | Design chain | Coordinates/SEQRES | Mapped identity | "
        "Chain breaks | Status |"
    )
    report = (
        f"""# Experimental Cas13 structure data card

Generated from locally hashed RCSB PDB/mmCIF files and RCSB Data API metadata.
This is real structural QC (Evidence Level 0), not evidence that a designed
sequence is an effective Cas13.

{table_header}
|---|---|---|---|---:|---|---|---|---:|---:|---:|---|
"""
        + "\n".join(rows)
        + """

The primary benchmark uses 6E9F and 5XWP. Their same-study binary states, 6E9E
and 5XWY, are retained for conformational comparison. RNA atoms are preserved
for contact annotation and LigandMPNN, but are never supplied as protein chains
to ESM-IF1. Coordinate gaps are recorded rather than imputed.

Related structures were queried and recorded in the manifest, but were not
silently added to the preregistered bootstrap benchmark.
"""
    )
    report_path = repo / "reports/experimental_structure_data_card.md"
    report_path.write_text(report, encoding="utf-8")
    print(json.dumps(funnel, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
