#!/usr/bin/env python
"""Build strict four-layer mappings and real structural annotations."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import yaml

from cas13_if.alignments.scaffold_mapping import build_scaffold_mapping
from cas13_if.evaluation.regions import build_structure_regions
from cas13_if.provenance import atomic_write_text
from cas13_if.structures.parser import parse_structure, structure_qc


def _config(path: Path) -> dict[str, Any]:
    value = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(value, dict):
        raise ValueError("Stage-0003A structure config must be a mapping")
    return value


def _split_chains(value: Any) -> set[str]:
    return {item for item in str(value).split(";") if item and item != "nan"}


def main() -> int:
    repo = Path(__file__).resolve().parents[1]
    config_path = repo / "configs/stage_0003a_structures.yaml"
    config = _config(config_path)
    root = repo / "reports/stage_0003a"
    scaffold_rows = pd.read_csv(root / "scaffolds.csv").set_index("scaffold_id")
    mapping_root = root / "residue_mapping"
    review_root = root / "manual_review"
    mapping_root.mkdir(parents=True, exist_ok=True)
    review_root.mkdir(parents=True, exist_ok=True)
    msa = config["msa"]
    summaries: list[dict[str, Any]] = []
    for scaffold in config["scaffolds"]:
        scaffold_id = str(scaffold["scaffold_id"])
        full_sequence = str(scaffold_rows.loc[scaffold_id, "full_natural_sequence"])
        hepn_positions = [int(item) for item in scaffold["hepn_residues"]]
        hepn_labels = {
            hepn_positions[0]: "HEPN1_R",
            hepn_positions[1]: "HEPN1_H",
            hepn_positions[2]: "HEPN2_R",
            hepn_positions[3]: "HEPN2_H",
        }
        for state in scaffold["states"]:
            pdb_id = str(state["pdb_id"]).upper()
            lower = pdb_id.lower()
            output = mapping_root / lower
            structure = repo / f"data/experimental_structures/{lower}.cif"
            entity = repo / (
                f"data/experimental_structures/{lower}.entity_"
                f"{state['protein_entity']}.json"
            )
            accepted = {int(item) for item in state.get("catalytic_restorations", {})}
            summary = build_scaffold_mapping(
                structure_path=structure,
                entity_path=entity,
                msa_path=repo / str(msa["path"]),
                conservation_path=repo / str(msa["conservation"]),
                output_dir=output,
                chain_id=str(state["protein_chain"]),
                subtype="VI-D",
                mafft_executable=str(repo / str(msa["mafft_executable"])),
                threads=int(msa["threads"]),
                minimum_conservation_coverage=float(
                    msa["minimum_conservation_coverage"]
                ),
                pdb_id=pdb_id,
                state=str(state["state"]),
                full_sequence=full_sequence,
                query_identifier=f"cas13_if__{lower}__{scaffold_id.lower()}__full",
                accepted_substitution_positions_1=accepted,
            )
            mapping = pd.read_csv(output / "mapping.csv")
            coordinate_hepn = set(
                mapping.loc[
                    mapping["biological_index_1"].isin(hepn_positions)
                    & mapping["coordinate_index_0"].notna(),
                    "coordinate_index_0",
                ].astype(int)
            )
            regions, region_rows = build_structure_regions(
                structure_path=structure,
                protein_chain=str(state["protein_chain"]),
                crrna_chains=set(state.get("crrna_chains", [])),
                target_rna_chains=set(state.get("target_rna_chains", [])),
                hepn_positions=coordinate_hepn,
                direct_cutoff=5.0,
                second_shell_cutoff=8.0,
                buried_rsa_threshold=0.2,
            )
            del regions
            region_frame = pd.DataFrame.from_records(region_rows).rename(
                columns={
                    "rna_interface": "RNA_contact",
                    "rna_second_shell": "RNA_second_shell",
                    "burial": "structural_burial",
                }
            )
            mapping = mapping.merge(
                region_frame[
                    [
                        "coordinate_index_0",
                        "relative_sasa",
                        "structural_burial",
                        "RNA_contact",
                        "RNA_second_shell",
                        "crrna_interface",
                        "target_rna_interface",
                        "minimum_rna_distance",
                    ]
                ],
                on="coordinate_index_0",
                how="left",
                # Multiple unresolved full-sequence rows carry a null key;
                # resolved coordinate indices remain unique on both sides.
                validate="many_to_one",
            )
            mapping["coordinate_available"] = mapping["coordinate_index_0"].notna()
            mapping["protein_core"] = mapping["structural_burial"].eq("buried_core")
            mapping["domain"] = "unassigned_manual_review"
            for position, label in hepn_labels.items():
                mapping.loc[mapping["biological_index_1"] == position, "domain"] = (
                    "HEPN1" if "HEPN1" in label else "HEPN2"
                )
            mapping["domain_interface"] = pd.NA
            mapping["domain_interface_status"] = "not_annotated_no_domain_boundaries"
            mapping["HEPN_annotation"] = mapping["biological_index_1"].map(hepn_labels)
            mapping["state"] = str(state["state"])
            mapping["mapping_gate_passed"] = mapping[
                "conservation_constraint_eligible"
            ].astype(bool)
            qc = structure_qc(parse_structure(structure))
            break_residues = {
                (key.residue_number, key.insertion_code)
                for left, right, _ in qc.chain_breaks
                for key in (left, right)
                if key.chain_id == str(state["protein_chain"])
            }
            mapping["chain_break_adjacent"] = [
                (int(number), str(insertion) if pd.notna(insertion) else "")
                in break_residues
                if pd.notna(number)
                else False
                for number, insertion in zip(
                    mapping["pdb_residue_number"],
                    mapping["pdb_insertion_code"],
                    strict=True,
                )
            ]
            atomic_write_text(output / "mapping.csv", mapping.to_csv(index=False))
            review = mapping.loc[
                mapping["manual_review_required"].astype(bool)
                | ~mapping["coordinate_available"]
                | mapping["chain_break_adjacent"]
            ].copy()
            atomic_write_text(
                review_root / f"{lower}_manual_review.csv", review.to_csv(index=False)
            )
            summary.update(
                {
                    "scaffold_id": scaffold_id,
                    "RNA_contact_positions": int(
                        mapping["RNA_contact"]
                        .map(lambda value: str(value).lower() == "true")
                        .sum()
                    ),
                    "RNA_second_shell_positions": int(
                        mapping["RNA_second_shell"]
                        .map(lambda value: str(value).lower() == "true")
                        .sum()
                    ),
                    "protein_core_positions": int(mapping["protein_core"].sum()),
                    "manual_review_positions": len(review),
                }
            )
            atomic_write_text(
                output / "summary.json",
                json.dumps(summary, indent=2, sort_keys=True) + "\n",
            )
            # Replace the basic mapping table in the generated HTML with the
            # enriched table while preserving an explicit Level-0 boundary.
            html_text = (
                """<!doctype html><html><head><meta charset="utf-8">
<title>Stage 0003A residue mapping</title><style>body{font-family:system-ui;
margin:2rem}table{border-collapse:collapse;font-size:11px}th,td{border:1px solid
#ccd3db;padding:3px 5px}th{position:sticky;top:0;background:#eef2f6}</style>
</head><body><h1>"""
                + f"{pdb_id} / {scaffold_id} / {state['state']}"
                + """</h1>
<p><strong>Evidence Level 0.</strong> Conservation is eligible only where
<code>mapping_gate_passed=true</code>. Unassigned domains and interfaces are
not interpreted as negatives.</p><pre>"""
                + json.dumps(summary, indent=2, sort_keys=True)
                + "</pre>"
                + mapping.to_html(index=False, escape=True)
                + "</body></html>\n"
            )
            atomic_write_text(output / "mapping.html", html_text)
            summaries.append(summary)
    summary_frame = pd.DataFrame.from_records(summaries)
    atomic_write_text(root / "mapping_summary.csv", summary_frame.to_csv(index=False))
    total_full = int(summary_frame["full_scaffold_length"].sum())
    aggregate = {
        "state_units": len(summaries),
        "scaffolds": int(summary_frame["scaffold_id"].nunique()),
        "total_full_positions_across_states": total_full,
        "four_layer_exact_or_restored_positions": int(
            summary_frame["mapping_status_counts"]
            .map(
                lambda value: sum(
                    count
                    for status, count in value.items()
                    if status.startswith("four_layer_")
                )
            )
            .sum()
        ),
        "weighted_high_confidence_coverage": float(
            (
                summary_frame["high_confidence_coverage"]
                * summary_frame["full_scaffold_length"]
            ).sum()
            / total_full
        ),
        "is_mock": False,
        "evidence_level": 0,
    }
    atomic_write_text(
        root / "mapping_aggregate.json",
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n",
    )
    print(json.dumps(aggregate, indent=2, sort_keys=True))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
