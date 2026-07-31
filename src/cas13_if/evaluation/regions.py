"""Reusable experimental-structure region annotations for fair comparisons."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from cas13_if.structures.contacts import annotate_rna_contacts
from cas13_if.structures.parser import (
    Atom,
    parse_structure,
    protein_chain_sequence,
    residue_polymer_type,
)
from cas13_if.structures.sasa import relative_solvent_accessibility


def build_structure_regions(
    *,
    structure_path: Path,
    protein_chain: str,
    crrna_chains: set[str],
    target_rna_chains: set[str],
    hepn_positions: set[int],
    direct_cutoff: float,
    second_shell_cutoff: float,
    buried_rsa_threshold: float,
) -> tuple[dict[str, set[int]], list[dict[str, Any]]]:
    """Return coordinate-index regions and auditable residue annotations."""
    atoms = parse_structure(structure_path)
    _, residue_keys = protein_chain_sequence(atoms, protein_chain)
    rna_chains = crrna_chains.union(target_rna_chains)
    selected_atoms: list[Atom] = [
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
    contacts = annotate_rna_contacts(
        selected_atoms,
        direct_cutoff=direct_cutoff,
        second_shell_cutoff=second_shell_cutoff,
    )
    rsa = relative_solvent_accessibility(structure_path, chain_id=protein_chain)
    index_by_key = {key: index for index, key in enumerate(residue_keys)}
    regions: dict[str, set[int]] = {
        "buried_core": set(),
        "surface": set(),
        "rna_interface": set(),
        "crrna_interface": set(),
        "target_rna_interface": set(),
        "rna_second_shell": set(),
        "hepn_region": set(hepn_positions),
    }
    rows: list[dict[str, Any]] = []
    for annotation in contacts:
        index = index_by_key[annotation.protein_residue]
        accessibility = rsa.get(annotation.protein_residue)
        if accessibility is None:
            raise ValueError(f"missing RSA for {annotation.protein_residue}")
        if accessibility < buried_rsa_threshold:
            regions["buried_core"].add(index)
        else:
            regions["surface"].add(index)
        contacted = set(annotation.contacted_rna_chains)
        if annotation.direct_rna_contact:
            regions["rna_interface"].add(index)
        if contacted.intersection(crrna_chains):
            regions["crrna_interface"].add(index)
        if contacted.intersection(target_rna_chains):
            regions["target_rna_interface"].add(index)
        if annotation.second_shell:
            regions["rna_second_shell"].add(index)
        rows.append(
            {
                "coordinate_index_0": index,
                "pdb_residue_number": annotation.protein_residue.residue_number,
                "pdb_insertion_code": annotation.protein_residue.insertion_code,
                "amino_acid": residue_keys[index].residue_name,
                "relative_sasa": accessibility,
                "burial": (
                    "buried_core" if accessibility < buried_rsa_threshold else "surface"
                ),
                "minimum_rna_distance": annotation.minimum_rna_distance,
                "rna_interface": annotation.direct_rna_contact,
                "rna_second_shell": annotation.second_shell,
                "crrna_interface": bool(contacted.intersection(crrna_chains)),
                "target_rna_interface": bool(contacted.intersection(target_rna_chains)),
                "hepn_region": index in hepn_positions,
            }
        )
    if len(rows) != len(residue_keys):
        raise ValueError(
            f"region rows {len(rows)} do not match coordinates {len(residue_keys)}"
        )
    return regions, rows
