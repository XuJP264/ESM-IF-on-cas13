"""Protein-to-RNA direct-contact and second-shell annotations."""

from __future__ import annotations

from dataclasses import dataclass

from cas13_if.structures.parser import (
    Atom,
    ResidueKey,
    atom_distance,
    group_residues,
    residue_polymer_type,
)


@dataclass(frozen=True)
class ContactAnnotation:
    protein_residue: ResidueKey
    minimum_rna_distance: float | None
    direct_rna_contact: bool
    second_shell: bool
    contacted_rna_chains: tuple[str, ...]


def annotate_rna_contacts(
    atoms: list[Atom],
    *,
    direct_cutoff: float = 5.0,
    second_shell_cutoff: float = 8.0,
) -> list[ContactAnnotation]:
    if direct_cutoff <= 0 or second_shell_cutoff < direct_cutoff:
        raise ValueError("contact cutoffs must satisfy 0 < direct <= second shell")
    residues = group_residues(atoms)
    rna_atoms = [
        atom
        for atom in atoms
        if residue_polymer_type(atom.residue.residue_name) == "rna"
        and atom.element != "H"
    ]
    output: list[ContactAnnotation] = []
    for key, residue_atoms in residues.items():
        if residue_polymer_type(key.residue_name) != "protein":
            continue
        heavy_atoms = [atom for atom in residue_atoms if atom.element != "H"]
        pairs = [
            (atom_distance(protein_atom, rna_atom), rna_atom.residue.chain_id)
            for protein_atom in heavy_atoms
            for rna_atom in rna_atoms
        ]
        minimum = min((distance for distance, _ in pairs), default=None)
        contacted_chains = tuple(
            sorted({chain for distance, chain in pairs if distance <= direct_cutoff})
        )
        output.append(
            ContactAnnotation(
                protein_residue=key,
                minimum_rna_distance=minimum,
                direct_rna_contact=minimum is not None and minimum <= direct_cutoff,
                second_shell=(
                    minimum is not None
                    and direct_cutoff < minimum <= second_shell_cutoff
                ),
                contacted_rna_chains=contacted_chains,
            )
        )
    return sorted(
        output,
        key=lambda item: (
            item.protein_residue.chain_id,
            item.protein_residue.residue_number,
            item.protein_residue.insertion_code,
        ),
    )
