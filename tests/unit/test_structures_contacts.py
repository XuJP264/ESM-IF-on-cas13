from pathlib import Path

import pytest

from cas13_if.structures.contacts import annotate_rna_contacts
from cas13_if.structures.parser import (
    parse_pdb,
    protein_chain_sequence,
    structure_qc,
)

PDB = Path("tests/fixtures/minimal_complex.pdb")


def test_pdb_insertion_missing_atoms_and_chain_types() -> None:
    atoms = parse_pdb(PDB)
    qc = structure_qc(atoms)
    assert qc.protein_chains == ("A",)
    assert qc.rna_chains == ("R",)
    assert qc.protein_residue_count == 2
    assert qc.rna_residue_count == 1
    assert len(qc.missing_backbone) == 1
    key, missing = qc.missing_backbone[0]
    assert key.insertion_code == "A"
    assert missing == ("C",)
    assert qc.coordinates_are_finite
    sequence, keys = protein_chain_sequence(atoms, "A")
    assert sequence == "AG"
    assert keys[1].residue_number == 2


def test_rna_contact_annotation() -> None:
    contacts = annotate_rna_contacts(
        parse_pdb(PDB), direct_cutoff=5, second_shell_cutoff=8
    )
    assert len(contacts) == 2
    assert contacts[0].direct_rna_contact
    assert contacts[0].contacted_rna_chains == ("R",)
    assert contacts[0].minimum_rna_distance == pytest.approx(1.526433752)
