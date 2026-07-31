from pathlib import Path

from cas13_if.evaluation.regions import build_structure_regions
from cas13_if.structures.parser import parse_structure, protein_chain_sequence


def test_region_builder_maps_contacts_shell_and_burial(monkeypatch) -> None:
    path = Path("tests/fixtures/minimal_complex.pdb")
    _, keys = protein_chain_sequence(parse_structure(path), "A")
    monkeypatch.setattr(
        "cas13_if.evaluation.regions.relative_solvent_accessibility",
        lambda *_args, **_kwargs: {keys[0]: 0.1, keys[1]: 0.6},
    )
    regions, rows = build_structure_regions(
        structure_path=path,
        protein_chain="A",
        crrna_chains={"R"},
        target_rna_chains=set(),
        hepn_positions={0},
        direct_cutoff=5.0,
        second_shell_cutoff=8.0,
        buried_rsa_threshold=0.2,
    )
    assert len(rows) == 2
    assert regions["buried_core"] == {0}
    assert regions["surface"] == {1}
    assert regions["hepn_region"] == {0}
    assert regions["rna_interface"]
    assert regions["crrna_interface"] == regions["rna_interface"]
