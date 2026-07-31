"""Relative solvent accessibility from the standard Shrake-Rupley algorithm."""

from __future__ import annotations

from pathlib import Path
from typing import Any

from Bio.PDB import MMCIFParser, PDBParser, ShrakeRupley  # type: ignore[attr-defined]

from cas13_if.structures.parser import PROTEIN_RESIDUES, ResidueKey

# Maximum solvent-accessible areas (Å²) for the extended Gly-X-Gly reference
# state from Tien et al., PLoS ONE 2013, doi:10.1371/journal.pone.0080635.
MAXIMUM_ASA = {
    "A": 129.0,
    "R": 274.0,
    "N": 195.0,
    "D": 193.0,
    "C": 167.0,
    "Q": 225.0,
    "E": 223.0,
    "G": 104.0,
    "H": 224.0,
    "I": 197.0,
    "L": 201.0,
    "K": 236.0,
    "M": 224.0,
    "F": 240.0,
    "P": 159.0,
    "S": 155.0,
    "T": 172.0,
    "W": 285.0,
    "Y": 263.0,
    "V": 174.0,
}


def relative_solvent_accessibility(
    path: Path,
    *,
    chain_id: str,
    probe_radius: float = 1.4,
    n_points: int = 100,
) -> dict[ResidueKey, float]:
    """Return per-residue relative SASA, retaining author numbering."""
    suffix = path.suffix.lower()
    parser: Any
    if suffix in {".cif", ".mmcif"}:
        parser = MMCIFParser(  # type: ignore[no-untyped-call]
            QUIET=True, auth_chains=True, auth_residues=True
        )
    elif suffix in {".pdb", ".ent"}:
        parser = PDBParser(QUIET=True)  # type: ignore[no-untyped-call]
    else:
        raise ValueError(f"unsupported structure format: {path}")
    structure = parser.get_structure(path.stem, str(path))
    model = next(structure.get_models())
    if chain_id not in model:
        raise ValueError(f"chain {chain_id!r} is missing from {path}")
    calculator = ShrakeRupley(  # type: ignore[no-untyped-call]
        probe_radius=probe_radius, n_points=n_points
    )
    calculator.compute(model, level="R")  # type: ignore[no-untyped-call]
    output: dict[ResidueKey, float] = {}
    for residue in model[chain_id]:
        name = str(residue.resname).strip().upper()
        amino_acid = PROTEIN_RESIDUES.get(name)
        if amino_acid is None:
            continue
        _, residue_number, insertion_code = residue.id
        absolute = float(residue.sasa)
        key = ResidueKey(
            chain_id=chain_id,
            residue_number=int(residue_number),
            insertion_code=str(insertion_code).strip(),
            residue_name=name,
        )
        output[key] = min(absolute / MAXIMUM_ASA[amino_acid], 1.0)
    if not output:
        raise ValueError(f"no protein residue SASA values for chain {chain_id!r}")
    return output
