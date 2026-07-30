"""Strict PDB/mmCIF atom parsing and protein/RNA residue classification."""

from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path
from typing import Literal

from Bio.PDB import MMCIFParser  # type: ignore[attr-defined]

PROTEIN_RESIDUES = {
    "ALA": "A",
    "ARG": "R",
    "ASN": "N",
    "ASP": "D",
    "CYS": "C",
    "GLN": "Q",
    "GLU": "E",
    "GLY": "G",
    "HIS": "H",
    "ILE": "I",
    "LEU": "L",
    "LYS": "K",
    "MET": "M",
    "PHE": "F",
    "PRO": "P",
    "SER": "S",
    "THR": "T",
    "TRP": "W",
    "TYR": "Y",
    "VAL": "V",
    "MSE": "M",
    "SEC": "C",
    "PYL": "K",
}
RNA_RESIDUES = {
    "A",
    "C",
    "G",
    "U",
    "I",
    "RA",
    "RC",
    "RG",
    "RU",
    "ADE",
    "CYT",
    "GUA",
    "URA",
    "PSU",
}


@dataclass(frozen=True)
class ResidueKey:
    chain_id: str
    residue_number: int
    insertion_code: str
    residue_name: str


@dataclass(frozen=True)
class Atom:
    record_type: str
    serial: int
    name: str
    alternate_location: str
    residue: ResidueKey
    x: float
    y: float
    z: float
    occupancy: float | None
    element: str

    @property
    def coordinate(self) -> tuple[float, float, float]:
        return self.x, self.y, self.z


@dataclass(frozen=True)
class StructureQC:
    protein_chains: tuple[str, ...]
    rna_chains: tuple[str, ...]
    other_chains: tuple[str, ...]
    protein_residue_count: int
    rna_residue_count: int
    missing_backbone: tuple[tuple[ResidueKey, tuple[str, ...]], ...]
    chain_breaks: tuple[tuple[ResidueKey, ResidueKey, float], ...]
    modified_residues: tuple[ResidueKey, ...]
    alternate_locations_discarded: int
    coordinates_are_finite: bool


def residue_polymer_type(
    residue_name: str,
) -> Literal["protein", "rna", "other"]:
    name = residue_name.strip().upper()
    if name in PROTEIN_RESIDUES:
        return "protein"
    if name in RNA_RESIDUES:
        return "rna"
    return "other"


def parse_structure(path: Path) -> list[Atom]:
    suffix = path.suffix.lower()
    if suffix in {".cif", ".mmcif"}:
        return parse_mmcif(path)
    if suffix in {".pdb", ".ent"}:
        return parse_pdb(path)
    raise ValueError(f"unsupported structure format: {path}")


def parse_pdb(path: Path) -> list[Atom]:
    atoms: list[Atom] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, line in enumerate(handle, start=1):
            record = line[0:6].strip()
            if record not in {"ATOM", "HETATM"}:
                continue
            try:
                alternate = line[16:17].strip()
                if alternate not in {"", "A"}:
                    continue
                atom = Atom(
                    record_type=record,
                    serial=int(line[6:11]),
                    name=line[12:16].strip(),
                    alternate_location=alternate,
                    residue=ResidueKey(
                        chain_id=line[21:22].strip() or "_",
                        residue_number=int(line[22:26]),
                        insertion_code=line[26:27].strip(),
                        residue_name=line[17:20].strip().upper(),
                    ),
                    x=float(line[30:38]),
                    y=float(line[38:46]),
                    z=float(line[46:54]),
                    occupancy=_optional_float(line[54:60]),
                    element=line[76:78].strip().upper()
                    or line[12:16].strip()[0].upper(),
                )
            except (ValueError, IndexError) as exc:
                raise ValueError(f"malformed PDB atom at line {line_number}") from exc
            atoms.append(atom)
    if not atoms:
        raise ValueError(f"structure contains no atoms: {path}")
    return atoms


def parse_mmcif(path: Path) -> list[Atom]:
    parser = MMCIFParser(  # type: ignore[no-untyped-call]
        QUIET=True, auth_chains=True, auth_residues=True
    )
    structure = parser.get_structure(  # type: ignore[no-untyped-call]
        path.stem, str(path)
    )
    atoms: list[Atom] = []
    serial = 0
    model = next(structure.get_models())
    for chain in model:
        for residue in chain:
            hetero_flag, residue_number, insertion_code = residue.id
            residue_key = ResidueKey(
                chain_id=str(chain.id).strip() or "_",
                residue_number=int(residue_number),
                insertion_code=str(insertion_code).strip(),
                residue_name=str(residue.resname).strip().upper(),
            )
            for bio_atom in residue:
                alternate = str(bio_atom.get_altloc()).strip()
                if alternate not in {"", "A"}:
                    continue
                serial += 1
                coordinate = bio_atom.get_coord()
                occupancy = bio_atom.get_occupancy()
                atoms.append(
                    Atom(
                        record_type="HETATM" if str(hetero_flag).strip() else "ATOM",
                        serial=int(bio_atom.get_serial_number() or serial),
                        name=str(bio_atom.get_name()).strip(),
                        alternate_location=alternate,
                        residue=residue_key,
                        x=float(coordinate[0]),
                        y=float(coordinate[1]),
                        z=float(coordinate[2]),
                        occupancy=float(occupancy) if occupancy is not None else None,
                        element=str(bio_atom.element).strip().upper(),
                    )
                )
    if not atoms:
        raise ValueError(f"structure contains no atoms: {path}")
    return atoms


def _optional_float(value: str) -> float | None:
    stripped = value.strip()
    return float(stripped) if stripped else None


def group_residues(atoms: list[Atom]) -> dict[ResidueKey, list[Atom]]:
    residues: dict[ResidueKey, list[Atom]] = {}
    for atom in atoms:
        residues.setdefault(atom.residue, []).append(atom)
    return residues


def structure_qc(atoms: list[Atom]) -> StructureQC:
    residues = group_residues(atoms)
    chain_types: dict[str, set[str]] = {}
    missing: list[tuple[ResidueKey, tuple[str, ...]]] = []
    modified: list[ResidueKey] = []
    protein_count = 0
    rna_count = 0
    for key, residue_atoms in residues.items():
        polymer_type = residue_polymer_type(key.residue_name)
        chain_types.setdefault(key.chain_id, set()).add(polymer_type)
        if polymer_type == "protein":
            protein_count += 1
            atom_names = {atom.name for atom in residue_atoms}
            missing_names = tuple(
                atom_name
                for atom_name in ("N", "CA", "C")
                if atom_name not in atom_names
            )
            if missing_names:
                missing.append((key, missing_names))
            if key.residue_name in {"MSE", "SEC", "PYL"}:
                modified.append(key)
        elif polymer_type == "rna":
            rna_count += 1
            if key.residue_name not in {"A", "C", "G", "U", "RA", "RC", "RG", "RU"}:
                modified.append(key)
    protein_chains = tuple(
        sorted(chain for chain, types in chain_types.items() if "protein" in types)
    )
    rna_chains = tuple(
        sorted(chain for chain, types in chain_types.items() if "rna" in types)
    )
    other_chains = tuple(
        sorted(
            chain
            for chain, types in chain_types.items()
            if not types.intersection({"protein", "rna"})
        )
    )
    breaks: list[tuple[ResidueKey, ResidueKey, float]] = []
    for chain in protein_chains:
        chain_residues = sorted(
            (
                (key, residue_atoms)
                for key, residue_atoms in residues.items()
                if key.chain_id == chain
                and residue_polymer_type(key.residue_name) == "protein"
            ),
            key=lambda item: (item[0].residue_number, item[0].insertion_code),
        )
        previous: tuple[ResidueKey, Atom] | None = None
        for key, residue_atoms in chain_residues:
            ca = next((atom for atom in residue_atoms if atom.name == "CA"), None)
            if ca is None:
                continue
            if previous is not None:
                distance = atom_distance(previous[1], ca)
                if distance > 4.5:
                    breaks.append((previous[0], key, distance))
            previous = key, ca
    return StructureQC(
        protein_chains=protein_chains,
        rna_chains=rna_chains,
        other_chains=other_chains,
        protein_residue_count=protein_count,
        rna_residue_count=rna_count,
        missing_backbone=tuple(missing),
        chain_breaks=tuple(breaks),
        modified_residues=tuple(modified),
        alternate_locations_discarded=0,
        coordinates_are_finite=all(
            math.isfinite(value) for atom in atoms for value in (atom.x, atom.y, atom.z)
        ),
    )


def atom_distance(left: Atom, right: Atom) -> float:
    return math.dist(left.coordinate, right.coordinate)


def protein_chain_sequence(
    atoms: list[Atom], chain_id: str
) -> tuple[str, list[ResidueKey]]:
    residues = group_residues(atoms)
    keys = sorted(
        (
            key
            for key in residues
            if key.chain_id == chain_id
            and residue_polymer_type(key.residue_name) == "protein"
        ),
        key=lambda key: (key.residue_number, key.insertion_code),
    )
    sequence = "".join(PROTEIN_RESIDUES[key.residue_name] for key in keys)
    return sequence, keys
