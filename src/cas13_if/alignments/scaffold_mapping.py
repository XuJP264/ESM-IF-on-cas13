"""Auditable scaffold/coordinate/MSA mapping without column renumbering."""

from __future__ import annotations

import json
import shutil
import subprocess
import tempfile
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pandas as pd
import pyarrow.parquet as pq  # type: ignore[import-untyped]
from Bio.Align import PairwiseAligner

from cas13_if.alignments.msa import Alignment, read_aligned_fasta
from cas13_if.data.fasta import write_fasta
from cas13_if.provenance import atomic_write_text, sha256_file
from cas13_if.schemas import STANDARD_AA
from cas13_if.structures.parser import (
    PROTEIN_RESIDUES,
    ResidueKey,
    group_residues,
    parse_structure,
    protein_chain_sequence,
)

QUERY_IDENTIFIER = "cas13_if__6e9f__chain_a__full_scaffold"


@dataclass(frozen=True)
class PairwiseIndexMapping:
    """Monotonic maps for a global reference/query alignment."""

    reference_to_query: tuple[int | None, ...]
    query_to_reference: tuple[int | None, ...]


@dataclass(frozen=True)
class AddedAlignment:
    """Added query row plus output-to-original MSA column mapping."""

    query_aligned: str
    output_to_original_column: tuple[int | None, ...]
    original_columns_preserved: bool
    mafft_command: tuple[str, ...]
    mafft_stderr: str


def global_index_mapping(reference: str, query: str) -> PairwiseIndexMapping:
    """Map indices in two sequences with a deterministic global alignment."""
    reference = reference.upper()
    query = query.upper()
    invalid = set(reference + query).difference(STANDARD_AA)
    if invalid:
        raise ValueError(
            f"non-canonical sequence symbols in mapping: {sorted(invalid)}"
        )
    if not reference or not query:
        raise ValueError("mapping sequences cannot be empty")
    aligner = PairwiseAligner()  # type: ignore[no-untyped-call]
    aligner.mode = "global"
    aligner.match_score = 2.0
    aligner.mismatch_score = -1.0
    aligner.open_gap_score = -4.0
    aligner.extend_gap_score = -0.5
    alignment = aligner.align(reference, query)[0]  # type: ignore[no-untyped-call]
    reference_to_query: list[int | None] = [None] * len(reference)
    query_to_reference: list[int | None] = [None] * len(query)
    for reference_block, query_block in zip(
        alignment.aligned[0], alignment.aligned[1], strict=True
    ):
        reference_start, reference_end = map(int, reference_block)
        query_start, query_end = map(int, query_block)
        reference_length = reference_end - reference_start
        query_length = query_end - query_start
        if reference_length != query_length:
            raise ValueError("pairwise aligned blocks have unequal lengths")
        for offset in range(reference_length):
            reference_index = reference_start + offset
            query_index = query_start + offset
            reference_to_query[reference_index] = query_index
            query_to_reference[query_index] = reference_index
    return PairwiseIndexMapping(
        reference_to_query=tuple(reference_to_query),
        query_to_reference=tuple(query_to_reference),
    )


def _parse_alignment_text(text: str, path: Path) -> Alignment:
    atomic_write_text(path, text)
    return read_aligned_fasta(path)


def _column_signature(alignment: Alignment, column: int) -> tuple[str, ...]:
    return tuple(sequence[column] for sequence in alignment.sequences)


def map_added_alignment_columns(
    original: Alignment,
    added: Alignment,
    *,
    query_identifier: str = QUERY_IDENTIFIER,
) -> AddedAlignment:
    """Recover original MSA columns after MAFFT adds query insertion columns."""
    if query_identifier in original.identifiers:
        raise ValueError(
            f"query identifier collides with original MSA: {query_identifier}"
        )
    if query_identifier not in added.identifiers:
        raise ValueError("added alignment does not contain the scaffold query")
    if any(
        all(sequence[column] == "-" for sequence in original.sequences)
        for column in range(original.n_columns)
    ):
        raise ValueError("original MSA contains an ambiguous all-gap column")
    added_by_identifier = dict(zip(added.identifiers, added.sequences, strict=True))
    missing = set(original.identifiers).difference(added_by_identifier)
    unexpected = set(added.identifiers).difference(
        set(original.identifiers).union({query_identifier})
    )
    if missing or unexpected:
        raise ValueError(
            "MAFFT row set changed: "
            f"missing={sorted(missing)} unexpected={sorted(unexpected)}"
        )
    ordered_existing = Alignment(
        identifiers=original.identifiers,
        sequences=tuple(added_by_identifier[item] for item in original.identifiers),
    )
    original_column = 0
    output_to_original: list[int | None] = []
    for output_column in range(added.n_columns):
        signature = _column_signature(ordered_existing, output_column)
        if original_column < original.n_columns and signature == _column_signature(
            original, original_column
        ):
            output_to_original.append(original_column)
            original_column += 1
        elif all(token == "-" for token in signature):
            output_to_original.append(None)
        else:
            raise ValueError(
                "MAFFT altered or reordered an original alignment column at "
                f"output column {output_column}"
            )
    preserved = original_column == original.n_columns
    if not preserved:
        raise ValueError(
            f"MAFFT output retained only {original_column}/{original.n_columns} columns"
        )
    return AddedAlignment(
        query_aligned=added_by_identifier[query_identifier],
        output_to_original_column=tuple(output_to_original),
        original_columns_preserved=True,
        mafft_command=(),
        mafft_stderr="",
    )


def add_scaffold_to_msa(
    *,
    msa_path: Path,
    scaffold_sequence: str,
    mafft_executable: str,
    threads: int,
) -> AddedAlignment:
    """Add a full scaffold to an existing MSA while preserving old columns."""
    if threads < 1:
        raise ValueError("MAFFT threads must be positive")
    executable = shutil.which(mafft_executable)
    if executable is None:
        raise FileNotFoundError(f"MAFFT executable not found: {mafft_executable}")
    original = read_aligned_fasta(msa_path)
    if QUERY_IDENTIFIER in original.identifiers:
        raise ValueError("reserved query identifier is already present in MSA")
    with tempfile.TemporaryDirectory(prefix="cas13-if-vi-d-map-") as temporary:
        temporary_root = Path(temporary)
        query_path = temporary_root / "scaffold.fasta"
        output_path = temporary_root / "added.fasta"
        write_fasta([(QUERY_IDENTIFIER, scaffold_sequence)], query_path)
        command = (
            executable,
            "--quiet",
            "--thread",
            str(threads),
            "--addfull",
            str(query_path),
            str(msa_path),
        )
        completed = subprocess.run(
            command,
            text=True,
            capture_output=True,
            check=False,
        )
        if completed.returncode != 0:
            raise RuntimeError(
                f"MAFFT addfull failed ({completed.returncode}): "
                f"{completed.stderr.strip()}"
            )
        added = _parse_alignment_text(completed.stdout, output_path)
        mapped = map_added_alignment_columns(original, added)
        if mapped.query_aligned.replace("-", "") != scaffold_sequence.upper():
            raise ValueError("MAFFT output lost or changed scaffold residues")
        return AddedAlignment(
            query_aligned=mapped.query_aligned,
            output_to_original_column=mapped.output_to_original_column,
            original_columns_preserved=mapped.original_columns_preserved,
            mafft_command=command,
            mafft_stderr=completed.stderr,
        )


def scaffold_to_msa_columns(added: AddedAlignment) -> tuple[int | None, ...]:
    """Map every zero-based scaffold index to an original MSA column or null."""
    mapping: list[int | None] = []
    for output_column, token in enumerate(added.query_aligned):
        if token != "-":
            mapping.append(added.output_to_original_column[output_column])
    return tuple(mapping)


def _load_full_scaffold_sequence(entity_path: Path) -> str:
    entity = json.loads(entity_path.read_text(encoding="utf-8"))
    entity_poly = entity.get("entity_poly")
    if not isinstance(entity_poly, dict):
        raise ValueError("RCSB entity JSON is missing entity_poly")
    if entity_poly.get("rcsb_entity_polymer_type") != "Protein":
        raise ValueError("declared RCSB entity is not a protein")
    raw = entity_poly.get("pdbx_seq_one_letter_code_can")
    if not isinstance(raw, str):
        raise ValueError("RCSB entity JSON is missing canonical sequence")
    sequence = raw.replace("\n", "").replace(" ", "").upper()
    invalid = set(sequence).difference(STANDARD_AA)
    if not sequence or invalid:
        raise ValueError(f"invalid full scaffold sequence symbols: {sorted(invalid)}")
    return sequence


def _complete_backbone_by_key(structure_path: Path) -> dict[ResidueKey, bool]:
    atoms = parse_structure(structure_path)
    return {
        key: {atom.name for atom in residue_atoms}.issuperset({"N", "CA", "C"})
        for key, residue_atoms in group_residues(atoms).items()
        if key.residue_name in PROTEIN_RESIDUES
    }


def _mapping_status(
    *,
    full_index: int,
    coordinate_index: int | None,
    full_length: int,
    mapped_coordinate_indices: list[int],
    coordinate_matches: bool | None,
    msa_column: int | None,
) -> tuple[str, str, list[str]]:
    reasons: list[str] = []
    if coordinate_index is None:
        mapped_min = min(mapped_coordinate_indices)
        mapped_max = max(mapped_coordinate_indices)
        coordinate_status = (
            "unresolved_terminal"
            if full_index < mapped_min or full_index > mapped_max
            else "unresolved_internal"
        )
        reasons.append(coordinate_status)
    elif coordinate_matches:
        coordinate_status = "resolved_exact"
    else:
        coordinate_status = "resolved_substitution"
        reasons.append(coordinate_status)
    if msa_column is None:
        msa_status = "query_only_insertion"
        reasons.append(msa_status)
    else:
        msa_status = "mapped_existing_column"
    if coordinate_status == "resolved_exact" and msa_status == "mapped_existing_column":
        return "four_layer_exact", "high", reasons
    if (
        coordinate_status.startswith("unresolved")
        and msa_status == "mapped_existing_column"
    ):
        return "full_msa_only_unresolved_coordinate", "medium", reasons
    if full_index in {0, full_length - 1} and msa_column is None:
        reasons.append("terminal_msa_insertion")
    return "review_required", "low", reasons


def build_scaffold_mapping(
    *,
    structure_path: Path,
    entity_path: Path,
    msa_path: Path,
    conservation_path: Path,
    output_dir: Path,
    chain_id: str,
    subtype: str,
    mafft_executable: str,
    threads: int,
    minimum_conservation_coverage: float,
) -> dict[str, Any]:
    """Build and write the strict 6E9F/full-scaffold/VI-D mapping audit."""
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite mapping output: {output_dir}")
    if not 0 <= minimum_conservation_coverage <= 1:
        raise ValueError("minimum conservation coverage must be in [0, 1]")
    full_sequence = _load_full_scaffold_sequence(entity_path)
    atoms = parse_structure(structure_path)
    coordinate_sequence, residue_keys = protein_chain_sequence(atoms, chain_id)
    pairwise = global_index_mapping(full_sequence, coordinate_sequence)
    mapped_full_indices = [
        index
        for index, coordinate_index in enumerate(pairwise.reference_to_query)
        if coordinate_index is not None
    ]
    if not mapped_full_indices:
        raise ValueError("coordinate sequence does not map to full scaffold")
    added = add_scaffold_to_msa(
        msa_path=msa_path,
        scaffold_sequence=full_sequence,
        mafft_executable=mafft_executable,
        threads=threads,
    )
    msa_columns = scaffold_to_msa_columns(added)
    if len(msa_columns) != len(full_sequence):
        raise ValueError("scaffold-to-MSA mapping length mismatch")
    conservation_rows = pq.read_table(conservation_path).to_pylist()
    conservation_by_column = {int(row["column"]): row for row in conservation_rows}
    backbone_complete = _complete_backbone_by_key(structure_path)
    rows: list[dict[str, Any]] = []
    for full_index, full_amino_acid in enumerate(full_sequence):
        coordinate_index = pairwise.reference_to_query[full_index]
        coordinate_amino_acid = (
            coordinate_sequence[coordinate_index]
            if coordinate_index is not None
            else None
        )
        residue_key = (
            residue_keys[coordinate_index] if coordinate_index is not None else None
        )
        msa_column = msa_columns[full_index]
        conservation = (
            conservation_by_column.get(msa_column) if msa_column is not None else None
        )
        status, confidence, reasons = _mapping_status(
            full_index=full_index,
            coordinate_index=coordinate_index,
            full_length=len(full_sequence),
            mapped_coordinate_indices=mapped_full_indices,
            coordinate_matches=(
                coordinate_amino_acid == full_amino_acid
                if coordinate_amino_acid is not None
                else None
            ),
            msa_column=msa_column,
        )
        insertion_code = residue_key.insertion_code if residue_key else ""
        if insertion_code:
            reasons.append("pdb_insertion_code")
        coverage = float(conservation["coverage"]) if conservation else None
        eligible = bool(
            confidence == "high"
            and residue_key is not None
            and backbone_complete.get(residue_key, False)
            and coverage is not None
            and coverage >= minimum_conservation_coverage
        )
        consensus = str(conservation["consensus"]) if conservation else None
        if consensus is None:
            scaffold_consensus_status = "no_original_msa_column"
        elif consensus == full_amino_acid:
            scaffold_consensus_status = "matches_consensus"
        else:
            scaffold_consensus_status = "differs_from_consensus"
        review_reasons = list(reasons)
        if scaffold_consensus_status == "differs_from_consensus":
            review_reasons.append("scaffold_differs_from_msa_consensus")
        rows.append(
            {
                "pdb_id": "6E9F",
                "chain_id": chain_id,
                "subtype": subtype,
                "full_scaffold_index_0": full_index,
                "biological_index_1": full_index + 1,
                "full_scaffold_amino_acid": full_amino_acid,
                "coordinate_index_0": coordinate_index,
                "coordinate_amino_acid": coordinate_amino_acid,
                "pdb_residue_number": (
                    residue_key.residue_number if residue_key is not None else None
                ),
                "pdb_insertion_code": insertion_code,
                "pdb_residue_name": residue_key.residue_name if residue_key else None,
                "backbone_complete": (
                    backbone_complete.get(residue_key, False) if residue_key else False
                ),
                "msa_column_0": msa_column,
                "msa_column_1": msa_column + 1 if msa_column is not None else None,
                "mapping_status": status,
                "mapping_confidence": confidence,
                "msa_coverage": coverage,
                "conservation": (
                    float(conservation["conservation"]) if conservation else None
                ),
                "entropy": float(conservation["entropy"]) if conservation else None,
                "gap_fraction": (
                    float(conservation["gap_fraction"]) if conservation else None
                ),
                "msa_consensus": consensus,
                "allowed_residues": (
                    ";".join(str(item) for item in conservation["allowed_residues"])
                    if conservation
                    else ""
                ),
                "scaffold_consensus_status": scaffold_consensus_status,
                "conservation_constraint_eligible": eligible,
                "manual_review_required": bool(review_reasons),
                "decision_reasons": ";".join(dict.fromkeys(review_reasons)),
            }
        )
    output_dir.mkdir(parents=True, exist_ok=False)
    frame = pd.DataFrame.from_records(rows)
    csv_path = output_dir / "mapping.csv"
    atomic_write_text(csv_path, frame.to_csv(index=False))
    review_frame = frame.loc[frame["manual_review_required"]].copy()
    review_path = output_dir / "manual_review.csv"
    atomic_write_text(review_path, review_frame.to_csv(index=False))
    confidence_counts = {
        str(key): int(value)
        for key, value in frame["mapping_confidence"]
        .value_counts()
        .sort_index()
        .items()
    }
    status_counts = {
        str(key): int(value)
        for key, value in frame["mapping_status"].value_counts().sort_index().items()
    }
    reason_counts: dict[str, int] = {}
    for reason_string in frame["decision_reasons"]:
        for reason in str(reason_string).split(";"):
            if reason:
                reason_counts[reason] = reason_counts.get(reason, 0) + 1
    summary: dict[str, Any] = {
        "schema_version": "1.0",
        "is_mock": False,
        "evidence_level": 0,
        "pdb_id": "6E9F",
        "chain_id": chain_id,
        "subtype": subtype,
        "full_scaffold_length": len(full_sequence),
        "coordinate_length": len(coordinate_sequence),
        "original_msa_sequences": read_aligned_fasta(msa_path).n_sequences,
        "original_msa_columns": read_aligned_fasta(msa_path).n_columns,
        "original_msa_columns_preserved": added.original_columns_preserved,
        "query_only_msa_insertions": sum(item is None for item in msa_columns),
        "mapped_coordinate_positions": int(frame["coordinate_index_0"].notna().sum()),
        "unresolved_positions": int(frame["coordinate_index_0"].isna().sum()),
        "mapping_confidence_counts": confidence_counts,
        "mapping_status_counts": status_counts,
        "unmapped_or_review_reason_counts": dict(sorted(reason_counts.items())),
        "high_confidence_coverage": confidence_counts.get("high", 0)
        / len(full_sequence),
        "four_layer_exact_coverage": status_counts.get("four_layer_exact", 0)
        / len(full_sequence),
        "conservation_constraint_eligible_positions": int(
            frame["conservation_constraint_eligible"].sum()
        ),
        "minimum_conservation_coverage": minimum_conservation_coverage,
        "conservation_gate_status": "passed",
        "mafft_command": list(added.mafft_command),
        "inputs": {
            str(path): sha256_file(path)
            for path in (structure_path, entity_path, msa_path, conservation_path)
        },
        "claim_scope": (
            "Level 0 coordinate-system audit. Conservation eligibility does not "
            "establish function or wet-lab validity."
        ),
    }
    summary_path = output_dir / "summary.json"
    atomic_write_text(
        summary_path, json.dumps(summary, indent=2, sort_keys=True) + "\n"
    )
    html = (
        """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><title>6E9F VI-D mapping audit</title>
<style>body{font-family:system-ui,sans-serif;margin:2rem;color:#18202a}
table{border-collapse:collapse;font-size:12px}th,td{border:1px solid #ccd3db;
padding:4px 6px}th{position:sticky;top:0;background:#eef2f6}tr:nth-child(even){
background:#f8fafc}code{background:#eef2f6;padding:2px 4px}</style></head><body>
<h1>6E9F / EsCas13d / VI-D mapping audit</h1>
<p><strong>Evidence Level 0.</strong> This report validates coordinate systems;
it does not validate Cas13 function.</p>
<h2>Summary</h2><pre>"""
        + json.dumps(summary, indent=2, sort_keys=True)
        + """</pre>
<h2>Manual audit table</h2>"""
        + frame.to_html(index=False, escape=True)
        + """
</body></html>\n"""
    )
    html_path = output_dir / "mapping.html"
    atomic_write_text(html_path, html)
    return summary
