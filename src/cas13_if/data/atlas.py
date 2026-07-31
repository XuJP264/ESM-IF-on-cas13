"""Streaming CRISPR-Cas Atlas normalization and conservative Cas13 pairing."""

from __future__ import annotations

import hashlib
import json
import os
import re
import sqlite3
from collections import Counter
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, TextIO

import pyarrow as pa  # type: ignore[import-untyped]
import pyarrow.parquet as pq  # type: ignore[import-untyped]

CAS13_NOMENCLATURE_SUBTYPES = {letter: f"VI-{letter}" for letter in "ABCDEFGHIJ"}


class AtlasParseError(ValueError):
    """Raised when the top-level Atlas JSON stream is malformed."""


@dataclass(frozen=True)
class Cas13Record:
    operon_id: str
    subtype: str
    subtype_raw: str
    subtype_source: str
    subtype_conflict: bool
    effector_name: str
    hmm_name: str | None
    protein_sequence: str
    sequence_sha256: str
    protein_length: int
    source_length_field: int | None
    evalue: float | None
    score: float | None
    truncated: str | None
    source: str | None
    taxonomy: str | None
    metadata_json: str


@dataclass(frozen=True)
class PairingRecord:
    operon_id: str
    subtype: str
    subtype_raw: str
    subtype_source: str
    subtype_conflict: bool
    effector_name: str
    protein_sequence: str
    direct_repeat: str
    direct_repeat_raw: str
    spacers_json: str
    orientation: str
    orientation_source: str
    pairing_confidence: str
    ambiguity_reason: str | None


@dataclass(frozen=True)
class CasEffectorRecord:
    operon_id: str
    subtype: str
    gene_name: str
    hmm_name: str | None
    protein_sequence: str
    sequence_sha256: str
    protein_length: int
    evalue: float | None
    score: float | None
    truncated: str | None
    source: str | None
    taxonomy: str | None


@dataclass(frozen=True)
class CrisprArrayRecord:
    operon_id: str
    subtype: str
    array_index: int
    direct_repeat_raw: str
    direct_repeat: str
    spacers_json: str
    spacer_count: int
    orientation: str
    orientation_source: str


def iter_json_array(path: Path, chunk_size: int = 1024 * 1024) -> Iterator[Any]:
    """Yield elements of one top-level JSON array without loading it all."""
    if chunk_size < 128:
        raise ValueError("chunk_size must be at least 128 bytes")
    decoder = json.JSONDecoder()
    with path.open(encoding="utf-8") as handle:
        yield from _iter_json_array_handle(handle, decoder, chunk_size)


def _iter_json_array_handle(
    handle: TextIO,
    decoder: json.JSONDecoder,
    chunk_size: int,
) -> Iterator[Any]:
    buffer = ""
    position = 0
    started = False
    ended = False
    eof = False
    while not ended:
        if not eof and (position >= len(buffer) or len(buffer) - position < 64):
            chunk = handle.read(chunk_size)
            buffer = buffer[position:] + chunk
            position = 0
            eof = chunk == ""
        while position < len(buffer) and buffer[position].isspace():
            position += 1
        if not started:
            if position >= len(buffer):
                if eof:
                    raise AtlasParseError("empty JSON document")
                continue
            if buffer[position] != "[":
                raise AtlasParseError("Atlas top level must be a JSON array")
            position += 1
            started = True
            continue
        while position < len(buffer) and (
            buffer[position].isspace() or buffer[position] == ","
        ):
            position += 1
        if position < len(buffer) and buffer[position] == "]":
            position += 1
            ended = True
            break
        if position >= len(buffer):
            if eof:
                raise AtlasParseError("unterminated JSON array")
            continue
        try:
            value, end = decoder.raw_decode(buffer, position)
        except json.JSONDecodeError as exc:
            if eof:
                raise AtlasParseError(
                    f"invalid JSON array near byte {exc.pos}"
                ) from exc
            chunk = handle.read(chunk_size)
            buffer = buffer[position:] + chunk
            position = 0
            eof = chunk == ""
            continue
        yield value
        position = end
    trailing = buffer[position:] + handle.read()
    if trailing.strip():
        raise AtlasParseError("unexpected content after top-level JSON array")


def is_cas13_effector(cas: dict[str, Any]) -> bool:
    values = " ".join(
        str(cas.get(key, "")) for key in ("gene_name", "hmm_name", "annotation")
    ).lower()
    return "cas13" in values or "c2c2" in values


def resolve_cas13_subtype(
    raw_subtype: str,
    cas: dict[str, Any],
) -> tuple[str, str, bool]:
    """Resolve precise Type VI subtype while retaining source disagreements."""
    raw = raw_subtype.strip().upper()
    hmm_name = str(cas.get("hmm_name") or "")
    effector_name = str(cas.get("gene_name") or cas.get("annotation") or "")
    hmm_match = re.search(r"CAS[-_]VI[-_]([A-Z][A-Z0-9]*)", hmm_name.upper())
    inferred: str | None = None
    source = "atlas_summary"
    if hmm_match:
        inferred = f"VI-{hmm_match.group(1)}"
        source = "cas_hmm_explicit"
    else:
        name_match = re.search(r"\bCAS13([A-Z])(?:\b|[^A-Z])", effector_name.upper())
        if name_match and name_match.group(1) in CAS13_NOMENCLATURE_SUBTYPES:
            inferred = CAS13_NOMENCLATURE_SUBTYPES[name_match.group(1)]
            source = "effector_name_nomenclature"

    if raw.startswith("VI-"):
        conflict = (
            inferred is not None and inferred != raw and not raw.startswith(inferred)
        )
        subtype_source = "atlas_summary" if not conflict else "atlas_summary_conflict"
        return raw, subtype_source, conflict
    if raw == "VI" or not raw:
        subtype_source = source if inferred else "atlas_summary_unresolved"
        return inferred or raw, subtype_source, False
    if inferred is not None:
        return inferred, f"{source}_summary_conflict", True
    return raw, "atlas_summary_non_type_vi", True


def normalize_rna(sequence: str) -> str:
    return "".join(sequence.split()).upper().replace("T", "U")


def reverse_complement_rna(sequence: str) -> str:
    normalized = normalize_rna(sequence)
    complement = str.maketrans("ACGUN", "UGCAN")
    return normalized.translate(complement)[::-1]


def _optional_float(value: Any) -> float | None:
    if value is None:
        return None
    parsed = float(value)
    return parsed if parsed == parsed else None


def extract_cas_effectors(operon: dict[str, Any]) -> list[CasEffectorRecord]:
    operon_id = str(operon.get("operon_id", "")).strip()
    if not operon_id:
        raise ValueError("operon is missing operon_id")
    summary = operon.get("summary") or {}
    metadata = operon.get("metadata") or {}
    subtype = str(summary.get("subtype") or "").strip()
    records: list[CasEffectorRecord] = []
    for cas in operon.get("cas") or []:
        if not isinstance(cas, dict):
            continue
        sequence = "".join(str(cas.get("protein") or "").split()).upper()
        if not sequence:
            continue
        gene_name = str(cas.get("gene_name") or cas.get("hmm_name") or "unknown")
        records.append(
            CasEffectorRecord(
                operon_id=operon_id,
                subtype=subtype,
                gene_name=gene_name,
                hmm_name=(
                    str(cas["hmm_name"]) if cas.get("hmm_name") is not None else None
                ),
                protein_sequence=sequence,
                sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                protein_length=len(sequence),
                evalue=_optional_float(cas.get("evalue")),
                score=_optional_float(cas.get("score")),
                truncated=(
                    str(cas["truncated"]) if cas.get("truncated") is not None else None
                ),
                source=(
                    str(metadata["source_db"])
                    if metadata.get("source_db") is not None
                    else None
                ),
                taxonomy=(
                    str(metadata["taxonomy"])
                    if metadata.get("taxonomy") is not None
                    else None
                ),
            )
        )
    return records


def _array_orientation(array: dict[str, Any]) -> tuple[str, str]:
    orientation = str(array.get("orientation") or "unknown").strip().lower()
    source = str(array.get("orientation_source") or "").strip()
    if orientation in {"forward", "reverse", "recovered"}:
        return orientation, source or "declared_in_source_record"
    return "unknown", source or "not_provided_by_atlas_v1.0"


def extract_crispr_arrays(operon: dict[str, Any]) -> list[CrisprArrayRecord]:
    operon_id = str(operon.get("operon_id", "")).strip()
    if not operon_id:
        raise ValueError("operon is missing operon_id")
    subtype = str((operon.get("summary") or {}).get("subtype") or "").strip()
    records: list[CrisprArrayRecord] = []
    for index, array in enumerate(operon.get("crispr") or []):
        if not isinstance(array, dict):
            continue
        orientation, source = _array_orientation(array)
        repeat_raw = normalize_rna(str(array.get("crispr_repeat") or ""))
        spacers_raw = [
            normalize_rna(str(spacer)) for spacer in (array.get("crispr_spacers") or [])
        ]
        if orientation == "reverse":
            repeat = reverse_complement_rna(repeat_raw)
            spacers = [reverse_complement_rna(spacer) for spacer in spacers_raw]
        else:
            repeat = repeat_raw
            spacers = spacers_raw
        records.append(
            CrisprArrayRecord(
                operon_id=operon_id,
                subtype=subtype,
                array_index=index,
                direct_repeat_raw=repeat_raw,
                direct_repeat=repeat,
                spacers_json=json.dumps(spacers),
                spacer_count=len(spacers),
                orientation=orientation,
                orientation_source=source,
            )
        )
    return records


def extract_cas13_records(operon: dict[str, Any]) -> list[Cas13Record]:
    operon_id = str(operon.get("operon_id", "")).strip()
    if not operon_id:
        raise ValueError("operon is missing operon_id")
    summary = operon.get("summary") or {}
    metadata = operon.get("metadata") or {}
    subtype = str(summary.get("subtype") or "").strip()
    records: list[Cas13Record] = []
    for cas in operon.get("cas") or []:
        if not isinstance(cas, dict) or not is_cas13_effector(cas):
            continue
        sequence = "".join(str(cas.get("protein") or "").split()).upper()
        if not sequence:
            continue
        name = str(cas.get("gene_name") or cas.get("hmm_name") or "Cas13")
        resolved_subtype, subtype_source, subtype_conflict = resolve_cas13_subtype(
            subtype, cas
        )
        records.append(
            Cas13Record(
                operon_id=operon_id,
                subtype=resolved_subtype,
                subtype_raw=subtype,
                subtype_source=subtype_source,
                subtype_conflict=subtype_conflict,
                effector_name=name,
                hmm_name=(
                    str(cas["hmm_name"]) if cas.get("hmm_name") is not None else None
                ),
                protein_sequence=sequence,
                sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                protein_length=len(sequence),
                source_length_field=(
                    int(cas["length"]) if cas.get("length") is not None else None
                ),
                evalue=_optional_float(cas.get("evalue")),
                score=_optional_float(cas.get("score")),
                truncated=(
                    str(cas["truncated"]) if cas.get("truncated") is not None else None
                ),
                source=(
                    str(metadata["source_db"])
                    if metadata.get("source_db") is not None
                    else None
                ),
                taxonomy=(
                    str(metadata["taxonomy"])
                    if metadata.get("taxonomy") is not None
                    else None
                ),
                metadata_json=json.dumps(metadata, sort_keys=True),
            )
        )
    return records


def pair_cas13_direct_repeat(operon: dict[str, Any]) -> PairingRecord | None:
    effectors = extract_cas13_records(operon)
    arrays = extract_crispr_arrays(operon)
    if not effectors:
        return None
    reasons: list[str] = []
    if len(effectors) != 1:
        reasons.append(f"effector_count={len(effectors)}")
    if len(arrays) != 1:
        reasons.append(f"array_count={len(arrays)}")
    repeat = arrays[0].direct_repeat if arrays else ""
    if not repeat:
        reasons.append("empty_direct_repeat")
    subtype = effectors[0].subtype
    if not subtype or not subtype.upper().startswith("VI-"):
        reasons.append("ambiguous_subtype")
    if effectors[0].subtype_conflict:
        reasons.append("subtype_conflict_with_operon_summary")
    orientation = arrays[0].orientation if arrays else "unknown"
    orientation_source = (
        arrays[0].orientation_source if arrays else "array_not_available"
    )
    if orientation == "unknown":
        reasons.append("orientation_not_recovered")
    confidence = "high" if not reasons else "ambiguous"
    spacers_json = arrays[0].spacers_json if arrays else "[]"
    return PairingRecord(
        operon_id=effectors[0].operon_id,
        subtype=subtype,
        subtype_raw=effectors[0].subtype_raw,
        subtype_source=effectors[0].subtype_source,
        subtype_conflict=effectors[0].subtype_conflict,
        effector_name=effectors[0].effector_name,
        protein_sequence=effectors[0].protein_sequence,
        direct_repeat=repeat,
        direct_repeat_raw=arrays[0].direct_repeat_raw if arrays else "",
        spacers_json=spacers_json,
        orientation=orientation,
        orientation_source=orientation_source,
        pairing_confidence=confidence,
        ambiguity_reason=";".join(reasons) or None,
    )


def exact_deduplicate(records: list[Cas13Record]) -> list[dict[str, Any]]:
    grouped: dict[str, list[Cas13Record]] = {}
    for record in records:
        grouped.setdefault(record.sequence_sha256, []).append(record)
    output: list[dict[str, Any]] = []
    for digest, members in sorted(grouped.items()):
        representative = min(members, key=lambda item: item.operon_id)
        evalues = [item.evalue for item in members if item.evalue is not None]
        scores = [item.score for item in members if item.score is not None]
        output.append(
            {
                "sequence_sha256": digest,
                "protein_sequence": representative.protein_sequence,
                "protein_length": representative.protein_length,
                "representative_operon_id": representative.operon_id,
                "record_count": len(members),
                "nonconflicting_record_count": sum(
                    not item.subtype_conflict for item in members
                ),
                "complete_record_count": sum(
                    item.truncated == "00" for item in members
                ),
                "minimum_evalue": min(evalues) if evalues else None,
                "maximum_score": max(scores) if scores else None,
                "truncated_flags": sorted(
                    {item.truncated for item in members if item.truncated}
                ),
                "subtypes": sorted({item.subtype for item in members}),
                "operon_ids": sorted(item.operon_id for item in members),
            }
        )
    return output


class _ParquetBatchWriter:
    def __init__(
        self,
        path: Path,
        schema: pa.Schema,
        *,
        batch_size: int = 10_000,
    ) -> None:
        self.path = path
        self.schema = schema
        self.batch_size = batch_size
        self.buffer: list[dict[str, Any]] = []
        self.writer = pq.ParquetWriter(path, schema, compression="zstd")
        self.row_count = 0

    def append(self, row: dict[str, Any]) -> None:
        self.buffer.append(row)
        if len(self.buffer) >= self.batch_size:
            self.flush()

    def flush(self) -> None:
        if not self.buffer:
            return
        table = pa.Table.from_pylist(self.buffer, schema=self.schema)
        self.writer.write_table(table)
        self.row_count += len(self.buffer)
        self.buffer.clear()

    def close(self) -> None:
        self.flush()
        self.writer.close()


def _schemas() -> dict[str, pa.Schema]:
    string = pa.string()
    integer = pa.int64()
    floating = pa.float64()
    boolean = pa.bool_()
    operon = pa.schema(
        [
            ("operon_id", string),
            ("subtype", string),
            ("is_type_vi", boolean),
            ("metadata_json", string),
            ("source", string),
            ("taxonomy", string),
            ("biome", string),
            ("assembly_type", string),
            ("n_cas", integer),
            ("n_crispr", integer),
            ("n_spacers", integer),
            ("crispr_arrays_json", string),
            ("cas_annotations_json", string),
            ("genomic_context_json", string),
            ("genomic_context_available", boolean),
        ]
    )
    cas_effector = pa.schema(
        [
            ("operon_id", string),
            ("subtype", string),
            ("gene_name", string),
            ("hmm_name", string),
            ("protein_sequence", string),
            ("sequence_sha256", string),
            ("protein_length", integer),
            ("evalue", floating),
            ("score", floating),
            ("truncated", string),
            ("source", string),
            ("taxonomy", string),
        ]
    )
    cas13 = pa.schema(
        [
            ("operon_id", string),
            ("subtype", string),
            ("subtype_raw", string),
            ("subtype_source", string),
            ("subtype_conflict", boolean),
            ("effector_name", string),
            ("hmm_name", string),
            ("protein_sequence", string),
            ("sequence_sha256", string),
            ("protein_length", integer),
            ("source_length_field", integer),
            ("evalue", floating),
            ("score", floating),
            ("truncated", string),
            ("source", string),
            ("taxonomy", string),
            ("metadata_json", string),
        ]
    )
    pair = pa.schema(
        [
            ("operon_id", string),
            ("subtype", string),
            ("subtype_raw", string),
            ("subtype_source", string),
            ("subtype_conflict", boolean),
            ("effector_name", string),
            ("protein_sequence", string),
            ("direct_repeat", string),
            ("direct_repeat_raw", string),
            ("spacers_json", string),
            ("orientation", string),
            ("orientation_source", string),
            ("pairing_confidence", string),
            ("ambiguity_reason", string),
        ]
    )
    return {
        "atlas_operons": operon,
        "type_vi_operons": operon,
        "cas_effectors": cas_effector,
        "cas13_records": cas13,
        "crispr_arrays": pa.schema(
            [
                ("operon_id", string),
                ("subtype", string),
                ("array_index", integer),
                ("direct_repeat_raw", string),
                ("direct_repeat", string),
                ("spacers_json", string),
                ("spacer_count", integer),
                ("orientation", string),
                ("orientation_source", string),
            ]
        ),
        "cas13_direct_repeat_pairs": pair,
        "ambiguous_pairs": pair,
        "processing_failures": pa.schema(
            [
                ("record_index", integer),
                ("operon_id", string),
                ("error_type", string),
                ("error", string),
            ]
        ),
        "cas13_exact_unique": pa.schema(
            [
                ("sequence_sha256", string),
                ("protein_sequence", string),
                ("protein_length", integer),
                ("representative_operon_id", string),
                ("record_count", integer),
                ("nonconflicting_record_count", integer),
                ("complete_record_count", integer),
                ("minimum_evalue", floating),
                ("maximum_score", floating),
                ("truncated_flags", pa.list_(string)),
                ("subtypes", pa.list_(string)),
                ("operon_ids", pa.list_(string)),
            ]
        ),
    }


def _sha256_path(path: Path, chunk_size: int = 8 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def _operon_row(raw: dict[str, Any]) -> dict[str, Any]:
    operon_id = str(raw.get("operon_id") or "").strip()
    if not operon_id:
        raise ValueError("operon is missing operon_id")
    summary = raw.get("summary") or {}
    metadata = raw.get("metadata") or {}
    subtype = str(summary.get("subtype") or "").strip()
    arrays = raw.get("crispr") or []
    cas_entries = raw.get("cas") or []
    inferred_type_vi = not subtype and any(
        resolve_cas13_subtype(subtype, cas)[0].startswith("VI-")
        for cas in cas_entries
        if isinstance(cas, dict) and is_cas13_effector(cas)
    )
    genomic_context = {
        key: raw[key]
        for key in (
            "contig",
            "contig_id",
            "start",
            "end",
            "strand",
            "genes",
            "genomic_context",
        )
        if key in raw
    }
    cas_annotations = [
        {key: value for key, value in cas.items() if key != "protein"}
        for cas in cas_entries
        if isinstance(cas, dict)
    ]
    return {
        "operon_id": operon_id,
        "subtype": subtype,
        "is_type_vi": (
            subtype.upper() == "VI"
            or subtype.upper().startswith("VI-")
            or inferred_type_vi
        ),
        "metadata_json": json.dumps(metadata, sort_keys=True),
        "source": (
            str(metadata["source_db"])
            if metadata.get("source_db") is not None
            else None
        ),
        "taxonomy": (
            str(metadata["taxonomy"]) if metadata.get("taxonomy") is not None else None
        ),
        "biome": str(metadata["biome"]) if metadata.get("biome") is not None else None,
        "assembly_type": (
            str(metadata["assembly_type"])
            if metadata.get("assembly_type") is not None
            else None
        ),
        "n_cas": len(cas_entries),
        "n_crispr": len(arrays),
        "n_spacers": sum(
            len(array.get("crispr_spacers") or [])
            for array in arrays
            if isinstance(array, dict)
        ),
        "crispr_arrays_json": json.dumps(arrays, sort_keys=True),
        "cas_annotations_json": json.dumps(cas_annotations, sort_keys=True),
        "genomic_context_json": json.dumps(genomic_context, sort_keys=True),
        "genomic_context_available": bool(genomic_context),
    }


def _write_exact_unique(
    connection: sqlite3.Connection,
    writer: _ParquetBatchWriter,
) -> int:
    cursor = connection.execute(
        """
        SELECT sequence_sha256, protein_sequence, protein_length, operon_id, subtype,
               subtype_conflict, truncated, evalue, score
        FROM cas13_dedup
        ORDER BY sequence_sha256, operon_id
        """
    )
    current_digest: str | None = None
    sequence = ""
    protein_length = 0
    operon_ids: list[str] = []
    subtypes: set[str] = set()
    nonconflicting_record_count = 0
    complete_record_count = 0
    evalues: list[float] = []
    scores: list[float] = []
    truncated_flags: set[str] = set()
    unique_count = 0

    def flush_group() -> None:
        nonlocal unique_count
        if current_digest is None:
            return
        writer.append(
            {
                "sequence_sha256": current_digest,
                "protein_sequence": sequence,
                "protein_length": protein_length,
                "representative_operon_id": operon_ids[0],
                "record_count": len(operon_ids),
                "nonconflicting_record_count": nonconflicting_record_count,
                "complete_record_count": complete_record_count,
                "minimum_evalue": min(evalues) if evalues else None,
                "maximum_score": max(scores) if scores else None,
                "truncated_flags": sorted(truncated_flags),
                "subtypes": sorted(subtypes),
                "operon_ids": operon_ids,
            }
        )
        unique_count += 1

    for (
        digest,
        row_sequence,
        row_length,
        operon_id,
        subtype,
        subtype_conflict,
        truncated,
        evalue,
        score,
    ) in cursor:
        digest_string = str(digest)
        if current_digest is not None and digest_string != current_digest:
            flush_group()
            operon_ids = []
            subtypes = set()
            nonconflicting_record_count = 0
            complete_record_count = 0
            evalues = []
            scores = []
            truncated_flags = set()
        if digest_string != current_digest:
            current_digest = digest_string
            sequence = str(row_sequence)
            protein_length = int(row_length)
        operon_ids.append(str(operon_id))
        subtypes.add(str(subtype))
        nonconflicting_record_count += int(not bool(subtype_conflict))
        complete_record_count += int(str(truncated) == "00")
        if truncated is not None:
            truncated_flags.add(str(truncated))
        if evalue is not None:
            evalues.append(float(evalue))
        if score is not None:
            scores.append(float(score))
    flush_group()
    return unique_count


def process_atlas(
    path: Path,
    output_dir: Path,
    *,
    batch_size: int = 10_000,
) -> dict[str, Any]:
    """Stream Atlas into atomic, audited Parquet tables with bounded memory."""
    if not path.is_file():
        raise FileNotFoundError(f"Atlas input is missing: {path}")
    if output_dir.exists():
        raise FileExistsError(f"refusing to overwrite Atlas output: {output_dir}")
    output_dir.parent.mkdir(parents=True, exist_ok=True)
    staging = output_dir.with_name(f".{output_dir.name}.building-{os.getpid()}")
    staging.mkdir(parents=False, exist_ok=False)
    schemas = _schemas()
    writers = {
        name: _ParquetBatchWriter(
            staging / f"{name}.parquet",
            schema,
            batch_size=batch_size,
        )
        for name, schema in schemas.items()
    }
    sqlite_path = staging / "cas13_exact_dedup.sqlite"
    connection = sqlite3.connect(sqlite_path)
    connection.execute(
        """
        CREATE TABLE cas13_dedup (
            sequence_sha256 TEXT NOT NULL,
            protein_sequence TEXT NOT NULL,
            protein_length INTEGER NOT NULL,
            operon_id TEXT NOT NULL,
            subtype TEXT NOT NULL,
            subtype_conflict INTEGER NOT NULL,
            truncated TEXT,
            evalue REAL,
            score REAL
        )
        """
    )
    connection.execute(
        "CREATE INDEX cas13_digest_index ON cas13_dedup(sequence_sha256, operon_id)"
    )
    subtype_counts: Counter[str] = Counter()
    cas13_subtype_counts: Counter[str] = Counter()
    cas13_subtype_sources: Counter[str] = Counter()
    total = 0
    type_vi = 0
    cas_effector_count = 0
    cas13_count = 0
    high_pair_count = 0
    ambiguous_pair_count = 0
    failure_count = 0
    cas13_subtype_conflicts = 0
    cas13_nonconflicting_records = 0
    cas13_complete_records = 0
    try:
        for index, raw in enumerate(iter_json_array(path)):
            total += 1
            try:
                if not isinstance(raw, dict):
                    raise ValueError("operon record is not a mapping")
                operon = _operon_row(raw)
                cas_effectors = extract_cas_effectors(raw)
                cas13_records = extract_cas13_records(raw)
                arrays = extract_crispr_arrays(raw)
                pair = pair_cas13_direct_repeat(raw)

                subtype = str(operon["subtype"])
                subtype_counts[subtype] += 1
                is_type_vi = bool(operon["is_type_vi"])
                type_vi += int(is_type_vi)
                writers["atlas_operons"].append(operon)
                if is_type_vi:
                    writers["type_vi_operons"].append(operon)
                for cas_effector in cas_effectors:
                    writers["cas_effectors"].append(asdict(cas_effector))
                    cas_effector_count += 1
                for crispr_array in arrays:
                    writers["crispr_arrays"].append(asdict(crispr_array))
                for cas13_record in cas13_records:
                    writers["cas13_records"].append(asdict(cas13_record))
                    cas13_subtype_counts[cas13_record.subtype] += 1
                    cas13_subtype_sources[cas13_record.subtype_source] += 1
                    cas13_subtype_conflicts += int(cas13_record.subtype_conflict)
                    cas13_nonconflicting_records += int(
                        not cas13_record.subtype_conflict
                    )
                    cas13_complete_records += int(cas13_record.truncated == "00")
                    connection.execute(
                        "INSERT INTO cas13_dedup VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?)",
                        (
                            cas13_record.sequence_sha256,
                            cas13_record.protein_sequence,
                            cas13_record.protein_length,
                            cas13_record.operon_id,
                            cas13_record.subtype,
                            int(cas13_record.subtype_conflict),
                            cas13_record.truncated,
                            cas13_record.evalue,
                            cas13_record.score,
                        ),
                    )
                    cas13_count += 1
                if pair is not None:
                    if pair.pairing_confidence == "high":
                        writers["cas13_direct_repeat_pairs"].append(asdict(pair))
                        high_pair_count += 1
                    else:
                        writers["ambiguous_pairs"].append(asdict(pair))
                        ambiguous_pair_count += 1
                if total % batch_size == 0:
                    connection.commit()
            except (TypeError, ValueError, KeyError, UnicodeError) as exc:
                failure_count += 1
                writers["processing_failures"].append(
                    {
                        "record_index": index,
                        "operon_id": (
                            raw.get("operon_id") if isinstance(raw, dict) else None
                        ),
                        "error_type": type(exc).__name__,
                        "error": str(exc),
                    }
                )
        connection.commit()
        exact_unique_count = _write_exact_unique(
            connection, writers["cas13_exact_unique"]
        )
        evolution_eligible_exact_unique = int(
            connection.execute(
                """
                SELECT COUNT(DISTINCT sequence_sha256)
                FROM cas13_dedup
                WHERE subtype_conflict = 0 AND truncated = '00'
                """
            ).fetchone()[0]
        )
        for writer in writers.values():
            writer.close()
        connection.close()
    except BaseException:
        connection.close()
        for writer in writers.values():
            try:
                writer.close()
            except Exception:
                pass
        raise
    funnel = {
        "is_mock": False,
        "evidence_level": 0,
        "atlas_version": "1.0",
        "input_path": str(path),
        "input_size_bytes": path.stat().st_size,
        "input_sha256": _sha256_path(path),
        "atlas_operons": total,
        "type_vi_operons": type_vi,
        "cas_effectors": cas_effector_count,
        "cas13_records": cas13_count,
        "cas13_exact_unique": exact_unique_count,
        "high_confidence_pairs": high_pair_count,
        "ambiguous_pairs": ambiguous_pair_count,
        "processing_failures": failure_count,
        "subtype_counts": dict(sorted(subtype_counts.items())),
        "cas13_subtype_counts": dict(sorted(cas13_subtype_counts.items())),
        "cas13_subtype_sources": dict(sorted(cas13_subtype_sources.items())),
        "cas13_subtype_conflicts": cas13_subtype_conflicts,
        "cas13_nonconflicting_records": cas13_nonconflicting_records,
        "cas13_complete_records": cas13_complete_records,
        "cas13_evolution_eligible_exact_unique": evolution_eligible_exact_unique,
        "pair_orientation_policy": (
            "unknown orientation is ambiguous; reverse is reverse-complemented; "
            "only declared/recovered orientation may be high confidence"
        ),
    }
    (staging / "data_funnel.json").write_text(
        json.dumps(funnel, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (staging / "data_card.md").write_text(_render_data_card(funnel), encoding="utf-8")
    staging.replace(output_dir)
    return funnel


def _render_data_card(funnel: dict[str, Any]) -> str:
    lines = [
        "# Processed Atlas data card",
        "",
        "This card was generated from a real declared Atlas input (`is_mock=false`).",
        "Counts are records unless explicitly labeled exact-unique; records are not",
        "assumed independent.",
        "",
        "```json",
        json.dumps(funnel, indent=2, sort_keys=True),
        "```",
        "",
        "Only one-effector/one-array/nonempty-repeat/unambiguous-subtype records",
        "with declared or reliably recovered orientation enter the high-confidence",
        "paired table. Atlas v1.0 records without orientation are retained in",
        "`ambiguous_pairs.parquet`; they are never silently used for coevolution.",
        "",
        "The raw Cas13 and exact-unique tables retain subtype conflicts and",
        "truncation flags. Evolutionary MSA eligibility requires at least one",
        "nonconflicting record explicitly marked `truncated=00`; this eligibility",
        "does not imply functional validation.",
    ]
    return "\n".join(lines) + "\n"
