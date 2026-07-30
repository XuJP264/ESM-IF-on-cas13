"""Streaming CRISPR-Cas Atlas normalization and conservative Cas13 pairing."""

from __future__ import annotations

import hashlib
import json
from collections import Counter
from collections.abc import Iterator
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, TextIO

import pandas as pd


class AtlasParseError(ValueError):
    """Raised when the top-level Atlas JSON stream is malformed."""


@dataclass(frozen=True)
class Cas13Record:
    operon_id: str
    subtype: str
    effector_name: str
    hmm_name: str | None
    protein_sequence: str
    sequence_sha256: str
    protein_length: int
    source: str | None
    taxonomy: str | None
    metadata_json: str


@dataclass(frozen=True)
class PairingRecord:
    operon_id: str
    subtype: str
    effector_name: str
    protein_sequence: str
    direct_repeat: str
    spacers_json: str
    orientation: str
    pairing_confidence: str
    ambiguity_reason: str | None


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
            buffer = buffer[position:] + handle.read(chunk_size)
            position = 0
            eof = len(buffer) == 0 or len(buffer) < chunk_size
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
            buffer = buffer[position:] + handle.read(chunk_size)
            position = 0
            eof = len(buffer) < chunk_size
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


def normalize_rna(sequence: str) -> str:
    return "".join(sequence.split()).upper().replace("T", "U")


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
        records.append(
            Cas13Record(
                operon_id=operon_id,
                subtype=subtype,
                effector_name=name,
                hmm_name=(
                    str(cas["hmm_name"]) if cas.get("hmm_name") is not None else None
                ),
                protein_sequence=sequence,
                sequence_sha256=hashlib.sha256(sequence.encode("ascii")).hexdigest(),
                protein_length=len(sequence),
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
    arrays = [item for item in (operon.get("crispr") or []) if isinstance(item, dict)]
    if not effectors:
        return None
    reasons: list[str] = []
    if len(effectors) != 1:
        reasons.append(f"effector_count={len(effectors)}")
    if len(arrays) != 1:
        reasons.append(f"array_count={len(arrays)}")
    repeat = normalize_rna(str(arrays[0].get("crispr_repeat") or "")) if arrays else ""
    if not repeat:
        reasons.append("empty_direct_repeat")
    subtype = effectors[0].subtype
    if not subtype or not subtype.upper().startswith("VI-"):
        reasons.append("ambiguous_subtype")
    orientation = (
        str(arrays[0].get("orientation") or "unknown") if arrays else "unknown"
    )
    if orientation not in {"forward", "reverse", "recovered", "unknown"}:
        reasons.append("invalid_orientation")
    confidence = "high" if not reasons else "ambiguous"
    spacers = arrays[0].get("crispr_spacers") or [] if arrays else []
    return PairingRecord(
        operon_id=effectors[0].operon_id,
        subtype=subtype,
        effector_name=effectors[0].effector_name,
        protein_sequence=effectors[0].protein_sequence,
        direct_repeat=repeat,
        spacers_json=json.dumps(spacers),
        orientation=orientation,
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
        output.append(
            {
                "sequence_sha256": digest,
                "protein_sequence": representative.protein_sequence,
                "protein_length": representative.protein_length,
                "representative_operon_id": representative.operon_id,
                "record_count": len(members),
                "subtypes": sorted({item.subtype for item in members}),
                "operon_ids": sorted(item.operon_id for item in members),
            }
        )
    return output


def process_atlas(path: Path, output_dir: Path) -> dict[str, Any]:
    """Stream Atlas into audited tables; Parquet output is mandatory."""
    output_dir.mkdir(parents=True, exist_ok=True)
    operons: list[dict[str, Any]] = []
    cas13_records: list[Cas13Record] = []
    high_pairs: list[PairingRecord] = []
    ambiguous_pairs: list[PairingRecord] = []
    failures: list[dict[str, Any]] = []
    subtype_counts: Counter[str] = Counter()
    total = 0
    type_vi = 0
    for index, raw in enumerate(iter_json_array(path)):
        total += 1
        try:
            if not isinstance(raw, dict):
                raise ValueError("operon record is not a mapping")
            summary = raw.get("summary") or {}
            subtype = str(summary.get("subtype") or "")
            subtype_counts[subtype] += 1
            is_type_vi = subtype.upper().startswith("VI-")
            type_vi += int(is_type_vi)
            operons.append(
                {
                    "operon_id": raw.get("operon_id"),
                    "subtype": subtype,
                    "is_type_vi": is_type_vi,
                    "metadata_json": json.dumps(
                        raw.get("metadata") or {}, sort_keys=True
                    ),
                    "n_cas": len(raw.get("cas") or []),
                    "n_crispr": len(raw.get("crispr") or []),
                }
            )
            extracted = extract_cas13_records(raw)
            cas13_records.extend(extracted)
            pair = pair_cas13_direct_repeat(raw)
            if pair is not None:
                target_pairs = (
                    high_pairs if pair.pairing_confidence == "high" else ambiguous_pairs
                )
                target_pairs.append(pair)
        except (TypeError, ValueError, KeyError) as exc:
            failures.append(
                {
                    "record_index": index,
                    "operon_id": (
                        raw.get("operon_id") if isinstance(raw, dict) else None
                    ),
                    "error": str(exc),
                }
            )
    exact_unique = exact_deduplicate(cas13_records)
    tables = {
        "atlas_operons": operons,
        "cas13_records": [asdict(item) for item in cas13_records],
        "cas13_exact_unique": exact_unique,
        "cas13_direct_repeat_pairs": [asdict(item) for item in high_pairs],
        "ambiguous_pairs": [asdict(item) for item in ambiguous_pairs],
        "processing_failures": failures,
    }
    for name, rows in tables.items():
        pd.DataFrame(rows).to_parquet(output_dir / f"{name}.parquet", index=False)
    pd.DataFrame([row for row in operons if row["is_type_vi"]]).to_parquet(
        output_dir / "type_vi_operons.parquet", index=False
    )
    funnel = {
        "is_mock": False,
        "atlas_operons": total,
        "type_vi_operons": type_vi,
        "cas13_records": len(cas13_records),
        "cas13_exact_unique": len(exact_unique),
        "high_confidence_pairs": len(high_pairs),
        "ambiguous_pairs": len(ambiguous_pairs),
        "processing_failures": len(failures),
        "subtype_counts": dict(sorted(subtype_counts.items())),
    }
    (output_dir / "data_funnel.json").write_text(
        json.dumps(funnel, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    (output_dir / "data_card.md").write_text(
        _render_data_card(funnel), encoding="utf-8"
    )
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
        "enter the high-confidence paired table. Ambiguity is retained separately.",
    ]
    return "\n".join(lines) + "\n"
