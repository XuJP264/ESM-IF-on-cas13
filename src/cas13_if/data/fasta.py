"""Strict FASTA I/O with deterministic identifiers."""

from __future__ import annotations

from collections.abc import Iterable, Iterator
from pathlib import Path


class FastaError(ValueError):
    """Raised for malformed FASTA input."""


def iter_fasta(path: Path) -> Iterator[tuple[str, str]]:
    identifier: str | None = None
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if identifier is not None:
                    yield identifier, "".join(chunks).upper()
                identifier = line[1:].split(maxsplit=1)[0]
                if not identifier:
                    raise FastaError(f"empty FASTA identifier at line {line_number}")
                chunks = []
            elif identifier is None:
                raise FastaError(f"sequence before header at line {line_number}")
            else:
                chunks.append("".join(line.split()))
    if identifier is not None:
        yield identifier, "".join(chunks).upper()
    elif not chunks:
        raise FastaError(f"no FASTA records in {path}")


def write_fasta(
    records: Iterable[tuple[str, str]],
    path: Path,
    *,
    line_width: int = 80,
) -> None:
    if line_width < 1:
        raise ValueError("line_width must be positive")
    lines: list[str] = []
    seen: set[str] = set()
    for identifier, sequence in records:
        clean_id = identifier.strip().split(maxsplit=1)[0]
        clean_sequence = "".join(sequence.split()).upper()
        if not clean_id or clean_id in seen:
            raise FastaError(f"empty or duplicate FASTA identifier: {clean_id!r}")
        if not clean_sequence:
            raise FastaError(f"empty FASTA sequence: {clean_id}")
        seen.add(clean_id)
        lines.append(f">{clean_id}")
        lines.extend(
            clean_sequence[index : index + line_width]
            for index in range(0, len(clean_sequence), line_width)
        )
    if not lines:
        raise FastaError("cannot write an empty FASTA")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
