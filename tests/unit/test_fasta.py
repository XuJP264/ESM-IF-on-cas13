from pathlib import Path

import pytest

from cas13_if.data.fasta import FastaError, iter_fasta, write_fasta


def test_fasta_roundtrip(tmp_path: Path) -> None:
    path = tmp_path / "sequences.fasta"
    write_fasta([("one", "ACDE"), ("two", "FGHI")], path, line_width=2)
    assert list(iter_fasta(path)) == [("one", "ACDE"), ("two", "FGHI")]


def test_fasta_rejects_malformed(tmp_path: Path) -> None:
    path = tmp_path / "bad.fasta"
    path.write_text("ACDE\n", encoding="utf-8")
    with pytest.raises(FastaError, match="before header"):
        list(iter_fasta(path))
    with pytest.raises(FastaError, match="duplicate"):
        write_fasta([("same", "AC"), ("same", "DE")], tmp_path / "dup.fa")
