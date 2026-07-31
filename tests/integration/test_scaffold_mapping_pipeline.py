import json
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from cas13_if.alignments.scaffold_mapping import (
    AddedAlignment,
    build_scaffold_mapping,
)


def test_fixture_scaffold_mapping_writes_confidence_audit(
    monkeypatch, tmp_path: Path
) -> None:
    entity = tmp_path / "entity.json"
    entity.write_text(
        json.dumps(
            {
                "entity_poly": {
                    "rcsb_entity_polymer_type": "Protein",
                    "pdbx_seq_one_letter_code_can": "MAG",
                }
            }
        ),
        encoding="utf-8",
    )
    msa = tmp_path / "msa.fasta"
    msa.write_text(">one\nMAG\n>two\nM-G\n", encoding="utf-8")
    conservation = tmp_path / "conservation.parquet"
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "column": index,
                    "coverage": coverage,
                    "conservation": 0.9,
                    "entropy": 0.1,
                    "gap_fraction": 1.0 - coverage,
                    "consensus": token,
                    "allowed_residues": [token],
                }
                for index, (coverage, token) in enumerate(
                    [(1.0, "M"), (0.5, "A"), (1.0, "G")]
                )
            ]
        ),
        conservation,
    )
    monkeypatch.setattr(
        "cas13_if.alignments.scaffold_mapping.add_scaffold_to_msa",
        lambda **_kwargs: AddedAlignment(
            query_aligned="MAG",
            output_to_original_column=(0, 1, 2),
            original_columns_preserved=True,
            mafft_command=("mafft", "--addfull"),
            mafft_stderr="",
        ),
    )
    output = tmp_path / "audit"
    summary = build_scaffold_mapping(
        structure_path=Path("tests/fixtures/minimal_complex.pdb"),
        entity_path=entity,
        msa_path=msa,
        conservation_path=conservation,
        output_dir=output,
        chain_id="A",
        subtype="VI-D",
        mafft_executable="mafft",
        threads=1,
        minimum_conservation_coverage=0.8,
    )
    assert summary["full_scaffold_length"] == 3
    assert summary["coordinate_length"] == 2
    assert summary["unresolved_positions"] == 1
    assert summary["original_msa_columns_preserved"]
    assert summary["conservation_constraint_eligible_positions"] == 0
    assert (output / "mapping.csv").is_file()
    assert (output / "manual_review.csv").is_file()
    assert "Evidence Level 0" in (output / "mapping.html").read_text()
