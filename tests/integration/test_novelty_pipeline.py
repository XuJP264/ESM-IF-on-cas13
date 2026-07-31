import json
from pathlib import Path

import pyarrow.parquet as pq

from cas13_if.novelty.pipeline import (
    NoveltyThresholds,
    run_candidate_novelty_pipeline,
)


def test_candidate_jsonl_to_mmseqs_novelty_report(tmp_path: Path) -> None:
    candidates = tmp_path / "candidates.jsonl"
    payloads = [
        {
            "pdb_id": "fixture",
            "method": "fixture",
            "candidate": {
                "candidate_id": "candidate-a",
                "scaffold_id": "fixture-A",
                "sequence": "ACDEFGHIK",
                "parent_sequence": "ACDEFGHIL",
                "temperature": 1.0,
                "seed": 7,
                "fixed_positions": {"0": "A"},
                "is_mock": False,
                "evidence_level": 2,
                "traces": [{"large": "field ignored by compact reader"}],
            },
            "recovery": {"fixed_position_violations": 0},
        },
        {
            "pdb_id": "fixture",
            "method": "fixture",
            "candidate": {
                "candidate_id": "candidate-b",
                "scaffold_id": "fixture-A",
                "sequence": "ACDEFGHIL",
                "parent_sequence": "ACDEFGHIK",
                "temperature": 1.0,
                "seed": 8,
                "fixed_positions": {"0": "A"},
                "is_mock": False,
                "evidence_level": 2,
            },
            "recovery": {"fixed_position_violations": 0},
        },
    ]
    candidates.write_text(
        "".join(json.dumps(payload) + "\n" for payload in payloads),
        encoding="utf-8",
    )
    atlas = tmp_path / "atlas.fasta"
    atlas.write_text(">atlas-a\nACDEFGHIK\n>atlas-b\nACDEFGHIL\n", encoding="utf-8")
    fake_mmseqs = tmp_path / "mmseqs"
    fake_mmseqs.write_text(
        "#!/bin/sh\n"
        "set -eu\n"
        'test "$1" = "easy-search"\n'
        "printf '%s\\n' "
        "'candidate-a\tatlas-a\t0.5\t9\t1.0\t1.0\t1e-5\t50' "
        "'candidate-b\tatlas-b\t0.6\t9\t1.0\t1.0\t1e-6\t60' "
        '> "$4"\n',
        encoding="utf-8",
    )
    fake_mmseqs.chmod(0o755)

    output = tmp_path / "novelty"
    summary = run_candidate_novelty_pipeline(
        candidate_jsonl=candidates,
        atlas_fasta=atlas,
        output_dir=output,
        executable=fake_mmseqs,
        threads=1,
        sensitivity=7.5,
        minimum_query_coverage=0.8,
        maximum_evalue=1000.0,
        maximum_sequences=10,
        thresholds=NoveltyThresholds(
            maximum_parent_identity=0.95,
            maximum_atlas_identity=0.8,
            maximum_homopolymer_length=4,
            maximum_low_complexity_windows=0,
            minimum_designed_position_entropy=0.0,
            low_complexity_window=4,
            low_complexity_maximum_fraction=0.5,
        ),
    )

    assert summary["candidate_count"] == 2
    assert summary["passes_level1_novelty"] == 2
    assert summary["evidence_level_max"] == 1
    assert summary["atlas_fasta_sha256"]
    assert (output / "summary.json").is_file()
    assert (output / "atlas_alignments.tsv").is_file()
    table = pq.read_table(output / "candidate_novelty.parquet")
    assert table.num_rows == 2
    assert table.column("maximum_atlas_identity").to_pylist() == [0.5, 0.6]
