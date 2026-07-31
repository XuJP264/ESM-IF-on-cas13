from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from cas13_if.alignments.pipeline import build_subtype_msas
from cas13_if.evolution.pipeline import compute_subtype_conservation


def test_subtype_msa_to_conservation_pipeline(tmp_path: Path) -> None:
    exact = tmp_path / "exact.parquet"
    mapping = tmp_path / "mapping.parquet"
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "sequence_sha256": "seq-a",
                    "protein_sequence": "ACDE",
                    "protein_length": 4,
                    "subtypes": ["VI-D"],
                },
                {
                    "sequence_sha256": "seq-b",
                    "protein_sequence": "ACDF",
                    "protein_length": 4,
                    "subtypes": ["VI-D"],
                },
                {
                    "sequence_sha256": "seq-invalid",
                    "protein_sequence": "ACDX",
                    "protein_length": 4,
                    "subtypes": ["VI-D"],
                },
            ]
        ),
        exact,
    )
    pq.write_table(
        pa.Table.from_pylist(
            [
                {
                    "sequence_sha256": name,
                    "representative_sha256": name,
                }
                for name in ("seq-a", "seq-b", "seq-invalid")
            ]
        ),
        mapping,
    )
    fake_mafft = tmp_path / "mafft"
    fake_mafft.write_text(
        "#!/bin/sh\n"
        'if [ "$1" = "--version" ]; then echo "v7.fixture" >&2; exit 0; fi\n'
        'for argument in "$@"; do input="$argument"; done\n'
        'sed "s/^/>/" /dev/null >/dev/null\n'
        'exec /bin/cp "$input" /dev/stdout\n',
        encoding="utf-8",
    )
    fake_mafft.chmod(0o755)
    msa_manifest = build_subtype_msas(
        exact_unique_path=exact,
        cluster_mapping_path=mapping,
        output_dir=tmp_path / "msa",
        executable=str(fake_mafft),
        threads=1,
    )
    assert msa_manifest["subtypes"]["VI-D"]["status"] == "success"
    assert msa_manifest["excluded_sequence_count"] == 1

    conservation = compute_subtype_conservation(
        msa_root=tmp_path / "msa",
        output_dir=tmp_path / "conservation",
        identity_threshold=0.8,
        allowed_frequency=0.05,
    )
    assert conservation["subtypes"]["VI-D"]["columns"] == 4
    assert (tmp_path / "conservation/vi-d.parquet").is_file()
