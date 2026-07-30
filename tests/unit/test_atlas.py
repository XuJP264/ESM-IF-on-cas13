import json
from pathlib import Path

from cas13_if.data.atlas import (
    exact_deduplicate,
    extract_cas13_records,
    iter_json_array,
    pair_cas13_direct_repeat,
    process_atlas,
)

FIXTURE = Path("data/fixtures/atlas_operons.json")


def test_streaming_parser_with_small_chunks() -> None:
    records = list(iter_json_array(FIXTURE, chunk_size=128))
    assert len(records) == 2
    assert records[0]["operon_id"] == "fixture-VI-D-001"


def test_pairing_routes_ambiguity() -> None:
    records = list(iter_json_array(FIXTURE, chunk_size=128))
    high = pair_cas13_direct_repeat(records[0])
    ambiguous = pair_cas13_direct_repeat(records[1])
    assert high is not None and high.pairing_confidence == "high"
    assert high.direct_repeat == "GAAACACCGUUGAAAGUG"
    assert ambiguous is not None
    assert ambiguous.pairing_confidence == "ambiguous"
    assert ambiguous.ambiguity_reason == "array_count=2"


def test_exact_dedup_and_fixture_pipeline(tmp_path: Path) -> None:
    raw = list(iter_json_array(FIXTURE, chunk_size=128))
    record = extract_cas13_records(raw[0])[0]
    unique = exact_deduplicate([record, record])
    assert unique[0]["record_count"] == 2
    funnel = process_atlas(FIXTURE, tmp_path / "processed")
    assert funnel["atlas_operons"] == 2
    assert funnel["cas13_exact_unique"] == 2
    assert funnel["high_confidence_pairs"] == 1
    assert funnel["ambiguous_pairs"] == 1
    saved = json.loads(
        (tmp_path / "processed/data_funnel.json").read_text(encoding="utf-8")
    )
    assert saved == funnel
