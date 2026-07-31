import json
from pathlib import Path

from cas13_if.data.atlas import (
    exact_deduplicate,
    extract_cas13_records,
    extract_crispr_arrays,
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


def test_unknown_orientation_is_ambiguous_and_reverse_is_oriented() -> None:
    records = list(iter_json_array(FIXTURE, chunk_size=128))
    unknown = json.loads(json.dumps(records[0]))
    unknown["crispr"][0].pop("orientation")
    unknown["crispr"][0].pop("orientation_source")
    pair = pair_cas13_direct_repeat(unknown)
    assert pair is not None
    assert pair.pairing_confidence == "ambiguous"
    assert pair.ambiguity_reason == "orientation_not_recovered"
    assert pair.orientation_source == "not_provided_by_atlas_v1.0"

    reverse = json.loads(json.dumps(records[0]))
    reverse["crispr"][0]["orientation"] = "reverse"
    reverse["crispr"][0]["crispr_repeat"] = "AAGT"
    oriented = extract_crispr_arrays(reverse)[0]
    assert oriented.direct_repeat_raw == "AAGU"
    assert oriented.direct_repeat == "ACUU"


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
