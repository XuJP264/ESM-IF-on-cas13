import json
from pathlib import Path

from cas13_if.data.atlas import (
    exact_deduplicate,
    extract_cas13_records,
    extract_crispr_arrays,
    iter_json_array,
    pair_cas13_direct_repeat,
    process_atlas,
    resolve_cas13_subtype,
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


def test_cas13_subtype_recovery_and_conflict_are_explicit(tmp_path: Path) -> None:
    explicit_hmm = {
        "gene_name": "Cas13d",
        "hmm_name": "Cas13d_0_CAS-VI-D",
    }
    assert resolve_cas13_subtype("VI", explicit_hmm) == (
        "VI-D",
        "cas_hmm_explicit",
        False,
    )
    assert resolve_cas13_subtype(
        "VI", {"gene_name": "Cas13f", "hmm_name": "Cas13f_c126"}
    ) == ("VI-F", "effector_name_nomenclature", False)
    assert resolve_cas13_subtype("I-C", explicit_hmm) == (
        "VI-D",
        "cas_hmm_explicit_summary_conflict",
        True,
    )

    record = next(iter(iter_json_array(FIXTURE)))
    record["summary"]["subtype"] = "I-C"
    record["cas"][0]["hmm_name"] = "Cas13d_0_CAS-VI-D"
    cas13 = extract_cas13_records(record)[0]
    assert cas13.subtype == "VI-D"
    assert cas13.subtype_raw == "I-C"
    assert cas13.subtype_conflict
    pair = pair_cas13_direct_repeat(record)
    assert pair is not None
    assert pair.pairing_confidence == "ambiguous"
    assert pair.ambiguity_reason == "subtype_conflict_with_operon_summary"

    generic = next(iter(iter_json_array(FIXTURE)))
    generic["summary"]["subtype"] = "VI"
    generic["cas"][0]["hmm_name"] = "Cas13d_0_CAS-VI-D"
    generic_path = tmp_path / "generic-vi.json"
    generic_path.write_text(json.dumps([generic]), encoding="utf-8")
    funnel = process_atlas(generic_path, tmp_path / "generic-processed")
    assert funnel["type_vi_operons"] == 1
    assert funnel["cas13_subtype_counts"] == {"VI-D": 1}
    assert funnel["cas13_subtype_sources"] == {"cas_hmm_explicit": 1}


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
    assert funnel["cas13_subtype_counts"] == {"VI-A": 1, "VI-D": 1}
    assert funnel["cas13_subtype_conflicts"] == 0
    saved = json.loads(
        (tmp_path / "processed/data_funnel.json").read_text(encoding="utf-8")
    )
    assert saved == funnel
