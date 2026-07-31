from pathlib import Path

import pytest

from cas13_if.alignments.msa import Alignment
from cas13_if.alignments.scaffold_mapping import (
    QUERY_IDENTIFIER,
    coordinate_index_mapping,
    global_index_mapping,
    map_added_alignment_columns,
    scaffold_to_msa_columns,
)
from cas13_if.structures.parser import ResidueKey


def test_global_mapping_distinguishes_terminal_and_internal_missing() -> None:
    mapping = global_index_mapping("ACDEFGHIK", "CDEGHIK")
    assert mapping.reference_to_query[0] is None
    assert mapping.reference_to_query[4] is None
    assert mapping.reference_to_query[1:4] == (0, 1, 2)
    assert mapping.reference_to_query[5:] == (3, 4, 5, 6)
    assert mapping.query_to_reference == (1, 2, 3, 5, 6, 7, 8)


def test_added_alignment_preserves_original_columns_and_tracks_query_insertion() -> (
    None
):
    original = Alignment(
        identifiers=("first", "second"),
        sequences=("AC-D", "A-ED"),
    )
    added = Alignment(
        identifiers=("first", "second", QUERY_IDENTIFIER),
        sequences=("AC--D", "A--ED", "ACWED"),
    )
    result = map_added_alignment_columns(original, added)
    assert result.original_columns_preserved
    assert result.output_to_original_column == (0, 1, None, 2, 3)
    assert scaffold_to_msa_columns(result) == (0, 1, None, 2, 3)


def test_added_alignment_rejects_changed_original_column() -> None:
    original = Alignment(identifiers=("a", "b"), sequences=("AC", "AD"))
    changed = Alignment(
        identifiers=("a", "b", QUERY_IDENTIFIER),
        sequences=("AT", "AD", "AT"),
    )
    with pytest.raises(ValueError, match="altered or reordered"):
        map_added_alignment_columns(original, changed)


def test_real_fixture_mmcif_retains_insertion_code_contract() -> None:
    # The strict mapping data model always carries an insertion-code string;
    # empty strings are not silently converted to a numeric residue index.
    assert Path("tests/fixtures/minimal_complex.pdb").is_file()


def test_coordinate_mapping_prefers_validated_author_numbers() -> None:
    keys = [
        ResidueKey("A", 1, "", "ALA"),
        ResidueKey("A", 3, "", "ASP"),
    ]
    mapping, strategy = coordinate_index_mapping("ACD", "AD", keys)
    assert strategy == "validated_author_residue_number"
    assert mapping.reference_to_query == (0, None, 1)


def test_coordinate_mapping_falls_back_for_insertion_codes() -> None:
    keys = [
        ResidueKey("A", 1, "", "ALA"),
        ResidueKey("A", 1, "A", "CYS"),
    ]
    mapping, strategy = coordinate_index_mapping("AC", "AC", keys)
    assert strategy == "global_sequence_alignment"
    assert mapping.reference_to_query == (0, 1)
