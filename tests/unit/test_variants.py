import pytest

from cas13_if.data.variants import activity_label, apply_mutation


def test_apply_substitution_deletion_and_insertion() -> None:
    assert apply_mutation("ACDE", "A1V,E4F") == "VCDF"
    assert apply_mutation("ACDEFG", "del2-3+del5-5") == "AEG"
    assert apply_mutation("ACDE", "ins2_3:GG") == "ACGGDE"


def test_apply_mutation_rejects_wrong_wild_type() -> None:
    with pytest.raises(ValueError, match="WT token mismatch"):
        apply_mutation("ACDE", "G2A")


def test_activity_labels_do_not_pool_endpoints() -> None:
    assert (
        activity_label(
            0.6,
            0.1,
            active_minimum=0.8,
            inactive_maximum=0.2,
            cis_retained_minimum=0.5,
            trans_reduced_maximum=0.2,
        )
        == "cis-retained/trans-reduced"
    )
    assert (
        activity_label(
            0.1,
            0.1,
            active_minimum=0.8,
            inactive_maximum=0.2,
            cis_retained_minimum=0.5,
            trans_reduced_maximum=0.2,
        )
        == "inactive"
    )
