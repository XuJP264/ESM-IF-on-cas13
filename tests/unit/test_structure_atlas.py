from cas13_if.structures.atlas import aligned_identity


def test_aligned_identity_handles_indels() -> None:
    assert aligned_identity("ACDE", "ACE") == 1.0
    assert aligned_identity("ACDE", "ACNE") == 0.75
