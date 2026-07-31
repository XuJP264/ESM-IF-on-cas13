import json
from pathlib import Path

import numpy as np
import pytest

from cas13_if.refold.level3 import (
    contact_recovery,
    domain_rmsd,
    hepn_geometry,
    interface_confidence,
    interface_pae,
    load_pae,
    pareto_front,
    validate_level3_result,
)


def test_level3_pae_geometry_contacts_and_confidence(tmp_path: Path) -> None:
    pae_path = tmp_path / "pae.json"
    pae_path.write_text(
        json.dumps({"predicted_aligned_error": [[0.0, 4.0], [6.0, 0.0]]}),
        encoding="utf-8",
    )
    pae = load_pae(pae_path)
    assert interface_pae(pae, {0}, {1}) == pytest.approx(5.0)
    reference = np.asarray([[0, 0, 0], [1, 0, 0], [0, 1, 0]], dtype=float)
    mobile = reference + np.asarray([4, 3, 2])
    assert domain_rmsd(mobile, reference, {"HEPN": [0, 1, 2]})["HEPN"] < 1e-10
    geometry = hepn_geometry({0: reference[0], 1: reference[1]}, [(0, 1)])
    assert geometry["0-1"] == pytest.approx(1.0)
    assert contact_recovery({1, 2}, {2, 3})["recall"] == pytest.approx(0.5)
    assert interface_confidence([70.0, 90.0], {0, 1}) == pytest.approx(80.0)


def test_level3_mock_gate_and_pareto() -> None:
    result = {
        "candidate_id": "a",
        "provider": "boltz",
        "mean_plddt": 80,
        "structure_path": "prediction.cif",
        "pae_path": "pae.json",
        "seed": 1,
        "is_mock": True,
    }
    validate_level3_result(result, expected_mock=True)
    with pytest.raises(ValueError, match="mock state"):
        validate_level3_result(result, expected_mock=False)
    rows = [
        {"candidate_id": "a", "novelty": 0.8, "rmsd": 2.0},
        {"candidate_id": "b", "novelty": 0.7, "rmsd": 3.0},
        {"candidate_id": "c", "novelty": 0.9, "rmsd": 4.0},
    ]
    assert pareto_front(rows, maximize=["novelty"], minimize=["rmsd"]) == ["a", "c"]
