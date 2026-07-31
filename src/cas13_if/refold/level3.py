"""Level-3 result validation, structural metrics, and Pareto utilities."""

from __future__ import annotations

import json
import math
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np
from numpy.typing import NDArray


@dataclass(frozen=True)
class AlignmentComparison:
    status: str
    tm_score: float | None
    rmsd: float | None
    executable: str
    is_mock: bool
    failure_reason: str | None = None


def validate_level3_result(result: dict[str, Any], *, expected_mock: bool) -> None:
    required = {
        "candidate_id",
        "provider",
        "mean_plddt",
        "structure_path",
        "pae_path",
        "seed",
        "is_mock",
    }
    missing = sorted(required.difference(result))
    if missing:
        raise ValueError(f"Level-3 result fields missing: {missing}")
    if bool(result["is_mock"]) != expected_mock:
        raise ValueError("Level-3 result mock state does not match ingest mode")
    plddt = float(result["mean_plddt"])
    if not 0 <= plddt <= 100:
        raise ValueError("mean pLDDT must be in [0, 100]")


def load_pae(path: Path) -> NDArray[np.float64]:
    value = json.loads(path.read_text(encoding="utf-8"))
    if isinstance(value, list) and value and isinstance(value[0], dict):
        value = value[0]
    if not isinstance(value, dict):
        raise ValueError("PAE file root must be a mapping")
    matrix = value.get("predicted_aligned_error", value.get("pae"))
    array = np.asarray(matrix, dtype=float)
    if array.ndim != 2 or array.shape[0] != array.shape[1] or array.size == 0:
        raise ValueError("PAE must be a non-empty square matrix")
    if not np.isfinite(array).all() or (array < 0).any():
        raise ValueError("PAE values must be finite and non-negative")
    return array


def interface_pae(pae: NDArray[np.float64], left: set[int], right: set[int]) -> float:
    if not left or not right:
        raise ValueError("interface PAE requires two non-empty position sets")
    if min(left | right) < 0 or max(left | right) >= pae.shape[0]:
        raise ValueError("interface position exceeds PAE dimensions")
    values = [float(pae[i, j]) for i in sorted(left) for j in sorted(right)]
    values.extend(float(pae[j, i]) for i in sorted(left) for j in sorted(right))
    return float(np.mean(values))


def kabsch_rmsd(mobile: NDArray[np.float64], reference: NDArray[np.float64]) -> float:
    if mobile.shape != reference.shape or mobile.ndim != 2 or mobile.shape[1] != 3:
        raise ValueError("coordinate arrays must share shape (n, 3)")
    if (
        mobile.shape[0] < 3
        or not np.isfinite(mobile).all()
        or not np.isfinite(reference).all()
    ):
        raise ValueError("RMSD requires at least three finite coordinate pairs")
    centered_mobile = mobile - mobile.mean(axis=0)
    centered_reference = reference - reference.mean(axis=0)
    covariance = centered_mobile.T @ centered_reference
    left, _, right = np.linalg.svd(covariance)
    sign = np.sign(np.linalg.det(left @ right))
    rotation = left @ np.diag([1.0, 1.0, sign]) @ right
    aligned = centered_mobile @ rotation
    return float(np.sqrt(np.mean(np.sum((aligned - centered_reference) ** 2, axis=1))))


def domain_rmsd(
    mobile: NDArray[np.float64],
    reference: NDArray[np.float64],
    domains: dict[str, list[int]],
) -> dict[str, float]:
    output: dict[str, float] = {}
    for name, indices in domains.items():
        if len(indices) < 3:
            raise ValueError(f"domain {name} has fewer than three mapped residues")
        selected = np.asarray(indices, dtype=int)
        output[name] = kabsch_rmsd(mobile[selected], reference[selected])
    return output


def hepn_geometry(
    coordinates: dict[int, NDArray[np.float64]], pairs: list[tuple[int, int]]
) -> dict[str, float]:
    output: dict[str, float] = {}
    for left, right in pairs:
        if left not in coordinates or right not in coordinates:
            raise ValueError(f"HEPN coordinate missing for pair {left}-{right}")
        distance = float(np.linalg.norm(coordinates[left] - coordinates[right]))
        if not math.isfinite(distance):
            raise ValueError("HEPN geometry produced a non-finite distance")
        output[f"{left}-{right}"] = distance
    return output


def contact_recovery(
    reference_contacts: set[Any], predicted_contacts: set[Any]
) -> dict[str, float]:
    if not reference_contacts:
        raise ValueError("reference RNA-contact set cannot be empty")
    overlap = len(reference_contacts.intersection(predicted_contacts))
    return {
        "recall": overlap / len(reference_contacts),
        "precision": overlap / len(predicted_contacts) if predicted_contacts else 0.0,
        "reference_count": float(len(reference_contacts)),
        "predicted_count": float(len(predicted_contacts)),
    }


def interface_confidence(plddt: list[float], positions: set[int]) -> float:
    if not positions or min(positions) < 0 or max(positions) >= len(plddt):
        raise ValueError("interface confidence positions are empty or invalid")
    values = [float(plddt[index]) for index in sorted(positions)]
    if any(not 0 <= value <= 100 for value in values):
        raise ValueError("pLDDT values must be in [0, 100]")
    return float(np.mean(values))


def run_usalign(
    prediction: Path,
    scaffold: Path,
    *,
    executable: Path,
    is_mock: bool,
) -> AlignmentComparison:
    if not executable.is_file():
        return AlignmentComparison(
            status="not_run",
            tm_score=None,
            rmsd=None,
            executable=str(executable),
            is_mock=is_mock,
            failure_reason="US-align/TM-align executable is missing",
        )
    completed = subprocess.run(
        [str(executable), str(prediction), str(scaffold)],
        text=True,
        capture_output=True,
        check=False,
    )
    if completed.returncode != 0:
        return AlignmentComparison(
            status="failed",
            tm_score=None,
            rmsd=None,
            executable=str(executable),
            is_mock=is_mock,
            failure_reason=completed.stderr[-1000:],
        )
    tm_scores: list[float] = []
    rmsd: float | None = None
    for line in completed.stdout.splitlines():
        if "TM-score=" in line:
            try:
                tm_scores.append(
                    float(line.split("TM-score=", 1)[1].strip().split()[0])
                )
            except ValueError:
                pass
        if "RMSD=" in line:
            try:
                rmsd = float(line.split("RMSD=", 1)[1].strip().split(",")[0])
            except ValueError:
                pass
    if not tm_scores or rmsd is None:
        raise ValueError("alignment output lacks a parseable TM-score or RMSD")
    return AlignmentComparison(
        status="success",
        tm_score=min(tm_scores),
        rmsd=rmsd,
        executable=str(executable),
        is_mock=is_mock,
    )


def consistency_summary(comparisons: list[AlignmentComparison]) -> dict[str, Any]:
    successful = [item for item in comparisons if item.status == "success"]
    values = [float(item.tm_score) for item in successful if item.tm_score is not None]
    return {
        "comparison_count": len(comparisons),
        "successful_count": len(successful),
        "minimum_tm_score": min(values) if values else None,
        "mean_tm_score": float(np.mean(values)) if values else None,
        "all_mock": all(item.is_mock for item in comparisons),
    }


def pareto_front(
    rows: list[dict[str, Any]],
    *,
    maximize: list[str],
    minimize: list[str],
) -> list[str]:
    if not rows or not maximize + minimize:
        raise ValueError("Pareto ranking needs rows and dimensions")
    output: list[str] = []
    for index, row in enumerate(rows):
        dominated = False
        for other_index, other in enumerate(rows):
            if index == other_index:
                continue
            no_worse = all(float(other[key]) >= float(row[key]) for key in maximize)
            no_worse &= all(float(other[key]) <= float(row[key]) for key in minimize)
            strictly_better = any(
                float(other[key]) > float(row[key]) for key in maximize
            ) or any(float(other[key]) < float(row[key]) for key in minimize)
            if no_worse and strictly_better:
                dominated = True
                break
        if not dominated:
            output.append(str(row["candidate_id"]))
    return sorted(output)
