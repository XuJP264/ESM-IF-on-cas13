"""Provider-neutral structure-prediction job exchange and QC."""

from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
from abc import ABC, abstractmethod
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from cas13_if.data.fasta import write_fasta
from cas13_if.schemas import Candidate

ProviderName = Literal["alphafold2", "colabfold", "alphafold3", "protenix", "boltz"]


@dataclass(frozen=True)
class PredictionJob:
    candidate_id: str
    sequence: str
    provider: str
    seed: int
    shard: int
    is_mock: bool


@dataclass(frozen=True)
class IngestedPrediction:
    candidate_id: str
    provider: str
    structure_path: str | None
    mean_plddt: float | None
    pae_path: str | None
    status: str
    is_mock: bool
    failure_reason: str | None = None


class StructurePredictionProvider(ABC):
    @abstractmethod
    def export_jobs(
        self, candidates: list[Candidate], output_dir: Path, *, shards: int
    ) -> list[PredictionJob]:
        """Export FASTA, JSONL, and deterministic shards."""

    @abstractmethod
    def shard_jobs(
        self, jobs: list[PredictionJob], output_dir: Path, *, shards: int
    ) -> None:
        """Write deterministic shard manifests."""

    @abstractmethod
    def validate_inputs(self, jobs: list[PredictionJob]) -> None:
        """Reject malformed or duplicate prediction jobs."""

    @abstractmethod
    def expected_outputs(self, job: PredictionJob) -> dict[str, Any]:
        """Describe exact expected outputs for one job."""

    @abstractmethod
    def ingest_outputs(
        self, jobs: list[PredictionJob], result_root: Path
    ) -> list[IngestedPrediction]:
        """Ingest results without silently ignoring missing jobs."""

    @abstractmethod
    def qc_outputs(self, predictions: list[IngestedPrediction]) -> dict[str, Any]:
        """Audit success, failures, mock state, pLDDT, and PAE."""

    @abstractmethod
    def compare_to_scaffold(self, prediction: Path, scaffold: Path) -> dict[str, Any]:
        """Use genuine US-align when installed; never approximate TM-score."""


class ManifestPredictionProvider(StructurePredictionProvider):
    """Job exchange implementation shared by declared production providers."""

    def __init__(self, provider: ProviderName, *, is_mock: bool = False) -> None:
        self.provider = provider
        self.is_mock = is_mock

    def export_jobs(
        self, candidates: list[Candidate], output_dir: Path, *, shards: int
    ) -> list[PredictionJob]:
        if shards < 1:
            raise ValueError("shards must be positive")
        output_dir.mkdir(parents=True, exist_ok=False)
        sorted_candidates = sorted(candidates, key=lambda item: item.candidate_id)
        jobs = [
            PredictionJob(
                candidate_id=candidate.candidate_id,
                sequence=candidate.sequence,
                provider=self.provider,
                seed=candidate.seed,
                shard=_stable_shard(candidate.candidate_id, shards),
                is_mock=self.is_mock,
            )
            for candidate in sorted_candidates
        ]
        self.validate_inputs(jobs)
        write_fasta(
            [(job.candidate_id, job.sequence) for job in jobs],
            output_dir / "candidates.fasta",
        )
        _write_jsonl(output_dir / "jobs.jsonl", [asdict(job) for job in jobs])
        self.shard_jobs(jobs, output_dir, shards=shards)
        expected = [
            {"job": asdict(job), "outputs": self.expected_outputs(job)} for job in jobs
        ]
        (output_dir / "expected_outputs.json").write_text(
            json.dumps(expected, indent=2, sort_keys=True) + "\n", encoding="utf-8"
        )
        (output_dir / "failed_job_retry.jsonl").write_text("", encoding="utf-8")
        return jobs

    def shard_jobs(
        self, jobs: list[PredictionJob], output_dir: Path, *, shards: int
    ) -> None:
        shard_dir = output_dir / "shards"
        shard_dir.mkdir()
        for shard in range(shards):
            members = [asdict(job) for job in jobs if job.shard == shard]
            _write_jsonl(shard_dir / f"shard-{shard:04d}.jsonl", members)

    def validate_inputs(self, jobs: list[PredictionJob]) -> None:
        identifiers = [job.candidate_id for job in jobs]
        if len(identifiers) != len(set(identifiers)):
            raise ValueError("duplicate candidate IDs in prediction jobs")
        for job in jobs:
            if not job.sequence or not job.candidate_id:
                raise ValueError("prediction job has empty ID or sequence")
            if job.provider != self.provider or job.is_mock != self.is_mock:
                raise ValueError("job provider/mock metadata mismatch")

    def expected_outputs(self, job: PredictionJob) -> dict[str, Any]:
        base = f"{job.candidate_id}"
        return {
            "candidate_id": job.candidate_id,
            "provider": self.provider,
            "result_json": f"{base}/result.json",
            "structure_path": f"{base}/prediction.cif",
            "pae_path": f"{base}/pae.json",
            "required_result_fields": [
                "candidate_id",
                "provider",
                "mean_plddt",
                "structure_path",
                "is_mock",
            ],
        }

    def ingest_outputs(
        self, jobs: list[PredictionJob], result_root: Path
    ) -> list[IngestedPrediction]:
        predictions: list[IngestedPrediction] = []
        retry: list[dict[str, Any]] = []
        for job in jobs:
            reason: str | None
            expected = self.expected_outputs(job)
            result_path = result_root / str(expected["result_json"])
            if not result_path.is_file():
                reason = f"missing {result_path}"
                predictions.append(
                    IngestedPrediction(
                        candidate_id=job.candidate_id,
                        provider=self.provider,
                        structure_path=None,
                        mean_plddt=None,
                        pae_path=None,
                        status="missing",
                        is_mock=self.is_mock,
                        failure_reason=reason,
                    )
                )
                retry.append({"job": asdict(job), "reason": reason})
                continue
            result = json.loads(result_path.read_text(encoding="utf-8"))
            if bool(result.get("is_mock")) != self.is_mock:
                raise ValueError(f"mock-state mismatch in {result_path}")
            structure_value = result.get("structure_path")
            structure_path = (
                result_path.parent / str(structure_value)
                if structure_value is not None
                else None
            )
            pae_value = result.get("pae_path")
            pae_path = (
                result_path.parent / str(pae_value) if pae_value is not None else None
            )
            missing = (
                structure_path is None
                or not structure_path.is_file()
                or pae_path is None
                or not pae_path.is_file()
            )
            status = "failed_qc" if missing else "success"
            reason = "declared structure or PAE output missing" if missing else None
            predictions.append(
                IngestedPrediction(
                    candidate_id=job.candidate_id,
                    provider=self.provider,
                    structure_path=str(structure_path) if structure_path else None,
                    mean_plddt=(
                        float(result["mean_plddt"])
                        if result.get("mean_plddt") is not None
                        else None
                    ),
                    pae_path=str(pae_path) if pae_path else None,
                    status=status,
                    is_mock=self.is_mock,
                    failure_reason=reason,
                )
            )
            if reason:
                retry.append({"job": asdict(job), "reason": reason})
        _write_jsonl(result_root / "failed_job_retry.jsonl", retry)
        return predictions

    def qc_outputs(self, predictions: list[IngestedPrediction]) -> dict[str, Any]:
        successful = [
            prediction for prediction in predictions if prediction.status == "success"
        ]
        plddts = [
            prediction.mean_plddt
            for prediction in successful
            if prediction.mean_plddt is not None
        ]
        return {
            "provider": self.provider,
            "is_mock": self.is_mock,
            "jobs": len(predictions),
            "successful": len(successful),
            "missing_or_failed": len(predictions) - len(successful),
            "mean_plddt": sum(plddts) / len(plddts) if plddts else None,
        }

    def compare_to_scaffold(self, prediction: Path, scaffold: Path) -> dict[str, Any]:
        executable = shutil.which("USalign") or shutil.which("TMalign")
        if executable is None:
            return {
                "status": "not_run",
                "reason": "US-align/TM-align executable not found",
                "tm_score": None,
                "rmsd": None,
                "is_mock": self.is_mock,
            }
        result = subprocess.run(
            [executable, str(prediction), str(scaffold)],
            text=True,
            capture_output=True,
            check=False,
        )
        if result.returncode != 0:
            return {
                "status": "failed",
                "reason": result.stderr.strip(),
                "tm_score": None,
                "rmsd": None,
                "is_mock": self.is_mock,
            }
        return {
            "status": "success",
            "raw_output": result.stdout,
            "tm_score": _parse_first_value(result.stdout, "TM-score="),
            "rmsd": _parse_rmsd(result.stdout),
            "executable": executable,
            "is_mock": self.is_mock,
        }


def _stable_shard(identifier: str, shards: int) -> int:
    digest = hashlib.sha256(identifier.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big") % shards


def _write_jsonl(path: Path, rows: list[dict[str, Any]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
        encoding="utf-8",
    )


def _parse_first_value(output: str, marker: str) -> float | None:
    for line in output.splitlines():
        if marker in line:
            tail = line.split(marker, maxsplit=1)[1].strip().split()[0]
            try:
                return float(tail)
            except ValueError:
                continue
    return None


def _parse_rmsd(output: str) -> float | None:
    for line in output.splitlines():
        if "RMSD=" in line:
            tail = line.split("RMSD=", maxsplit=1)[1].strip().split(",")[0]
            try:
                return float(tail)
            except ValueError:
                continue
    return None
