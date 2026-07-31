"""Deterministic Stage-0003 monomer and RNA-complex job manifests."""

from __future__ import annotations

import hashlib
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any, Literal

from cas13_if.provenance import atomic_write_text
from cas13_if.schemas import STANDARD_AA

RefoldState = Literal["monomer", "binary", "ternary"]
Stage3Backend = Literal["colabfold", "alphafold2", "alphafold3", "boltz"]


@dataclass(frozen=True)
class Stage3RefoldJob:
    job_id: str
    candidate_id: str
    parent_scaffold: str
    sequence: str
    crrna: str | None
    target_rna: str | None
    state: RefoldState
    seed: int
    msa_policy: str
    template_policy: str
    recycles: int
    backend: Stage3Backend
    expected_output: str
    shard: int
    is_mock: bool = False


def stable_shard(identifier: str, shards: int) -> int:
    if shards < 1:
        raise ValueError("shards must be positive")
    digest = hashlib.sha256(identifier.encode("utf-8")).digest()
    return int.from_bytes(digest[:8], byteorder="big") % shards


def make_stage3_job(
    *,
    candidate_id: str,
    parent_scaffold: str,
    sequence: str,
    crrna: str | None,
    target_rna: str | None,
    state: RefoldState,
    seed: int,
    msa_policy: str,
    template_policy: str,
    recycles: int,
    backend: Stage3Backend,
    shards: int,
) -> Stage3RefoldJob:
    payload = {
        "candidate_id": candidate_id,
        "parent_scaffold": parent_scaffold,
        "state": state,
        "seed": seed,
        "backend": backend,
        "msa_policy": msa_policy,
        "template_policy": template_policy,
        "recycles": recycles,
    }
    digest = hashlib.sha256(
        json.dumps(payload, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()[:16]
    job_id = f"{state}-{backend}-{digest}"
    job = Stage3RefoldJob(
        job_id=job_id,
        candidate_id=candidate_id,
        parent_scaffold=parent_scaffold,
        sequence=sequence,
        crrna=crrna,
        target_rna=target_rna,
        state=state,
        seed=seed,
        msa_policy=msa_policy,
        template_policy=template_policy,
        recycles=recycles,
        backend=backend,
        expected_output=f"predictions/{job_id}/result.json",
        shard=stable_shard(job_id, shards),
    )
    validate_stage3_job(job)
    return job


def validate_stage3_job(job: Stage3RefoldJob) -> None:
    if not job.candidate_id or not job.parent_scaffold:
        raise ValueError("candidate and parent scaffold identifiers are required")
    invalid = sorted(set(job.sequence).difference(STANDARD_AA))
    if not job.sequence or invalid:
        raise ValueError(f"invalid protein sequence; invalid={invalid}")
    if job.recycles < 1 or job.seed < 0 or job.shard < 0:
        raise ValueError("seed, recycles, or shard is invalid")
    if "target_scaffold" in job.template_policy and "no_target" not in (
        job.template_policy
    ):
        raise ValueError("target scaffold cannot be forced as a template")
    if job.state == "monomer" and (job.crrna or job.target_rna):
        raise ValueError("monomer job cannot include RNA")
    if job.state == "binary" and (not job.crrna or job.target_rna):
        raise ValueError("binary job requires crRNA only")
    if job.state == "ternary" and (not job.crrna or not job.target_rna):
        raise ValueError("ternary job requires crRNA and target RNA")
    for name, sequence in (("crRNA", job.crrna), ("target RNA", job.target_rna)):
        if sequence is not None and (not sequence or set(sequence).difference("ACGU")):
            raise ValueError(f"{name} contains noncanonical RNA symbols")
    compatible: dict[str, set[str]] = {
        "monomer": {"colabfold", "alphafold2", "alphafold3", "boltz"},
        "binary": {"alphafold3", "boltz"},
        "ternary": {"alphafold3", "boltz"},
    }
    if job.backend not in compatible[job.state]:
        raise ValueError(f"{job.backend} is not enabled for {job.state}")


def backend_input(job: Stage3RefoldJob) -> tuple[str, str]:
    """Return the backend-neutral filename extension and exact input text."""
    validate_stage3_job(job)
    if job.backend in {"colabfold", "alphafold2"}:
        return "fasta", f">{job.candidate_id}\n{job.sequence}\n"
    chains: list[dict[str, Any]] = [
        {"protein": {"id": ["A"], "sequence": job.sequence}}
    ]
    if job.crrna:
        chains.append({"rna": {"id": ["B"], "sequence": job.crrna}})
    if job.target_rna:
        chains.append({"rna": {"id": ["C"], "sequence": job.target_rna}})
    if job.backend == "alphafold3":
        return (
            "json",
            json.dumps(
                {
                    "name": job.job_id,
                    "modelSeeds": [job.seed],
                    "sequences": chains,
                    "dialect": "alphafold3",
                    "version": 1,
                },
                indent=2,
                sort_keys=True,
            )
            + "\n",
        )
    yaml_lines = ["version: 1", "sequences:"]
    for chain in chains:
        kind, value = next(iter(chain.items()))
        yaml_lines.extend(
            [
                f"  - {kind}:",
                f"      id: {value['id'][0]}",
                f"      sequence: {value['sequence']}",
            ]
        )
    return "yaml", "\n".join(yaml_lines) + "\n"


def export_stage3_jobs(
    jobs: list[Stage3RefoldJob], output_root: Path, *, shards: int
) -> dict[str, Any]:
    identifiers = [job.job_id for job in jobs]
    if not jobs or len(identifiers) != len(set(identifiers)):
        raise ValueError("jobs are empty or contain duplicate IDs")
    if output_root.exists():
        raise FileExistsError(f"refusing to overwrite job root: {output_root}")
    for job in jobs:
        validate_stage3_job(job)
        if job.shard >= shards:
            raise ValueError("job shard exceeds declared shard count")
    for name in (
        "monomer",
        "binary",
        "ternary",
        "manifests",
        "shards",
        "expected_outputs",
        "retry_manifests",
    ):
        (output_root / name).mkdir(parents=True, exist_ok=False)
    sorted_jobs = sorted(jobs, key=lambda item: item.job_id)
    for job in sorted_jobs:
        extension, content = backend_input(job)
        path = output_root / job.state / job.backend / f"{job.job_id}.{extension}"
        atomic_write_text(path, content)
    for state in ("monomer", "binary", "ternary"):
        selected = [asdict(job) for job in sorted_jobs if job.state == state]
        _jsonl(output_root / f"manifests/{state}.jsonl", selected)
    _jsonl(output_root / "manifests/all_jobs.jsonl", map(asdict, sorted_jobs))
    for state in ("monomer", "binary", "ternary"):
        for backend in ("colabfold", "alphafold2", "alphafold3", "boltz"):
            selected_jobs = [
                job
                for job in sorted_jobs
                if job.state == state and job.backend == backend
            ]
            if not selected_jobs:
                continue
            for shard in range(shards):
                selected = [asdict(job) for job in selected_jobs if job.shard == shard]
                path = (
                    output_root
                    / "shards"
                    / state
                    / backend
                    / f"shard-{shard:04d}.jsonl"
                )
                _jsonl(path, selected)
    expected = [
        {
            "job_id": job.job_id,
            "candidate_id": job.candidate_id,
            "provider": job.backend,
            "result_json": job.expected_output,
            "required_files": [
                "result.json",
                "prediction.cif",
                "pae.json",
                "execution.json",
            ],
            "is_mock": False,
        }
        for job in sorted_jobs
    ]
    _jsonl(output_root / "expected_outputs/expected_outputs.jsonl", expected)
    _jsonl(output_root / "retry_manifests/failed_jobs.jsonl", [])
    counts: dict[str, int] = {}
    for job in sorted_jobs:
        key = f"{job.state}:{job.backend}"
        counts[key] = counts.get(key, 0) + 1
    return {
        "job_count": len(sorted_jobs),
        "candidate_count": len({job.candidate_id for job in sorted_jobs}),
        "counts": counts,
        "shards": shards,
        "is_mock": False,
        "evidence_level": 0,
        "status": "prepared_not_run",
    }


def _jsonl(path: Path, rows: Any) -> None:
    atomic_write_text(
        path,
        "".join(json.dumps(row, sort_keys=True) + "\n" for row in rows),
    )
