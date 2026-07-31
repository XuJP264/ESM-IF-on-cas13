from pathlib import Path

import pytest

from cas13_if.refold.stage3 import (
    backend_input,
    export_stage3_jobs,
    make_stage3_job,
)


def _job(state: str = "ternary", backend: str = "alphafold3"):
    return make_stage3_job(
        candidate_id="candidate-a",
        parent_scaffold="EsCas13d",
        sequence="ACDE",
        crrna=None if state == "monomer" else "ACGU",
        target_rna="UGCA" if state == "ternary" else None,
        state=state,  # type: ignore[arg-type]
        seed=7,
        msa_policy="backend_default_reproducible_cache",
        template_policy="no_target_scaffold_as_forced_template",
        recycles=3,
        backend=backend,  # type: ignore[arg-type]
        shards=3,
    )


def test_stage3_job_is_deterministic_and_has_backend_input() -> None:
    first = _job()
    second = _job()
    assert first == second
    extension, content = backend_input(first)
    assert extension == "json"
    assert '"modelSeeds": [' in content
    assert '"rna"' in content


def test_stage3_export_writes_all_contract_directories(tmp_path: Path) -> None:
    jobs = [_job(), _job("binary", "boltz"), _job("monomer", "colabfold")]
    summary = export_stage3_jobs(jobs, tmp_path / "jobs", shards=3)
    assert summary["job_count"] == 3
    for name in (
        "monomer",
        "binary",
        "ternary",
        "manifests",
        "shards",
        "expected_outputs",
        "retry_manifests",
    ):
        assert (tmp_path / "jobs" / name).is_dir()


def test_stage3_job_rejects_forced_target_template() -> None:
    with pytest.raises(ValueError, match="target scaffold"):
        make_stage3_job(
            candidate_id="candidate-a",
            parent_scaffold="EsCas13d",
            sequence="ACDE",
            crrna=None,
            target_rna=None,
            state="monomer",
            seed=7,
            msa_policy="none",
            template_policy="force_target_scaffold_template",
            recycles=3,
            backend="colabfold",
            shards=1,
        )
