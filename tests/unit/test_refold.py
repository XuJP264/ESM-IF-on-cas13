from pathlib import Path

from cas13_if.refold.mock import write_mock_predictions
from cas13_if.refold.providers import ManifestPredictionProvider
from cas13_if.schemas import Candidate, EvidenceLevel


def candidates() -> list[Candidate]:
    return [
        Candidate(
            candidate_id=f"candidate-{index}",
            scaffold_id="fixture",
            backend="mock",
            sequence="ACDE",
            seed=index,
            temperature=1,
            is_mock=True,
            evidence_level=EvidenceLevel.IO_VALIDATED,
        )
        for index in range(3)
    ]


def test_refold_export_ingest_and_mock_labels(tmp_path: Path) -> None:
    provider = ManifestPredictionProvider("colabfold", is_mock=True)
    export = tmp_path / "export"
    jobs = provider.export_jobs(candidates(), export, shards=2)
    assert len(jobs) == 3
    assert (export / "shards/shard-0000.jsonl").is_file()
    result_root = tmp_path / "predictions"
    write_mock_predictions(jobs, result_root)
    ingested = provider.ingest_outputs(jobs, result_root)
    qc = provider.qc_outputs(ingested)
    assert qc["successful"] == 3
    assert qc["is_mock"] is True
    comparison = provider.compare_to_scaffold(
        result_root / "candidate-0/prediction.cif",
        Path("tests/fixtures/minimal_complex.pdb"),
    )
    assert comparison["status"] in {"not_run", "success"}


def test_missing_refold_output_enters_retry(tmp_path: Path) -> None:
    provider = ManifestPredictionProvider("boltz", is_mock=True)
    jobs = provider.export_jobs(candidates()[:1], tmp_path / "export", shards=1)
    result_root = tmp_path / "missing"
    result_root.mkdir()
    predictions = provider.ingest_outputs(jobs, result_root)
    assert predictions[0].status == "missing"
    assert (result_root / "failed_job_retry.jsonl").read_text(encoding="utf-8")
