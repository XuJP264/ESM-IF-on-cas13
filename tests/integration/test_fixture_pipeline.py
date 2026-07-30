from pathlib import Path

from cas13_if.backends.mock import MockBackend
from cas13_if.data.atlas import process_atlas
from cas13_if.refold.mock import write_mock_predictions
from cas13_if.refold.providers import ManifestPredictionProvider
from cas13_if.reporting.reports import render_run_report
from cas13_if.schemas import SampleRequest


def test_fixture_atlas_to_candidate_to_refold_report(tmp_path: Path) -> None:
    funnel = process_atlas(
        Path("data/fixtures/atlas_operons.json"), tmp_path / "processed"
    )
    backend = MockBackend()
    backend.load()
    generated = backend.sample(
        SampleRequest(
            scaffold_id="fixture",
            structure_path="tests/fixtures/minimal_complex.pdb",
            parent_sequence="ACDEFGHIK",
            count=2,
            fixed_positions={0: "A", 8: "K"},
            seed=20260731,
        )
    )
    provider = ManifestPredictionProvider("colabfold", is_mock=True)
    jobs = provider.export_jobs(generated, tmp_path / "jobs", shards=2)
    write_mock_predictions(jobs, tmp_path / "predictions")
    ingested = provider.ingest_outputs(jobs, tmp_path / "predictions")
    qc = provider.qc_outputs(ingested)
    render_run_report(
        title="Fixture E2E",
        evidence_level=0,
        is_mock=True,
        metrics={"funnel": funnel, "refold": qc},
        failures=[],
        markdown_path=tmp_path / "report.md",
        html_path=tmp_path / "report.html",
    )
    assert qc["successful"] == 2
    assert "MOCK" in (tmp_path / "report.md").read_text(encoding="utf-8")
