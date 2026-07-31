from datetime import datetime, timezone
from pathlib import Path

from cas13_if.provenance import RunRecorder, make_run_id


def test_run_id_is_deterministic_for_fixed_inputs() -> None:
    now = datetime(2026, 7, 31, tzinfo=timezone.utc)
    first = make_run_id("Example Run", {"a": 1}, "abc1234", now=now)
    second = make_run_id("Example Run", {"a": 1}, "abc1234", now=now)
    assert first == second == "20260731-example-run-015abd7f5c-abc1234"


def test_run_recorder_uses_immutable_retry_suffix(
    tmp_path: Path, monkeypatch: object
) -> None:
    import cas13_if.provenance as provenance

    monkeypatch.setattr(  # type: ignore[attr-defined]
        provenance,
        "make_run_id",
        lambda *args, **kwargs: "20260731-benchmark-beced1fe69-abc1234",
    )
    arguments = {
        "root": tmp_path,
        "experiment": "benchmark",
        "resolved_config": {"seed": 7},
        "command": ["cas13-if", "benchmark"],
        "repo_root": Path(__file__).resolve().parents[2],
        "is_mock": False,
    }
    first = RunRecorder(**arguments)  # type: ignore[arg-type]
    second = RunRecorder(**arguments)  # type: ignore[arg-type]
    assert first.run_dir.name == "20260731-benchmark-beced1fe69-abc1234"
    assert second.run_dir.name == "20260731-benchmark-beced1fe69-abc1234-r001"
    assert first.run_dir.is_dir()
    assert second.run_dir.is_dir()
