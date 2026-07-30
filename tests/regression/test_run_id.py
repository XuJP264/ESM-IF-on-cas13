from datetime import datetime, timezone

from cas13_if.provenance import make_run_id


def test_run_id_is_deterministic_for_fixed_inputs() -> None:
    now = datetime(2026, 7, 31, tzinfo=timezone.utc)
    first = make_run_id("Example Run", {"a": 1}, "abc1234", now=now)
    second = make_run_id("Example Run", {"a": 1}, "abc1234", now=now)
    assert first == second == "20260731-example-run-015abd7f5c-abc1234"
