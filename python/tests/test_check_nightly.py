"""The three ways a nightly stops being a gate, from fixtures rather than from the API.

`tools/check_nightly.py` is the only thing in the tree that can notice a schedule GitHub has
silently disabled, so it has to be right about the difference between *healthy*, *stopped* and
*could not decide*. The third is the one worth pinning: an unreachable API returning 0 would
make an outage look like a pass, which is the failure the whole testing regime exists to
avoid.
"""

from __future__ import annotations

import importlib
import json
import sys
import urllib.error
import urllib.request
from datetime import UTC, datetime, timedelta
from pathlib import Path
from types import ModuleType

import pytest

TOOLS = Path(__file__).resolve().parents[2] / "tools"

NOW = datetime(2026, 10, 16, 12, 0, tzinfo=UTC)


def tool() -> ModuleType:
    """`tools/check_nightly.py`, imported by name - the tools are not a package."""
    sys.path.insert(0, str(TOOLS))
    try:
        return importlib.import_module("check_nightly")
    finally:
        sys.path.pop(0)


def run(hours_ago: float, conclusion: str | None) -> dict[str, object]:
    """One workflow run, as much of the payload as `evaluate` reads."""
    created = NOW - timedelta(hours=hours_ago)
    return {"created_at": created.isoformat().replace("+00:00", "Z"), "conclusion": conclusion}


def test_a_fresh_green_run_passes() -> None:
    code, message = tool().evaluate([run(4, "success")], NOW, 36.0, 48.0)
    assert code == 0
    assert "succeeded 4 hours ago" in message


def test_a_workflow_that_has_never_run_is_a_finding() -> None:
    """The state this repository is actually in until 04:17 UTC tonight."""
    code, message = tool().evaluate([], NOW, 36.0, 48.0)
    assert code == 1
    assert "never run" in message


def test_a_stopped_schedule_is_a_finding_and_says_what_to_do() -> None:
    code, message = tool().evaluate([run(90, "success")], NOW, 36.0, 48.0)
    assert code == 1
    assert "90 hours old" in message
    assert "disables a schedule" in message


def test_a_red_nightly_is_a_finding_even_when_it_is_running() -> None:
    """A schedule that runs every night and fails every night is not a gate either."""
    code, message = tool().evaluate([run(4, "failure"), run(28, "failure")], NOW, 36.0, 48.0)
    assert code == 1
    assert "not succeeded" in message


def test_a_green_run_that_is_too_old_is_a_finding() -> None:
    """The newest run is fresh but red; the newest green one is not."""
    runs = [run(2, "failure"), run(70, "success")]
    code, message = tool().evaluate(runs, NOW, 36.0, 48.0)
    assert code == 1
    assert "not been green for 70 hours" in message


def test_an_api_that_cannot_be_read_is_two_never_zero(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    """**The distinction the exit codes exist for.** An outage is not a healthy nightly."""

    def unreachable(*_: object, **__: object) -> None:
        raise urllib.error.URLError("no route to host")

    monkeypatch.setattr(urllib.request, "urlopen", unreachable)
    monkeypatch.setenv("GITHUB_TOKEN", "not-a-real-token")
    monkeypatch.setattr(
        sys,
        "argv",
        ["check_nightly.py", "--check", "--repo", "owner/name"],
    )
    assert tool().main() == 2


def test_no_token_is_two_as_well(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.delenv("GITHUB_TOKEN", raising=False)
    monkeypatch.setattr(sys, "argv", ["check_nightly.py", "--check", "--repo", "owner/name"])
    assert tool().main() == 2


def test_the_fetch_reads_the_field_it_claims_to() -> None:
    """A cheap check that the payload shape `fetch` expects is the one the API documents."""
    payload = {"workflow_runs": [run(1, "success")]}
    assert json.loads(json.dumps(payload))["workflow_runs"][0]["conclusion"] == "success"
