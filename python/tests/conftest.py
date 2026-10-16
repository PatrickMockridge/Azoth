"""Shared test configuration.

The one job here is to make "the extension is missing" behave differently
depending on why it is missing.

Ordinarily a test needing the Rust extension skips when it has not been built -
that is correct for a developer working only on Python. But in CI, where the
cross-language agreement is a promise the project makes, a silent skip means the
promise is checked by nothing and the suite stays green anyway. Setting
``AZOTH_REQUIRE_RUST=1`` turns those skips into failures.
"""

from __future__ import annotations

import os

import pytest

from azoth import available

#: Above this, a test that is **not parametrised** and not marked `slow_sweep` is a finding.
#: A parametrised test is selectable by `-k`, so being slow is a property of the case, not a
#: sweep; an unparametrised one that runs for minutes is a loop over every registered id, which
#: is exactly what the marker means. Measured 2026-10-16: the two marked sweeps run 836 s and
#: 462 s, and the slowest *parametrised* test is 138 s, so four minutes separates them cleanly.
SLOW_SECONDS = 240.0

_durations: dict[str, tuple[float, bool]] = {}


def pytest_runtest_logreport(report: pytest.TestReport) -> None:
    if report.when == "call":
        _durations[report.nodeid] = (report.duration, "slow_sweep" in report.keywords)


def pytest_sessionfinish(session: pytest.Session, exitstatus: int) -> None:
    """**A slow test that is not marked is a finding**, reported every night.

    The marker is what the iterate lane deselects, and a marker that stops being applied is
    invisible - which is exactly how a 462-second sweep went unnamed while the tool's own list
    claimed to have them all. So the suite times itself: under `AZOTH_SLOW_SWEEP_GATE=1` (the
    nightly and the release preflight) an unmarked sweep fails the session; without it the
    offenders are printed, because a developer running one file should not be stopped by it.
    """
    unmarked = sorted(
        (node, seconds)
        for node, (seconds, marked) in _durations.items()
        if seconds > SLOW_SECONDS and not marked and "[" not in node
    )
    if not unmarked:
        return
    for node, seconds in unmarked:
        print(f"\nslow and unmarked: {node} ran for {seconds:.0f}s")
    if os.environ.get("AZOTH_SLOW_SWEEP_GATE") == "1":
        print(
            "Mark it `@pytest.mark.slow_sweep` (or `pytestmark` for the whole module): an "
            "unparametrised test this slow is a sweep no `-k` can narrow, and the iterate lane "
            "needs to know."
        )
        session.exitstatus = 1


def pytest_collection_modifyitems(items: list[pytest.Item]) -> None:
    """Skip Rust-only tests, or fail them when Rust is required."""
    if "rust" in available():
        return
    require_rust = os.environ.get("AZOTH_REQUIRE_RUST") == "1"
    for item in items:
        if "requires_rust" not in item.keywords:
            continue
        if require_rust:
            item.add_marker(
                pytest.mark.xfail(
                    reason=(
                        "AZOTH_REQUIRE_RUST=1 but azoth._core is not built. The "
                        "cross-implementation guarantee must not quietly degrade into "
                        "a skip - build the extension with `maturin develop`, or unset "
                        "the variable if you are working on Python only."
                    ),
                    strict=True,
                    run=True,
                )
            )
        else:
            item.add_marker(
                pytest.mark.skip(reason="azoth._core is not built; run `maturin develop`")
            )
