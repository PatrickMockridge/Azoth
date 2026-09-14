"""Spec-driven tests for ``eos.srk_z_factor``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth.core.result import SrkZFactorResult
from azoth.eos import srk_z_factor

CALC_ID = "eos.srk_z_factor"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> SrkZFactorResult:
    return srk_z_factor(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    for field in ("z_min", "z_max"):
        h.assert_close(
            getattr(result, field),
            h.expected(case, field),
            case.get("tolerance", 1e-12),
            f"{case['id']} ({field})",
        )
    h.assert_consistent(result, case["id"])


def test_every_active_case_ran() -> None:
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 2, f"expected several active cases, found {total}"


def test_property_tests_are_covered() -> None:
    declared = {str(c["property"]) for c in PROPERTIES if c.get("property")}
    missing = declared - set(_PROPERTY_IMPLEMENTATIONS)
    assert not missing, (
        f"{CALC_ID}: spec declares property tests {sorted(missing)} that no test implements"
    )


def test_the_skipped_case_says_why_it_is_skipped() -> None:
    h.assert_skips_are_explained(SPEC)
    declared = {str(c["property"]) for c in SKIPPED if c.get("property")}
    assert "unit_round_trip" in declared, (
        "this spec skips the unit round trip on purpose; if that changed, remove this test"
    )
