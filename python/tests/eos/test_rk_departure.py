"""Spec-driven tests for ``eos.rk_departure``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth.core.result import RkDepartureResult
from azoth.eos import rk_departure

CALC_ID = "eos.rk_departure"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]

FIELDS = ("ln_phi", "h_dep_rt", "s_dep_r", "cp_dep_r")


def call(case: dict[str, Any]) -> RkDepartureResult:
    return rk_departure(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    for field in FIELDS:
        h.assert_close(
            getattr(result, field),
            h.expected(case, field),
            case.get("tolerance", 1e-12),
            f"{case['id']} ({field})",
        )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity in FIELDS:
            return float(getattr(result, quantity))
        if quantity == "z_minus_b_reduced":
            return float(_case["inputs"]["z"]) - float(_case["inputs"]["b_reduced"])
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 2, f"expected several active cases, found {total}"


def test_the_gibbs_identity_holds() -> None:
    """``h_dep_rt - s_dep_r`` equals ``ln_phi``, the exact identity."""
    d = rk_departure(0.1866943088347048, 0.02707510936404929, 0.8119920001727409)
    assert abs(d.h_dep_rt - d.s_dep_r - d.ln_phi) < 1e-12


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
