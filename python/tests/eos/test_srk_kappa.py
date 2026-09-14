"""Spec-driven tests for ``eos.srk_kappa``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth.core.result import SrkKappaResult
from azoth.eos import srk_kappa

CALC_ID = "eos.srk_kappa"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> SrkKappaResult:
    return srk_kappa(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.kappa,
        h.expected(case, "kappa"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (kappa)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "kappa":
            return result.kappa
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 2, f"expected several active cases, found {total}"


def test_the_worked_example_is_exact_in_binary() -> None:
    assert srk_kappa(0.152).kappa == 0.715181696


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
