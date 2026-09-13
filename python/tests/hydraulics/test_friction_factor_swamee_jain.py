"""Spec-driven tests for ``hydraulics.friction_factor_swamee_jain``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import SwameeJainResult
from azoth.core.warnings import WarningCode
from azoth.hydraulics import friction_factor_swamee_jain, reynolds_number

CALC_ID = "hydraulics.friction_factor_swamee_jain"
Q = ureg.Quantity

SPEC = h.spec(CALC_ID)
#: Property test name -> the test function implementing it. Name-based rather
#: than a direct reference because the functions are defined further down.
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {"unit_round_trip": "test_unit_round_trip"}

ALL_CASES = h.all_tests(SPEC)
# Property tests assert invariants rather than comparing against a spec value, so
# they have no `inputs` and are exercised by their own named test functions rather
# than by the parametrized runner below.
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> SwameeJainResult:
    return friction_factor_swamee_jain(h.input_(case, "re"), h.input_(case, "relative_roughness"))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.f, h.expected(case, "f"), case.get("tolerance", 1e-12), f"{case['id']} (f)"
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "f":
            return result.f
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_unit_round_trip() -> None:
    """The same physical state in SI and US customary units gives the same f."""
    rr = 4.6e-4
    re_si = reynolds_number(Q(998.0, "kg/m**3"), Q(1.5, "m/s"), Q(0.1, "m"), Q(1.002e-3, "Pa*s")).re
    re_us = reynolds_number(
        Q(998.0, "kg/m**3"), Q(1.5 / 0.3048, "ft/s"), Q(0.1 / 0.3048, "ft"), Q(1.002e-3, "Pa*s")
    ).re
    h.assert_close(
        friction_factor_swamee_jain(re_us, rr).f,
        friction_factor_swamee_jain(re_si, rr).f,
        1e-12,
        "f from SI vs US customary state",
    )


def test_outside_the_published_range_warns_but_computes() -> None:
    """Below Re = 5000 the equation is outside the range the paper fitted it to.

    It is still well defined, so a warning rather than an error - but the
    accuracy claim no longer covers the answer.
    """
    low = friction_factor_swamee_jain(1000.0, 4.6e-4)
    assert low.has_warning(WarningCode.OUT_OF_VALID_RANGE)

    high = friction_factor_swamee_jain(1e9, 4.6e-4)
    assert high.has_warning(WarningCode.OUT_OF_VALID_RANGE)

    rough = friction_factor_swamee_jain(1e5, 0.5)
    assert rough.has_warning(WarningCode.OUT_OF_VALID_RANGE)


def test_inside_the_published_range_is_clean() -> None:
    """A mid-range case must produce no warnings at all.

    If this ever fires, the bounds or their inclusive flags have drifted.
    """
    result = friction_factor_swamee_jain(1e5, 4.6e-4)
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_zero_reynolds_number_is_an_error() -> None:
    """``Re**0.9`` is zero at the origin, so ``5.74/Re**0.9`` is singular."""
    with pytest.raises(OutOfRangeError) as excinfo:
        friction_factor_swamee_jain(0.0, 4.6e-4)
    assert excinfo.value.field() == "re"


def test_the_explicit_form_needs_no_iteration() -> None:
    """Deliberate contrast with Colebrook: the result type has no iteration count,
    because there is no iteration. The two result shapes differ for that reason."""
    first = friction_factor_swamee_jain(1e5, 4.6e-4)
    second = friction_factor_swamee_jain(1e5, 4.6e-4)
    assert first.f == second.f
    assert not hasattr(first, "iterations")


def test_property_tests_are_covered() -> None:
    """Every property test the spec declares must be one we actually run.

    Not a name check: for each property the spec asks for, the function that
    implements it is looked up and called here. A spec asking for an invariant
    that nothing checks is the failure mode this exists to catch, and a mapping
    that pointed at a renamed or deleted function would otherwise pass silently.
    """
    declared = {str(c["property"]) for c in PROPERTIES if c.get("property")}
    missing = declared - set(_PROPERTY_IMPLEMENTATIONS)
    assert not missing, (
        f"{CALC_ID}: spec declares property tests {sorted(missing)} that no test implements"
    )
    for property_name in sorted(declared):
        function_name = _PROPERTY_IMPLEMENTATIONS[property_name]
        function = globals().get(function_name)
        assert callable(function), (
            f"{CALC_ID}: property {property_name!r} claims to be implemented by "
            f"{function_name}, which does not exist"
        )
        function()
