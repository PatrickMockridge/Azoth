"""Spec-driven tests for ``hydraulics.friction_factor_haaland``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import HaalandResult
from azoth.core.warnings import WarningCode
from azoth.hydraulics import (
    friction_factor_colebrook,
    friction_factor_haaland,
    reynolds_number,
)

CALC_ID = "hydraulics.friction_factor_haaland"
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


def call(case: dict[str, Any]) -> HaalandResult:
    return friction_factor_haaland(h.input_(case, "re"), h.input_(case, "relative_roughness"))


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
        friction_factor_haaland(re_us, rr).f,
        friction_factor_haaland(re_si, rr).f,
        1e-12,
        "f from SI vs US customary state",
    )


def test_the_two_explicit_forms_are_the_same_cloth() -> None:
    """Haaland and Swamee-Jain both approximate Colebrook, and are close to it.

    Neither is checked against the other - each is checked against Colebrook by its
    own spec - so this asserts only that they sit within a few percent of one
    another at a mid-range point. It is a sanity check on the registry rather than
    on either equation: if these two ever disagree wildly, one of the three
    friction factor calcs has been changed in a way its own tests did not catch.
    """
    re, rr = 1.0e5, 4.6e-4
    haaland = friction_factor_haaland(re, rr).f
    colebrook = friction_factor_colebrook(re, rr).f
    assert abs(haaland / colebrook - 1.0) < 0.05, (
        f"Haaland is {100 * (haaland / colebrook - 1):.2f}% from Colebrook, outside the "
        f"published accuracy of about 1-2% by a wide margin"
    )


def test_outside_the_published_range_warns_but_computes() -> None:
    """Below Re ~ 4000 the equation is outside the turbulent range it describes.

    It is still well defined, so a warning rather than an error - but the accuracy
    claim no longer covers the answer.
    """
    low = friction_factor_haaland(1000.0, 4.6e-4)
    assert low.has_warning(WarningCode.OUT_OF_VALID_RANGE)

    rough = friction_factor_haaland(1e5, 0.5)
    assert rough.has_warning(WarningCode.OUT_OF_VALID_RANGE)


def test_inside_the_published_range_is_clean() -> None:
    """A mid-range case must produce no warnings at all.

    If this ever fires, the bounds or their inclusive flags have drifted.
    """
    result = friction_factor_haaland(1e5, 4.6e-4)
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_zero_reynolds_number_is_an_error() -> None:
    """``6.9/Re`` is singular at the origin."""
    with pytest.raises(OutOfRangeError) as excinfo:
        friction_factor_haaland(0.0, 4.6e-4)
    assert excinfo.value.field() == "re"


def test_negative_roughness_is_an_error() -> None:
    """A negative roughness has no physical meaning and is rejected."""
    with pytest.raises(OutOfRangeError) as excinfo:
        friction_factor_haaland(1e5, -1.0e-4)
    assert excinfo.value.field() == "relative_roughness"


def test_the_explicit_form_needs_no_iteration() -> None:
    """Deliberate contrast with Colebrook: the result type has no iteration count,
    because there is no iteration."""
    first = friction_factor_haaland(1e5, 4.6e-4)
    second = friction_factor_haaland(1e5, 4.6e-4)
    assert first.f == second.f
    assert not hasattr(first, "iterations")


def test_property_tests_are_covered() -> None:
    """Every property test the spec declares must be one we actually run.

    Not a name check: for each property the spec asks for, the function that
    implements it is looked up and called here. A spec asking for an invariant
    that nothing checks is the failure mode this exists to catch.
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
