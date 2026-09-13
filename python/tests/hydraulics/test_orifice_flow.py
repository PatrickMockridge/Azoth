"""Spec-driven tests for ``hydraulics.orifice_flow``."""

from __future__ import annotations

import math
from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import OrificeFlowResult
from azoth.hydraulics import orifice_flow

CALC_ID = "hydraulics.orifice_flow"
Q = ureg.Quantity

SPEC = h.spec(CALC_ID)
#: Property test name -> the test function implementing it. Name-based rather than
#: a direct reference because the functions are defined further down.
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "monotonic": "test_monotonic",
    "unit_round_trip": "test_unit_round_trip",
}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> OrificeFlowResult:
    # `kwargs_for` reads each input's declared unit from the spec, so `d` arrives as
    # `Q(value, "mm")` and `Cd` as a bare float, with no unit handling here.
    return orifice_flow(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.q.magnitude,
        h.expected(case, "q"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (q)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "q":
            return result.q.magnitude
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_monotonic() -> None:
    """Flow rises with the bore, the differential and the coefficient, and is zero
    when the differential is.

    The directions locate each input in the expression: a wrong power on ``d``, a
    missing square root or an inverted ratio all show up as a direction changing.
    The zero case is the one a form that divided by ``dP`` would get wrong.
    """

    def at(d: float, dp: float, rho: float, cd: float) -> float:
        return orifice_flow(Q(d, "mm"), Q(dp, "Pa"), Q(rho, "kg/m**3"), cd).q.magnitude

    base = at(50.0, 25000.0, 998.0, 0.62)
    assert at(100.0, 25000.0, 998.0, 0.62) > base, "q must rise with d"
    assert at(50.0, 50000.0, 998.0, 0.62) > base, "q must rise with dP"
    assert at(50.0, 25000.0, 500.0, 0.62) > base, "q must rise as rho falls"
    assert at(50.0, 25000.0, 998.0, 0.80) > base, "q must rise with Cd"
    assert at(50.0, 0.0, 998.0, 0.62) == 0.0, "no difference, no flow"


def test_unit_round_trip() -> None:
    """The same orifice and operating point described in SI and US customary units.

    This is the test the ``mm`` fix was waiting for: ``d`` is the one input in the
    registry whose spec unit is not its own SI base unit, so this is the first calc
    through which the millimetre conversion runs end to end.

    The US quantities come from pint's own conversion rather than restated factors -
    writing them out by hand is how the pump_power round trip first failed, at 7e-8
    against a 1e-12 tolerance.
    """
    si = orifice_flow(Q(50.0, "mm"), Q(25000.0, "Pa"), Q(998.0, "kg/m**3"), 0.62).q
    us = orifice_flow(
        Q(50.0, "mm").to("in"),
        Q(25000.0, "Pa").to("psi"),
        Q(998.0, "kg/m**3").to("lb/ft**3"),
        0.62,
    ).q

    h.assert_close(us.magnitude, si.magnitude, 1e-12, "q from SI vs US customary state")


def test_a_millimetre_bore_matches_the_same_bore_in_metres() -> None:
    """The millimetre conversion stated directly, in the calc that finally uses it.

    A 50 mm bore and a 0.05 m bore are the same orifice and must give the same flow.
    If the millimetre path handed the calculation millimetres instead of metres - the
    bug M1 fixed - this is the assertion that reports it, and it is the only test in
    the suite that exercises it through an actual calculation rather than through the
    units module.
    """
    in_mm = orifice_flow(Q(50.0, "mm"), Q(25000.0, "Pa"), Q(998.0, "kg/m**3"), 0.62)
    in_m = orifice_flow(Q(0.05, "m"), Q(25000.0, "Pa"), Q(998.0, "kg/m**3"), 0.62)
    assert in_mm.q.magnitude == in_m.q.magnitude


def test_an_ideal_coefficient_gives_the_ideal_flow() -> None:
    """The limiting case, which separates the coefficient from the rest.

    If the worked example and this both fail, the fault is in the area or the root;
    if only the worked example fails, it is in the coefficient.
    """
    ideal = math.pi * 0.05**2 / 4.0 * math.sqrt(2.0 * 25000.0 / 998.0)
    result = orifice_flow(Q(50.0, "mm"), Q(25000.0, "Pa"), Q(998.0, "kg/m**3"), 1.0)
    h.assert_close(result.q.magnitude, ideal, 1e-15, "ideal orifice")
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_a_zero_bore_is_an_error() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        orifice_flow(Q(0.0, "mm"), Q(25000.0, "Pa"), Q(998.0, "kg/m**3"), 0.62)
    assert excinfo.value.field() == "d"


def test_a_negative_differential_is_an_error_rather_than_reversed_flow() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        orifice_flow(Q(50.0, "mm"), Q(-25000.0, "Pa"), Q(998.0, "kg/m**3"), 0.62)
    assert excinfo.value.field() == "dP"


def test_a_coefficient_above_one_is_an_error() -> None:
    """More flow than an obstruction of that area can pass at that difference is not
    an uncertainty to warn about but a violation of the relation."""
    with pytest.raises(OutOfRangeError) as excinfo:
        orifice_flow(Q(50.0, "mm"), Q(25000.0, "Pa"), Q(998.0, "kg/m**3"), 1.01)
    assert excinfo.value.field() == "Cd"


def test_a_zero_coefficient_is_an_error() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        orifice_flow(Q(50.0, "mm"), Q(25000.0, "Pa"), Q(998.0, "kg/m**3"), 0.0)
    assert excinfo.value.field() == "Cd"


def test_property_tests_are_covered() -> None:
    """Every property test the spec declares must be one we actually run."""
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
