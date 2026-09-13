"""Spec-driven tests for ``hydraulics.pump_power``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PumpPowerResult
from azoth.hydraulics import pump_power
from azoth.hydraulics.reference.pump_power import STANDARD_GRAVITY_M_S2

CALC_ID = "hydraulics.pump_power"
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


def call(case: dict[str, Any]) -> PumpPowerResult:
    # `kwargs_for` reads each input's declared unit from the spec, so `eta` arrives
    # as a bare float (dimensionless) and the rest as quantities, with no unit
    # handling here.
    return pump_power(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.power.magnitude,
        h.expected(case, "power"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (power)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "power":
            return result.power.magnitude
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_monotonic() -> None:
    """Power rises with density, flow and head, and falls with efficiency.

    The four partial derivatives. A point-value test can pass with a wrong formula
    that happens to be right at one operating point; an inverted ratio or a sign
    error cannot survive all four directions at once.
    """

    def at(rho: float, q: float, head: float, eta: float) -> float:
        return pump_power(Q(rho, "kg/m**3"), Q(q, "m**3/s"), Q(head, "m"), eta).power.magnitude

    base = at(998.0, 0.01, 30.0, 0.75)
    assert at(1200.0, 0.01, 30.0, 0.75) > base, "power must rise with rho"
    assert at(998.0, 0.02, 30.0, 0.75) > base, "power must rise with q"
    assert at(998.0, 0.01, 60.0, 0.75) > base, "power must rise with H"
    assert at(998.0, 0.01, 30.0, 0.50) > base, "power must rise as eta falls"


def test_unit_round_trip() -> None:
    """The same operating point described in SI and in US customary units.

    All three dimensioned inputs are converted, so an accidental factor between them
    shows up as a different power rather than cancelling out.

    The US quantities are produced by pint's own conversion rather than by restating
    the factors here. The Rust half of this test learned that the hard way: it wrote
    the factors out, they were accurate to about 7e-8, and the test failed at the
    1e-12 tolerance the rest of the suite uses. Restating what the units library
    already knows is the mistake that produced the `mm` bug, and it hid here as a
    test that looked like it was checking the calculation.
    """
    si = pump_power(Q(998.0, "kg/m**3"), Q(0.01, "m**3/s"), Q(30.0, "m"), 0.75).power
    us = pump_power(
        Q(998.0, "kg/m**3").to("lb/ft**3"),
        Q(0.01, "m**3/s").to("ft**3/s"),
        Q(30.0, "m").to("ft"),
        0.75,
    ).power

    h.assert_close(us.magnitude, si.magnitude, 1e-12, "power from SI vs US customary state")


def test_standard_gravity_is_the_defined_value() -> None:
    """Pinned in both languages, so a change to one of them fails a test rather
    than putting a 0.5% error between the two implementations that neither
    reports."""
    assert STANDARD_GRAVITY_M_S2 == 9.80665


def test_an_efficiency_of_one_reduces_to_the_hydraulic_power() -> None:
    """The limiting case: a lossless pump delivers exactly ``rho * g * q * H``.

    This separates the two halves of the equation, so a failure says which half
    broke rather than only that something did.
    """
    hydraulic = 998.0 * STANDARD_GRAVITY_M_S2 * 0.01 * 30.0
    result = pump_power(Q(998.0, "kg/m**3"), Q(0.01, "m**3/s"), Q(30.0, "m"), 1.0)
    h.assert_close(result.power.magnitude, hydraulic, 1e-12, "ideal pump")
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_a_negative_flow_is_an_error_rather_than_reversed_flow() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        pump_power(Q(998.0, "kg/m**3"), Q(-0.01, "m**3/s"), Q(30.0, "m"), 0.75)
    assert excinfo.value.field() == "q"


def test_an_efficiency_above_one_is_an_error() -> None:
    """Not an engineering uncertainty to warn about: more hydraulic power out than
    shaft power in is a violation of the first law."""
    with pytest.raises(OutOfRangeError) as excinfo:
        pump_power(Q(998.0, "kg/m**3"), Q(0.01, "m**3/s"), Q(30.0, "m"), 1.01)
    assert excinfo.value.field() == "eta"


def test_a_zero_efficiency_is_an_error_rather_than_infinite_power() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        pump_power(Q(998.0, "kg/m**3"), Q(0.01, "m**3/s"), Q(30.0, "m"), 0.0)
    assert excinfo.value.field() == "eta"


def test_zero_flow_needs_zero_power_and_is_not_a_warning() -> None:
    """Allowed, and correct rather than merely tolerable: no flow moves no fluid."""
    result = pump_power(Q(998.0, "kg/m**3"), Q(0.0, "m**3/s"), Q(30.0, "m"), 0.75)
    assert result.power.magnitude == 0.0
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


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
