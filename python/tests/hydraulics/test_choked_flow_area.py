"""Spec-driven tests for ``hydraulics.choked_flow_area``."""

from __future__ import annotations

import math
from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import ChokedFlowAreaResult
from azoth.core.warnings import WarningCode
from azoth.hydraulics import choked_flow_area

CALC_ID = "hydraulics.choked_flow_area"
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


def call(case: dict[str, Any]) -> ChokedFlowAreaResult:
    # `kwargs_for` reads each input's declared unit from the spec, so the dimensioned
    # inputs arrive as quantities and `k` as a bare float.
    return choked_flow_area(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.a.magnitude,
        h.expected(case, "a"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (a)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "a":
            return result.a.magnitude
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_monotonic() -> None:
    """The area rises with mass flow and falls with pressure, density and ``k``.

    The last is the least obvious and the most worth having: a higher isentropic
    exponent means a larger geometric factor and so a higher critical flux, so the
    same mass flow needs *less* area. A form that got the exponent's sign or its
    numerator wrong would move that direction.
    """

    def at(m: float, p: float, r: float, k: float) -> float:
        return choked_flow_area(Q(m, "kg/s"), Q(p, "Pa"), Q(r, "kg/m**3"), k).a.magnitude

    base = at(1.0, 1.0e6, 10.0, 1.4)
    assert at(2.0, 1.0e6, 10.0, 1.4) > base, "a must rise with m_dot"
    assert at(1.0, 5.0e5, 10.0, 1.4) > base, "a must rise as P0 falls"
    assert at(1.0, 1.0e6, 5.0, 1.4) > base, "a must rise as rho0 falls"
    assert at(1.0, 1.0e6, 10.0, 5.0 / 3.0) < base, "a must fall as k rises"
    assert at(0.0, 1.0e6, 10.0, 1.4) == 0.0, "no flow, no area"


def test_unit_round_trip() -> None:
    """The same throat described in SI and in US customary units.

    All three dimensioned inputs are converted, through pint's own conversion rather
    than through factors restated here. Writing them out is how the pump_power round
    trip first failed - the factors were good to 7e-8 against a 1e-12 tolerance.
    """
    si = choked_flow_area(Q(1.0, "kg/s"), Q(1.0e6, "Pa"), Q(10.0, "kg/m**3"), 1.4).a
    us = choked_flow_area(
        Q(1.0, "kg/s").to("lb/s"),
        Q(1.0e6, "Pa").to("psi"),
        Q(10.0, "kg/m**3").to("lb/ft**3"),
        1.4,
    ).a

    h.assert_close(us.magnitude, si.magnitude, 1e-12, "a from SI vs US customary units")


def test_one_square_metre_comes_back_for_the_critical_flux() -> None:
    """The physics, checked without restating the implementation.

    At k = 1.4 the exponent is exactly 3 and the geometric factor is exactly
    (5/6)**3 = 125/216, so with unit density and unit pressure the critical mass flux
    is exactly ``sqrt(1.4) * 125/216``. Feeding that in as the mass flow must give
    exactly one square metre of throat.

    Stated as "the area must be 1" rather than by recomputing the expression: a test
    that restates the implementation cannot catch the implementation being wrong, and
    this one moves off 1 for a wrong exponent, a misplaced root, or a unit that did
    not convert.
    """
    critical_flux = math.sqrt(1.4) * 125.0 / 216.0
    result = choked_flow_area(Q(critical_flux, "kg/s"), Q(1.0, "Pa"), Q(1.0, "kg/m**3"), 1.4)
    h.assert_close(result.a.magnitude, 1.0, 1e-15, "the critical flux needs unit area")


def test_agrees_with_the_stagnation_temperature_form() -> None:
    """The other algebraic form of the same physics, which shares no arithmetic.

    The relation is usually written with stagnation temperature and the gas constant:
    ``G* = P0 sqrt(k/(R T0)) (2/(k+1))**((k+1)/(2(k-1)))``. That is this calc's
    expression rewritten through ``rho0 = P0/(R T0)``, so feeding in that flux must
    give one square metre here.

    A cross-check against an independent *form* rather than the same expression
    evaluated twice - which is the difference between a test that can fail and one
    that cannot.
    """
    p0, t0, r_gas, k = 1.0e6, 300.0, 287.0, 1.4
    textbook_flux = (
        p0 * math.sqrt(k / (r_gas * t0)) * (2.0 / (k + 1.0)) ** ((k + 1.0) / (2.0 * (k - 1.0)))
    )
    result = choked_flow_area(
        Q(textbook_flux, "kg/s"), Q(p0, "Pa"), Q(p0 / (r_gas * t0), "kg/m**3"), k
    )
    h.assert_close(result.a.magnitude, 1.0, 1e-9, "the textbook flux needs unit area")


def test_an_exponent_above_the_monatomic_limit_warns_but_computes() -> None:
    """The first warning bound in five calcs, and why it exists here.

    ``k`` is a property of a substance rather than a coefficient in a convention, and
    substances have a known range: above 5/3 is above every real gas. The value is
    still returned, because the relation evaluates perfectly well and refusing would
    be less useful than saying so.
    """
    result = choked_flow_area(Q(1.0, "kg/s"), Q(1.0e6, "Pa"), Q(10.0, "kg/m**3"), 1.8)
    assert result.has_warning(WarningCode.OUT_OF_VALID_RANGE)
    assert result.a.magnitude > 0.0


def test_the_monatomic_limit_itself_is_clean() -> None:
    """5/3 is the boundary and is inclusive: argon and helium, not a mistake."""
    result = choked_flow_area(Q(1.0, "kg/s"), Q(1.0e6, "Pa"), Q(10.0, "kg/m**3"), 5.0 / 3.0)
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_an_exponent_of_one_is_an_error_rather_than_a_division_by_zero() -> None:
    """``k - 1`` is a denominator in the exponent, and no gas has equal specific
    heats."""
    with pytest.raises(OutOfRangeError) as excinfo:
        choked_flow_area(Q(1.0, "kg/s"), Q(1.0e6, "Pa"), Q(10.0, "kg/m**3"), 1.0)
    assert excinfo.value.field() == "k"


def test_non_positive_pressure_or_density_is_an_error() -> None:
    for p, r in [(0.0, 10.0), (1.0e6, 0.0)]:
        with pytest.raises(OutOfRangeError):
            choked_flow_area(Q(1.0, "kg/s"), Q(p, "Pa"), Q(r, "kg/m**3"), 1.4)


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
