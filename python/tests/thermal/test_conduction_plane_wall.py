"""Spec-driven tests for ``thermal.conduction_plane_wall``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError, UnitMismatchError
from azoth.core.result import ConductionPlaneWallResult
from azoth.thermal import conduction_plane_wall

CALC_ID = "thermal.conduction_plane_wall"
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


def call(case: dict[str, Any]) -> ConductionPlaneWallResult:
    # `h.kwargs_for` reads each input's declared unit from the spec, so this needs
    # no unit handling of its own and a spec change needs no change here.
    return conduction_plane_wall(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    expected = h.expected(case, "q")
    h.assert_close(result.q.magnitude, expected, case.get("tolerance", 1e-12), f"{case['id']} (q)")
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
    """``q`` rises with ``k``, ``A`` and ``dT``, and falls with ``L``.

    The four partial derivatives of the equation. A point-value test can pass with
    a wrong formula that happens to be right at one set of inputs; an inverted ratio
    or a sign error cannot survive all four directions at once.
    """

    def at(k: float, a: float, dt: float, length: float) -> float:
        return conduction_plane_wall(
            Q(k, "W/(m*K)"), Q(a, "m**2"), Q(dt, "K"), Q(length, "m")
        ).q.magnitude

    base = at(45.0, 2.0, 30.0, 0.05)
    assert at(90.0, 2.0, 30.0, 0.05) > base, "q must increase with k"
    assert at(45.0, 4.0, 30.0, 0.05) > base, "q must increase with A"
    assert at(45.0, 2.0, 60.0, 0.05) > base, "q must increase with dT"
    assert at(45.0, 2.0, 30.0, 0.10) < base, "q must decrease with L"


def test_unit_round_trip() -> None:
    """The same physical wall described in SI and in US customary units.

    Three of the four inputs are converted, mirroring the Rust half of this test:
    the area to square feet, the thickness to feet, and the temperature difference
    to a Fahrenheit interval. That last one is the interesting one - a *difference*
    of 30 K is 54 degF, while an absolute 30 K is -243.15 degC, so a test that
    reached for plain ``degF`` here would be off by hundreds of degrees and would
    say so loudly.

    The conductivity stays in SI because ``uom`` carries no US customary thermal
    conductivity, and the two languages' tests are kept identical rather than
    converting one more quantity in one of them.
    """
    si = conduction_plane_wall(Q(45.0, "W/(m*K)"), Q(2.0, "m**2"), Q(30.0, "K"), Q(0.05, "m")).q
    us = conduction_plane_wall(
        Q(45.0, "W/(m*K)"),
        Q(2.0 / 0.3048 / 0.3048, "ft**2"),
        Q(30.0 * 9.0 / 5.0, "delta_degF"),
        Q(0.05 / 0.3048, "ft"),
    ).q

    # `assert_same_state` exists on the Rust side's shared helper and not here; on
    # this side the two values are simply two floats to compare.
    h.assert_close(si.magnitude, us.magnitude, 1e-12, "q from SI vs US customary state")


def test_a_non_positive_conductivity_is_an_error() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        conduction_plane_wall(Q(0.0, "W/(m*K)"), Q(2.0, "m**2"), Q(30.0, "K"), Q(0.05, "m"))
    assert excinfo.value.field() == "k"


def test_a_zero_thickness_is_an_error_rather_than_infinite_heat_flow() -> None:
    """``L`` is a divisor, so zero is the singular limit.

    Returning infinity would be arithmetic; refusing is the honest answer, because
    a zero-thickness wall is not a state this model describes.
    """
    with pytest.raises(OutOfRangeError) as excinfo:
        conduction_plane_wall(Q(45.0, "W/(m*K)"), Q(2.0, "m**2"), Q(30.0, "K"), Q(0.0, "m"))
    assert excinfo.value.field() == "L"


def test_a_negative_temperature_difference_reverses_the_flow_instead_of_failing() -> None:
    """The one modelling decision this calc makes, and the reason ``dT`` is
    unbounded below: it models a difference across a slab, not a named hot face, so
    the signed answer is more useful than a refusal.

    A result carrying no warnings is the point - a negative heat flow is a
    direction, not a problem.
    """
    result = conduction_plane_wall(Q(45.0, "W/(m*K)"), Q(2.0, "m**2"), Q(-30.0, "K"), Q(0.05, "m"))
    h.assert_close(result.q.magnitude, -54000.0, 1e-12, "reversed dT")
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_an_absolute_temperature_is_refused_rather_than_offset() -> None:
    """A difference and an absolute temperature are the same dimension and are not
    interchangeable, so the absolute one is rejected.

    ``dT`` is declared ``interval: true`` in the spec. Without that, ``to_si`` would
    convert an *absolute* ``Q(30, "degC")`` to 303.15 K rather than to 30, and this
    calc would return 545670 W where the caller meant 54000 - a plausible number
    wrong by a factor of ten, from a unit that looks right.

    Rust never had this problem, because ``dT`` takes a ``TemperatureInterval`` and a
    ``ThermodynamicTemperature`` does not fit. This test is what holds the Python
    side to the same promise, on both backends, since the refusal must not depend on
    which one answered.
    """
    right = conduction_plane_wall(
        Q(45.0, "W/(m*K)"), Q(2.0, "m**2"), Q(30.0, "delta_degC"), Q(0.05, "m")
    ).q.magnitude
    assert right == pytest.approx(54000.0)

    for absolute in (Q(30.0, "degC"), Q(86.0, "degF")):
        with pytest.raises(UnitMismatchError) as excinfo:
            conduction_plane_wall(Q(45.0, "W/(m*K)"), Q(2.0, "m**2"), absolute, Q(0.05, "m"))
        assert excinfo.value.field() == "dT"
        # The message has to carry the fix, not just the refusal: the caller's
        # mistake is a plausible one and "wrong dimensions" would not explain it.
        assert "delta_degC" in str(excinfo.value)


def test_a_scaled_but_unoffset_temperature_unit_still_works() -> None:
    """`delta_degF` is not an offset unit and must not be caught by the interval check.

    This is the case that makes the check's *implementation* matter. An offset unit
    is one where zero of it is not zero of its base - `degC` is 273.15 K at zero.
    Testing "one of this unit is not one kelvin" instead would also reject
    `delta_degF`, whose one is 5/9 K, and that is the unit a US-customary caller is
    supposed to reach for.
    """
    celsius = conduction_plane_wall(
        Q(45.0, "W/(m*K)"), Q(2.0, "m**2"), Q(30.0, "delta_degC"), Q(0.05, "m")
    ).q.magnitude
    fahrenheit = conduction_plane_wall(
        Q(45.0, "W/(m*K)"), Q(2.0, "m**2"), Q(54.0, "delta_degF"), Q(0.05, "m")
    ).q.magnitude
    h.assert_close(fahrenheit, celsius, 1e-12, "54 degF interval is 30 degC interval")


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
