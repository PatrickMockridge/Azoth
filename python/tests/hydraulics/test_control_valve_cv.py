"""Spec-driven tests for ``hydraulics.control_valve_cv``."""

from __future__ import annotations

import math
from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import ControlValveCvResult
from azoth.hydraulics import control_valve_cv
from azoth.hydraulics.reference.control_valve_cv import CV_TO_SI, GALLON_M3, PSI_PA

CALC_ID = "hydraulics.control_valve_cv"
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


def call(case: dict[str, Any]) -> ControlValveCvResult:
    # `kwargs_for` reads each input's declared unit from the spec, so `dP` arrives as
    # a quantity and `Cv`/`SG` as bare floats, with no unit handling here.
    return control_valve_cv(**h.kwargs_for(SPEC, case["inputs"]))


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
    """Flow rises with the coefficient and the drop, falls with gravity, and is zero
    when the drop is.

    The directions place each input in the expression: a wrong power, a misplaced
    square root or an inverted ratio all show up as a direction changing, and the
    zero case is the one a form that divided by ``dP`` would get wrong.
    """

    def at(cv: float, dp: float, sg: float) -> float:
        return control_valve_cv(cv, Q(dp, "Pa"), sg).q.magnitude

    base = at(10.0, 100000.0, 1.0)
    assert at(20.0, 100000.0, 1.0) > base, "q must rise with Cv"
    assert at(10.0, 200000.0, 1.0) > base, "q must rise with dP"
    assert at(10.0, 100000.0, 1.5) < base, "q must fall as SG rises"
    assert at(10.0, 0.0, 1.0) == 0.0, "no drop, no flow"


def test_unit_round_trip() -> None:
    """The same valve described with the drop in psi rather than pascals.

    ``Cv`` and ``SG`` are dimensionless and cannot round-trip, so ``dP`` is the only
    input here carrying a unit. Converting it is a real check rather than a
    restatement of the coefficient convention: it runs the pressure through pint's
    conversion and back into the relation.
    """
    si = control_valve_cv(10.0, Q(100000.0, "Pa"), 1.0).q
    us = control_valve_cv(10.0, Q(100000.0, "Pa").to("psi"), 1.0).q
    h.assert_close(us.magnitude, si.magnitude, 1e-12, "q from Pa vs psi")


def test_the_conversion_constant_is_derived_from_the_definitions() -> None:
    """The constant is arithmetic, not a fitted value.

    Python computes it at import from the definitions of the gallon and the
    pound-force rather than storing a literal, so this re-derives it independently
    and checks the module's value agrees. If either the definitions or the constant
    were edited, one of the two would move and this would report it.
    """
    independent = (3.785411784e-3 / 60.0) / math.sqrt(4.4482216152605 / 0.0254**2)
    assert pytest.approx(independent, rel=0, abs=1e-21) == CV_TO_SI
    assert GALLON_M3 == 3.785411784e-3
    assert pytest.approx(6894.757293168361, rel=0, abs=1e-9) == PSI_PA


def test_the_conversion_constant_matches_the_rust_literal() -> None:
    """Both languages must use the same number, and results alone would not prove it.

    Python derives the constant and Rust stores it as a literal, because `const`
    cannot call `sqrt`. A disagreement smaller than a test's tolerance would pass
    every cross-language comparison while being a real difference - so the literal
    is pinned here as well as on the Rust side.
    """
    assert pytest.approx(7.59805421208337e-07, rel=0, abs=1e-21) == CV_TO_SI


def test_a_coefficient_of_one_passes_one_gallon_per_minute() -> None:
    """The definition of ``Cv``, stated directly: 1 gpm of water at 1 psi.

    A check on the conversion constant rather than on the arithmetic - if it were
    wrong by any factor, the flow at unit coefficient and unit drop would not be one
    gallon per minute.
    """
    result = control_valve_cv(1.0, Q(PSI_PA, "Pa"), 1.0)
    h.assert_close(result.q.magnitude, GALLON_M3 / 60.0, 1e-15, "one Cv at one psi")
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


def test_a_non_positive_coefficient_is_an_error() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        control_valve_cv(0.0, Q(100000.0, "Pa"), 1.0)
    assert excinfo.value.field() == "Cv"


def test_a_negative_drop_is_an_error_rather_than_reversed_flow() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        control_valve_cv(10.0, Q(-1000.0, "Pa"), 1.0)
    assert excinfo.value.field() == "dP"


def test_a_non_positive_specific_gravity_is_an_error() -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        control_valve_cv(10.0, Q(100000.0, "Pa"), 0.0)
    assert excinfo.value.field() == "SG"


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
