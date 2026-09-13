"""Spec-driven tests for ``hydraulics.friction_factor_colebrook``."""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pytest

import _helpers as h
from chemeng import ureg
from chemeng.core.errors import OutOfRangeError
from chemeng.core.result import ColebrookResult
from chemeng.core.warnings import WarningCode
from chemeng.hydraulics import (
    friction_factor_colebrook,
    friction_factor_swamee_jain,
    reynolds_number,
)
from chemeng.hydraulics.reference.friction_factor_colebrook import fully_rough_limit

CALC_ID = "hydraulics.friction_factor_colebrook"
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


def call(case: dict[str, Any]) -> ColebrookResult:
    return friction_factor_colebrook(h.input_(case, "re"), h.input_(case, "relative_roughness"))


def resolve(case: dict[str, Any], f: float) -> Callable[[str], float | None]:
    def lookup(quantity: str) -> float | None:
        if quantity == "f":
            return f
        return h.input_(case, quantity) if quantity in case["inputs"] else None

    return lookup


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    assert result.converged, (
        f"{case['id']}: the solver reported non-convergence, which should have been "
        f"raised as an error rather than returned"
    )
    h.assert_close(
        result.f, h.expected(case, "f"), case.get("tolerance", 1e-12), f"{case['id']} (f)"
    )
    h.assert_consistent(result, case["id"])
    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve(case, result.f), case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_unit_round_trip() -> None:
    """The same physical state in SI and US customary units gives the same f.

    Both Colebrook inputs are dimensionless, so this runs the whole chain -
    velocity and diameter in feet, through Reynolds number, into the friction
    factor - rather than converting a ratio that has no units to convert.
    """
    rr = 4.6e-4
    re_si = reynolds_number(Q(998.0, "kg/m**3"), Q(1.5, "m/s"), Q(0.1, "m"), Q(1.002e-3, "Pa*s")).re
    re_us = reynolds_number(
        Q(998.0, "kg/m**3"), Q(1.5 / 0.3048, "ft/s"), Q(0.1 / 0.3048, "ft"), Q(1.002e-3, "Pa*s")
    ).re
    f_si = friction_factor_colebrook(re_si, rr).f
    f_us = friction_factor_colebrook(re_us, rr).f
    h.assert_close(f_us, f_si, 1e-12, "f from SI vs US customary state")


def test_converges_to_the_fully_rough_asymptote() -> None:
    """Confirms the solver lands on the right fixed point.

    The fully-rough asymptote is explicit, so this checks the iteration against an
    independent value rather than against its own stopping rule.
    """
    rr = 0.01
    solved = friction_factor_colebrook(1e12, rr).f
    h.assert_close(solved, fully_rough_limit(rr), 1e-8, "fully rough limit at Re = 1e12")


def test_smooth_pipe_agrees_with_the_textbook_value() -> None:
    """At zero roughness the accepted value at Re = 1e5 is about 0.0180.

    This pins the implementation to an independent value rather than to itself.
    """
    f = friction_factor_colebrook(1e5, 0.0).f
    assert 0.0178 <= f <= 0.0182, f"smooth-pipe f at Re=1e5 was {f}, outside the textbook range"


def test_iteration_count_is_deterministic() -> None:
    """Deterministic iteration count is what makes the Python and Rust results
    comparable rather than merely close."""
    first = friction_factor_colebrook(1e5, 4.6e-4)
    second = friction_factor_colebrook(1e5, 4.6e-4)
    assert first.iterations == second.iterations
    assert first.f == second.f, "f must be bit-identical across runs"


def test_below_transition_warns_but_still_returns_a_number() -> None:
    """Out of the fitted range, not undefined: a warning, not an error."""
    result = friction_factor_colebrook(3000.0, 4.6e-4)
    assert result.has_warning(WarningCode.OUT_OF_VALID_RANGE)
    assert result.converged


@pytest.mark.parametrize(
    ("re_value", "roughness", "field"),
    [(0.0, 4.6e-4, "re"), (1e5, -1e-4, "relative_roughness")],
)
def test_undefined_inputs_raise(re_value: float, roughness: float, field: str) -> None:
    """These are undefined, not merely out of range.

    ``2.51/(Re*sqrt(f))`` is singular at Re = 0, and negative roughness is
    unphysical.
    """
    with pytest.raises(OutOfRangeError) as excinfo:
        friction_factor_colebrook(re_value, roughness)
    assert excinfo.value.field() == field


def test_comparable_to_swamee_jain_within_the_published_claim() -> None:
    """The two correlations in this slice must be consistent with each other.

    Swamee and Jain claim 1% agreement with Colebrook over 5000 < Re < 1e8, so
    that published figure is the tolerance rather than a number chosen here. The
    points are strictly inside the range - see the boundary test below.
    """
    for re_value, rr in [(1e4, 4.6e-4), (1e5, 4.6e-4), (1e6, 1e-3), (1e7, 5e-3)]:
        f_colebrook = friction_factor_colebrook(re_value, rr).f
        f_swamee = friction_factor_swamee_jain(re_value, rr).f
        h.assert_close(
            f_swamee, f_colebrook, 1e-2, f"swamee-jain vs colebrook at Re={re_value}, rr={rr}"
        )


def test_swamee_jain_exceeds_its_claim_at_the_lower_bound() -> None:
    """Measured, not assumed: the explicit approximation degrades near its bound.

    At Re = 5000 with epsilon/D = 4.6e-4 it is 1.39% from Colebrook, against a
    headline claim of 1%; by Re = 1e4 it is back to 0.59%.

    Recorded rather than hidden, so the published 1% is not quietly treated as
    uniform across the claimed range, and so a future change to the low-Re
    behaviour is noticed rather than absorbed by a loosened tolerance.
    """
    re_value, rr = 5000.0, 4.6e-4
    f_colebrook = friction_factor_colebrook(re_value, rr).f
    f_swamee = friction_factor_swamee_jain(re_value, rr).f
    relative = abs(f_swamee - f_colebrook) / f_colebrook
    assert 0.013 <= relative <= 0.015, (
        f"swamee-jain is {relative:.4f} from colebrook at Re=5000; expected about "
        f"0.0139. If this moved, the low-Re behaviour changed - check whether that "
        f"was intended before widening this range."
    )


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
