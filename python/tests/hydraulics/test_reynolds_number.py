"""Spec-driven tests for ``hydraulics.reynolds_number``.

Every case comes from ``specs/calcs/hydraulics/reynolds_number.toml``. Add a test
to that file and it runs here and in Rust with no new test code.
"""

from __future__ import annotations

from collections.abc import Callable
from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import FlowRegime, ReynoldsNumberResult
from azoth.core.warnings import WarningCode
from azoth.hydraulics import reynolds_number

CALC_ID = "hydraulics.reynolds_number"
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


def call(case: dict[str, Any]) -> ReynoldsNumberResult:
    return reynolds_number(
        Q(h.input_(case, "rho"), "kg/m**3"),
        Q(h.input_(case, "v"), "m/s"),
        Q(h.input_(case, "D"), "m"),
        Q(h.input_(case, "mu"), "Pa*s"),
    )


def resolve(case: dict[str, Any], result: ReynoldsNumberResult) -> Callable[[str], float | None]:
    def lookup(quantity: str) -> float | None:
        if quantity == "re":
            return result.re
        return h.input_(case, quantity) if quantity in case["inputs"] else None

    return lookup


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.re,
        h.expected(case, "re"),
        case.get("tolerance", 1e-12),
        f"{CALC_ID}::{case['id']} (re)",
    )
    h.assert_consistent(result, f"{CALC_ID}::{case['id']}")
    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve(case, result), case["id"])
    if case["type"] == "worked_example":
        assert result.regime is FlowRegime.from_reynolds_number(result.re), (
            "regime must follow from the Reynolds number"
        )


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 3, f"expected several active cases, found {total}"


def test_skipped_cases_say_why() -> None:
    """A skip with no reason is indistinguishable from an oversight."""
    h.assert_skips_are_explained(SPEC)


def test_unit_round_trip() -> None:
    """The same physical state in SI and US customary units gives the same Re.

    Reynolds number is dimensionless, so this is the sharpest check available
    that unit handling is right: any stray conversion factor that survived into
    the ratio would show up as a discrepancy rather than cancelling.
    """
    si = reynolds_number(Q(998.0, "kg/m**3"), Q(1.5, "m/s"), Q(0.1, "m"), Q(1.002e-3, "Pa*s"))
    # The identical state with length and velocity in feet.
    us = reynolds_number(
        Q(998.0, "kg/m**3"), Q(1.5 / 0.3048, "ft/s"), Q(0.1 / 0.3048, "ft"), Q(1.002e-3, "Pa*s")
    )
    h.assert_close(us.re, si.re, 1e-12, "Re from SI vs US customary")
    assert si.regime is us.regime


@pytest.mark.parametrize(
    ("field", "kwargs"),
    [
        ("rho", {"rho": Q(0.0, "kg/m**3")}),
        ("v", {"v": Q(-1.0, "m/s")}),
        ("D", {"D": Q(0.0, "m")}),
        ("mu", {"mu": Q(0.0, "Pa*s")}),
    ],
)
def test_hard_bounds_raise_rather_than_warn(field: str, kwargs: dict[str, Any]) -> None:
    """Each of these makes the Reynolds number undefined or meaningless.

    A warning here would return ``inf`` dressed as a result.
    """
    args: dict[str, Any] = {
        "rho": Q(998.0, "kg/m**3"),
        "v": Q(1.5, "m/s"),
        "D": Q(0.1, "m"),
        "mu": Q(1.002e-3, "Pa*s"),
    }
    args.update(kwargs)
    with pytest.raises(OutOfRangeError) as excinfo:
        reynolds_number(**args)
    assert excinfo.value.field() == field


def test_laminar_flow_gets_no_regime_warning() -> None:
    """The transitional band is the only regime that warns.

    Laminar and turbulent are both determinate, so warning about them would be
    noise - and noise is how a caller learns to ignore the warnings field.
    """
    result = reynolds_number(Q(1000.0, "kg/m**3"), Q(0.01, "m/s"), Q(0.05, "m"), Q(1e-2, "Pa*s"))
    assert result.regime is FlowRegime.LAMINAR
    assert result.is_clean, result.warnings
    assert not result.has_warning(WarningCode.TRANSITIONAL_FLOW)


def test_transitional_flow_warns_with_the_specific_code() -> None:
    """The spec names TRANSITIONAL_FLOW for this band.

    A caller can then react to indeterminate friction without string-matching a
    message.
    """
    result = reynolds_number(Q(1000.0, "kg/m**3"), Q(0.1, "m/s"), Q(0.05, "m"), Q(1.5e-3, "Pa*s"))
    assert result.regime is FlowRegime.TRANSITIONAL
    assert result.has_warning(WarningCode.TRANSITIONAL_FLOW)
    assert not result.has_warning(WarningCode.OUT_OF_VALID_RANGE), (
        "the band check names its own code; the generic one must not also fire"
    )


def test_nan_input_is_rejected() -> None:
    """NaN comparisons are all false, so without an explicit check a NaN would
    sail through every bound and produce a confident wrong answer."""
    with pytest.raises(OutOfRangeError):
        reynolds_number(Q(float("nan"), "kg/m**3"), Q(1.5, "m/s"), Q(0.1, "m"), Q(1e-3, "Pa*s"))


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
