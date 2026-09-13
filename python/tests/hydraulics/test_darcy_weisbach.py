"""Spec-driven tests for ``hydraulics.darcy_weisbach``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from chemeng import ureg
from chemeng.core.errors import OutOfRangeError
from chemeng.core.result import DarcyWeisbachResult, FlowRegime
from chemeng.core.warnings import WarningCode
from chemeng.hydraulics import darcy_weisbach
from chemeng.hydraulics.reference.darcy_weisbach import add_fitting_loss

CALC_ID = "hydraulics.darcy_weisbach"
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


def call(case: dict[str, Any]) -> DarcyWeisbachResult:
    """The spec declares ``mu`` optional, so its absence is a legitimate case
    rather than a missing argument."""
    inputs = case["inputs"]
    mu = Q(float(inputs["mu"]), "Pa*s") if "mu" in inputs else None
    return darcy_weisbach(
        h.input_(case, "f"),
        Q(h.input_(case, "L"), "m"),
        Q(h.input_(case, "D"), "m"),
        Q(h.input_(case, "rho"), "kg/m**3"),
        Q(h.input_(case, "v"), "m/s"),
        mu,
    )


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.dp.magnitude,
        h.expected(case, "dp"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (dp)",
    )
    if "re" in case["expected"]:
        assert result.re is not None, f"{case['id']}: re expected but missing"
        h.assert_close(
            result.re, h.expected(case, "re"), case.get("tolerance", 1e-12), f"{case['id']} (re)"
        )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "re":
            # None when mu was omitted, which is exactly the case the spec's range
            # check has to report as unchecked.
            return result.re
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_the_unverifiable_sources_are_skipped_not_deleted() -> None:
    """The two source_needed references must be visibly skipped, with reasons.

    A skipped test that vanishes silently is indistinguishable from one that was
    never written, and the whole point of the convention is that a gap stays
    visible.
    """
    h.assert_skips_are_explained(SPEC)
    assert len(SKIPPED) == 2, f"expected the two source_needed references, found {len(SKIPPED)}"
    for case in SKIPPED:
        assert "source_needed" in case["skip_reason"], (
            f"{case['id']} is skipped but not marked source_needed"
        )


def test_unit_round_trip() -> None:
    """The same physical state with length and velocity in feet gives the same dp."""
    si = darcy_weisbach(0.02, Q(100.0, "m"), Q(0.1, "m"), Q(998.0, "kg/m**3"), Q(1.5, "m/s"))
    us = darcy_weisbach(
        0.02,
        Q(100.0 / 0.3048, "ft"),
        Q(0.1 / 0.3048, "ft"),
        Q(998.0, "kg/m**3"),
        Q(1.5 / 0.3048, "ft/s"),
    )
    h.assert_close(us.dp.magnitude, si.dp.magnitude, 1e-12, "dp from SI vs US customary")


def test_omitting_viscosity_marks_the_regime_unchecked() -> None:
    """The heart of the optional-input design: "checked and fine" and "never
    checked" must not look the same."""
    unchecked = darcy_weisbach(0.02, Q(100.0, "m"), Q(0.1, "m"), Q(998.0, "kg/m**3"), Q(1.5, "m/s"))
    assert unchecked.re is None
    assert unchecked.regime is None
    assert unchecked.has_warning(WarningCode.RANGE_CHECK_SKIPPED)

    checked = darcy_weisbach(
        0.02,
        Q(100.0, "m"),
        Q(0.1, "m"),
        Q(998.0, "kg/m**3"),
        Q(1.5, "m/s"),
        Q(1.002e-3, "Pa*s"),
    )
    # Supplying mu must not change the pressure drop at all - it only enables the
    # check.
    assert unchecked.dp.magnitude == checked.dp.magnitude, "supplying mu must not change dp"
    assert checked.re is not None
    assert checked.regime is FlowRegime.TURBULENT
    assert not checked.has_warning(WarningCode.RANGE_CHECK_SKIPPED), (
        "with mu supplied the regime check did run, so it must not report as skipped"
    )


def test_transitional_flow_warns_through_this_calc_too() -> None:
    """Chosen so Re lands in the 2000-4000 band."""
    result = darcy_weisbach(
        0.032,
        Q(10.0, "m"),
        Q(0.05, "m"),
        Q(1000.0, "kg/m**3"),
        Q(0.1, "m/s"),
        Q(1.5e-3, "Pa*s"),
    )
    assert result.regime is FlowRegime.TRANSITIONAL
    assert result.has_warning(WarningCode.TRANSITIONAL_FLOW)


def test_zero_diameter_is_an_error() -> None:
    """L/D is singular at D = 0."""
    with pytest.raises(OutOfRangeError) as excinfo:
        darcy_weisbach(0.02, Q(100.0, "m"), Q(0.0, "m"), Q(998.0, "kg/m**3"), Q(1.5, "m/s"))
    assert excinfo.value.field() == "D"


def test_zero_length_is_allowed_and_gives_zero_drop() -> None:
    """Zero length is physical - a degenerate pipe - so it is not a hard bound."""
    result = darcy_weisbach(0.02, Q(0.0, "m"), Q(0.1, "m"), Q(998.0, "kg/m**3"), Q(1.5, "m/s"))
    assert result.dp.magnitude == 0.0


def test_pressure_drop_scales_with_density_and_velocity_squared() -> None:
    """The equation is linear in rho and quadratic in v.

    Checking the scaling catches a transcription error that a single point would
    not.
    """
    base = darcy_weisbach(
        0.02, Q(100.0, "m"), Q(0.1, "m"), Q(998.0, "kg/m**3"), Q(1.5, "m/s")
    ).dp.magnitude
    denser = darcy_weisbach(
        0.02, Q(100.0, "m"), Q(0.1, "m"), Q(1996.0, "kg/m**3"), Q(1.5, "m/s")
    ).dp.magnitude
    h.assert_close(denser, 2.0 * base, 1e-12, "dp is linear in rho")

    faster = darcy_weisbach(
        0.02, Q(100.0, "m"), Q(0.1, "m"), Q(998.0, "kg/m**3"), Q(3.0, "m/s")
    ).dp.magnitude
    h.assert_close(faster, 4.0 * base, 1e-12, "dp is quadratic in v")


def test_fitting_loss_adds_the_velocity_head_term() -> None:
    """The CLI's composition, exercised here so the CLI inherits a tested path.

    ``dP_total = f*(L/D)*(rho v**2/2) + K*(rho v**2/2)``
    """
    f, length, diameter, rho, velocity = 0.018, 100.0, 0.1, 998.0, 1.5
    straight = darcy_weisbach(
        f,
        Q(length, "m"),
        Q(diameter, "m"),
        Q(rho, "kg/m**3"),
        Q(velocity, "m/s"),
    ).dp.magnitude

    k_total = 0.684
    total = add_fitting_loss(straight, k_total, rho, velocity)
    velocity_head = rho * velocity**2 / 2.0

    h.assert_close(
        total, straight + k_total * velocity_head, 1e-15, "fitting loss is K times the head"
    )

    # Pin the actual magnitude rather than asserting a threshold chosen by feel.
    # For these inputs: velocity head = 998*1.5**2/2 = 1122.75 Pa, the fitting term
    # is 0.684 * 1122.75 = 767.961 Pa, and the straight-pipe term is
    # 0.018 * 1000 * 1122.75 = 20209.5 Pa, so the fittings add 3.8%.
    fitting_term = k_total * velocity_head
    h.assert_close(fitting_term, 767.961, 1e-9, "fitting term in Pa")
    ratio = fitting_term / straight
    assert 0.037 <= ratio <= 0.039, (
        f"fitting loss is {ratio:.4f} of the straight-pipe loss, expected about 0.038"
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
