"""Spec-driven tests for ``eos.pr_mass_density``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PrMassDensityResult
from azoth.eos import pr_mass_density

CALC_ID = "eos.pr_mass_density"
Q = ureg.Quantity

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "consistency_with": "test_consistency_with",
    "unit_round_trip": "test_unit_round_trip",
}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> PrMassDensityResult:
    return pr_mass_density(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.rho.magnitude,
        h.expected(case, "rho"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (rho)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "rho":
            return float(result.rho.magnitude)
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 3, f"expected several active cases, found {total}"


def test_consistency_with() -> None:
    """``rho`` doubles with ``M`` and halves with ``v``.

    The two relations that pin the formula, using no reference value at all. A
    point-value case cannot tell ``M/v`` from ``M*v`` or from ``v/M``; these can.
    """
    m, v = 0.0440956, 0.0018317107825229842
    base = pr_mass_density(Q(m, "kg/mol"), Q(v, "m**3/mol")).rho.magnitude

    doubled_m = pr_mass_density(Q(2.0 * m, "kg/mol"), Q(v, "m**3/mol")).rho.magnitude
    h.assert_close(doubled_m, 2.0 * base, 1e-12, "rho should double with M")

    doubled_v = pr_mass_density(Q(m, "kg/mol"), Q(2.0 * v, "m**3/mol")).rho.magnitude
    h.assert_close(doubled_v, base / 2.0, 1e-12, "rho should halve with v")


def test_unit_round_trip() -> None:
    """The same state in converted units gives the same density.

    This is the first test to exercise ``kg/mol`` through a real calculation. The
    conversion factors are large - 1e-6 for ``cm**3/mol``, and 1e3 between ``g/mol``
    and ``kg/mol`` - so a direction error fails by orders of magnitude rather than
    marginally.
    """
    si = pr_mass_density(Q(0.0440956, "kg/mol"), Q(0.002, "m**3/mol")).rho
    converted = pr_mass_density(Q(44.0956, "g/mol"), Q(2000.0, "cm**3/mol")).rho
    h.assert_close(si.magnitude, converted.magnitude, 1e-12, "kg/mol+cm**3 vs g/mol+m**3")

    # And the g/mol mistake is a factor of a thousand, asserted so the input's
    # warning is attached to something rather than merely written down.
    in_grams = pr_mass_density(Q(44.0956, "kg/mol"), Q(0.002, "m**3/mol")).rho.magnitude
    h.assert_close(in_grams, 1000.0 * si.magnitude, 1e-6, "the g/mol/ kg/mol factor")


@pytest.mark.parametrize(("m", "v"), [(0.0, 1e-3), (-0.044, 1e-3), (0.044, 0.0), (0.044, -1e-3)])
def test_a_non_positive_input_is_an_error(m: float, v: float) -> None:
    with pytest.raises(OutOfRangeError):
        pr_mass_density(Q(m, "kg/mol"), Q(v, "m**3/mol"))


def test_property_tests_are_covered() -> None:
    """Every property test the spec declares must be one we actually run."""
    declared = {str(c["property"]) for c in PROPERTIES if c.get("property")}
    missing = declared - set(_PROPERTY_IMPLEMENTATIONS)
    assert not missing, (
        f"{CALC_ID}: spec declares property tests {sorted(missing)} that no test implements"
    )
    for property_name in sorted(declared):
        function = globals().get(_PROPERTY_IMPLEMENTATIONS[property_name])
        assert callable(function), f"{CALC_ID}: property {property_name!r} has no implementation"
        function()


def test_no_skipped_cases_are_inexplicable() -> None:
    h.assert_skips_are_explained(SPEC)
    assert not SKIPPED, "this spec has no skipped tests"
