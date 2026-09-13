"""Spec-driven tests for ``eos.pr_molar_volume``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PrMolarVolumeResult
from azoth.eos import pr_molar_volume
from azoth.eos.reference import MOLAR_GAS_CONSTANT

CALC_ID = "eos.pr_molar_volume"
Q = ureg.Quantity

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "monotonic": "test_monotonic",
    "unit_round_trip": "test_unit_round_trip",
}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> PrMolarVolumeResult:
    return pr_molar_volume(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.v.magnitude,
        h.expected(case, "v"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (v)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "v":
            return float(result.v.magnitude)
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_the_gas_constant_is_the_exact_product() -> None:
    """``R`` is the exact product, not the truncated literal everyone quotes.

    The difference is invisible at a glance - `8.314462618` and `8.31446261815324` are
    both "the gas constant" to a reader - so it is pinned here rather than left to a
    point-value case with a loose tolerance to notice.
    """
    assert MOLAR_GAS_CONSTANT == 8.31446261815324
    assert MOLAR_GAS_CONSTANT > 8.314462618


def test_monotonic() -> None:
    """``v`` rises with ``z`` and ``T``, and falls with ``P``.

    The three partial derivatives. A point-value case cannot tell a correct formula
    from one with ``T`` and ``P`` transposed; three directions at once can.
    """

    def at(z: float, t: float, p: float) -> float:
        return pr_molar_volume(z, Q(t, "K"), Q(p, "Pa")).v.magnitude

    base = at(0.8, 300.0, 1.0e5)
    assert at(0.9, 300.0, 1.0e5) > base
    assert at(0.8, 400.0, 1.0e5) > base
    assert at(0.8, 300.0, 2.0e5) < base
    # And it is exactly linear in z, which a bound-check or a plot would not say.
    assert at(1.6, 300.0, 1.0e5) == 2.0 * base


def test_unit_round_trip() -> None:
    """The same physical state in converted units gives the same molar volume.

    Unlike the rest of this namespace, this calc is dimensional, so the round trip has
    something to convert - and it is the first test to exercise ``m**3/mol`` through a
    real calculation rather than only in the units module's own tests.
    """
    z, T, P = 0.7907789662973796, 295.864, 1_062_000.0
    si = pr_molar_volume(z, Q(T, "K"), Q(P, "Pa")).v
    converted = pr_molar_volume(z, Q(T - 273.15, "degC"), Q(P / 1e5, "bar")).v
    h.assert_close(si.magnitude, converted.magnitude, 1e-12, "Pa vs bar, K vs degC")


@pytest.mark.parametrize(
    ("z", "t", "p"), [(0.0, 300.0, 1e5), (-0.5, 300.0, 1e5), (0.8, 0.0, 1e5), (0.8, 300.0, 0.0)]
)
def test_a_non_positive_input_is_an_error(z: float, t: float, p: float) -> None:
    with pytest.raises(OutOfRangeError):
        pr_molar_volume(z, Q(t, "K"), Q(p, "Pa"))


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
