"""Spec-driven tests for ``eos.rachford_rice_binary``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth.core.errors import OutOfRangeError
from azoth.core.result import RachfordRiceBinaryResult
from azoth.eos import rachford_rice_binary

CALC_ID = "eos.rachford_rice_binary"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {"consistency_with": "test_consistency_with"}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def rachford_rice_sum(z1: float, k1: float, k2: float, beta: float) -> float:
    """The Rachford-Rice sum, written out here rather than imported.

    The calc returns the closed form; this is the equation the closed form solves.
    Keeping them separate is the point - the residual test compares the answer against
    the *equation*, not against the algebra that produced it, so an error in the
    derivation cannot hide by being present in both.
    """
    return z1 * (k1 - 1.0) / (1.0 + beta * (k1 - 1.0)) + (1.0 - z1) * (k2 - 1.0) / (
        1.0 + beta * (k2 - 1.0)
    )


def call(case: dict[str, Any]) -> RachfordRiceBinaryResult:
    return rachford_rice_binary(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.beta,
        h.expected(case, "beta"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (beta)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "beta":
            return result.beta
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_consistency_with() -> None:
    """The returned ``beta`` solves the Rachford-Rice equation, over a broad sweep.

    The check that does not depend on knowing the right answer, which is what matters
    for a closed form with several roundings in it - the worked example's expected
    value is one ulp below the exact 2/3 and cannot be otherwise.

    Swept **including pairs whose beta falls outside [0, 1]**: the equation still has
    a solution there and the residual must still vanish. A solver that clamped beta to
    the physical interval, or that bailed out when the root left it, fails here while
    passing every two-phase point-value case.
    """
    checked = 0
    outside = 0
    for z1_step in range(1, 10):
        for k1 in (1.5, 2.0, 4.0, 8.0, 30.0):
            for k2 in (0.02, 0.2, 0.5, 0.8, 0.95):
                z1 = z1_step / 10.0
                result = rachford_rice_binary(z1, k1, k2)
                residual = rachford_rice_sum(z1, k1, k2, result.beta)
                assert abs(residual) < 1e-12, (
                    f"z1={z1}, K1={k1}, K2={k2}: beta={result.beta} leaves {residual:e}"
                )
                if not 0.0 <= result.beta <= 1.0:
                    outside += 1
                checked += 1

    assert checked > 100, f"the sweep should be broad, checked {checked}"
    assert outside > 0, "the sweep should reach single-phase feeds; none did"


def test_a_feed_that_does_not_split_carries_a_warning() -> None:
    """A feed whose K-values all exceed 1 does not split, and says so.

    The warning is on the *output* rather than on the inputs because whether the feed
    splits depends on all three. Clamping ``beta`` to ``[0, 1]`` would pass a value
    check and fail here, and would also discard the information that the feed is
    subcooled liquid rather than merely not-split.
    """
    result = rachford_rice_binary(0.5, 2.0, 1.5)
    assert result.beta < 0.0
    assert not result.is_clean

    vapour = rachford_rice_binary(0.5, 0.8, 0.5)
    assert vapour.beta > 1.0
    assert not vapour.is_clean


def test_a_two_phase_feed_carries_no_warning() -> None:
    for z1, k1, k2 in ((0.6, 4.0, 0.25), (0.3, 5.0, 0.2), (0.5, 3.0, 0.3)):
        result = rachford_rice_binary(z1, k1, k2)
        assert 0.0 <= result.beta <= 1.0, f"expected two-phase for {z1}, {k1}, {k2}"
        assert result.is_clean, f"unexpected warnings: {result.warnings}"


@pytest.mark.parametrize(("k1", "k2"), [(1.0, 0.5), (2.0, 1.0)])
def test_a_k_value_of_one_is_an_error(k1: float, k2: float) -> None:
    """``K - 1`` is a divisor, and the component degenerates besides.

    This is the case the `equals` bound exists for, and it is worth noting that the
    bound did not work when this spec was first written: the schema permitted
    ``equals``, the range-check machinery implemented only min/max, and the check
    silently never fired - so the calc returned ``-inf`` with a warning instead of
    raising. See `crates/azoth-core/src/range.rs`.
    """
    with pytest.raises(OutOfRangeError):
        rachford_rice_binary(0.5, k1, k2)


@pytest.mark.parametrize(("k1", "k2"), [(0.0, 0.5), (-1.0, 0.5), (2.0, 0.0)])
def test_a_non_positive_k_value_is_an_error(k1: float, k2: float) -> None:
    with pytest.raises(OutOfRangeError):
        rachford_rice_binary(0.5, k1, k2)


@pytest.mark.parametrize("z1", [-0.01, 1.01])
def test_a_mole_fraction_outside_zero_to_one_is_an_error(z1: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        rachford_rice_binary(z1, 4.0, 0.25)
    assert excinfo.value.field() == "z1"


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


def test_the_skipped_case_says_why_it_is_skipped() -> None:
    h.assert_skips_are_explained(SPEC)
    declared = {str(c["property"]) for c in SKIPPED if c.get("property")}
    assert "unit_round_trip" in declared, "this spec skips the unit round trip on purpose"
