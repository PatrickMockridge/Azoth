"""Spec-driven tests for ``eos.vdw1f_mix_binary``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth.core.errors import OutOfRangeError
from azoth.core.result import Vdw1fMixBinaryResult
from azoth.eos import vdw1f_mix_binary

CALC_ID = "eos.vdw1f_mix_binary"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {"monotonic": "test_monotonic"}

#: A deliberately arbitrary `k12`. This library ships no fitted binary parameters,
#: so the tests pick a round number and say so rather than borrowing one that looks
#: like data.
SOME_K12 = 0.05

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]
FIELDS = ("a_mix", "b_mix")


def call(case: dict[str, Any]) -> Vdw1fMixBinaryResult:
    return vdw1f_mix_binary(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    for field in FIELDS:
        h.assert_close(
            getattr(result, field),
            h.expected(case, field),
            case.get("tolerance", 1e-12),
            f"{case['id']} ({field})",
        )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity in FIELDS:
            return float(getattr(result, quantity))
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 3, f"expected several active cases, found {total}"


def test_a_pure_component_reproduces_its_own_parameters() -> None:
    """The check that a mixing rule is a mixing rule.

    At ``z1 = 1`` every term carrying ``z2`` vanishes and the mixture parameters must
    be component 1's own; at ``z1 = 0``, component 2's. Swept over parameters rather
    than fixed at the spec's case, because it is the *structure* being tested: an
    implementation with a transposed weight passes a mid-range composition and fails
    here.
    """
    cases = (
        (
            0.20206500174625697,
            0.08448417260831159,
            0.02431127309496514,
            0.025932024634629486,
            SOME_K12,
        ),
        (0.115, 0.31, 0.018, 0.041, 0.4),
        (1.0, 2.0, 3.0, 4.0, 0.0),
    )
    for a1, a2, b1, b2, k12 in cases:
        pure1 = vdw1f_mix_binary(1.0, a1, a2, b1, b2, k12)
        assert pure1.a_mix == a1 and pure1.b_mix == b1, "z1 = 1 should give component 1"
        pure2 = vdw1f_mix_binary(0.0, a1, a2, b1, b2, k12)
        assert pure2.a_mix == a2 and pure2.b_mix == b2, "z1 = 0 should give component 2"


def test_b_mix_is_linear_in_composition() -> None:
    """``b_mix`` is the weighted mean, exactly, with no ``k12`` in it.

    Asserted bit-for-bit rather than approximately: it is linear in ``z1`` and has no
    square root, so a ``k12`` that leaked into the covolume would show up here and
    nowhere else. vdW1f has no ``l12``.
    """
    b1, b2 = 0.02431127309496514, 0.025932024634629486
    for z1 in (0.0, 0.25, 0.5, 0.6, 0.75, 1.0):
        result = vdw1f_mix_binary(z1, 0.2, 0.1, b1, b2, SOME_K12)
        assert result.b_mix == z1 * b1 + (1.0 - z1) * b2


def test_monotonic() -> None:
    """Each parameter moves in the direction its own pure values order.

    ``a_mix`` rises with ``z1`` when ``a1 > a2`` and falls when it does not, and the
    two need not agree - nothing relates one component's attraction to its covolume.
    The propane-like and methane-like pair has ``a1 > a2`` but ``b1 < b2``, so the two
    parameters move in *opposite* directions. An implementation with a transposed
    weight would move both the same way.
    """
    cases = (
        (0.20206500174625697, 0.08448417260831159, 0.02431127309496514, 0.025932024634629486),
        (0.31, 0.115, 0.041, 0.018),
        (0.05, 0.4, 0.01, 0.09),
    )
    for a1, a2, b1, b2 in cases:
        previous: tuple[float, float] | None = None
        for step in range(11):
            z1 = step / 10.0
            result = vdw1f_mix_binary(z1, a1, a2, b1, b2, SOME_K12)
            if previous is not None:
                prev_a, prev_b = previous
                assert (result.a_mix > prev_a) == (a1 > a2), f"a_mix direction at z1={z1}"
                assert (result.b_mix > prev_b) == (b1 > b2), f"b_mix direction at z1={z1}"
            previous = (result.a_mix, result.b_mix)


def test_an_unphysical_mixture_attraction_is_refused() -> None:
    """A ``k12`` outside ``[0, 2]`` can drive ``a_mix`` negative, and that is refused.

    The bound is on ``a_mix`` rather than on ``k12`` because a ``k12`` outside the
    band is not wrong in itself - most compositions still give a physical ``a_mix`` -
    and refusing the parameter would reject calls that are fine.
    """
    with pytest.raises(OutOfRangeError) as excinfo:
        vdw1f_mix_binary(0.5, 1.0, 1.0, 0.1, 0.1, 3.0)
    assert excinfo.value.field() == "a_mix"

    # The same k12 at a composition that stays physical is fine.
    assert vdw1f_mix_binary(0.99, 1.0, 1.0, 0.1, 0.1, 3.0).a_mix > 0.0


@pytest.mark.parametrize("z1", [-0.01, 1.01, 2.0])
def test_a_mole_fraction_outside_zero_to_one_is_an_error(z1: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        vdw1f_mix_binary(z1, 0.2, 0.1, 0.02, 0.03, SOME_K12)
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
