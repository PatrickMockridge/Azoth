"""Spec-driven tests for ``eos.pr_alpha_ab``."""

from __future__ import annotations

from decimal import Decimal, getcontext
from itertools import pairwise
from typing import Any

import pytest

import _helpers as h
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PrAlphaAbResult
from azoth.eos import pr_alpha_ab
from azoth.eos.reference import OMEGA_A, OMEGA_B

CALC_ID = "eos.pr_alpha_ab"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "monotonic": "test_monotonic",
}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> PrAlphaAbResult:
    return pr_alpha_ab(**h.kwargs_for(SPEC, case["inputs"]))


#: Every quantity this calc touches, in the order the spec declares them.
FIELDS = ("alpha", "a_reduced", "b_reduced")


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


def test_the_specs_expected_values_are_arithmetically_correct() -> None:
    """Re-derive every expected value in exact decimal arithmetic.

    This exists because the ordinary test above cannot catch the failure it is
    aimed at. `test_spec_case` compares the implementation against the spec's
    `expected` block, and *both were written by the same person from the same
    working* - so an arithmetic slip in the derivation appears in both and the
    comparison passes. The only thing that catches it is a third route to the
    answer, which is what `decimal` at 50 digits is: not another implementation of
    the equation of state, but the same equation evaluated without any rounding to
    argue about.

    It is not hypothetical. Every intermediate in this spec's derivation was
    written by hand first and was wrong in the sixth significant figure; the
    decimals below found it. The spec makes a claim - that each value is the
    correctly rounded exact result - and this is that claim, asserted.
    """
    getcontext().prec = 50
    omega_a, omega_b = Decimal(repr(OMEGA_A)), Decimal(repr(OMEGA_B))

    checked = 0
    cases = [("worked_example", SPEC["worked_example"])] + [
        (case["id"], case) for case in ACTIVE if "expected" in case
    ]
    for case_id, case in cases:
        inputs = case["inputs"]
        kappa = Decimal(str(inputs["kappa"]))
        tr = Decimal(str(inputs["Tr"]))
        pr = Decimal(str(inputs["Pr"]))

        attraction = 1 + kappa * (1 - tr.sqrt())
        alpha = attraction**2
        exact = {
            "alpha": alpha,
            "a_reduced": omega_a * alpha * pr / (tr**2),
            "b_reduced": omega_b * pr / tr,
        }

        for field in FIELDS:
            expected = Decimal(repr(h.expected(case, field)))
            deviation = abs(expected - exact[field]) / abs(exact[field])
            assert deviation < Decimal("1e-15"), (
                f"{case_id} ({field}): the spec says {expected}, exact arithmetic gives "
                f"{exact[field]} - a relative deviation of {deviation:.2e}. "
                f"The derivation in the spec needs correcting, not the test."
            )
            checked += 1

    assert checked >= 6, f"expected several values to check, checked {checked}"


def test_monotonic() -> None:
    """The four relations the spec declares, sampled rather than spot-checked."""

    def at(tr: float, pr: float) -> PrAlphaAbResult:
        return pr_alpha_ab(0.60282728832, tr, pr)

    # Falls with Tr, all three of them.
    for low, high in pairwise([0.6, 0.7, 0.8, 0.9, 1.0]):
        for field in FIELDS:
            assert getattr(at(high, 0.25), field) < getattr(at(low, 0.25), field), (
                f"{field} must fall as Tr rises: {low} -> {high}"
            )

    # Rises with Pr, and linearly for the two reduced parameters - which is what
    # `* Pr` appearing once and not as an exponent means.
    base = at(0.8, 0.25)
    for pr in [0.5, 1.0, 4.0]:
        scaled = at(0.8, pr)
        assert scaled.alpha == base.alpha, "alpha does not depend on Pr at all"
        for field in ("a_reduced", "b_reduced"):
            expected = getattr(base, field) * (pr / 0.25)
            h.assert_close(getattr(scaled, field), expected, 1e-12, f"{field} at Pr={pr}")


def depressed(omega_a: float, omega_b: float) -> tuple[float, float]:
    """`(p, q)` of the depressed reduced cubic at `A = omega_a`, `B = omega_b`.

    `z**3 - (1-B)z**2 + (A - 3B**2 - 2B)z - (A*B - B**2 - B**3)` written as
    `z**3 + c2*z**2 + c1*z + c0` and depressed by `z = w - c2/3`. A triple root
    leaves both zero and nothing else does.
    """
    c2 = -(1.0 - omega_b)
    c1 = omega_a - 3.0 * omega_b**2 - 2.0 * omega_b
    c0 = -(omega_a * omega_b - omega_b**2 - omega_b**3)
    return c1 - c2 * c2 / 3.0, 2.0 * c2**3 / 27.0 - c2 * c1 / 3.0 + c0


def test_the_exact_omegas_are_the_triple_root_solution() -> None:
    """The pair the cubic's condition gives, which this library does *not* ship.

    Asserted because it is where the constants come from and what the departure
    below is measured against - not because it is what the library uses.
    """
    p, q = depressed(0.4572355289213822, 0.07779607390388846)
    assert abs(p) < 1e-12, f"the exact pair should be a triple root: p = {p:e}"
    assert abs(q) < 1e-12, f"the exact pair should be a triple root: q = {q:e}"


def test_the_shipped_omegas_are_neqsims_and_are_not_the_triple_root_solution() -> None:
    """NeqSim's literals, deliberately, and the departure is on record.

    The port rule settles it: azoth is a port of NeqSim 3.20.0, NeqSim's
    `ComponentPR` constructor sets these two, and carrying them is what makes the
    two agree - see
    `validation/eos/methane_butane_flash_against_neqsim.json`. Asserted rather than
    commented because "correcting" them back to the exact pair is the most natural
    edit a reader could make, and the one that would end the agreement.
    """
    assert (OMEGA_A, OMEGA_B) == (0.45724333333, 0.077803333)
    p, q = depressed(OMEGA_A, OMEGA_B)
    assert abs(p) > 1e-9 and abs(q) > 1e-9, (
        f"the shipped pair was expected to depart from the triple root; p = {p:e}, q = {q:e}"
    )
    assert abs(p) < 1e-5 and abs(q) < 1e-5, (
        f"the departure is larger than NeqSim's literals explain; p = {p:e}, q = {q:e}"
    )


def test_the_printed_omegas_depart_from_the_triple_root_too() -> None:
    """The paper's five digits are a third pair, and they break it as well.

    Not the test that keeps the constants honest - the equality above is. This
    records why the printed pair is not an option either, and it is the same reason
    in kind: a triple root is cubically ill-conditioned.
    """
    p, q = depressed(0.45724, 0.07780)
    assert abs(p) > 1e-9 and abs(q) > 1e-9, (
        f"the printed pair was expected to break it; p = {p:e}, q = {q:e}"
    )


def test_alpha_is_exactly_one_at_the_critical_temperature() -> None:
    """``1 - sqrt(1)`` is exactly zero, so the square is exactly 1 for any kappa."""
    for kappa in (-0.5, 0.0, 0.60282728832, 3.0):
        assert pr_alpha_ab(kappa, 1.0, 1.0).alpha == 1.0, f"for kappa = {kappa}"


@pytest.mark.parametrize("tr", [0.0, -0.5])
def test_a_non_positive_reduced_temperature_is_an_error(tr: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        pr_alpha_ab(0.60282728832, tr, 0.25)
    assert excinfo.value.field() == "Tr"


@pytest.mark.parametrize("pr", [0.0, -0.1])
def test_a_non_positive_reduced_pressure_is_an_error(pr: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        pr_alpha_ab(0.60282728832, 0.8, pr)
    assert excinfo.value.field() == "Pr"


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
    """A skipped test is still reported, and has to carry its reason."""
    h.assert_skips_are_explained(SPEC)
    declared = {str(c["property"]) for c in SKIPPED if c.get("property")}
    assert "unit_round_trip" in declared, "this spec skips the unit round trip on purpose"
