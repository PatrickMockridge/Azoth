"""Spec-driven tests for ``eos.prsv_kappa``."""

from __future__ import annotations

from decimal import Decimal, getcontext
from itertools import pairwise
from typing import Any

import pytest

import _helpers as h
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PrsvKappaResult
from azoth.eos import pr_alpha_ab, prsv_kappa

CALC_ID = "eos.prsv_kappa"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "consistency_with": "test_consistency_with",
    "monotonic": "test_monotonic",
}

#: A deliberately arbitrary, non-zero ``kappa1``.
#:
#: This library ships no fitted values for it - a table of them is the databank it
#: does not have - so the tests pick a round number and say so rather than borrowing
#: one that looks like data.
SOME_KAPPA1 = 0.05

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> PrsvKappaResult:
    return prsv_kappa(**h.kwargs_for(SPEC, case["inputs"]))


def f64_chain(omega: float, Tr: float, kappa1: float) -> float:
    """The arithmetic as written, so a spec value can be compared to it exactly."""
    acentric_only: float = (
        0.378893
        + 1.4897153 * omega
        - 0.17131848 * (omega * omega)
        + 0.0196554 * (omega * omega * omega)
    )
    temperature: float = (1.0 + Tr**0.5) * (0.7 - Tr)
    return acentric_only + kappa1 * temperature


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    h.assert_close(
        result.kappa,
        h.expected(case, "kappa"),
        case.get("tolerance", 1e-12),
        f"{case['id']} (kappa)",
    )
    h.assert_consistent(result, case["id"])

    def resolve(quantity: str, _case: dict[str, Any] = case) -> float | None:
        if quantity == "kappa":
            return result.kappa
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 3, f"expected several active cases, found {total}"


def test_the_expected_values_are_exactly_what_the_arithmetic_produces() -> None:
    """Bit-exact, not merely within tolerance.

    The tolerance hides a one-ulp disagreement, and this spec had one: its expected
    values were the correctly-rounded *exact decimals* while the arithmetic produces
    a double one ulp below, so the case passed against a number the implementation
    never returns. A reader checking the spec by hand would compute the exact value
    and conclude the code was wrong; a caller comparing against the spec's number
    would find they never match.

    `test_spec_case` cannot catch that - it is comparing within 1e-12, and one ulp
    is 1.7e-16. This can, and does.
    """
    for case in ACTIVE:
        inputs = case["inputs"]
        produced = f64_chain(float(inputs["omega"]), float(inputs["Tr"]), float(inputs["kappa1"]))
        stated = h.expected(case, "kappa")
        assert produced == stated, (
            f"{case['id']}: the spec states {stated!r} but the arithmetic produces "
            f"{produced!r}. The expected value should be what the code returns."
        )


def test_the_spec_says_how_far_its_value_is_from_the_exact_decimal() -> None:
    """The derivation's claim of "one ulp below exact" is the claim, so check it.

    The spec states that the four acentric terms sum to 0.6014406094290432 as exact
    decimals and 0.6014406094290431 in double precision. That is a specific,
    checkable assertion about *why* the worked example is not bit-exact, and it is
    the kind of sentence that rots silently if the equation is ever edited.
    """
    getcontext().prec = 50
    omega = Decimal("0.152")
    exact = (
        Decimal("0.378893")
        + Decimal("1.4897153") * omega
        - Decimal("0.17131848") * omega**2
        + Decimal("0.0196554") * omega**3
    )
    assert exact == Decimal("0.6014406094290432"), f"exact decimal sum is {exact}"

    in_floats = f64_chain(0.152, 0.7, 0.0)
    assert in_floats == 0.6014406094290431, f"the f64 sum is {in_floats!r}"

    # One ulp apart, which is the claim. `math.ulp` is 3.11+; the comparison is
    # written as an equality of the difference against one ulp so a future Python
    # without `ulp` still reads clearly.
    import math

    gap = abs(Decimal(repr(in_floats)) - exact)
    assert gap <= Decimal(math.ulp(in_floats)), (
        f"the spec claims a one-ulp gap; it is {gap} against an ulp of {math.ulp(in_floats)}"
    )


def test_consistency_with() -> None:
    """The PRSV coefficient feeds Peng-Robinson's alpha function unchanged.

    The whole argument for splitting the coefficients out of the equation of state
    rests on this: PRSV modifies the coefficient's temperature dependence and
    nothing else, so `pr_alpha_ab` served it without one line changing. If that were
    wrong, `pr_alpha_ab` would have needed a variant and the decomposition would
    have been a fork after all.
    """
    for omega, tr, kappa1 in ((0.152, 0.8, SOME_KAPPA1), (0.01142, 1.2, -0.03)):
        kappa = prsv_kappa(omega, tr, kappa1).kappa
        ab = pr_alpha_ab(kappa, tr, 0.25)

        # The alpha function restated here, so the test checks it rather than
        # trusting the calc that computes it.
        expected_alpha = (1.0 + kappa * (1.0 - tr**0.5)) ** 2
        h.assert_close(ab.alpha, expected_alpha, 1e-12, "alpha from the PRSV coefficient")
        # And the reduced parameters must be built from *that* alpha - a calc that
        # took the coefficient and ignored it would pass the line above.
        h.assert_close(
            ab.a_reduced, 0.4572355289213822 * expected_alpha * 0.25 / tr**2, 1e-12, "a_reduced"
        )


def test_monotonic() -> None:
    """``kappa`` rises with ``omega``, for any ``kappa1``.

    The derivative of the acentric-only polynomial has a negative discriminant and so
    is positive everywhere, and the ``kappa1`` term does not involve ``omega`` at all
    - so monotonicity holds whatever the fitted parameter is, which is what makes it
    safe to assert without one.
    """
    for kappa1 in (0.0, SOME_KAPPA1, -0.5):
        for low, high in pairwise((0.0, 0.1, 0.2, 0.35, 0.6, 1.0)):
            assert prsv_kappa(high, 0.8, kappa1).kappa > prsv_kappa(low, 0.8, kappa1).kappa, (
                f"kappa must rise with omega at kappa1={kappa1}: {low} -> {high}"
            )


def test_the_kappa1_term_vanishes_at_the_anchor() -> None:
    """At ``Tr = 0.7`` the coefficient does not depend on ``kappa1`` at all.

    ``(1 + sqrt(Tr))`` is never zero, so the whole term is zero exactly when
    ``0.7 - Tr`` is. The coefficient collapses to its acentric-only part at the
    temperature where the acentric factor is defined - the anchor PRSV is built
    around - and it is the sharpest cheap test of the temperature term.
    """
    reference = prsv_kappa(0.152, 0.7, 0.0).kappa
    for kappa1 in (0.0, SOME_KAPPA1, -3.2, 100.0, 1e6):
        assert prsv_kappa(0.152, 0.7, kappa1).kappa == reference, (
            f"at Tr = 0.7 the coefficient should not depend on kappa1, but "
            f"kappa1={kappa1} changed it"
        )


def test_a_zero_kappa1_is_temperature_independent() -> None:
    """With ``kappa1 = 0`` the coefficient does not move with ``Tr`` at all."""
    reference = prsv_kappa(0.152, 0.7, 0.0).kappa
    for tr in (0.3, 0.5, 0.7, 0.9, 1.0, 1.5):
        assert prsv_kappa(0.152, tr, 0.0).kappa == reference, f"moved at Tr = {tr}"


@pytest.mark.parametrize("tr", [0.0, -0.5])
def test_a_non_positive_reduced_temperature_is_an_error(tr: float) -> None:
    with pytest.raises(OutOfRangeError) as excinfo:
        prsv_kappa(0.152, tr, SOME_KAPPA1)
    assert excinfo.value.field() == "Tr"


def test_kappa1_is_unconstrained_in_sign() -> None:
    """Both signs occur in published fits, so neither may be refused.

    Asserting it keeps a future bound from being added on the assumption that the
    parameter behaves like a coefficient that has to be positive.
    """
    for kappa1 in (-1.0, -0.03, 0.0, 0.05, 12.0):
        result = prsv_kappa(0.152, 0.8, kappa1)
        assert result.is_clean, (
            f"kappa1 = {kappa1} produced {result.warnings}; the spec declares no bound "
            f"on the parameter because both signs occur in published fits"
        )


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
