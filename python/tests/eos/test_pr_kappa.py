"""Spec-driven tests for ``eos.pr_kappa``."""

from __future__ import annotations

from itertools import pairwise
from typing import Any

import pytest

import _helpers as h
from azoth.core.result import PrKappaResult
from azoth.eos import pr_kappa

CALC_ID = "eos.pr_kappa"

SPEC = h.spec(CALC_ID)
#: Property test name -> the test function implementing it. Name-based rather than
#: a direct reference because the functions are defined further down.
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "monotonic": "test_monotonic",
}

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> PrKappaResult:
    # `h.kwargs_for` reads each input's declared unit from the spec, so this needs
    # no unit handling of its own and a spec change needs no change here. For a
    # `dimensionless` input it hands over a plain float, which is what this calc
    # takes - see the note on `test_the_acentric_factor_is_a_plain_float`.
    return pr_kappa(**h.kwargs_for(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", ACTIVE, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case declared in the spec."""
    result = call(case)
    expected = h.expected(case, "kappa")
    h.assert_close(result.kappa, expected, case.get("tolerance", 1e-12), f"{case['id']} (kappa)")
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


def test_the_acentric_factor_is_a_plain_float() -> None:
    """A dimensionless input takes a plain float, not a `pint` quantity.

    This is the project's convention for genuinely dimensionless quantities - the
    same rule that makes ``Re``, ``f`` and ``epsilon/D`` plain floats in the
    hydraulics calcs - rather than a decision made here. Asserting it explicitly
    because the opposite is tempting for a namespace whose *every* quantity is
    dimensionless: wrapping `omega` in a quantity would be ceremony with no unit
    behind it, and would then have to be unwrapped in the batch path, which takes
    plain numbers by construction.
    """
    import inspect

    parameter = inspect.signature(pr_kappa).parameters["omega"]
    assert parameter.annotation in (float, "float"), (
        f"`omega` should be annotated `float`, got {parameter.annotation!r}"
    )


def test_monotonic() -> None:
    """``kappa`` rises with ``omega`` across the range real fluids occupy.

    The derivative of the polynomial, ``1.54226 - 0.53984*omega``, is positive for
    every ``omega`` below 2.856883521043272, so the coefficient is strictly
    increasing in the acentric factor over the whole range a substance occupies.

    Sampled rather than checked at one adjacent pair, because a coefficient wrong in
    its curvature can be locally increasing and globally wrong - which is exactly
    the failure a mis-signed ``omega**2`` term produces.
    """
    samples = [-0.4 + 1.5 * i / 40.0 for i in range(41)]
    for low, high in pairwise(samples):
        assert pr_kappa(high).kappa > pr_kappa(low).kappa, (
            f"kappa must increase with omega: kappa({high}) = {pr_kappa(high).kappa} "
            f"is not greater than kappa({low}) = {pr_kappa(low).kappa}"
        )


def test_the_worked_example_is_exact_in_binary() -> None:
    """The spec's derivation claims the result is the double nearest the exact decimal.

    That claim is what justifies a 1e-12 tolerance on arithmetic whose inputs are
    decimal, so it is asserted rather than left as prose. If this fails, the spec's
    derivation has become false and the tolerance is no longer justified by it -
    which is a documentation failure, not a wrong number.
    """
    assert pr_kappa(0.152).kappa == 0.60282728832


def test_a_negative_kappa_is_returned_with_a_warning_rather_than_refused() -> None:
    """Helium's acentric factor is below the polynomial's root, so kappa is negative.

    A refusal would be wrong: "what does this correlation give for helium" is a real
    question, and the answer is a real number that happens to be out of range. The
    warning is what tells the caller not to use it as an attraction parameter.
    """
    result = pr_kappa(-0.385)
    h.assert_close(result.kappa, -0.259138992, 1e-12, "helium")
    assert not result.is_clean, "a negative kappa must be visible without checking the sign"


def test_a_positive_kappa_at_the_same_magnitude_is_silent() -> None:
    """The bound discriminates rather than merely firing near zero.

    Hydrogen sits just the other side of the root at -0.23338349942403008: its
    coefficient is small but positive, and the result must be clean.
    """
    result = pr_kappa(-0.216)
    h.assert_close(result.kappa, 0.02891845248, 1e-12, "hydrogen")
    assert result.is_clean, f"unexpected warnings: {result.warnings}"


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


def test_the_skipped_case_says_why_it_is_skipped() -> None:
    """A skipped test is still reported, and has to carry its reason."""
    h.assert_skips_are_explained(SPEC)
    declared = {str(c["property"]) for c in SKIPPED if c.get("property")}
    assert "unit_round_trip" in declared, (
        "this spec skips the unit round trip on purpose; if that changed, remove this test"
    )
