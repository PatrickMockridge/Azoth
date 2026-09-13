"""Spec-driven tests for ``eos.pr_departure``."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth.core.errors import OutOfRangeError
from azoth.core.result import PrDepartureResult
from azoth.eos import pr_departure

CALC_ID = "eos.pr_departure"

SPEC = h.spec(CALC_ID)
_PROPERTY_IMPLEMENTATIONS: dict[str, str] = {
    "consistency_with": "test_consistency_with",
}

#: The three outputs the spec declares, in order.
FIELDS: tuple[str, ...] = ("ln_phi", "h_dep_rt", "s_dep_r")

ALL_CASES = h.all_tests(SPEC)
ACTIVE = [c for c in ALL_CASES if c["status"] == "active" and c["type"] != "property"]
PROPERTIES = [c for c in ALL_CASES if c["status"] == "active" and c["type"] == "property"]
SKIPPED = [c for c in ALL_CASES if c["status"] != "active"]


def call(case: dict[str, Any]) -> PrDepartureResult:
    return pr_departure(**h.kwargs_for(SPEC, case["inputs"]))


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
        if quantity == "z_minus_b_reduced":
            return float(_case["inputs"]["z"]) - float(_case["inputs"]["b_reduced"])
        return h.input_(_case, quantity) if quantity in _case["inputs"] else None

    h.assert_warnings_agree_with_spec(SPEC, result.warnings, resolve, case["id"])


def test_every_active_case_ran() -> None:
    """Guard against a spec edit that silently removes all runnable cases."""
    total = len(ACTIVE) + len(PROPERTIES)
    assert total >= 4, f"expected several active cases, found {total}"


def test_consistency_with() -> None:
    """The two `consistency_with` properties the spec declares.

    The first is the Gibbs identity: ``h_dep_rt - s_dep_r`` equals ``ln_phi``. It is
    exact in real arithmetic - the two ``psi`` terms cancel, as the spec's
    ``verification`` shows in three lines - so any disagreement is rounding, and it
    is the only check that catches a sign error in either departure function.

    The second is the low-pressure limit: as ``A`` and ``B`` vanish at ``z = 1``, all
    three outputs vanish linearly. That one uses no reference values at all.
    """
    # The identity, swept rather than spot-checked.
    for a in (0.05, 0.20206500174625697, 0.4572355289213822):
        for b in (0.005, 0.02431127309496514, 0.07779607390388846):
            for z in (0.05, 0.3, 0.79, 0.95):
                if z <= b:
                    continue
                result = pr_departure(a, b, z, 0.60282728832, 0.8)
                identity = result.h_dep_rt - result.s_dep_r
                assert abs(identity - result.ln_phi) < 1e-12, (
                    f"A={a}, B={b}, z={z}: h - s = {identity} against ln_phi = {result.ln_phi}"
                )

    # The low-pressure limit, and its linearity in B.
    previous: list[float] | None = None
    for exponent in (3, 5, 7, 9):
        b = 10.0**-exponent
        result = pr_departure(8.0 * b, b, 1.0, 0.6, 0.8)
        magnitudes = [abs(getattr(result, field)) for field in FIELDS]
        if previous is not None:
            for field, before, after in zip(FIELDS, previous, magnitudes, strict=True):
                shrink = before / after
                assert 50.0 <= shrink <= 200.0, (
                    f"{field} should shrink roughly 100-fold per decade of B, shrank {shrink}"
                )
        previous = magnitudes

    # And negligible at the smallest B. Threshold from the measurement - 1.2e-8 at
    # B = 1e-9 - rather than picked to look small.
    b = 1e-9
    result = pr_departure(8.0 * b, b, 1.0, 0.6, 0.8)
    for field in FIELDS:
        assert abs(getattr(result, field)) < 1e-7, f"{field} should vanish at B = 1e-9"


def test_kappa_and_tr_do_not_reach_the_fugacity_coefficient() -> None:
    """The quiet failure the spec warns about, asserted so it stays visible.

    ``kappa`` and ``Tr`` enter only through ``psi``, which the ``ln_phi`` expression
    does not use. So a caller who mixes states gets an *unchanged* fugacity
    coefficient and shifted departures, rather than something obviously wrong. This
    test pins that behaviour: if a future revision made ``ln_phi`` depend on them,
    the spec's `assumptions` entry would be wrong and this would say so.
    """
    base = pr_departure(0.20206500174625697, 0.02431127309496514, 0.79, 0.60282728832, 0.8)
    for kappa, tr in ((0.2, 0.8), (0.60282728832, 1.4), (0.9, 0.5)):
        other = pr_departure(0.20206500174625697, 0.02431127309496514, 0.79, kappa, tr)
        assert other.ln_phi == base.ln_phi, (
            f"ln_phi should not depend on kappa or Tr, but ({kappa}, {tr}) moved it"
        )
        assert abs(other.h_dep_rt - base.h_dep_rt) > 1e-9, (
            f"h_dep_rt should move with ({kappa}, {tr})"
        )


def test_a_non_positive_b_reduced_is_an_error() -> None:
    for b in (0.0, -0.01):
        with pytest.raises(OutOfRangeError) as excinfo:
            pr_departure(0.2, b, 0.79, 0.6, 0.8)
        assert excinfo.value.field() == "b_reduced"


@pytest.mark.parametrize("z", [0.02431127309496514, 0.01, 0.0, -1.0])
def test_a_z_at_or_below_b_reduced_is_an_error_naming_the_difference(z: float) -> None:
    """The bound is on ``z - B``, not on ``z``.

    That is the quantity ``ln(z - B)`` needs to be positive, and nothing about ``z``
    alone says where the boundary is. ``eos.pr_z_factor`` already discards roots at or
    below ``B``, so this is for the caller who supplied ``z`` some other way - and it
    names the failure rather than letting a domain error do it.
    """
    with pytest.raises(OutOfRangeError) as excinfo:
        pr_departure(0.2, 0.02431127309496514, z, 0.6, 0.8)
    assert excinfo.value.field() == "z_minus_b_reduced"


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
