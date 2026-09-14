"""Tests for the ``eos.ps_flash`` model.

The same shape as ``test_ph_flash.py``, and deliberately so - the two models differ in
which property they invert and agree on everything else, so a reader who has read one has
read both. What is particular to it is ``test_the_two_phase_entropy_carries_the_entropy_of_mixing``:
the two-phase entropy is a weighted sum at each phase's **own** composition, which
carries the entropy of mixing, and a model that summed at the feed composition would
still pass every round-trip test here.

That check is the reason this file is not merely ``test_ph_flash`` with the letters
changed.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.errors import InvalidInputError, OutOfRangeError, SolverNotConvergedError
from azoth.core.result import PsFlashResult
from azoth.eos import IdealGasModel, Mixture, component, mixture, ps_flash
from azoth.eos.reference.ps_flash import entropy_at

MODEL_ID = "eos.ps_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]
TOLERANCE = SPEC["algorithm"]["tolerance"]

# The substances the cases and the identities below use, resolved through the
# databank rather than typed here. A `Component` written out longhand is a second
# copy of NeqSim's table, and this one had drifted: its methane was 0.01142 and
# 4 599 200 Pa, where COMP.csv says 0.0115 and 4 599 000.
METHANE = component("methane")
BUTANE = component("n-butane")

#: The mixture the spec's cases describe, at the pressure they use.
P_BAR = 20.0


def a_mixture() -> Mixture:
    return mixture([METHANE, BUTANE], kij={(0, 1): 0.01289789})


def an_ideal_gas() -> IdealGasModel:
    """A deliberately trivial datum, as in the enthalpy model's tests."""
    return IdealGasModel(
        cp_a=(3.0, 5.0),
        cp_b=(0.0, 0.0),
        cp_c=(0.0, 0.0),
        cp_d=(0.0, 0.0),
        cp_e=(0.0, 0.0),
    )


def call(case: dict[str, Any]) -> PsFlashResult:
    """Run one case declared in the model spec.

    Through :func:`_helpers.model_kwargs`, which is the one place a case's
    declared inputs become arguments: it resolves `components` against the
    databank and hands over the mixture and the ideal-gas model the function
    takes. A hand-built mixture here would be a second fluid, described by the
    case file rather than by NeqSim's tables.
    """
    return ps_flash(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case and compare against the values the spec records."""
    result = call(case)
    tolerance = case["tolerance"]

    for name, expected_value in case["expected"].items():
        got = getattr(result, name)
        if name == "T":
            got = got.to("K").magnitude
        h.assert_close(got, expected_value, tolerance, f"{case['id']} ({name})")


def test_every_case_ran() -> None:
    """A suite that silently collected nothing passes every other test here."""
    assert CASES, "the spec declares no cases"
    for case in CASES:
        call(case)


@pytest.mark.parametrize("temperature", [220.0, 300.0, 360.0, 450.0])
def test_the_round_trip_returns_the_temperature_that_went_in(temperature: float) -> None:
    """``S(T)`` inverted gives ``T`` back, across both phases."""
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    entropy, _ = entropy_at(fluid, ideal_gas, temperature, P_BAR * 1.0e5, [0.6, 0.4])
    result = azoth.eos.ps_flash(
        fluid, ideal_gas, Q(P_BAR, "bar"), Q(entropy, "J/(mol*K)"), [0.6, 0.4]
    )
    h.assert_close(
        result.T.to("K").magnitude, temperature, TOLERANCE, f"round trip at {temperature} K"
    )


def test_the_two_phase_entropy_carries_the_entropy_of_mixing() -> None:
    """The check that distinguishes this model from a naive weighted sum.

    Below the dew point the entropy of the feed is **not** the phase-weighted sum of the
    two phases' entropies evaluated at the *feed* composition. Each phase's entropy is
    computed at its own composition, and doing otherwise omits the entropy of mixing -
    which for this mixture changes the answer by far more than the solver's tolerance.

    The test computes both and requires them to differ. It is deliberately an assertion
    of *inequality*: if someone "fixed" the model to sum at the feed, every round-trip
    test in this file would still pass, and this is the one that would not.
    """
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    P = Q(P_BAR, "bar")
    P_si = P_BAR * 1.0e5

    flash = azoth.eos.pt_flash(fluid, Q(300.0, "K"), P, [0.6, 0.4])
    assert flash.phase == "two_phase", "this test needs a split to mean anything"

    correct, _ = entropy_at(fluid, ideal_gas, 300.0, P_si, [0.6, 0.4])

    # The wrong assembly: both phases' entropies computed at the feed composition.
    naive_l = (
        azoth.eos.molar_enthalpy_entropy(
            fluid, ideal_gas, Q(300.0, "K"), P, [0.6, 0.4], flash.z_liquid
        )
        .s.to("J/(mol*K)")
        .magnitude
    )
    naive_v = (
        azoth.eos.molar_enthalpy_entropy(
            fluid, ideal_gas, Q(300.0, "K"), P, [0.6, 0.4], flash.z_vapour
        )
        .s.to("J/(mol*K)")
        .magnitude
    )
    beta = flash.beta
    assert beta is not None, "a two-phase flash reports a vapour fraction"
    naive = (1.0 - beta) * naive_l + beta * naive_v

    assert abs(correct - naive) > 1.0e-6, (
        "the phase-weighted sum at each phase's own composition agrees with the sum at "
        "the feed composition, which means the entropy of mixing is missing"
    )


def test_a_single_phase_feed_reports_no_vapour_fraction() -> None:
    """Above the dew point there is no vapour fraction, and no extrapolation is used."""
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    entropy, _ = entropy_at(fluid, ideal_gas, 450.0, P_BAR * 1.0e5, [0.6, 0.4])
    result = azoth.eos.ps_flash(
        fluid, ideal_gas, Q(P_BAR, "bar"), Q(entropy, "J/(mol*K)"), [0.6, 0.4]
    )

    assert result.phase == "all_vapour"
    assert result.beta is None


def test_an_entropy_no_temperature_produces_is_a_solver_failure() -> None:
    """Outside the bracket is a failure, not a number."""
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    with pytest.raises(SolverNotConvergedError):
        azoth.eos.ps_flash(fluid, ideal_gas, Q(P_BAR, "bar"), Q(-500.0, "J/(mol*K)"), [0.6, 0.4])


def test_a_non_positive_pressure_is_refused() -> None:
    """A pressure is positive by definition."""
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    with pytest.raises(OutOfRangeError):
        azoth.eos.ps_flash(fluid, ideal_gas, Q(0.0, "Pa"), Q(0.0, "J/(mol*K)"), [0.6, 0.4])


def test_a_composition_that_is_not_a_composition_is_refused() -> None:
    """Checked rather than renormalised, as everywhere else in this library."""
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    with pytest.raises(InvalidInputError):
        azoth.eos.ps_flash(fluid, ideal_gas, Q(P_BAR, "bar"), Q(-39.0, "J/(mol*K)"), [0.6, 0.6])


@pytest.mark.requires_rust
def test_the_two_implementations_agree_with_matching_iteration_counts() -> None:
    """Both implementations reach the same temperature, in the same number of steps."""
    import azoth
    from azoth._dispatch import use_backend

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    P = Q(P_BAR, "bar")

    for temperature in (240.0, 300.0, 360.0, 450.0):
        entropy, _ = entropy_at(fluid, ideal_gas, temperature, P_BAR * 1.0e5, [0.6, 0.4])
        S = Q(entropy, "J/(mol*K)")

        with use_backend("python"):
            py = azoth.eos.ps_flash(fluid, ideal_gas, P, S, [0.6, 0.4])
        with use_backend("rust"):
            rs = azoth.eos.ps_flash(fluid, ideal_gas, P, S, [0.6, 0.4])

        context = f"ps_flash at {temperature} K"
        h.assert_close(rs.T.to("K").magnitude, py.T.to("K").magnitude, TOLERANCE, context)
        assert rs.iterations == py.iterations, (
            f"{context}: python took {py.iterations} steps and Rust {rs.iterations}. "
            f"The bracket and the stopping rule have to be the same on both sides."
        )
        assert rs.phase == py.phase, f"{context}: phase disagrees"
