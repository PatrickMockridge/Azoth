"""Tests for the ``eos.ph_flash`` model.

The spec's cases pin the answers, and they are the weakest of the checks here: a
recorded number only says the model still does what it did when the number was
written. What the model actually claims is that it **inverts** the enthalpy, and that
is testable without knowing the enthalpy - compute ``H(T)`` at a temperature the test
chose, ask for it back, and require the temperature that comes out to be the one that
went in.

Two further properties matter and each has its own test below. The answer must not
depend on the *bracket* - bisection converges to the root, so a scan over 2000 points
and a scan over 40 must agree - and the single-phase branch must be reachable and must
report no ``beta``, because the flash's own value there is an extrapolation rather
than a vapour fraction.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import PhFlashResult
from azoth.eos import Component, IdealGasModel, Mixture, mixture
from azoth.eos.reference.ph_flash import enthalpy_at

MODEL_ID = "eos.ph_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]
TOLERANCE = SPEC["algorithm"]["tolerance"]

METHANE = Component(Q(190.56, "K"), Q(4_599_000.0, "Pa"), 0.0115)
BUTANE = Component(Q(425.12, "K"), Q(3_796_000.0, "Pa"), 0.2002)


def a_mixture() -> Mixture:
    """The mixture every test here uses - the one the spec's cases describe."""
    return mixture([METHANE, BUTANE], kij={(0, 1): 0.01289789})


def an_ideal_gas() -> IdealGasModel:
    """A deliberately trivial ideal-gas model: five zero coefficients.

    Nothing about this model depends on the coefficients being physical - a unit
    operation's arithmetic is the same whatever they are - so they are chosen to make
    a failure readable rather than to describe a substance.
    """
    return IdealGasModel(
        cp_a=(3.0, 5.0),
        cp_b=(0.0, 0.0),
        cp_c=(0.0, 0.0),
        cp_d=(0.0, 0.0),
        cp_e=(0.0, 0.0),
    )


def call(case: dict[str, Any]) -> PhFlashResult:
    """Run one case declared in the model spec."""
    inputs = case["inputs"]
    fluid = mixture(
        [
            Component(Q(tc, "K"), Q(pc, "Pa"), omega)
            for tc, pc, omega in zip(inputs["Tc"], inputs["Pc"], inputs["omega"], strict=True)
        ],
        kij={(0, 1): inputs["kij"][0][1]},
    )
    ideal_gas = IdealGasModel(
        cp_a=inputs["cp_a"],
        cp_b=inputs["cp_b"],
        cp_c=inputs["cp_c"],
        cp_d=inputs["cp_d"],
        cp_e=(0.0, 0.0),
    )
    import azoth

    return azoth.eos.ph_flash(
        fluid, ideal_gas, Q(inputs["P"], "Pa"), Q(inputs["H"], "J/mol"), inputs["z"]
    )


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

    assert result.warnings is not None


def test_every_case_ran() -> None:
    """A suite that silently collected nothing passes every other test here."""
    assert CASES, "the spec declares no cases"
    for case in CASES:
        call(case)


def test_the_answer_does_not_depend_on_the_bracket() -> None:
    """A narrow scan and a wide one converge to the same root.

    This is the check that is not self-referential. The round trip below proves the
    inversion is exact at whatever temperature was asked for; this proves the answer
    is a property of the *state* rather than of the search - two brackets that do not
    share an endpoint still bracket the same root, and bisection finds it in both.
    """
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    P, H = Q(20.0, "bar"), Q(-6723.003844102389, "J/mol")

    wide = azoth.eos.ph_flash(fluid, ideal_gas, P, H, [0.6, 0.4])

    # A bracket a hundred times narrower, centred away from the wide one's midpoint,
    # so an answer that was really the scan's resolution would move.
    narrow = _with_bracket(
        {"lower": 250.0, "upper": 350.0, "steps": 40},
        fluid,
        ideal_gas,
        P,
        H,
    )

    h.assert_close(
        narrow.T.to("K").magnitude,
        wide.T.to("K").magnitude,
        TOLERANCE,
        "narrow vs wide bracket",
    )


def _with_bracket(
    bracket: dict[str, Any],
    fluid: Mixture,
    ideal_gas: IdealGasModel,
    P: Any,
    H: Any,
) -> PhFlashResult:
    """Solve with a substituted bracket, through the reference implementation.

    The public entry point reads its bracket from the spec, which is the point of the
    spec - so a test that needs a *different* bracket calls the procedure directly
    rather than pretending the spec said something else.
    """
    from azoth.core.units import input_to_si
    from azoth.eos.reference.ph_flash import solve_temperature

    algorithm = dict(SPEC["algorithm"])
    algorithm["bracket"] = bracket
    solved = solve_temperature(
        fluid,
        ideal_gas,
        input_to_si(SPEC, "P", P),
        input_to_si(SPEC, "H", H),
        [0.6, 0.4],
        algorithm,
    )
    from azoth.core.units import from_si

    flash = solved["state"]["flash"]
    return PhFlashResult(
        T=from_si(solved["T"], "K"),
        beta=flash.beta,
        x=tuple(flash.x),
        y=tuple(flash.y),
        k=tuple(flash.k),
        phase=flash.phase,
        z_liquid=flash.z_liquid,
        z_vapour=flash.z_vapour,
        iterations=solved["iterations"],
        residual=solved["residual"],
        warnings=(),
    )


@pytest.mark.parametrize("temperature", [220.0, 300.0, 360.0, 450.0])
def test_the_round_trip_returns_the_temperature_that_went_in(temperature: float) -> None:
    """``H(T)`` inverted gives ``T`` back, across both phases.

    The property the model actually claims, and the only check here that does not
    depend on a number somebody recorded. Two of these temperatures are in the
    two-phase region and two are not, so a branch that was wrong would fail on one
    side of the bubble point rather than being flattered by the other.
    """
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    P = Q(20.0, "bar")

    enthalpy, _ = enthalpy_at(fluid, ideal_gas, temperature, 20.0e5, [0.6, 0.4])
    result = azoth.eos.ph_flash(fluid, ideal_gas, P, Q(enthalpy, "J/mol"), [0.6, 0.4])

    h.assert_close(
        result.T.to("K").magnitude, temperature, TOLERANCE, f"round trip at {temperature} K"
    )


def test_a_single_phase_feed_reports_no_vapour_fraction() -> None:
    """`beta` is absent above the dew point, and the answer is still the temperature.

    The branch this pins is the one that is easy to get wrong and impossible to see:
    the flash *does* return a split for a single-phase feed, and it is an
    extrapolation - values like 1.9 are ordinary. A model that used it would produce a
    wrong enthalpy shaped exactly like a right one.
    """
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    enthalpy, _ = enthalpy_at(fluid, ideal_gas, 450.0, 20.0e5, [0.6, 0.4])

    result = azoth.eos.ph_flash(fluid, ideal_gas, Q(20.0, "bar"), Q(enthalpy, "J/mol"), [0.6, 0.4])

    assert result.phase == "all_vapour"
    assert result.beta is None, "a single-phase feed has no vapour fraction to report"
    h.assert_close(result.T.to("K").magnitude, 450.0, TOLERANCE, "single-phase round trip")


def test_the_inversion_holds_across_the_whole_range() -> None:
    """Swept, not spot-checked, and that is the whole point of the test.

    The four temperatures in the round trip above are `[220, 300, 360, 450]`, and every
    one of them passed while this model was wrong for **105 of 111** enthalpy targets
    between 800 and 1900 J/mol. They sat either side of the damaged band rather than in
    it: at 20 bar the flash switched between a negative flash and no root at about
    386.7 K, the enthalpy of the feed jumped across that switch, and the bisection
    collapsed onto the jump and returned 386.6944 K for every target inside it - with
    residuals up to 1.275 against a declared tolerance of 1e-8, and no error.

    So what is asserted is the claim itself, at every point rather than at four: the
    temperature that comes back has the enthalpy that was asked for. A target this
    model genuinely cannot reach is allowed to raise; returning a wrong number is not.
    """
    import azoth
    from azoth.core.errors import SolverNotConvergedError

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    z = [0.6, 0.4]
    P = Q(20.0, "bar")

    unreachable = 0
    for step in range(111):
        target = 800.0 + 10.0 * step
        try:
            result = azoth.eos.ph_flash(fluid, ideal_gas, P, Q(target, "J/mol"), z)
        except SolverNotConvergedError:
            unreachable += 1
            continue
        achieved, _ = enthalpy_at(fluid, ideal_gas, result.T.to("K").magnitude, 20.0e5, z)
        residual = abs(achieved - target) / max(abs(target), 1.0)
        assert residual <= TOLERANCE, (
            f"H = {target}: the flash returned {result.T.to('K').magnitude} K, whose "
            f"enthalpy is {achieved}, a residual of {residual:.3e} against the declared "
            f"tolerance {TOLERANCE:.0e}"
        )

    assert unreachable < 111, "every target raised - the sweep proved nothing"


def test_an_enthalpy_no_temperature_produces_is_a_solver_failure() -> None:
    """Outside the bracket is a failure, not a number.

    The cost of a fixed scan instead of an adaptive expansion, and reported rather
    than hidden: an enthalpy below anything 100 K produces has no answer this model
    can give, and returning the bracket's end would be a plausible wrong number.
    """
    import azoth
    from azoth.core.errors import SolverNotConvergedError

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    # Far below the enthalpy at 100 K, which is the bottom of the scan.
    with pytest.raises(SolverNotConvergedError):
        azoth.eos.ph_flash(fluid, ideal_gas, Q(20.0, "bar"), Q(-1.0e7, "J/mol"), [0.6, 0.4])


def test_a_non_positive_pressure_is_refused() -> None:
    """A pressure is positive by definition; the range check says so."""
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    with pytest.raises(OutOfRangeError):
        azoth.eos.ph_flash(fluid, ideal_gas, Q(0.0, "Pa"), Q(0.0, "J/mol"), [0.6, 0.4])


def test_a_composition_that_is_not_a_composition_is_refused() -> None:
    """Checked rather than renormalised, as everywhere else in this library."""
    import azoth

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    with pytest.raises(InvalidInputError):
        azoth.eos.ph_flash(fluid, ideal_gas, Q(20.0, "bar"), Q(-6723.0, "J/mol"), [0.6, 0.6])


@pytest.mark.requires_rust
def test_the_two_implementations_agree_with_matching_iteration_counts() -> None:
    """Both implementations reach the same temperature, in the same number of steps.

    The iteration count is the load-bearing half. Two implementations that agree on a
    temperature reached by different paths have not been shown to agree on anything
    durable - and the count is also what pins the bracket, since a scan that stopped
    somewhere else would bisect from a different interval.
    """
    import azoth
    from azoth._dispatch import use_backend

    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    P = Q(20.0, "bar")

    for temperature in (240.0, 300.0, 360.0, 450.0):
        enthalpy, _ = enthalpy_at(fluid, ideal_gas, temperature, 20.0e5, [0.6, 0.4])
        H = Q(enthalpy, "J/mol")

        with use_backend("python"):
            py = azoth.eos.ph_flash(fluid, ideal_gas, P, H, [0.6, 0.4])
        with use_backend("rust"):
            rs = azoth.eos.ph_flash(fluid, ideal_gas, P, H, [0.6, 0.4])

        context = f"ph_flash at {temperature} K"
        h.assert_close(rs.T.to("K").magnitude, py.T.to("K").magnitude, TOLERANCE, context)
        assert rs.iterations == py.iterations, (
            f"{context}: python took {py.iterations} steps and Rust {rs.iterations}. "
            f"The bracket and the stopping rule have to be the same on both sides."
        )
        assert rs.phase == py.phase, f"{context}: phase disagrees"
