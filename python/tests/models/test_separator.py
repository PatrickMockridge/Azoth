"""Tests for the ``process.separator`` unit operation.

The spec's cases pin the answers, and they are the weakest of the checks here: a
recorded number only says the model still does what it did when the number was written.

What the model actually claims is two things, and both are testable without a recorded
number. **Mass is conserved**: the two outlet flows add back to the feed, at every state
rather than at one, which is the analogue of a worked example that a flowsheet has and a
single calculation does not. And **the split is the flash's split**: with no duty this
model is ``eos.pt_flash`` plus arithmetic, so its vapour fraction and both compositions
must be that model's, computed independently in the same test. A separator whose split
agreed with nothing would still pass a recorded-number test.

The rest pin the branches that are easy to get wrong and impossible to see - a
single-phase feed leaving one outlet empty, a duty moving the temperature, and the two
range checks whose whole purpose is to stop a division by zero and a pressure *rise*
being modelled as a separator.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.errors import OutOfRangeError
from azoth.core.result import Phase, SeparatorResult
from azoth.eos import IdealGasModel, Mixture, component, mixture

MODEL_ID = "process.separator"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

# The substances the cases and the identities below use, resolved through the
# databank rather than typed here. A `Component` written out longhand is a second
# copy of NeqSim's table, and this one had drifted: its methane was 0.01142 and
# 4 599 200 Pa, where COMP.csv says 0.0115 and 4 599 000.
METHANE = component("methane")
BUTANE = component("n-butane")


def a_mixture() -> Mixture:
    """The mixture every test here uses - the one the spec's cases describe."""
    return mixture([METHANE, BUTANE], kij={(0, 1): 0.01289789})


def an_ideal_gas() -> IdealGasModel:
    """A deliberately trivial datum: four zero coefficients and zero reference values.

    Nothing about this model depends on the datum being physical - the duty branch
    measures the feed's enthalpy from whatever it is and adds to it - so the
    coefficients are chosen to make a failure readable rather than to describe a
    substance.
    """
    return IdealGasModel(
        cp_a=(3.0, 5.0),
        cp_b=(0.0, 0.0),
        cp_c=(0.0, 0.0),
        cp_d=(0.0, 0.0),
        cp_e=(0.0, 0.0),
    )


def run(
    *,
    T: float = 300.0,
    P: float = 20.0e5,
    n: float = 10.0,
    z: list[float] | None = None,
    pressure_drop: float = 0.0,
    heat_duty: float = 0.0,
) -> SeparatorResult:
    """Call the model at a state, with the spec's own mixture and datum."""
    import azoth

    return azoth.process.separator(
        a_mixture(),
        an_ideal_gas(),
        T=Q(T, "K"),
        P=Q(P, "Pa"),
        n=Q(n, "mol/s"),
        z=[0.6, 0.4] if z is None else z,
        pressure_drop=Q(pressure_drop, "Pa"),
        heat_duty=Q(heat_duty, "W"),
    )


def call(case: dict[str, Any]) -> SeparatorResult:
    """Run one case declared in the model spec, through the public API."""
    import azoth

    inputs = case["inputs"]
    return azoth.process.separator(
        a_mixture(),
        an_ideal_gas(),
        T=Q(inputs["T"], "K"),
        P=Q(inputs["P"], "Pa"),
        n=Q(inputs["n"], "mol/s"),
        z=inputs["z"],
        pressure_drop=Q(inputs["pressure_drop"], "Pa"),
        heat_duty=Q(inputs["heat_duty"], "W"),
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    """Run one case and compare against the values the spec records."""
    result = call(case)

    for name, expected_value in case["expected"].items():
        got = getattr(result, name)
        if name in ("T", "P"):
            got = got.to("K" if name == "T" else "Pa").magnitude
        h.assert_close(got, expected_value, case["tolerance"], f"{case['id']} ({name})")


def test_every_case_ran() -> None:
    """A suite that silently collected nothing passes every other test here."""
    assert CASES, "the spec declares no cases"
    for case in CASES:
        call(case)


@pytest.mark.parametrize("pressure_drop", [0.0, 5.0e5, 1.5e6])
@pytest.mark.parametrize("T", [300.0, 450.0])
def test_the_two_outlet_flows_add_back_to_the_feed(T: float, pressure_drop: float) -> None:
    """Mass is conserved at every state, not at one recorded one.

    The check a flowsheet carries that a single calculation cannot, and the reason this
    model is worth more than its worked examples: a split that lost or invented moles
    would still match a recorded number if the number were recorded from the same wrong
    code.
    """
    result = run(T=T, P=20.0e5, n=10.0, pressure_drop=pressure_drop)
    h.assert_close(
        result.gas_flow + result.liquid_flow,
        10.0,
        1.0e-9,
        f"mole balance at {T} K with a {pressure_drop} Pa drop",
    )


def test_the_split_is_the_flash_split() -> None:
    """The separator's split equals ``eos.pt_flash``'s answer, computed separately.

    The cross-model check, and the one here that is not self-referential. With no duty
    and no pressure drop this model *is* that model plus arithmetic, so the vapour
    fraction and both compositions have to come out identical - and this computes them
    from the other model rather than recording them, so the two cannot drift together.
    """
    import azoth

    flash = azoth.eos.pt_flash(a_mixture(), T=Q(300.0, "K"), P=Q(20.0e5, "Pa"), z=[0.6, 0.4])
    result = run()

    assert result.phase == flash.phase
    assert result.beta is not None
    assert flash.beta is not None
    h.assert_close(result.beta, flash.beta, 1.0e-12, "vapour fraction against the flash")

    # The gas outlet is the vapour phase and the liquid outlet the liquid one, so the
    # two vectors are the flash's `y` and `x` - not swapped, and not renormalised.
    for got, expected in zip(result.gas_z, flash.y, strict=True):
        h.assert_close(got, expected, 1.0e-12, "gas composition")
    for got, expected in zip(result.liquid_z, flash.x, strict=True):
        h.assert_close(got, expected, 1.0e-12, "liquid composition")


def test_a_single_phase_feed_leaves_one_outlet_empty() -> None:
    """Above the dew point the whole feed leaves in one stream, with no ``beta``.

    The branch that is easy to get wrong and impossible to see: the flash *does* return
    a vapour fraction for a single-phase feed, and it is an extrapolation - values like
    1.9 are ordinary. A model that split on it would move moles into an outlet that
    should be empty.

    The assertion is on the totals rather than on which outlet, because `all_vapour` and
    `trivial` are both reachable at this state and are routed the same way.
    """
    result = run(T=450.0)

    assert result.phase in (Phase.ALL_VAPOUR, Phase.TRIVIAL)
    assert result.beta is None, "a single-phase feed has no vapour fraction to report"
    assert result.liquid_flow == 0.0
    h.assert_close(result.gas_flow, 10.0, 1.0e-12, "the whole feed in one outlet")


def test_a_duty_moves_the_temperature_and_a_drop_does_not() -> None:
    """The two branches of the procedure, told apart by what they do to ``T``.

    A separator with no duty is isothermal at its feed temperature - that is what the
    ``TPflash`` branch means and it is why the spec's second case can assert 300 K
    exactly. A duty makes the flash isenthalpic instead, and the temperature becomes an
    answer rather than an input. Both feed the same mixture, so the only difference is
    the branch.
    """
    heated = run(T=300.0, heat_duty=5.0e4)
    isothermal = run(T=300.0, pressure_drop=5.0e5)

    h.assert_close(isothermal.T.to("K").magnitude, 300.0, 1.0e-12, "a drop leaves T alone")
    assert heated.T.to("K").magnitude > 300.0, (
        "a positive duty must raise the temperature; the isenthalpic branch is not being "
        "taken, or the duty is not reaching it"
    )
    h.assert_close(
        heated.gas_flow + heated.liquid_flow, 10.0, 1.0e-9, "mole balance on the duty branch"
    )


def test_a_negative_pressure_drop_is_refused() -> None:
    """A separator drops pressure. A rise is a compressor, and is not this model."""
    with pytest.raises(OutOfRangeError):
        run(pressure_drop=-1.0e5)


def test_a_zero_flow_is_refused() -> None:
    """Refused loudly rather than reaching the duty branch's division by zero."""
    with pytest.raises(OutOfRangeError):
        run(n=0.0)


def test_a_non_positive_pressure_is_refused() -> None:
    """A pressure is positive by definition; the range check says so."""
    with pytest.raises(OutOfRangeError):
        run(P=0.0)


@pytest.mark.requires_rust
@pytest.mark.parametrize("T", [300.0, 450.0])
def test_the_two_implementations_agree_with_matching_iteration_counts(T: float) -> None:
    """Both implementations reach the same split, in the same number of steps.

    The iteration count is the load-bearing half. Two implementations that agree on a
    split reached by different paths have not been shown to agree on anything durable -
    and here the count is the *flash's*, so it also pins that both sides delegate to the
    same model rather than one of them reimplementing the iteration.
    """
    from azoth._dispatch import use_backend

    with use_backend("python"):
        py = run(T=T)
    with use_backend("rust"):
        rs = run(T=T)

    context = f"separator at {T} K"
    h.assert_close(rs.T.to("K").magnitude, py.T.to("K").magnitude, 1.0e-12, f"{context}: T")
    h.assert_close(rs.gas_flow, py.gas_flow, 1.0e-12, f"{context}: gas flow")
    h.assert_close(rs.liquid_flow, py.liquid_flow, 1.0e-12, f"{context}: liquid flow")
    assert rs.phase == py.phase, f"{context}: phase disagrees"
    assert rs.iterations == py.iterations, (
        f"{context}: python took {py.iterations} steps and Rust {rs.iterations}. "
        f"Both sides have to delegate to the same flash."
    )
    for a, b in zip(rs.gas_z, py.gas_z, strict=True):
        h.assert_close(a, b, 1.0e-12, f"{context}: gas composition")
