"""Tests for the three machines: ``process.compressor``, ``process.pump`` and ``process.expander``.

**One file for three models, and the reason is that they are one procedure.** NeqSim has
three classes and the same eleven lines in each; here the procedure lives in
``crates/azoth-process/src/isentropic.rs`` and ``azoth/process/reference/_isentropic.py``,
and the only difference between the three is which way the efficiency scales the ideal
enthalpy change:

    compressor, pump   h_out = h_in + (h_ideal - h_in) / eta     Compressor.java:1738
    expander           h_out = h_in + (h_ideal - h_in) * eta     Expander.java:653

What is worth testing, then, is not each model's arithmetic in isolation - there is only
one arithmetic - but **the difference between them**, and that is what most of this file
does. Three separate files would test the same procedure three times and would each miss
the contrast.

The spec cases still pin each model's own recorded numbers, parametrised below.

# What is checked, and what is not

Two identities, computed in the test from `eos.ps_flash` and
`eos.molar_enthalpy_entropy` independently of any of the three models:

* **the isentropic state is isentropic**: ``S(T_is, P_out) == S(T_in, P_in)``;
* **the efficiency is honoured**: the real outlet enthalpy is the inlet's plus the ideal
  change scaled by the efficiency, in the direction that model claims.

The second is what catches the division and the multiplication being the wrong way round,
and the tests below run a compressor and an expander on the *same* input to make that
concrete: they are the two directions, and a sign error swaps them.
"""

from __future__ import annotations

from typing import Any

import pytest
from _process import Q, a_mixture, an_ideal_gas, assert_case, call, h

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.eos.reference.ph_flash import enthalpy_at
from azoth.eos.reference.ps_flash import entropy_at, ps_flash

MODELS = {name: _models_gen.model(f"process.{name}") for name in ("compressor", "pump", "expander")}

INLET = (300.0, 20.0e5)
COMPOSITION = [0.6, 0.4]
FLOW = 10.0

#: Which way each machine moves the pressure, and which way the efficiency therefore
#: scales. The whole of the difference between the three.
DIRECTIONS = {"compressor": "consuming", "pump": "consuming", "expander": "producing"}


def run(
    name: str,
    *,
    outlet_pressure: float,
    efficiency: float = 0.8,
    T: float = INLET[0],
    P: float = INLET[1],
    n: float = FLOW,
) -> Any:
    """Call one machine with the spec's own mixture and datum."""
    import azoth

    return getattr(azoth.process, name)(
        a_mixture(),
        an_ideal_gas(),
        T=Q(T, "K"),
        P=Q(P, "Pa"),
        n=Q(n, "mol/s"),
        z=COMPOSITION,
        outlet_pressure=Q(outlet_pressure, "Pa"),
        efficiency=efficiency,
    )


#: Every case of every one of the three, tagged with the model whose spec declares it, so
#: that one parametrised test covers what three per-model files would each do once.
CASES = [
    pytest.param(name, case, id=f"{name}-{case['id']}")
    for name, model in sorted(MODELS.items())
    for case in model["cases"]
]


@pytest.mark.parametrize(("name", "case"), CASES)
def test_spec_case(name: str, case: dict[str, Any]) -> None:
    """Run one case against the model whose spec declares it."""
    assert_case(MODELS[name], case)


@pytest.mark.parametrize("name", sorted(MODELS))
def test_every_model_declares_cases(name: str) -> None:
    assert MODELS[name]["cases"], f"{name} declares no cases"
    for case in MODELS[name]["cases"]:
        call(MODELS[name], case)


# ---------------------------------------------------------------------------------
# The two identities, per direction
# ---------------------------------------------------------------------------------


@pytest.mark.parametrize(
    ("name", "outlet_pressure"),
    [("compressor", 40.0e5), ("pump", 30.0e5), ("expander", 5.0e5)],
)
def test_the_ideal_state_is_isentropic(name: str, outlet_pressure: float) -> None:
    """``S(T_is, P_out) == S(T_in, P_in)``, computed from the equation of state.

    The ideal outlet that ``isentropic_temperature`` reports is where the fluid would get
    to at constant entropy, and this asks the entropy model directly rather than trusting
    the flash that produced it.
    """
    mixture, ideal_gas = a_mixture(), an_ideal_gas()
    result = run(name, outlet_pressure=outlet_pressure)

    s_in, _ = entropy_at(mixture, ideal_gas, *INLET, COMPOSITION)
    s_ideal, _ = entropy_at(
        mixture,
        ideal_gas,
        result.isentropic_temperature.to("K").magnitude,
        outlet_pressure,
        COMPOSITION,
    )
    h.assert_close(s_ideal, s_in, 1.0e-8, f"{name}: entropy at the ideal outlet")


@pytest.mark.parametrize(
    ("name", "outlet_pressure"),
    [("compressor", 40.0e5), ("pump", 30.0e5), ("expander", 5.0e5)],
)
@pytest.mark.parametrize("efficiency", [0.6, 0.8, 1.0])
def test_the_efficiency_is_honoured(name: str, outlet_pressure: float, efficiency: float) -> None:
    """The enthalpy change is the ideal one scaled the way that machine scales it.

    The check that catches the division and the multiplication being the wrong way round.
    It computes the ideal change from the model's own reported `isentropic_temperature`
    and applies the direction's scaling itself, so a model that did the arithmetic the
    other way fails here and nowhere else.
    """
    mixture, ideal_gas = a_mixture(), an_ideal_gas()
    result = run(name, outlet_pressure=outlet_pressure, efficiency=efficiency)

    h_in, _ = enthalpy_at(mixture, ideal_gas, *INLET, COMPOSITION)
    h_ideal, _ = enthalpy_at(
        mixture,
        ideal_gas,
        result.isentropic_temperature.to("K").magnitude,
        outlet_pressure,
        COMPOSITION,
    )
    h_real, _ = enthalpy_at(
        mixture,
        ideal_gas,
        result.T.to("K").magnitude,
        outlet_pressure,
        COMPOSITION,
    )

    if DIRECTIONS[name] == "consuming":
        expected = h_in + (h_ideal - h_in) / efficiency
    else:
        expected = h_in + (h_ideal - h_in) * efficiency

    h.assert_close(h_real, expected, 1.0e-6, f"{name} at eta={efficiency}: outlet enthalpy")


@pytest.mark.parametrize(
    ("name", "outlet_pressure"),
    [("compressor", 40.0e5), ("pump", 30.0e5), ("expander", 5.0e5)],
)
def test_the_power_is_the_enthalpy_change_times_the_flow(name: str, outlet_pressure: float) -> None:
    """``power == n * (H_out - H_in)``, with both enthalpies computed here.

    The sign convention is the same for all three - positive into the fluid - so this is
    also the check that an expander's power comes back negative.
    """
    mixture, ideal_gas = a_mixture(), an_ideal_gas()
    result = run(name, outlet_pressure=outlet_pressure)

    h_in, _ = enthalpy_at(mixture, ideal_gas, *INLET, COMPOSITION)
    h_out, _ = enthalpy_at(
        mixture, ideal_gas, result.T.to("K").magnitude, outlet_pressure, COMPOSITION
    )
    h.assert_close(result.power, FLOW * (h_out - h_in), 1.0e-6, f"{name}: power")


# ---------------------------------------------------------------------------------
# The difference between the three
# ---------------------------------------------------------------------------------


def test_the_power_sign_is_what_distinguishes_a_compressor_from_an_expander() -> None:
    """One convention, three machines, and the expander's power is negative.

    The same inlet, the same flow, and the same *magnitude* of pressure change in
    opposite directions. A model that got the efficiency the wrong way round would still
    produce a plausible number here and the sign is where it would show.
    """
    compressed = run("compressor", outlet_pressure=40.0e5, efficiency=0.75)
    expanded = run("expander", outlet_pressure=5.0e5, efficiency=0.75)

    assert compressed.power > 0.0, "a compressor's shaft power must be positive"
    assert expanded.power < 0.0, (
        "an expander's shaft power must be negative - the fluid is doing the work - and "
        "NeqSim reports the opposite sign on its own energy port, so a port that copied "
        "that convention would fail here"
    )
    assert expanded.T.to("K").magnitude < INLET[0] < compressed.T.to("K").magnitude, (
        "an expansion must cool the fluid and a compression must heat it"
    )


@pytest.mark.parametrize(
    ("name", "outlet_pressure"),
    [("compressor", 40.0e5), ("pump", 30.0e5), ("expander", 5.0e5)],
)
def test_a_real_machine_leaves_the_fluid_warmer_than_an_ideal_one(
    name: str, outlet_pressure: float
) -> None:
    """``T >= isentropic_temperature`` at every efficiency below one.

    An irreversibility adds entropy, and for a compressor that means the fluid leaves
    hotter than the ideal path would leave it; for an expander it means the fluid leaves
    *warmer than the ideal path* too, because less work was extracted. Both are the same
    inequality, and it is the one statement here that needs no reference to the source at
    all - it is a consequence of the second law.
    """
    # At unit efficiency the two temperatures are the same state, reached by two
    # different inversions - `eos.ps_flash` for the ideal one and `eos.ph_flash` for the
    # real one - so they agree to their own tolerances rather than to the last bit. The
    # slack is that tolerance carried into kelvin: an absolute enthalpy tolerance of
    # `1e-8 * |H|` over a heat capacity of order 100 J/(mol*K) is a temperature of order
    # `1e-6 K`. A tighter slack would be asserting that two solvers round alike.
    slack = 1.0e-6
    for efficiency in (0.5, 0.8, 1.0):
        result = run(name, outlet_pressure=outlet_pressure, efficiency=efficiency)
        assert (
            result.T.to("K").magnitude >= result.isentropic_temperature.to("K").magnitude - slack
        ), (
            f"{name} at eta={efficiency}: the real outlet is colder than the ideal one, "
            f"which no machine can produce"
        )


@pytest.mark.parametrize(("name", "outlet_pressure"), [("compressor", 40.0e5), ("expander", 5.0e5)])
def test_unit_efficiency_reaches_the_ideal_temperature(name: str, outlet_pressure: float) -> None:
    """At an efficiency of one the two reported temperatures are the same.

    They come from different flashes - one bisects an entropy, the other an enthalpy - so
    they agree to their own tolerances rather than to the last bit, and the bound below
    says so. At unity the two branches of the procedure meet, which is what pins the
    direction of the efficiency without a recorded number needing to be trusted.
    """
    result = run(name, outlet_pressure=outlet_pressure, efficiency=1.0)
    h.assert_close(
        result.T.to("K").magnitude,
        result.isentropic_temperature.to("K").magnitude,
        1.0e-6,
        f"{name} at unit efficiency",
    )


@pytest.mark.parametrize("efficiency", [0.5, 0.8])
def test_a_worse_efficiency_costs_more_and_runs_hotter(efficiency: float) -> None:
    """A smaller efficiency is a greater power and a hotter outlet, for a compressor.

    Asserted as a *direction* rather than as a number, so it holds for any state: this is
    what the efficiency means, and a model that had it inverted would reduce both.
    """
    better = run("compressor", outlet_pressure=40.0e5, efficiency=1.0)
    worse = run("compressor", outlet_pressure=40.0e5, efficiency=efficiency)

    assert worse.power > better.power, "a worse efficiency must cost more power"
    assert worse.T.to("K").magnitude > better.T.to("K").magnitude, (
        "a worse efficiency must leave the fluid hotter"
    )


@pytest.mark.parametrize("efficiency", [0.5, 0.8])
def test_a_worse_efficiency_delivers_less_and_runs_hotter(efficiency: float) -> None:
    """The same statement for an expander, and the two directions are opposite.

    A worse expander extracts **less** work - so its power is closer to zero - and leaves
    the fluid **warmer**. It is the mirror of the compressor's test above, and a model
    that used one rule for both would pass one of them and fail the other.
    """
    better = run("expander", outlet_pressure=5.0e5, efficiency=1.0)
    worse = run("expander", outlet_pressure=5.0e5, efficiency=efficiency)

    assert worse.power > better.power, (
        "a worse expander delivers less work, and its power is already negative - so it "
        "must be closer to zero, which is greater"
    )
    assert worse.T.to("K").magnitude > better.T.to("K").magnitude, (
        "a worse expander must leave the fluid hotter"
    )


# ---------------------------------------------------------------------------------
# Refusals
# ---------------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["compressor", "pump"])
def test_a_machine_that_would_lower_the_pressure_is_refused(name: str) -> None:
    """A compressor that lowers the pressure is an expander, and there is one."""
    with pytest.raises(InvalidInputError, match="above"):
        run(name, outlet_pressure=10.0e5)


def test_an_expander_that_would_raise_the_pressure_is_refused() -> None:
    """The mirror of the check above, and it is a separate model for that reason."""
    with pytest.raises(InvalidInputError, match="below"):
        run("expander", outlet_pressure=40.0e5)


@pytest.mark.parametrize("name", sorted(MODELS))
@pytest.mark.parametrize("efficiency", [0.0, 1.5])
def test_an_impossible_efficiency_is_refused(name: str, efficiency: float) -> None:
    """NeqSim clamps the efficiency to ``(0, 1]``; this refuses instead.

    A clamp turns a caller's mistake into a slightly different answer that looks
    deliberate, which is the opposite of what the check is for.
    """
    outlet = 40.0e5 if name != "expander" else 5.0e5
    with pytest.raises(OutOfRangeError):
        run(name, outlet_pressure=outlet, efficiency=efficiency)


@pytest.mark.parametrize(
    ("name", "outlet_pressure"),
    [("compressor", 40.0e5), ("pump", 30.0e5), ("expander", 5.0e5)],
)
def test_a_zero_flow_is_refused(name: str, outlet_pressure: float) -> None:
    """A flow of zero moles has no power to compute."""
    with pytest.raises(OutOfRangeError):
        run(name, outlet_pressure=outlet_pressure, n=0.0)


# ---------------------------------------------------------------------------------
# Cross-language
# ---------------------------------------------------------------------------------


@pytest.mark.requires_rust
@pytest.mark.parametrize(
    ("name", "outlet_pressure"),
    [("compressor", 40.0e5), ("pump", 30.0e5), ("expander", 5.0e5)],
)
def test_the_two_implementations_agree_with_matching_iteration_counts(
    name: str, outlet_pressure: float
) -> None:
    """Both implementations reach the same state, in the same number of steps.

    The iteration count is the load-bearing half, and for these models it is the sum over
    *two* flashes - an entropy bisection and an enthalpy one. Two implementations that
    agreed on a temperature reached by different paths have not been shown to agree on
    anything durable.
    """
    import azoth
    from azoth._dispatch import use_backend

    arguments: dict[str, Any] = {
        "mixture": a_mixture(),
        "ideal_gas": an_ideal_gas(),
        "T": Q(INLET[0], "K"),
        "P": Q(INLET[1], "Pa"),
        "n": Q(FLOW, "mol/s"),
        "z": COMPOSITION,
        "outlet_pressure": Q(outlet_pressure, "Pa"),
        "efficiency": 0.8,
    }
    call_it = getattr(azoth.process, name)
    with use_backend("python"):
        py = call_it(**arguments)
    with use_backend("rust"):
        rs = call_it(**arguments)

    context = f"{name}"
    h.assert_close(rs.T.to("K").magnitude, py.T.to("K").magnitude, 1.0e-8, f"{context}: T")
    h.assert_close(rs.power, py.power, 1.0e-6, f"{context}: power")
    h.assert_close(
        rs.isentropic_temperature.to("K").magnitude,
        py.isentropic_temperature.to("K").magnitude,
        1.0e-8,
        f"{context}: isentropic temperature",
    )
    assert rs.phase == py.phase, f"{context}: phase disagrees"
    assert rs.iterations == py.iterations, (
        f"{context}: python took {py.iterations} steps and Rust {rs.iterations}"
    )


@pytest.mark.parametrize(
    ("name", "outlet_pressure"),
    [("compressor", 40.0e5), ("pump", 30.0e5), ("expander", 5.0e5)],
)
def test_the_isentropic_outlet_comes_from_the_entropy_flash(
    name: str, outlet_pressure: float
) -> None:
    """``isentropic_temperature`` is ``eos.ps_flash``'s answer, computed here.

    The three models report where the fluid would get to at constant entropy, and that
    number is the whole basis of the efficiency. This computes it from the flash directly
    rather than trusting the model that reported it, which is the one cross-model check
    in this file that does not go through an enthalpy.
    """
    mixture, ideal_gas = a_mixture(), an_ideal_gas()
    result = run(name, outlet_pressure=outlet_pressure)

    s_in, _ = entropy_at(mixture, ideal_gas, *INLET, COMPOSITION)
    ideal = ps_flash(
        mixture, ideal_gas, Q(outlet_pressure, "Pa"), Q(s_in, "J/(mol*K)"), COMPOSITION
    )
    h.assert_close(
        result.isentropic_temperature.to("K").magnitude,
        ideal.T.to("K").magnitude,
        1.0e-8,
        f"{name}: isentropic temperature against eos.ps_flash",
    )


def test_this_file_covers_the_two_directions() -> None:
    """A check on this file rather than on the models.

    The contrast tests above are the reason the three share a file, and they only mean
    something if the three really do cover both directions of the one rule. If a later
    change gave one machine its own arithmetic, this is what would stop the file from
    quietly testing the same direction twice.
    """
    assert set(DIRECTIONS) == set(MODELS), "every machine must declare a direction"
    assert {direction for direction in DIRECTIONS.values()} == {
        "consuming",
        "producing",
    }, "the three machines must cover both directions of the efficiency rule"
