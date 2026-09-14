"""Tests for the ``process.mixer`` unit operation.

The spec's cases pin the answers, and they are the weakest of the checks here.

What the model claims is three things, and each is testable without a recorded number:
**moles are conserved**, **the composition is the flow-weighted blend**, and **the
energy is conserved**. The third is the one that matters most, because the port changed
its expression: NeqSim sums *total* enthalpies in joules, while `eos.ph_flash` takes a
*molar* one, so this model passes the flow-weighted **mean**. A model that passed the sum
would be wrong by a factor of the total flow and would still produce a temperature - and
no recorded number in the spec would necessarily catch it, because the number would have
been recorded from the same wrong code.
"""

from __future__ import annotations

from typing import Any

import pytest
from _process import Q, a_mixture, an_ideal_gas, assert_case, call, h

from azoth import _models_gen
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.eos.reference.ph_flash import enthalpy_at

MODEL = _models_gen.model("process.mixer")
CASES = MODEL["cases"]


def run(
    *,
    T: list[float],
    P: list[float],
    n: list[float],
    z: list[list[float]],
) -> Any:
    """Call the model with the spec's own mixture and datum."""
    import azoth

    return azoth.process.mixer(
        a_mixture(),
        an_ideal_gas(),
        T=[Q(value, "K") for value in T],
        P=[Q(value, "Pa") for value in P],
        n=[Q(value, "mol/s") for value in n],
        z=z,
    )


#: Two feeds that differ in temperature, pressure and composition, so every one of the
#: three identities below is exercised by something other than a coincidence.
INLETS: dict[str, Any] = {
    "T": [300.0, 350.0],
    "P": [20.0e5, 15.0e5],
    "n": [6.0, 4.0],
    "z": [[0.6, 0.4], [0.4, 0.6]],
}


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    assert_case(MODEL, case)


def test_every_case_ran() -> None:
    assert CASES, "the spec declares no cases"
    for case in CASES:
        call(MODEL, case)


def test_moles_are_conserved() -> None:
    """The outlet flow is the sum of the inlets', at every set of flows."""
    for flows in ([6.0, 4.0], [1.0, 99.0], [5.0, 5.0]):
        result = run(**{**INLETS, "n": flows})
        h.assert_close(result.flow, sum(flows), 1.0e-12, f"flow with inlets {flows}")


def test_the_composition_is_the_flow_weighted_blend() -> None:
    """``z_out == sum(n_s z_s) / sum(n_s)`` - arithmetic a reader can redo by hand.

    At 6 and 4 mol/s of ``[0.6, 0.4]`` and ``[0.4, 0.6]`` that is ``0.52`` and ``0.48``.
    Nothing about it needs the equation of state, so it checks the blend rather than the
    flash.
    """
    result = run(**INLETS)
    h.assert_close(result.z_out[0], (6.0 * 0.6 + 4.0 * 0.4) / 10.0, 1.0e-12, "z_out[0]")
    h.assert_close(result.z_out[1], (6.0 * 0.4 + 4.0 * 0.6) / 10.0, 1.0e-12, "z_out[1]")


def test_the_energy_is_conserved() -> None:
    """``flow * H_out == sum_s (n_s * H_s)``, with each inlet's enthalpy computed here.

    **This is the check that would catch the mean-versus-sum mistake**, and it is the
    reason this model's test file is worth more than its spec cases. Each inlet's enthalpy
    comes from `eos.molar_enthalpy_entropy` at that inlet's own state, and the outlet's
    from the same model at the outlet state, so neither side is read out of the result.
    """
    mixture, ideal_gas = a_mixture(), an_ideal_gas()
    result = run(**INLETS)

    expected = 0.0
    for temperature, pressure, flow, composition in zip(
        INLETS["T"], INLETS["P"], INLETS["n"], INLETS["z"], strict=True
    ):
        enthalpy, _ = enthalpy_at(mixture, ideal_gas, temperature, pressure, composition)
        expected += flow * enthalpy

    h_out, _ = enthalpy_at(
        mixture,
        ideal_gas,
        result.T.to("K").magnitude,
        result.P.to("Pa").magnitude,
        list(result.z_out),
    )
    h.assert_close(result.flow * h_out, expected, 1.0e-6, "energy balance")


def test_the_outlet_pressure_is_the_lowest_inlet_pressure() -> None:
    """Not the highest and not an average: a mixer is a vessel.

    A feed arriving at 20 bar into a vessel held at 15 bar flashes as it enters, so
    nothing in it can be above the lowest pressure any feed arrives at. The maximum, or
    the arithmetic mean, would give a number that is arithmetically fine and describes a
    pump.
    """
    for pressures in ([20.0e5, 15.0e5], [15.0e5, 20.0e5], [20.0e5, 20.0e5]):
        result = run(**{**INLETS, "P": pressures})
        h.assert_close(result.P.to("Pa").magnitude, min(pressures), 1.0e-12, "outlet P")


def test_two_identical_feeds_double_the_flow_and_change_nothing_else() -> None:
    """A mixture of a stream with itself is itself.

    The case that catches a mean-versus-sum mistake without relying on a flash: a model
    that passed the summed enthalpy would produce a temperature that is not 300 K, and
    the recorded case for this in the spec says 300 K exactly.
    """
    result = run(
        T=[300.0, 300.0],
        P=[20.0e5, 20.0e5],
        n=[6.0, 4.0],
        z=[[0.6, 0.4], [0.6, 0.4]],
    )
    # 1e-7, not the last bit: the mix is an isenthalpic flash, and that bisects the
    # temperature to 1e-8 relative.
    h.assert_close(result.T.to("K").magnitude, 300.0, 1.0e-7, "identical feeds: T")
    h.assert_close(result.P.to("Pa").magnitude, 20.0e5, 1.0e-12, "identical feeds: P")
    h.assert_close(result.flow, 10.0, 1.0e-12, "identical feeds: flow")
    h.assert_close(result.z_out[0], 0.6, 1.0e-12, "identical feeds: z_out")


def test_inlets_that_disagree_about_how_many_streams_there_are_are_refused() -> None:
    """Two temperatures and three pressures is not a set of streams."""
    with pytest.raises(InvalidInputError):
        run(T=[300.0, 350.0], P=[20.0e5, 15.0e5, 10.0e5], n=[6.0, 4.0], z=INLETS["z"])


def test_a_composition_matrix_of_the_wrong_shape_is_refused() -> None:
    """One row per inlet, or it is not a blend of these streams."""
    with pytest.raises(InvalidInputError):
        run(**{**INLETS, "z": [[0.6, 0.4]]})


def test_a_zero_flow_inlet_is_refused() -> None:
    """A zero-flow inlet contributes no moles, and its composition row means nothing."""
    with pytest.raises(OutOfRangeError):
        run(**{**INLETS, "n": [0.0, 10.0]})


@pytest.mark.requires_rust
def test_the_two_implementations_agree_with_matching_iteration_counts() -> None:
    """Both implementations reach the same blend, in the same number of steps."""
    import azoth
    from azoth._dispatch import use_backend

    arguments = {
        "mixture": a_mixture(),
        "ideal_gas": an_ideal_gas(),
        "T": [Q(value, "K") for value in INLETS["T"]],
        "P": [Q(value, "Pa") for value in INLETS["P"]],
        "n": [Q(value, "mol/s") for value in INLETS["n"]],
        "z": INLETS["z"],
    }
    with use_backend("python"):
        py = azoth.process.mixer(**arguments)
    with use_backend("rust"):
        rs = azoth.process.mixer(**arguments)

    assert rs.phase == py.phase, "phase disagrees"
    assert rs.iterations == py.iterations, (
        f"python took {py.iterations} steps and Rust {rs.iterations}"
    )
    for index in range(2):
        h.assert_close(rs.z_out[index], py.z_out[index], 1.0e-12, f"z_out[{index}]")
    h.assert_close(rs.T.to("K").magnitude, py.T.to("K").magnitude, 1.0e-8, "T")
    h.assert_close(rs.flow, py.flow, 1.0e-12, "flow")
