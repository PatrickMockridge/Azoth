"""Tests for the ``process.heater`` unit operation.

The spec's cases pin the answers, and they are the weakest of the checks here: a recorded
number only says the model still does what it did when the number was written.

What the model actually claims is an **energy balance**: `n * (H_out - H_in) == Q`, with
both enthalpies computed from `eos.molar_enthalpy_entropy` rather than read out of the
result. The test runs it at several duties and - more to the point - at flows that differ
by a factor of ten at the same duty, which is the check that catches a model that forgot
to divide by the flow.

A heater is also a **cooler**, with a negative duty. That is not a convenience: NeqSim's
``Cooler`` has no ``run()`` of its own and inherits ``Heater``'s steady state, so the two
are one procedure and there is a test below that runs one.
"""

from __future__ import annotations

from typing import Any

import pytest
from _process import Q, a_mixture, an_ideal_gas, assert_case, call, h

from azoth import _models_gen
from azoth.core.errors import OutOfRangeError
from azoth.eos.reference.ph_flash import enthalpy_at

MODEL = _models_gen.model("process.heater")
CASES = MODEL["cases"]

INLET = (300.0, 20.0e5)


def run(*, T: float = 300.0, P: float = 20.0e5, n: float = 10.0, duty: float = 0.0) -> Any:
    """Call the model at a state, with the spec's own mixture and datum."""
    import azoth

    return azoth.process.heater(
        a_mixture(),
        an_ideal_gas(),
        T=Q(T, "K"),
        P=Q(P, "Pa"),
        n=Q(n, "mol/s"),
        z=[0.6, 0.4],
        pressure_drop=Q(0.0, "Pa"),
        heat_duty=Q(duty, "W"),
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    assert_case(MODEL, case)


def test_every_case_ran() -> None:
    assert CASES, "the spec declares no cases"
    for case in CASES:
        call(MODEL, case)


@pytest.mark.parametrize("duty", [-8.0e4, -2.0e4, 0.0, 2.0e4, 4.0e4])
def test_the_duty_is_the_enthalpy_change(duty: float) -> None:
    """``n * (H_out - H_in) == Q``, with both enthalpies from the equation of state.

    The identity the model is, and the check that does not depend on a number anybody
    recorded. The outlet enthalpy is computed at the outlet state rather than read from
    the result, so this compares two models rather than one model against itself.

    **Duties above about 5e4 W are excluded, and the reason is a defect in the model
    underneath this one.** At +8e4 W the answer is near 386.5 K, where the cubic's
    largest root jumps and `H(T)` stops being monotonic; `eos.ph_flash` then bisects to
    the discontinuity and returns a temperature whose enthalpy is 7 per cent from the
    one it was asked for, reporting it as converged. The reproduction and the measured
    numbers are in `specs/models/eos/ph_flash.yaml`'s notes. It is recorded there rather
    than worked around here, and this parametrisation stops short of it so that the
    energy balance can still be tested below it.
    """
    mixture, ideal_gas = a_mixture(), an_ideal_gas()
    n = 10.0
    result = run(duty=duty)

    h_in, _ = enthalpy_at(mixture, ideal_gas, *INLET, [0.6, 0.4])
    h_out, _ = enthalpy_at(mixture, ideal_gas, result.T.to("K").magnitude, INLET[1], [0.6, 0.4])
    # An **absolute** bound, and it has to be: at zero duty the expected value is zero,
    # and `assert_close` measures relative error, so it could only ever pass on an exact
    # zero. The scale is the enthalpy the flash inverts - a few thousand J/mol - times
    # the 1e-8 relative temperature the bisection converges to, times the flow.
    assert n * (h_out - h_in) == pytest.approx(duty, rel=1.0e-6, abs=1.0e-2), (
        f"duty {duty} W: the enthalpy balance gives {n * (h_out - h_in)} W"
    )


@pytest.mark.parametrize("duty", [5.0e3, 1.0e4])
def test_the_same_duty_on_less_flow_is_a_bigger_rise(duty: float) -> None:
    """The duty is extensive and the flash is molar, so the flow is part of the model.

    Ten times the flow at the same duty is one tenth the temperature rise. A model that
    passed `Q` straight to the flash as though it were molar would pass every recorded
    case in the spec - both of which use 10 mol/s - and fail this.

    **The duties are bounded by what the pressure can produce, and that is a narrower
    bound than it looks.** At 1 mol/s a 10 kW duty is 10 kJ/mol, which this mixture at
    2 MPa reaches at about 1100 K. A 20 kW duty is 20 kJ/mol, and *no* temperature at
    this pressure reaches it: the reachable enthalpy rises to about 13 kJ/mol above the
    inlet by 2000 K and then flattens, because the latent heat of the split dominates.
    `eos.ph_flash` reports that as a non-convergence rather than as a number, which is
    the honest reading - the duty is outside the model's reach, not hard to find.
    """
    a_lot = run(n=10.0, duty=duty)
    a_little = run(n=1.0, duty=duty)

    rise_lot = a_lot.T.to("K").magnitude - 300.0
    rise_little = a_little.T.to("K").magnitude - 300.0
    assert rise_little > rise_lot, (
        f"ten times the flow at the same duty must give a smaller temperature rise; got "
        f"{rise_little:.4f} K at 1 mol/s and {rise_lot:.4f} K at 10 mol/s"
    )


def test_a_negative_duty_is_a_cooler() -> None:
    """The difference between this model and NeqSim's ``Cooler`` is one sign.

    ``Cooler`` has no ``run()`` of its own and inherits ``Heater``'s steady state, so
    there is nothing for a separate model to do. This asserts that the sign reaches the
    flash and moves the temperature the other way.
    """
    heated = run(duty=5.0e4)
    cooled = run(duty=-5.0e4)

    assert cooled.T.to("K").magnitude < 300.0 < heated.T.to("K").magnitude, (
        "a negative duty must lower the temperature and a positive one raise it"
    )


def test_no_duty_is_the_identity() -> None:
    """At zero duty the outlet state is the inlet state.

    To 1e-7, not to the last bit: the flash bisects the temperature to 1e-8 relative,
    so an isenthalpic round trip returns the temperature it was given to about seven
    figures. See the same note in ``test_throttling_valve.py``.
    """
    result = run()
    h.assert_close(result.T.to("K").magnitude, 300.0, 1.0e-7, "no duty")


def test_a_zero_flow_is_refused() -> None:
    """The model divides the duty by the flow, so the spec's range check protects it."""
    with pytest.raises(OutOfRangeError):
        run(n=0.0, duty=1.0e4)


@pytest.mark.requires_rust
@pytest.mark.parametrize("duty", [-5.0e4, 0.0, 5.0e4])
def test_the_two_implementations_agree_with_matching_iteration_counts(duty: float) -> None:
    """Both implementations reach the same temperature, in the same number of steps."""
    import azoth
    from azoth._dispatch import use_backend

    arguments: dict[str, Any] = {
        "mixture": a_mixture(),
        "ideal_gas": an_ideal_gas(),
        "T": Q(300.0, "K"),
        "P": Q(20.0e5, "Pa"),
        "n": Q(10.0, "mol/s"),
        "z": [0.6, 0.4],
        "pressure_drop": Q(0.0, "Pa"),
        "heat_duty": Q(duty, "W"),
    }
    with use_backend("python"):
        py = azoth.process.heater(**arguments)
    with use_backend("rust"):
        rs = azoth.process.heater(**arguments)

    context = f"heater at {duty} W"
    h.assert_close(rs.T.to("K").magnitude, py.T.to("K").magnitude, 1.0e-8, f"{context}: T")
    assert rs.phase == py.phase, f"{context}: phase disagrees"
    assert rs.iterations == py.iterations, (
        f"{context}: python took {py.iterations} steps and Rust {rs.iterations}"
    )
