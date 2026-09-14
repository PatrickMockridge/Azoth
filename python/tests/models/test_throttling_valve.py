"""Tests for the ``process.throttling_valve`` unit operation.

The spec's cases pin the answers, and they are the weakest of the checks here: a recorded
number only says the model still does what it did when the number was written.

What the model actually claims is that **the enthalpy does not change**. That is testable
without a recorded number, and it is the only thing a throttling valve does: compute the
feed's enthalpy, compute the outlet's from `eos.molar_enthalpy_entropy` at the outlet
state, and require them equal at several drops. A model that flashed at the *inlet*
pressure and relabelled it would pass every recorded number and fail that.

The second claim is that the temperature *moves*, because a real gas has a non-zero
Joule-Thomson coefficient and a model that returned the inlet temperature would be wrong
everywhere the pressure changes.
"""

from __future__ import annotations

from typing import Any

import pytest
from _process import Q, a_mixture, an_ideal_gas, assert_case, call, h

from azoth import _models_gen
from azoth.core.errors import OutOfRangeError
from azoth.eos.reference.ph_flash import enthalpy_at

MODEL = _models_gen.model("process.throttling_valve")
CASES = MODEL["cases"]

INLET = (300.0, 20.0e5)


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    assert_case(MODEL, case)


def test_every_case_ran() -> None:
    """A suite that silently collected nothing passes every other test here."""
    assert CASES, "the spec declares no cases"
    for case in CASES:
        call(MODEL, case)


@pytest.mark.parametrize("drop", [0.0, 1.0e5, 5.0e5, 1.0e6])
def test_the_enthalpy_is_unchanged(drop: float) -> None:
    """``H(T_out, P_out) == H(T_in, P_in)``, computed from the equation of state.

    The identity an isenthalpic flash *is*, and the check that does not depend on a
    number anybody recorded. Both enthalpies come from `eos.molar_enthalpy_entropy` at
    the two states - the outlet one is not read out of the result - so this compares the
    model against a different model rather than against itself.
    """
    import azoth

    mixture, ideal_gas = a_mixture(), an_ideal_gas()
    t_in, p_in = INLET

    result = azoth.process.throttling_valve(
        mixture,
        ideal_gas,
        T=Q(t_in, "K"),
        P=Q(p_in, "Pa"),
        z=[0.6, 0.4],
        pressure_drop=Q(drop, "Pa"),
    )

    h_in, _ = enthalpy_at(mixture, ideal_gas, t_in, p_in, [0.6, 0.4])
    h_out, _ = enthalpy_at(mixture, ideal_gas, result.T.to("K").magnitude, p_in - drop, [0.6, 0.4])
    # 1e-6 and not 1e-8: the flash converges the *temperature* to 1e-8 relative, so the
    # outlet enthalpy recomputed at the returned temperature carries that much error and
    # no less. A tighter bound would be asserting something about the flash's stopping
    # rule rather than about the enthalpy being conserved.
    h.assert_close(h_out, h_in, 1.0e-6, f"enthalpy across a {drop} Pa drop")


def test_a_drop_cools_a_real_gas() -> None:
    """The Joule-Thomson effect, which is the reason this is a model and not a copy.

    Methane/n-butane at 300 K falls to 293.91 K across 5 bar. That is a *small* change,
    which is why it is worth asserting separately: a model that returned the inlet
    temperature would be indistinguishable from a correct one at zero drop, and this is
    the case that tells them apart.
    """
    import azoth

    result = azoth.process.throttling_valve(
        a_mixture(),
        an_ideal_gas(),
        T=Q(300.0, "K"),
        P=Q(20.0e5, "Pa"),
        z=[0.6, 0.4],
        pressure_drop=Q(5.0e5, "Pa"),
    )

    cooled = result.T.to("K").magnitude
    assert cooled < 300.0, (
        f"a real gas cools across a throttling valve and this one reached {cooled} K. "
        f"The isenthalpic branch is not being taken."
    )


def test_no_drop_is_the_identity() -> None:
    """At zero drop the outlet state is the inlet state.

    A round trip rather than a comparison: an isenthalpic flash at the pressure the
    enthalpy was measured at must return the temperature it was measured at. **To
    1e-7, not to the last bit** - the flash bisects the temperature until its interval
    is within 1e-8 relative, so 300 K comes back as 300 K to about seven figures. The
    pressure is arithmetic and is exact.
    """
    result = call(MODEL, {"id": "zero", "inputs": dict(CASES[1]["inputs"]), "expected": {}})
    h.assert_close(result.T.to("K").magnitude, 300.0, 1.0e-7, "no drop")
    h.assert_close(result.P.to("Pa").magnitude, 20.0e5, 1.0e-12, "no drop")


def test_a_negative_drop_is_refused() -> None:
    """A valve that raises the pressure is a compressor, and there is one of those."""
    import azoth

    with pytest.raises(OutOfRangeError):
        azoth.process.throttling_valve(
            a_mixture(),
            an_ideal_gas(),
            T=Q(300.0, "K"),
            P=Q(20.0e5, "Pa"),
            z=[0.6, 0.4],
            pressure_drop=Q(-1.0e5, "Pa"),
        )


@pytest.mark.requires_rust
@pytest.mark.parametrize("drop", [0.0, 5.0e5, 1.0e6])
def test_the_two_implementations_agree_with_matching_iteration_counts(drop: float) -> None:
    """Both implementations reach the same state, in the same number of steps.

    The iteration count is the load-bearing half: the flash bisects an enthalpy, so two
    implementations that reached the same temperature by different paths have not been
    shown to agree on anything durable.
    """
    import azoth
    from azoth._dispatch import use_backend

    arguments: dict[str, Any] = {
        "mixture": a_mixture(),
        "ideal_gas": an_ideal_gas(),
        "T": Q(300.0, "K"),
        "P": Q(20.0e5, "Pa"),
        "z": [0.6, 0.4],
        "pressure_drop": Q(drop, "Pa"),
    }
    with use_backend("python"):
        py = azoth.process.throttling_valve(**arguments)
    with use_backend("rust"):
        rs = azoth.process.throttling_valve(**arguments)

    context = f"throttling_valve across {drop} Pa"
    h.assert_close(rs.T.to("K").magnitude, py.T.to("K").magnitude, 1.0e-8, f"{context}: T")
    assert rs.phase == py.phase, f"{context}: phase disagrees"
    assert rs.iterations == py.iterations, (
        f"{context}: python took {py.iterations} steps and Rust {rs.iterations}"
    )
