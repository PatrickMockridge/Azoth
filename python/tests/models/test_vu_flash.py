"""Tests for the ``eos.vu_flash`` model.

The spec's cases pin the answers; the round-trip below proves the model *inverts* the
volume and internal energy together - compute them at a state the test chose, ask for them
back, and require the pressure and temperature that come out to be the ones that went in.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.eos import IdealGasModel, Mixture, components, vu_flash
from azoth.eos.reference._flash_property import property_at

MODEL_ID = "eos.vu_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def a_mixture() -> tuple[Mixture, IdealGasModel]:
    return components.mixture_of(["methane", "n-butane"])


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = vu_flash(**h.model_kwargs(SPEC, case["inputs"]))
    for name, expected_value in case["expected"].items():
        got = getattr(result, name)
        if name in ("P", "T"):
            got = got.to("Pa" if name == "P" else "K").magnitude
        h.assert_close(got, expected_value, case["tolerance"], f"{case['id']} ({name})")


def test_the_round_trip_inverts_volume_and_energy() -> None:
    fluid, ideal_gas = a_mixture()
    z = [0.6, 0.4]
    volume, _ = property_at(fluid, ideal_gas, 400.0, 1.0e6, z, "v")
    energy, _ = property_at(fluid, ideal_gas, 400.0, 1.0e6, z, "u")
    result = vu_flash(fluid, ideal_gas, Q(volume, "m**3/mol"), Q(energy, "J/mol"), z)
    h.assert_close(result.P.to("Pa").magnitude, 1.0e6, 1e-3, "P")
    h.assert_close(result.T.to("K").magnitude, 400.0, 1e-3, "T")


def test_an_answer_the_iteration_never_settled_at_says_so() -> None:
    """The specification's acceptance is looser than the iteration's, and the gap is reported.

    A liquid's volume barely moves with pressure, so the relative volume error is a poor
    judge of a pressure: pure propane at 250 K, the 10 bar state comes back as 15.5 bar -
    55% out - with a relative volume error of 5.5e-4, under the 1e-3 the specification is
    accepted at. NeqSim 3.20.0 returns that number and sets `lastRunConverged = false`;
    returning it *without* a word is the one thing neither does.
    """
    fluid, ideal_gas = components.mixture_of(["propane"])
    z = [1.0]
    volume, _ = property_at(fluid, ideal_gas, 250.0, 1.0e6, z, "v")
    energy, _ = property_at(fluid, ideal_gas, 250.0, 1.0e6, z, "u")
    result = vu_flash(fluid, ideal_gas, Q(volume, "m**3/mol"), Q(energy, "J/mol"), z)

    assert result.P.to("bar").magnitude > 12.0, (
        "the pressure is 15.5 bar rather than the 10 asked for; if this ever converges, "
        "the warning below should stop firing too"
    )
    codes = {str(warning.code) for warning in result.warnings}
    assert any("NOT_CONVERGED" in code for code in codes), (
        f"an unconverged answer must say so; got {codes}"
    )
