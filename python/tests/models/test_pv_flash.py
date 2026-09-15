"""Tests for the ``eos.pv_flash`` model.

The spec's cases pin the answers; the round-trip below proves the model *inverts* the
volume without knowing the temperature - compute ``V(T)`` at a temperature the test
chose, ask for it back, and require the temperature that comes out to be the one that
went in.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.eos import IdealGasModel, Mixture, component, mixture, pv_flash
from azoth.eos.reference._flash_property import property_at

MODEL_ID = "eos.pv_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]

METHANE = component("methane")
BUTANE = component("n-butane")


def a_mixture() -> Mixture:
    return mixture([METHANE, BUTANE], kij={(0, 1): 0.01289789})


def an_ideal_gas() -> IdealGasModel:
    return IdealGasModel(
        cp_a=(3.0, 5.0),
        cp_b=(0.0, 0.0),
        cp_c=(0.0, 0.0),
        cp_d=(0.0, 0.0),
        cp_e=(0.0, 0.0),
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = pv_flash(**h.model_kwargs(SPEC, case["inputs"]))
    for name, expected_value in case["expected"].items():
        got = getattr(result, name)
        if name == "T":
            got = got.to("K").magnitude
        h.assert_close(got, expected_value, case["tolerance"], f"{case['id']} ({name})")


def test_the_round_trip_inverts_the_volume() -> None:
    fluid, ideal_gas = a_mixture(), an_ideal_gas()
    z = [0.6, 0.4]
    volume, _ = property_at(fluid, ideal_gas, 300.0, 1.0e6, z, "v")
    result = pv_flash(fluid, ideal_gas, Q(1.0e6, "Pa"), Q(volume, "m**3/mol"), z)
    h.assert_close(result.T.to("K").magnitude, 300.0, 1e-5, "T")
