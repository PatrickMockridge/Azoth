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
