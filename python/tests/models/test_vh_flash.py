"""Tests for the ``eos.vh_flash`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.eos import components, vh_flash
from azoth.eos.reference._flash_property import property_at

MODEL_ID = "eos.vh_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def call(case: dict[str, Any]) -> Any:
    return vh_flash(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    for name in ("P", "T"):
        got = getattr(result, name).to("Pa" if name == "P" else "K").magnitude
        h.assert_close(got, case["expected"][name], case["tolerance"], f"{case['id']} ({name})")
    h.assert_consistent(result, case["id"])


def test_the_round_trip_inverts_volume_and_enthalpy() -> None:
    """Compute the volume and enthalpy at a state the test chose, ask for them back.

    The check that does not depend on knowing the answer: whatever the model does
    internally, the pressure and temperature that come out must be the ones that went
    in, or the two specifications were not both met.
    """
    fluid, ideal_gas = components.mixture_of(["methane", "n-butane"])
    z = [0.6, 0.4]
    for t_c, p_pa in ((400.0, 1.0e6), (330.0, 2.5e6), (300.0, 5.0e6), (450.0, 3.0e6)):
        volume, _ = property_at(fluid, ideal_gas, t_c, p_pa, z, "v")
        enthalpy, _ = property_at(fluid, ideal_gas, t_c, p_pa, z, "h")
        result = vh_flash(fluid, ideal_gas, Q(volume, "m**3/mol"), Q(enthalpy, "J/mol"), z)
        h.assert_close(result.P.to("Pa").magnitude, p_pa, 1e-3, f"{t_c} K ({p_pa} Pa)")
        h.assert_close(result.T.to("K").magnitude, t_c, 1e-3, f"{t_c} K ({p_pa} Pa)")


def test_the_two_backends_agree_on_every_spec_case() -> None:
    from azoth._dispatch import use_backend

    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        h.assert_close(py.P.to("Pa").magnitude, rs.P.to("Pa").magnitude, 1e-9, case["id"])
        h.assert_close(py.T.to("K").magnitude, rs.T.to("K").magnitude, 1e-9, case["id"])
        assert py.phase is rs.phase, f"{case['id']}: phase"
        assert py.beta == rs.beta or (
            py.beta is not None and rs.beta is not None and abs(py.beta - rs.beta) < 1e-9
        ), f"{case['id']}: beta"
