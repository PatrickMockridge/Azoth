"""Tests for the ``eos.vs_flash`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.eos import components, vs_flash
from azoth.eos.reference._flash_property import property_at

MODEL_ID = "eos.vs_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def call(case: dict[str, Any]) -> Any:
    return vs_flash(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    for name in ("P", "T"):
        got = getattr(result, name).to("Pa" if name == "P" else "K").magnitude
        h.assert_close(got, case["expected"][name], case["tolerance"], f"{case['id']} ({name})")
    h.assert_consistent(result, case["id"])


def test_both_specifications_hold_at_the_answer() -> None:
    """The volume and the entropy at the answer are the two that were asked for.

    A model that satisfied one and not the other would pass every round-trip case whose
    two happened to be consistent, so this recomputes both from the answer.
    """
    fluid, ideal_gas = components.mixture_of(["methane", "propane", "n-butane"])
    z = [0.5, 0.3, 0.2]
    for t_c, p_pa in ((380.0, 1.5e6), (320.0, 4.0e6), (430.0, 8.0e6)):
        volume, _ = property_at(fluid, ideal_gas, t_c, p_pa, z, "v")
        entropy, _ = property_at(fluid, ideal_gas, t_c, p_pa, z, "s")
        result = vs_flash(fluid, ideal_gas, Q(volume, "m**3/mol"), Q(entropy, "J/mol/K"), z)
        found_v, _ = property_at(
            fluid,
            ideal_gas,
            result.T.to("K").magnitude,
            result.P.to("Pa").magnitude,
            z,
            "v",
        )
        found_s, _ = property_at(
            fluid,
            ideal_gas,
            result.T.to("K").magnitude,
            result.P.to("Pa").magnitude,
            z,
            "s",
        )
        h.assert_close(found_v, volume, 1e-5, f"{t_c} K ({p_pa} Pa) volume")
        h.assert_close(found_s, entropy, 1e-5, f"{t_c} K ({p_pa} Pa) entropy")


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
