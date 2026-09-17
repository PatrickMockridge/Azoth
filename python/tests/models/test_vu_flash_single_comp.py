"""Tests for the ``eos.vu_flash_single_comp`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.core.result import Phase, VuFlashSingleCompResult
from azoth.eos import components, vu_flash_single_comp

MODEL_ID = "eos.vu_flash_single_comp"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def call(case: dict[str, Any]) -> VuFlashSingleCompResult:
    return vu_flash_single_comp(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    for name, expected_value in case["expected"].items():
        got = getattr(result, name)
        if name in ("T", "V"):
            got = got.to("K" if name == "T" else "m**3/mol").magnitude
        h.assert_close(got, expected_value, case["tolerance"], f"{case['id']} ({name})")
    assert result.phase is Phase.TWO_PHASE
    assert result.is_clean, "the cases are states whose volume is consistent"


def test_the_split_is_the_fraction_of_the_span_the_energy_is() -> None:
    """The lever rule, at any internal energy between the two saturated ends.

    The case where ``beta`` is one half could pass on a model that always answers one
    half, so this walks the span. The ends are not inputs anyone knows - they are what
    the model computes - so they are recovered from two probes and used to predict a
    third.
    """
    fluid, ideal_gas = components.mixture_of(["propane"])
    p = Q(1.0e6, "Pa")

    def at(u: float) -> VuFlashSingleCompResult:
        return vu_flash_single_comp(fluid, ideal_gas, p, Q(0.001, "m**3/mol"), Q(u, "J/mol"))

    (u_a, u_b) = (-10_000.0, -5_000.0)
    a, b = at(u_a), at(u_b)
    slope = (u_b - u_a) / (b.beta - a.beta)
    u_liq = u_a - a.beta * slope
    span = slope

    for fraction in (0.1, 0.25, 0.5, 0.75, 0.9):
        result = at(u_liq + fraction * span)
        h.assert_close(result.beta, fraction, 1e-5, f"the split at {fraction} of the span")
        assert result.phase is Phase.TWO_PHASE


def test_a_volume_the_answer_does_not_imply_is_reported() -> None:
    """``V`` is checked, not used: the pressure and the energy fix it for a pure component.

    NeqSim's ``VUflashSingleComp`` takes a volume specification and never reads it. This
    model reads it, and says so when the caller's is not the one the split implies -
    silently ignoring an input is the one thing that would make a caller's mistake
    invisible.
    """
    fluid, ideal_gas = components.mixture_of(["propane"])
    case = CASES[0]["inputs"]
    good = call(CASES[0])
    assert good.is_clean

    wrong = vu_flash_single_comp(
        fluid,
        ideal_gas,
        Q(case["P"], "Pa"),
        Q(good.V.to("m**3/mol").magnitude * 2.0, "m**3/mol"),
        Q(case["U"], "J/mol"),
    )
    assert not wrong.is_clean, "a volume the split does not imply must warn"
    assert any("OUT_OF_VALID_RANGE" in str(w.code) for w in wrong.warnings)


def test_the_states_the_saturation_line_does_not_reach_are_refused() -> None:
    fluid, ideal_gas = components.mixture_of(["propane"])
    pc = fluid.components[0].Pc.to("Pa").magnitude
    with pytest.raises(OutOfRangeError):
        vu_flash_single_comp(
            fluid, ideal_gas, Q(pc * 1.1, "Pa"), Q(0.001, "m**3/mol"), Q(-5000.0, "J/mol")
        )
    with pytest.raises(OutOfRangeError):
        vu_flash_single_comp(
            fluid, ideal_gas, Q(1.0e6, "Pa"), Q(0.001, "m**3/mol"), Q(-1.0e6, "J/mol")
        )


def test_a_mixture_is_refused_by_name() -> None:
    fluid, ideal_gas = components.mixture_of(["methane", "n-butane"])
    with pytest.raises(InvalidInputError):
        vu_flash_single_comp(
            fluid, ideal_gas, Q(1.0e6, "Pa"), Q(0.001, "m**3/mol"), Q(-5000.0, "J/mol")
        )
