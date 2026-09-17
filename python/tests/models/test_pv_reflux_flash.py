"""Tests for the ``eos.pv_reflux_flash`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.result import PvRefluxFlashResult
from azoth.eos import components, pv_reflux_flash

MODEL_ID = "eos.pv_reflux_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def call(case: dict[str, Any]) -> PvRefluxFlashResult:
    fluid, _ = components.mixture_of(case["inputs"]["components"])
    return pv_reflux_flash(
        fluid,
        Q(case["inputs"]["P"], "Pa"),
        case["inputs"]["reflux"],
        case["inputs"]["phase"],
        Q(case["inputs"]["temperature"], "K"),
        case["inputs"]["z"],
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    h.assert_close(
        result.T.to("K").magnitude, case["expected"]["T"], case["tolerance"], f"{case['id']} (T)"
    )
    assert result.residual <= SPEC["algorithm"]["tolerance"]
    h.assert_consistent(result, case["id"])


def test_the_phase_named_decides_the_answer() -> None:
    """The two ratios are reciprocals, and a ratio of one is equal phase amounts.

    A model that ignored ``phase`` would answer the first half with one temperature for
    both requests, and the second with two.
    """
    fluid, _ = components.mixture_of(["methane", "n-butane"])
    z = [0.6, 0.4]
    p = Q(2.5e6, "Pa")

    vapour = pv_reflux_flash(fluid, p, 0.1873586001338361, "vapour", Q(330.0, "K"), z)
    liquid = pv_reflux_flash(fluid, p, 5.3373584094120519, "liquid", Q(330.0, "K"), z)
    h.assert_close(vapour.T.to("K").magnitude, liquid.T.to("K").magnitude, 1e-5, "reciprocals")

    half_vapour = pv_reflux_flash(fluid, p, 1.0, "vapour", Q(330.0, "K"), z)
    half_liquid = pv_reflux_flash(fluid, p, 1.0, "liquid", Q(330.0, "K"), z)
    h.assert_close(half_vapour.T.to("K").magnitude, half_liquid.T.to("K").magnitude, 1e-4, "half")
    assert half_vapour.beta is not None
    h.assert_close(half_vapour.beta, 0.5, 1e-6, "half beta")


def test_the_two_backends_agree_on_every_spec_case() -> None:
    from azoth._dispatch import use_backend

    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        h.assert_close(py.T.to("K").magnitude, rs.T.to("K").magnitude, 1e-9, case["id"])
        assert py.phase is rs.phase, f"{case['id']}: phase"
