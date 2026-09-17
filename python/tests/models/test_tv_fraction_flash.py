"""Tests for the ``eos.tv_fraction_flash`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.errors import InvalidInputError
from azoth.core.result import Phase, TvFractionFlashResult
from azoth.eos import components, tv_fraction_flash

MODEL_ID = "eos.tv_fraction_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def call(case: dict[str, Any]) -> TvFractionFlashResult:
    fluid, _ = components.mixture_of(case["inputs"]["components"])
    return tv_fraction_flash(
        fluid,
        Q(case["inputs"]["T"], "K"),
        case["inputs"]["fraction"],
        Q(case["inputs"]["P"], "Pa"),
        case["inputs"]["z"],
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    h.assert_close(
        result.P.to("Pa").magnitude, case["expected"]["P"], case["tolerance"], f"{case['id']} (P)"
    )
    h.assert_consistent(result, case["id"])


def test_the_answer_reproduces_the_fraction_it_was_asked_for() -> None:
    """The two cases sit at 0.5 and 0.9, so this walks the interior.

    The answer is defined by the fraction; a pressure is only the way there, so this
    asserts the *reported* fraction rather than the pressure.
    """
    fluid, _ = components.mixture_of(["methane", "n-butane"])
    for fraction in (0.05, 0.2, 0.5, 0.8, 0.95):
        result = tv_fraction_flash(fluid, Q(330.0, "K"), fraction, Q(2.5e6, "Pa"), [0.6, 0.4])
        h.assert_close(result.volume_fraction, fraction, 1e-6, f"the fraction at {fraction}")
        assert result.phase is Phase.TWO_PHASE


def test_the_pressure_falls_as_the_gas_fraction_rises() -> None:
    """The monotonicity the search relies on, asserted rather than assumed."""
    fluid, _ = components.mixture_of(["methane", "n-butane"])
    previous = float("inf")
    for fraction in (0.1, 0.3, 0.5, 0.7, 0.9):
        p = tv_fraction_flash(fluid, Q(330.0, "K"), fraction, Q(2.5e6, "Pa"), [0.6, 0.4])
        assert p.P.to("Pa").magnitude < previous, (
            f"the fraction {fraction} did not lower the pressure"
        )
        previous = p.P.to("Pa").magnitude


def test_a_fraction_outside_the_two_phase_region_is_refused() -> None:
    fluid, _ = components.mixture_of(["methane", "n-butane"])
    for fraction in (0.0, 1.0):
        with pytest.raises(InvalidInputError):
            tv_fraction_flash(fluid, Q(330.0, "K"), fraction, Q(2.5e6, "Pa"), [0.6, 0.4])


def test_the_two_backends_agree_on_every_spec_case() -> None:
    from azoth._dispatch import use_backend

    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        h.assert_close(py.P.to("Pa").magnitude, rs.P.to("Pa").magnitude, 1e-9, case["id"])
        h.assert_close(py.volume_fraction, rs.volume_fraction, 1e-9, case["id"])
        assert py.phase is rs.phase, f"{case['id']}: phase"
