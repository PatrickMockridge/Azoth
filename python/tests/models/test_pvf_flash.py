"""Tests for the ``eos.pvf_flash`` model."""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth.core.errors import InvalidInputError
from azoth.core.result import Phase, PvfFlashResult
from azoth.eos import components, pvf_flash

MODEL_ID = "eos.pvf_flash"
Q = ureg.Quantity

SPEC = _models_gen.model(MODEL_ID)
CASES = SPEC["cases"]


def call(case: dict[str, Any]) -> PvfFlashResult:
    fluid, _ = components.mixture_of(case["inputs"]["components"])
    return pvf_flash(
        fluid,
        Q(case["inputs"]["P"], "Pa"),
        case["inputs"]["beta"],
        Q(case["inputs"]["temperature"], "K"),
        case["inputs"]["z"],
    )


@pytest.mark.parametrize("case", CASES, ids=lambda c: c["id"])
def test_spec_case(case: dict[str, Any]) -> None:
    result = call(case)
    h.assert_close(
        result.T.to("K").magnitude, case["expected"]["T"], case["tolerance"], f"{case['id']} (T)"
    )
    h.assert_close(result.beta, case["inputs"]["beta"], case["tolerance"], f"{case['id']} (beta)")
    h.assert_consistent(result, case["id"])


def test_the_answer_reproduces_the_fraction_it_was_asked_for() -> None:
    """The spec cases sit at 0.84 and 0.99 - both near the dew end - so this walks the interior.

    A model that answered a temperature for *any* fraction would pass both.
    """
    fluid, _ = components.mixture_of(["methane", "n-butane"])
    for fraction in (0.05, 0.2, 0.5, 0.8, 0.95):
        result = pvf_flash(fluid, Q(2.5e6, "Pa"), fraction, Q(330.0, "K"), [0.6, 0.4])
        h.assert_close(result.beta, fraction, 1e-7, f"the fraction at {fraction}")
        assert result.phase is Phase.TWO_PHASE


def test_the_endpoints_are_refused_by_name() -> None:
    """Zero and one are the bubble and dew points, which are other models."""
    fluid, _ = components.mixture_of(["methane", "n-butane"])
    for beta, which in ((0.0, "bubble"), (1.0, "dew")):
        with pytest.raises(InvalidInputError) as caught:
            pvf_flash(fluid, Q(2.5e6, "Pa"), beta, Q(330.0, "K"), [0.6, 0.4])
        assert which in str(caught.value)


def test_the_two_backends_agree_on_every_spec_case() -> None:
    from azoth._dispatch import use_backend

    for case in CASES:
        with use_backend("python"):
            py = call(case)
        with use_backend("rust"):
            rs = call(case)
        h.assert_close(py.T.to("K").magnitude, rs.T.to("K").magnitude, 1e-9, case["id"])
        h.assert_close(py.beta, rs.beta, 1e-9, case["id"])
        assert py.phase is rs.phase, f"{case['id']}: phase"
