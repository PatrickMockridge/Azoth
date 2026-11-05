"""`process.expander`'s stated-power and outlet-temperature modes.

The spec's cases pin the two modes to the capture. These tests are for what a case cannot
say: the ordinary branch reports no efficiency, and the two modes are the isentropic
branch's, one at a time.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth._dispatch import available, resolve, use_backend
from azoth.core.errors import InvalidInputError
from azoth.process import expander as run_expander

MODEL_ID = "process.expander"
SPEC = _models_gen.model(MODEL_ID)

Q = ureg.Quantity


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


def test_the_public_wrapper_takes_every_declared_input() -> None:
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_expander).parameters)
    assert accepted == declared, sorted(declared - accepted)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_ordinary_branch_reports_no_efficiency(backend: str) -> None:
    result = _run(_case("the_efficiency_multiplies_the_step_and_not_divides_it"), backend)
    assert result.reported_isentropic_efficiency is None


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_stated_power_is_the_captures(backend: str) -> None:
    result = _run(_case("the_stated_power_reports_its_efficiency"), backend)
    assert result.outlet_t.to("K").magnitude == pytest.approx(246.50812370121898, rel=5e-4)
    assert result.reported_isentropic_efficiency == pytest.approx(0.8000000060644207, rel=5e-4)


@pytest.mark.parametrize("backend", sorted(available()))
def test_zero_power_is_the_isenthalpic_throttle(backend: str) -> None:
    result = _run(_case("zero_power_is_the_isenthalpic_throttle"), backend)
    assert result.reported_isentropic_efficiency == pytest.approx(0.0, abs=1e-12)
    assert result.outlet_t.to("K").magnitude == pytest.approx(263.30703321422186, rel=5e-4)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_outlet_temperature_solve_is_the_captures(backend: str) -> None:
    result = _run(_case("the_outlet_temperature_solve_finds_two_hundred_and_fifty"), backend)
    assert result.outlet_t.to("K").magnitude == pytest.approx(250.00000266550109, rel=5e-4)
    assert result.reported_isentropic_efficiency == pytest.approx(0.6354905789224377, rel=5e-4)


def test_both_modes_are_refused() -> None:
    base = h.model_kwargs(
        SPEC, _case("the_efficiency_multiplies_the_step_and_not_divides_it")["inputs"]
    )
    with pytest.raises(InvalidInputError, match="shaft_power"):
        resolve(MODEL_ID)(
            **base,
            shaft_power=Q(-1000.0, "W"),
            outlet_temperature=Q(280.0, "K"),
        )


def test_the_stepped_branch_refuses_the_modes() -> None:
    base = h.model_kwargs(SPEC, _case("the_stepped_polytropic_expander")["inputs"])
    with pytest.raises(InvalidInputError, match="shaft_power"):
        resolve(MODEL_ID)(**base, shaft_power=Q(-1000.0, "W"))
