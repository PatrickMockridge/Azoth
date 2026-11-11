"""`process.cooler`'s dynamic temperature control: the gate, the carried state and the refusals.

The spec's cases pin the arithmetic to `process_cooler_transient.tsv`. These tests are for
what a case cannot say: a transient input without a step, the active branch's refusals, and
the opening clamp.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth._dispatch import available, resolve, use_backend
from azoth.core.errors import InvalidInputError
from azoth.process import cooler as run_cooler

MODEL_ID = "process.cooler"
SPEC = _models_gen.model(MODEL_ID)

Q = ureg.Quantity


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


def test_the_public_wrapper_takes_every_declared_input() -> None:
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_cooler).parameters)
    assert accepted == declared, sorted(declared - accepted)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_dynamic_outputs_are_absent_on_the_fallback(backend: str) -> None:
    result = _run(_case("the_dynamic_gate_falls_back_to_the_steady_run"), backend)
    assert result.cooling_valve_opening is None
    assert result.equilibrium_outlet_temperature is None


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_second_step_carries_the_first_states(backend: str) -> None:
    first = _run(_case("the_manual_controller_opens_the_valve"), backend)
    assert first.cooling_valve_opening == pytest.approx(64.28571428571428)
    second = _run(_case("the_second_step_carries_the_opening"), backend)
    assert second.cooling_valve_opening == pytest.approx(74.48979591836734)
    assert second.outlet_t.to("K").magnitude < first.outlet_t.to("K").magnitude


def test_a_transient_input_without_time_step_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("the_design_point_is_the_equilibrium")["inputs"])
    base.pop("time_step")
    with pytest.raises(InvalidInputError, match="a transient input is stated without"):
        resolve(MODEL_ID)(**base)


def test_a_stated_outlet_temperature_is_refused_on_the_dynamic_branch() -> None:
    base = h.model_kwargs(SPEC, _case("the_design_point_is_the_equilibrium")["inputs"])
    with pytest.raises(InvalidInputError, match="outlet_temperature"):
        resolve(MODEL_ID)(**base, outlet_temperature=Q(320.0, "K"))


def test_the_design_point_is_required_on_the_dynamic_branch() -> None:
    base = h.model_kwargs(SPEC, _case("the_design_point_is_the_equilibrium")["inputs"])
    for name in (
        "design_mass_flow",
        "design_inlet_temperature",
        "design_outlet_temperature",
        "cooling_medium_temperature",
    ):
        base.pop(name)
    with pytest.raises(InvalidInputError, match="design_mass_flow"):
        resolve(MODEL_ID)(**base)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_opening_is_clamped(backend: Any) -> None:
    base = h.model_kwargs(SPEC, _case("the_manual_controller_opens_the_valve")["inputs"])
    base["cooling_valve_opening"] = 150.0
    base["requested_valve_opening"] = 150.0
    with use_backend(backend):
        result = resolve(MODEL_ID)(**base)
    assert result.cooling_valve_opening == pytest.approx(100.0)
