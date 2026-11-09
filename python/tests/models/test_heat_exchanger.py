"""`process.heat_exchanger`'s mechanical design: the estimate, its selection and its refusals.

The spec's cases pin the surface to the capture. These tests are for what a case cannot say:
the estimate is absent without its flag, and its inputs are refused without it.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth._dispatch import available, resolve, use_backend
from azoth.core.errors import InvalidInputError
from azoth.process import heat_exchanger as run_heat_exchanger

MODEL_ID = "process.heat_exchanger"
SPEC = _models_gen.model(MODEL_ID)

Q = ureg.Quantity

MECHANICAL_FIELDS = (
    "mechanical_selected_type",
    "mechanical_selected_area",
    "log_mean_temperature_difference",
    "mechanical_approach_temperature",
    "used_overall_heat_transfer_coefficient",
    "calculated_ua",
)


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


def test_the_public_wrapper_takes_every_declared_input() -> None:
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_heat_exchanger).parameters)
    assert accepted == declared, sorted(declared - accepted)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_mechanical_outputs_are_absent_without_their_flag(backend: str) -> None:
    result = _run(_case("the_effectiveness_ntu_rating"), backend)
    for name in MECHANICAL_FIELDS:
        assert getattr(result, name) is None, name
    assert result.mechanical_types == ()


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_rating_state_is_the_captures(backend: str) -> None:
    result = _run(_case("the_rating_state_is_the_twostream_estimate"), backend)
    assert result.mechanical_selected_type == "plate_and_frame"
    assert result.mechanical_selected_area.to("m**2").magnitude == pytest.approx(0.1)
    assert result.log_mean_temperature_difference.to("K").magnitude == pytest.approx(
        34.38376494949978, rel=5e-4
    )
    assert result.mechanical_approach_temperature.to("K").magnitude == pytest.approx(
        13.75707990040786, rel=5e-4
    )
    # The class's 500, not the exchanger's stated 100.
    assert result.used_overall_heat_transfer_coefficient.to(
        "W/(m**2*K)"
    ).magnitude == pytest.approx(500.0)
    assert result.calculated_ua.to("W/K").magnitude == pytest.approx(100.0)
    assert len(result.mechanical_types) == 5
    assert result.mechanical_types[0].mechanical_type == "shell_and_tube"


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_area_divides_by_the_types_own_coefficient(backend: str) -> None:
    result = _run(_case("the_rating_state_is_the_twostream_estimate"), backend)
    shell, plate, air = result.mechanical_types[:3]
    assert shell.required_area.to("m**2").magnitude == pytest.approx(100.0 / 580.0, rel=5e-4)
    # Plate-and-frame's 0.08696 is under the class's 0.1 floor.
    assert plate.required_area.to("m**2").magnitude == pytest.approx(0.1)
    assert air.required_area.to("m**2").magnitude == pytest.approx(100.0 / 90.0, rel=5e-4)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_selection_is_the_stated_one_or_the_criterion(backend: str) -> None:
    manual = _run(_case("the_manual_selection_wins_over_the_area_ranking"), backend)
    assert manual.mechanical_selected_type == "shell_and_tube"
    lightest = _run(_case("the_lightest_candidate_is_selected"), backend)
    assert lightest.mechanical_selected_type == "plate_fin"
    smallest = _run(_case("the_smallest_pressure_drop_is_selected"), backend)
    assert smallest.mechanical_selected_type == "plate_fin"


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_pinned_state_reads_the_500_default_or_the_stated_ua(backend: str) -> None:
    default = _run(_case("a_pinned_state_reads_the_500_default"), backend)
    assert default.used_overall_heat_transfer_coefficient.to(
        "W/(m**2*K)"
    ).magnitude == pytest.approx(500.0)
    assert default.calculated_ua.to("W/K").magnitude == pytest.approx(500.0)
    assert default.mechanical_selected_area.to("m**2").magnitude == pytest.approx(
        500.0 / 1150.0, rel=5e-4
    )
    stated = _run(_case("a_pinned_state_reads_the_stated_ua"), backend)
    assert stated.calculated_ua.to("W/K").magnitude == pytest.approx(200.0)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_delta_t_state_reads_its_computed_conductance(backend: str) -> None:
    result = _run(_case("the_delta_t_state_reads_its_computed_conductance"), backend)
    assert result.ua_value.to("W/K").magnitude == pytest.approx(157.4384599391579, rel=5e-4)
    assert result.calculated_ua.to("W/K").magnitude == pytest.approx(157.4384599391579, rel=5e-4)


def test_a_mechanical_input_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("the_effectiveness_ntu_rating")["inputs"])
    with pytest.raises(InvalidInputError, match="mechanical_type"):
        resolve(MODEL_ID)(**base, mechanical_type="shell_and_tube")
    with pytest.raises(InvalidInputError, match="mechanical_selection"):
        resolve(MODEL_ID)(**base, mechanical_selection="min_weight")


def test_an_unknown_mechanical_type_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("the_effectiveness_ntu_rating")["inputs"])
    with pytest.raises(InvalidInputError, match="mechanical_type"):
        resolve(MODEL_ID)(**base, mechanical_design=True, mechanical_type="coalescer")


def test_an_unknown_selection_criterion_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("the_effectiveness_ntu_rating")["inputs"])
    with pytest.raises(InvalidInputError, match="mechanical_selection"):
        resolve(MODEL_ID)(**base, mechanical_design=True, mechanical_selection="cheapest")


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_dynamic_state_is_carried_by_the_outputs(backend: str) -> None:
    result = _run(_case("the_first_dynamic_step_seeds_the_wall"), backend)
    assert result.wall_temperature is not None
    assert result.wall_temperature.to("K").magnitude == pytest.approx(321.35, rel=1e-9)
    # A holdup was not stated, so the two CSTR temperatures stay absent.
    assert result.shell_fluid_temperature is None
    assert result.tube_fluid_temperature is None


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_holdup_outputs_are_present_only_with_a_holdup(backend: str) -> None:
    holdup = _run(_case("the_holdup_step_carries_three_temperatures"), backend)
    assert holdup.shell_fluid_temperature is not None
    assert holdup.tube_fluid_temperature is not None
    control = _run(_case("the_zero_holdup_control_is_the_wall_only"), backend)
    assert control.shell_fluid_temperature is None
    assert control.tube_fluid_temperature is None


def test_a_transient_input_without_time_step_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("the_first_dynamic_step_seeds_the_wall")["inputs"])
    base.pop("time_step")
    with pytest.raises(InvalidInputError, match="a transient input is stated without"):
        resolve(MODEL_ID)(**base)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_fallback_echoes_a_stated_state(backend: Any) -> None:
    base = h.model_kwargs(SPEC, _case("the_dynamic_gate_falls_back_to_the_rating")["inputs"])
    with use_backend(backend):
        result = resolve(MODEL_ID)(**base, wall_temperature=Q(330.0, "K"))
    assert result.wall_temperature.to("K").magnitude == pytest.approx(330.0)
    assert result.shell_fluid_temperature is None


def test_a_two_phase_side_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("the_first_dynamic_step_seeds_the_wall")["inputs"])
    # The hot side is the shell here; a 50/50 methane/n-butane at 250 K and 5 bar is two
    # phases, and the transient step needs one density and one Cp.
    base["hot_components"] = ["methane", "n-butane"]
    base["hot_in_z"] = [0.5, 0.5]
    base["hot_in_p"] = Q(5.0, "bar")
    base["hot_in_t"] = Q(250.0, "K")
    with pytest.raises(InvalidInputError, match="stream"):
        resolve(MODEL_ID)(**base)
