"""`process.filter`'s loaded state: the booleans, the defaults and the refusals.

The spec's cases pin the two implementations to the captures. These tests are for what a
case cannot say: `replacement_required` has no slot in a generated case row, the absent
loading inputs are the class's zero state, and the two setters that throw are refused.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth._dispatch import available, resolve, use_backend
from azoth.core.errors import InvalidInputError, OutOfRangeError
from azoth.process import filter as run_filter

MODEL_ID = "process.filter"
SPEC = _models_gen.model(MODEL_ID)

Q = ureg.Quantity


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


def test_the_public_wrapper_takes_every_declared_input() -> None:
    """The harness drops a case input the wrapper does not accept; this cannot."""
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_filter).parameters)
    assert accepted == declared, sorted(declared - accepted)


@pytest.mark.parametrize("backend", sorted(available()))
def test_no_loading_inputs_is_the_class_s_zero_state(backend: str) -> None:
    result = _run(_case("pressure_drop_1_bar"), backend)
    assert result.solids_loading.to("kg").magnitude == 0.0
    assert result.loading_fraction == 0.0
    assert result.breakthrough_fraction == 0.0
    assert result.bypass_fraction == 0.0
    assert result.unrestricted_drop.to("Pa").magnitude == result.clean_drop.to("Pa").magnitude
    assert result.differential_pressure_utilization == pytest.approx(0.2)
    assert result.replacement_required is False


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_non_positive_capacity_is_the_sentinel(backend: str) -> None:
    result = _run(_case("a_non_positive_capacity_is_no_capacity"), backend)
    assert result.loading_fraction == 0.0
    assert result.breakthrough_fraction == 0.0
    assert result.unrestricted_drop.to("Pa").magnitude == result.clean_drop.to("Pa").magnitude


@pytest.mark.parametrize("backend", sorted(available()))
def test_replacement_required_reaches_both_branches(backend: str) -> None:
    terminal = _run(_case("the_terminal_drop_calls_for_replacement"), backend)
    assert terminal.differential_pressure_utilization == pytest.approx(4.0)
    assert terminal.replacement_required is True

    loading = _run(_case("loading_past_capacity_breaks_through"), backend)
    assert loading.loading_fraction == pytest.approx(1.1)
    assert loading.differential_pressure_utilization < 1.0
    assert loading.replacement_required is True

    boundary = _run(_case("breakthrough_starts_above_its_boundary"), backend)
    assert boundary.replacement_required is False


def test_a_non_positive_bypass_pressure_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(OutOfRangeError, match="bypass_cracking_pressure"):
        resolve(MODEL_ID)(**base, bypass_cracking_pressure=Q(0.0, "Pa"))


def test_a_non_positive_terminal_drop_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(OutOfRangeError, match="terminal_pressure_drop"):
        resolve(MODEL_ID)(**base, terminal_pressure_drop=Q(0.0, "Pa"))


def test_a_non_finite_loading_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(OutOfRangeError, match="solids_loading"):
        resolve(MODEL_ID)(**base, solids_loading=Q(float("nan"), "kg"))


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_step_accumulates_and_reports_its_residence(backend: str) -> None:
    result = _run(_case("the_first_step_accumulates_five_kilograms"), backend)
    assert result.solids_loading.to("kg").magnitude == pytest.approx(5.0)
    assert result.loading_fraction == pytest.approx(0.5)
    assert result.holdup_residence_time is not None
    assert result.holdup_residence_time.to("s").magnitude == pytest.approx(
        48.47531506255213, rel=1e-6
    )


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_pre_step_efficiency_is_in_the_accretion(backend: str) -> None:
    result = _run(_case("the_accretion_uses_the_pre_step_efficiency"), backend)
    assert result.solids_loading.to("kg").magnitude == pytest.approx(8.23264)
    assert result.current_removal_efficiency == pytest.approx(0.34993728)


def test_a_zero_flow_has_no_residence_time() -> None:
    base = h.model_kwargs(SPEC, _case("the_first_step_accumulates_five_kilograms")["inputs"])
    base["inlet_n"] = Q(0.0, "mol/s")
    for backend in sorted(available()):
        with use_backend(backend):
            result = resolve(MODEL_ID)(**base)
        assert result.holdup_residence_time is None, backend


def test_both_capture_sources_are_refused() -> None:
    base = h.model_kwargs(SPEC, _case("the_first_step_accumulates_five_kilograms")["inputs"])
    with pytest.raises(InvalidInputError, match="solids_loading_rate"):
        resolve(MODEL_ID)(**base, inlet_particle_concentration=1.0e-4)


DESIGN_FIELDS = (
    "design_pressure",
    "design_temperature",
    "required_elements",
    "vessel_inner_diameter",
    "vessel_outer_diameter",
    "vessel_length",
    "shell_thickness",
    "head_thickness",
    "effective_element_area",
    "calculated_face_velocity",
    "required_nozzle_diameter",
    "selected_nozzle_diameter",
    "calculated_nozzle_velocity",
    "empty_vessel_weight",
    "total_equipped_weight",
    "collapse_pressure_utilization",
    "differential_pressure_design_acceptable",
    "design_media_area",
)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_design_outputs_are_absent_without_their_flag(backend: str) -> None:
    result = _run(_case("pressure_drop_1_bar"), backend)
    for name in DESIGN_FIELDS:
        assert getattr(result, name) is None, name


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_coalescer_design_is_the_captures(backend: str) -> None:
    result = _run(_case("the_coalescer_design_is_the_captures"), backend)
    assert result.vessel_inner_diameter.to("m").magnitude == pytest.approx(0.5)
    assert result.vessel_length.to("m").magnitude == pytest.approx(1.5)
    assert result.shell_thickness.to("m").magnitude == pytest.approx(0.014832889206146602)
    assert result.selected_nozzle_diameter.to("m").magnitude == pytest.approx(0.08)
    assert result.differential_pressure_design_acceptable is True


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_basket_design_reports_four_warnings(backend: str) -> None:
    result = _run(_case("the_basket_strainer_reports_four_warnings"), backend)
    assert result.differential_pressure_design_acceptable is False
    design_warnings = [
        warning for warning in result.warnings if warning.field == "mechanical_design"
    ]
    assert len(design_warnings) == 4


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_granular_design_publishes_its_media_area(backend: str) -> None:
    result = _run(_case("the_activated_carbon_sizes_the_vessel_cross_section"), backend)
    assert result.required_elements == pytest.approx(1.0)
    assert result.design_media_area is not None
    assert result.design_media_area.to("m**2").magnitude == pytest.approx(
        3.9964880977495576, rel=5e-4
    )


def test_a_design_input_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(InvalidInputError, match="element_area"):
        resolve(MODEL_ID)(**base, element_area=Q(0.5, "m**2"))


def test_an_unknown_filter_type_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(InvalidInputError, match="filter_type"):
        resolve(MODEL_ID)(**base, mechanical_design=True, filter_type="nozzle")


def test_an_inverted_ld_window_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("pressure_drop_1_bar")["inputs"])
    with pytest.raises(InvalidInputError, match="min_ld_ratio"):
        resolve(MODEL_ID)(**base, mechanical_design=True, min_ld_ratio=6.0, max_ld_ratio=5.0)
