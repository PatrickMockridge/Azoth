"""`process.heater`'s instrument design: the list, its counts and its refusals.

The spec's cases pin the four lists to the capture. These tests are for what a case cannot
say: the design is absent without its flag, and its inputs are refused without it.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
from azoth import _models_gen, ureg
from azoth._dispatch import available, resolve, use_backend
from azoth.core.errors import InvalidInputError
from azoth.process import heater as run_heater

MODEL_ID = "process.heater"
SPEC = _models_gen.model(MODEL_ID)

Q = ureg.Quantity

COUNT_FIELDS = (
    "instrument_count",
    "analog_input_count",
    "analog_output_count",
    "digital_input_count",
    "digital_output_count",
    "total_io_count",
    "safety_instrument_count",
)


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


def test_the_public_wrapper_takes_every_declared_input() -> None:
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_heater).parameters)
    assert accepted == declared, sorted(declared - accepted)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_design_is_absent_without_its_flag(backend: str) -> None:
    result = _run(_case("outlet_temperature_380"), backend)
    for name in COUNT_FIELDS:
        assert getattr(result, name) is None, name
    assert result.instruments == ()


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_electric_heater_list_is_the_captures(backend: str) -> None:
    result = _run(_case("the_electric_heater_instrument_list"), backend)
    assert result.instrument_count == pytest.approx(8.0)
    assert result.analog_input_count == pytest.approx(6.0)
    assert result.digital_input_count == pytest.approx(2.0)
    assert result.safety_instrument_count == pytest.approx(2.0)
    assert len(result.instruments) == 8
    assert result.instruments[0].tag == "TT-1001"
    assert result.instruments[5].tag == "TSHH-1006"
    assert result.instruments[5].sil_rating == pytest.approx(2.0)
    assert result.instruments[5].range_unit is None


def test_an_instrument_input_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="include_safety_instruments"):
        resolve(MODEL_ID)(**base, include_safety_instruments=False)


def test_an_unknown_heat_exchanger_type_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="heat_exchanger_type"):
        resolve(MODEL_ID)(**base, instrument_design=True, heat_exchanger_type="coalescer")


ELECTRICAL_FIELDS = (
    "shaft_power",
    "electrical_input_power",
    "apparent_power",
    "reactive_power",
    "power_factor",
    "rated_voltage",
    "frequency_hz",
    "full_load_current",
    "starting_current",
    "total_electrical_losses",
    "motor_rated_power",
    "motor_efficiency",
    "motor_power_factor",
    "motor_rated_speed",
    "motor_locked_rotor_multiplier",
    "total_auxiliary_power",
    "total_connected_power",
    "cable_length",
    "cable_cross_section",
    "cable_ampacity",
    "cable_voltage_drop",
    "cable_short_circuit_withstand",
    "switchgear_starter",
    "switchgear_type",
    "switchgear_circuit_breaker",
    "switchgear_fuse",
    "switchgear_rated_current",
    "switchgear_short_circuit",
    "switchgear_weight",
    "vfd_rated_power",
    "vfd_rated_current",
    "vfd_input_voltage",
    "vfd_topology",
    "vfd_pulse_configuration",
    "vfd_efficiency",
    "vfd_input_power_factor",
    "vfd_thd_current",
    "vfd_has_active_rectifier",
    "vfd_requires_input_filter",
    "vfd_heat_dissipation",
    "vfd_cooling_method",
    "vfd_weight",
    "vfd_enclosure_rating",
)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_electrical_outputs_are_absent_without_their_flag(backend: str) -> None:
    result = _run(_case("outlet_temperature_380"), backend)
    for name in ELECTRICAL_FIELDS:
        assert getattr(result, name) is None, name


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_five_kilowatt_electrical_design_is_the_captures(backend: str) -> None:
    result = _run(_case("the_five_kilowatt_heater_is_the_captures"), backend)
    assert result.motor_rated_power.to("W").magnitude == pytest.approx(5500.0)
    assert result.motor_efficiency == pytest.approx(0.91)
    assert result.rated_voltage == pytest.approx(400.0)
    assert result.full_load_current == pytest.approx(9.55498260672858, rel=5e-4)
    assert result.total_connected_power.to("W").magnitude == pytest.approx(
        6994.50549297655, rel=5e-4
    )


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_frequency_does_not_reach_the_motor(backend: str) -> None:
    """The base sizes the motor before it copies the design's frequency to it."""
    result = _run(_case("a_sixty_hertz_supply_moves_the_synchronous_speed"), backend)
    assert result.frequency_hz == pytest.approx(60.0)
    assert result.motor_rated_speed == pytest.approx(1470.0)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_hundred_kilowatt_cable_and_switchgear_are_the_captures(backend: str) -> None:
    result = _run(_case("the_hundred_kilowatt_heater_takes_the_690_volt_tap"), backend)
    assert result.cable_length.to("m").magnitude == pytest.approx(50.0)
    assert result.cable_cross_section == pytest.approx(35.0)
    assert result.cable_ampacity == pytest.approx(101.64817796366708, rel=5e-4)
    assert result.cable_voltage_drop == pytest.approx(0.5868336821889366, rel=5e-4)
    assert result.cable_short_circuit_withstand == pytest.approx(5.005)
    assert result.switchgear_starter == "Star-Delta"
    assert result.switchgear_type == "MCC"
    assert result.switchgear_circuit_breaker == pytest.approx(160.0)
    assert result.switchgear_fuse == pytest.approx(200.0)
    assert result.switchgear_rated_current == pytest.approx(160.0)
    assert result.switchgear_short_circuit == pytest.approx(25.0)
    assert result.switchgear_weight.to("kg").magnitude == pytest.approx(74.0)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_cable_upsizes_for_voltage_drop(backend: str) -> None:
    result = _run(_case("the_cable_upsizes_for_two_kilometres"), backend)
    assert result.cable_cross_section == pytest.approx(300.0)
    assert result.cable_voltage_drop == pytest.approx(4.609117155205686, rel=5e-4)
    # The ampacity stays at the base-selected 35 mm2 row; only the section upsizes.
    assert result.cable_ampacity == pytest.approx(101.64817796366708, rel=5e-4)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_shell_branch_carries_no_cable(backend: str) -> None:
    result = _run(_case("the_shell_and_tube_design_is_the_auxiliary_load"), backend)
    assert result.cable_cross_section == pytest.approx(0.0)
    assert result.cable_short_circuit_withstand == pytest.approx(0.0)
    assert result.switchgear_starter == "DOL"
    assert result.switchgear_circuit_breaker == pytest.approx(0.0)
    assert result.switchgear_short_circuit == pytest.approx(25.0)
    assert result.switchgear_weight.to("kg").magnitude == pytest.approx(0.0)


VFD_FIELDS = (
    "vfd_rated_power",
    "vfd_rated_current",
    "vfd_input_voltage",
    "vfd_topology",
    "vfd_pulse_configuration",
    "vfd_efficiency",
    "vfd_input_power_factor",
    "vfd_thd_current",
    "vfd_has_active_rectifier",
    "vfd_requires_input_filter",
    "vfd_heat_dissipation",
    "vfd_cooling_method",
    "vfd_weight",
    "vfd_enclosure_rating",
)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_vfd_outputs_are_absent_without_their_flag(backend: str) -> None:
    result = _run(_case("the_five_kilowatt_heater_is_the_captures"), backend)
    for name in VFD_FIELDS:
        assert getattr(result, name) is None, name


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_vfd_sets_the_starting_current_and_the_starter(backend: str) -> None:
    result = _run(_case("the_hundred_kilowatt_vfd_needs_an_input_filter"), backend)
    assert result.switchgear_starter == "VFD"
    assert result.starting_current == pytest.approx(95.09079001292898, rel=5e-4)
    assert result.starting_current == pytest.approx(result.full_load_current)
    assert result.vfd_topology == "2-level"
    assert result.vfd_pulse_configuration == "6-pulse"
    assert result.vfd_efficiency == pytest.approx(0.975)
    assert result.vfd_rated_current == pytest.approx(208.90963687782514, rel=5e-4)
    assert result.vfd_heat_dissipation.to("W").magnitude == pytest.approx(
        2750.0000000000027, rel=5e-4
    )


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_vfd_flags_are_the_captures(backend: str) -> None:
    filtered = _run(_case("the_hundred_kilowatt_vfd_needs_an_input_filter"), backend)
    assert filtered.vfd_requires_input_filter is True
    assert filtered.vfd_has_active_rectifier is False
    three_level = _run(_case("the_five_hundred_kilowatt_vfd_is_a_three_level_drive"), backend)
    assert three_level.vfd_topology == "3-level"
    assert three_level.vfd_requires_input_filter is False
    assert three_level.vfd_has_active_rectifier is False
    assert three_level.vfd_cooling_method == "Water"
    multi = _run(_case("the_multi_level_drive_uses_its_active_rectifier_table"), backend)
    assert multi.vfd_topology == "Multi-level"
    assert multi.vfd_has_active_rectifier is True
    assert multi.vfd_input_power_factor == pytest.approx(0.98)
    air = _run(_case("the_air_coolers_low_voltage_drive_can_use_an_active_rectifier"), backend)
    assert air.vfd_topology == "2-level"
    assert air.vfd_pulse_configuration == "AFE"
    assert air.vfd_has_active_rectifier is True
    assert air.vfd_thd_current == pytest.approx(5.0)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_shell_branch_ignores_the_vfd(backend: Any) -> None:
    case = _case("the_shell_and_tube_design_is_the_auxiliary_load")
    with use_backend(backend):
        result = resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]), use_vfd=True)
    for name in VFD_FIELDS:
        assert getattr(result, name) is None, name
    assert result.switchgear_starter == "DOL"


def test_a_vfd_ask_without_the_electrical_design_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="use_vfd"):
        resolve(MODEL_ID)(**base, use_vfd=True)


def test_a_cable_length_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="cable_length"):
        resolve(MODEL_ID)(**base, cable_length=Q(500.0, "m"))


def test_an_electrical_input_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="motor_poles"):
        resolve(MODEL_ID)(**base, motor_poles=2.0)


def test_an_unknown_efficiency_class_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="efficiency_class"):
        resolve(MODEL_ID)(**base, electrical_design=True, efficiency_class="ie5")


UTILITY_FIELDS = (
    "has_supply_temperature",
    "has_return_temperature",
    "has_approach_temperature",
    "has_heat_capacity_rate",
    "has_overall_heat_transfer_coefficient",
)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_utility_outputs_are_absent_without_their_flag(backend: str) -> None:
    result = _run(_case("outlet_temperature_380"), backend)
    for name in UTILITY_FIELDS:
        assert getattr(result, name) is None, name


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_utility_predicates_are_the_captures(backend: str) -> None:
    supply = _run(_case("a_utility_supply_temperature_is_specified"), backend)
    assert supply.has_supply_temperature is True
    assert supply.has_return_temperature is False
    zero = _run(_case("a_zero_heat_capacity_rate_is_not_specified"), backend)
    assert zero.has_heat_capacity_rate is False
    negative = _run(_case("a_negative_u_value_is_not_specified"), backend)
    assert negative.has_overall_heat_transfer_coefficient is False


def test_a_utility_input_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="utility_supply_temperature"):
        resolve(MODEL_ID)(**base, utility_supply_temperature=Q(450.0, "K"))


MECHANICAL_FIELDS = (
    "mechanical_selected_type",
    "mechanical_selected_area",
    "log_mean_temperature_difference",
    "mechanical_approach_temperature",
    "used_overall_heat_transfer_coefficient",
    "calculated_ua",
)


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_mechanical_outputs_are_absent_without_their_flag(backend: str) -> None:
    result = _run(_case("outlet_temperature_380"), backend)
    for name in MECHANICAL_FIELDS:
        assert getattr(result, name) is None, name
    assert result.mechanical_types == ()


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_no_utility_mechanical_design_is_the_captures(backend: str) -> None:
    result = _run(_case("mechanical_no_utility"), backend)
    assert result.mechanical_selected_type == "plate_and_frame"
    assert result.mechanical_selected_area.to("m**2").magnitude == pytest.approx(
        4.167901647148911, rel=5e-4
    )
    assert result.log_mean_temperature_difference.to("K").magnitude == pytest.approx(
        20.863381404017616, rel=5e-4
    )
    assert result.calculated_ua.to("W/K").magnitude == pytest.approx(4793.086894221247, rel=5e-4)
    assert len(result.mechanical_types) == 5
    assert result.mechanical_types[1].mechanical_type == "plate_and_frame"


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_stated_u_is_the_used_one(backend: str) -> None:
    result = _run(_case("mechanical_supply_and_u"), backend)
    assert result.used_overall_heat_transfer_coefficient.to(
        "W/(m**2*K)"
    ).magnitude == pytest.approx(300.0)


def test_a_mechanical_input_without_the_flag_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="mechanical_type"):
        resolve(MODEL_ID)(**base, mechanical_type="plate_and_frame")


def test_an_unknown_mechanical_type_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="mechanical_type"):
        resolve(MODEL_ID)(**base, mechanical_design=True, mechanical_type="coalescer")


def test_an_unknown_selection_criterion_is_refused() -> None:
    base = h.model_kwargs(SPEC, _case("outlet_temperature_380")["inputs"])
    with pytest.raises(InvalidInputError, match="mechanical_selection"):
        resolve(MODEL_ID)(**base, mechanical_design=True, mechanical_selection="cheapest")
