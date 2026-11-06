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
