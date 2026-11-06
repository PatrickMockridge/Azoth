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
