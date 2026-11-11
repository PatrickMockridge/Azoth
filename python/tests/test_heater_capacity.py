"""`process.heater`'s outlet-temperature constraints: the booleans and the absent branches.

The six capacity fields are `Option`s or booleans, and neither a `null` nor a boolean has a
`TestCase` slot, so this file reads the raw spec. The formulas and the hard-limit quirk are
held by the cases; what is here is the flag-off branch and the wrapper's signature, which
`_helpers.py` would silently drop an input through if it were missing.
"""

from __future__ import annotations

import inspect
from typing import Any

import pytest

import _helpers as h
import azoth
from azoth import _models_gen
from azoth._dispatch import available, resolve, use_backend
from azoth.process import heater as run_heater

MODEL_ID = "process.heater"
SPEC = _models_gen.model(MODEL_ID)

q = azoth.ureg.Quantity

FLAGGED = [
    case
    for case in SPEC["cases"]
    if "capacity_exceeded" in case["expected"] or "hard_limit_exceeded" in case["expected"]
]


def _case(case_id: str) -> dict[str, Any]:
    return next(case for case in SPEC["cases"] if case["id"] == case_id)


def _run(case: dict[str, Any], backend: Any) -> Any:
    with use_backend(backend):
        return resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))


@pytest.mark.parametrize("case", FLAGGED, ids=lambda case: case["id"])
def test_the_capacity_flags_are_the_recorded_ones(case: dict[str, Any]) -> None:
    """`isCapacityExceeded` and `isHardLimitExceeded`, read from the raw spec."""
    result = resolve(MODEL_ID)(**h.model_kwargs(SPEC, case["inputs"]))
    expected = case["expected"]
    if "capacity_exceeded" in expected:
        assert result.capacity_exceeded is expected["capacity_exceeded"]
    if "hard_limit_exceeded" in expected:
        assert result.hard_limit_exceeded is expected["hard_limit_exceeded"]


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_heater_without_limits_has_no_per_constraint_outputs(backend: str) -> None:
    result = _run(_case("the_outlet_temperature_limits_are_absent_without_setters"), backend)
    assert result.max_outlet_temperature_utilization is None
    assert result.min_outlet_temperature_utilization is None
    assert result.bottleneck_constraint is None
    assert result.max_utilization == 0.0
    assert result.capacity_exceeded is False
    assert result.hard_limit_exceeded is False


@pytest.mark.parametrize("backend", sorted(available()))
def test_a_zero_minimum_is_not_a_constraint(backend: str) -> None:
    result = _run(_case("a_minimum_of_zero_is_not_a_constraint"), backend)
    assert result.min_outlet_temperature_utilization is None
    assert result.bottleneck_constraint is None


@pytest.mark.parametrize("backend", sorted(available()))
def test_only_the_minimum_trips_the_hard_flag(backend: str) -> None:
    maximum = _run(_case("a_maximum_below_the_outlet_is_exceeded_but_not_hard"), backend)
    assert maximum.capacity_exceeded is True
    assert maximum.hard_limit_exceeded is False

    minimum = _run(_case("a_minimum_above_the_outlet_is_hard_exceeded"), backend)
    assert minimum.capacity_exceeded is True
    assert minimum.hard_limit_exceeded is True


@pytest.mark.parametrize("backend", sorted(available()))
def test_the_minimum_is_the_inverse_utilization(backend: str) -> None:
    case = _case("a_minimum_below_the_outlet_is_its_inverse_utilization")
    result = _run(case, backend)
    assert result.min_outlet_temperature_utilization == pytest.approx(
        case["expected"]["min_outlet_temperature_utilization"]
    )
    assert result.bottleneck_constraint == "minOutletTemperature"


def test_the_public_wrapper_takes_every_declared_input() -> None:
    """The harness drops a case input the wrapper does not accept; this cannot."""
    declared = set(SPEC["inputs"])
    accepted = set(inspect.signature(run_heater).parameters)
    assert accepted == declared, sorted(declared - accepted)
