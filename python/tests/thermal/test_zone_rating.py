"""Shape tests for ``thermal.zone_rating`` beyond the spec cases.

The spec cases and their cross-implementation agreement are covered by
``test_cross_impl.py``; this file holds the record-input boundary the spec cannot state: an
interval field, an empty zone list, and a zone the permissive path would leave unanswered.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import InvalidInputError
from azoth.core.record_outputs_gen import ZoneRatingRow
from azoth.core.result import ZoneRatingResult
from azoth.thermal import zone_rating

CALC_ID = "thermal.zone_rating"
Q = ureg.Quantity

SPEC = h.spec(CALC_ID)
ACTIVE = [c for c in h.all_tests(SPEC) if c["status"] == "active"]


def call(case: dict[str, Any]) -> ZoneRatingResult:
    return zone_rating(**h.kwargs_for(SPEC, case["inputs"]))


def test_the_result_carries_generated_rows() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    result = call(case)
    assert isinstance(result, ZoneRatingResult)
    assert result.zone_count == 1.0
    assert isinstance(result.zones[0], ZoneRatingRow)
    assert result.zones[0].zone_name == "desuperheating"
    assert result.zones[0].required_area.magnitude == pytest.approx(0.6829250133843835)


def kwargs(case: dict[str, Any]) -> dict[str, Any]:
    return h.kwargs_for(SPEC, case["inputs"])


def test_a_celsius_lmtd_interval_reaches_the_same_kelvin_magnitude() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["zones"][0]["lmtd"] = Q(40.0, "delta_degC")
    result = zone_rating(**arguments)
    assert result.zones[0].lmtd.magnitude == pytest.approx(40.0)


def test_an_absolute_temperature_is_refused_for_the_lmtd_field() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["zones"][0]["lmtd"] = Q(40.0, "degC")
    with pytest.raises(Exception) as excinfo:
        zone_rating(**arguments)
    assert "lmtd" in str(excinfo.value)


def test_an_empty_zone_list_answers_no_zones() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["zones"] = []
    result = zone_rating(**arguments)
    assert result.zone_count == 0.0
    assert result.zones == ()


def test_a_zone_with_a_zero_viscosity_is_refused() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["zones"][0]["tube_viscosity"] = Q(0.0, "Pa*s")
    with pytest.raises(InvalidInputError) as excinfo:
        zone_rating(**arguments)
    assert excinfo.value.field() == "tube_viscosity"
