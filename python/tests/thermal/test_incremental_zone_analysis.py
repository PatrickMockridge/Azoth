"""Shape tests for ``thermal.incremental_zone_analysis`` beyond the spec cases.

The spec cases and their cross-implementation agreement are covered by ``test_cross_impl.py``;
this file holds the record-input boundary the spec cannot state: an empty list, a zero flow,
and the class spellings the rows echo.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import InvalidInputError
from azoth.core.record_outputs_gen import IncrementalZoneRow
from azoth.core.result import IncrementalZoneAnalysisResult
from azoth.thermal import incremental_zone_analysis

CALC_ID = "thermal.incremental_zone_analysis"
Q = ureg.Quantity

SPEC = h.spec(CALC_ID)
ACTIVE = [c for c in h.all_tests(SPEC) if c["status"] == "active"]


def kwargs(case: dict[str, Any]) -> dict[str, Any]:
    return h.kwargs_for(SPEC, case["inputs"])


def test_the_result_carries_generated_rows() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    result = incremental_zone_analysis(**kwargs(case))
    assert isinstance(result, IncrementalZoneAnalysisResult)
    assert result.number_of_zones == 1.0
    assert isinstance(result.zones[0], IncrementalZoneRow)
    assert result.zones[0].tube_phase_regime == "LIQUID"
    assert result.total_required_area.magnitude == pytest.approx(2.1944646424361007)


def test_an_empty_zone_list_answers_the_sentinels() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["zones"] = []
    result = incremental_zone_analysis(**arguments)
    assert result.number_of_zones == 0.0
    assert result.zones == ()
    assert result.minimum_approach_temperature.magnitude == pytest.approx(1.7976931348623157e308)


def test_a_zero_flow_is_refused() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["zones"][0]["tube_mass_flow"] = Q(0.0, "kg/s")
    with pytest.raises(InvalidInputError) as excinfo:
        incremental_zone_analysis(**arguments)
    assert excinfo.value.field() == "tube_mass_flow"


def test_the_tube_regime_selects_the_correlation() -> None:
    case_cycle = next(c for c in ACTIVE if c["id"] == "the_condensing_evaporating_vapour_cycle")
    result = incremental_zone_analysis(**kwargs(case_cycle))
    assert [row.tube_phase_regime for row in result.zones] == [
        "CONDENSING",
        "EVAPORATING",
        "VAPOR",
    ]
    assert result.zones[0].tube_side_heat_transfer_coefficient.magnitude == pytest.approx(
        7667.176474980685
    )
    assert result.zones[2].tube_side_heat_transfer_coefficient.magnitude == pytest.approx(
        3397.728059324884
    )
