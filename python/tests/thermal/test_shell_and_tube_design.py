"""Shape tests for ``thermal.shell_and_tube_design`` beyond the spec cases.

The spec cases and their cross-implementation agreement are covered by ``test_cross_impl.py``;
this file holds the defects and refusals the class's permissive path does not make.
"""

from __future__ import annotations

from typing import Any

import pytest

import _helpers as h
from azoth import ureg
from azoth.core.errors import InvalidInputError
from azoth.core.record_outputs_gen import AppliedStandardRow
from azoth.core.result import ShellAndTubeDesignResult
from azoth.thermal import shell_and_tube_design

CALC_ID = "thermal.shell_and_tube_design"
Q = ureg.Quantity

SPEC = h.spec(CALC_ID)
ACTIVE = [c for c in h.all_tests(SPEC) if c["status"] == "active"]


def kwargs(case: dict[str, Any]) -> dict[str, Any]:
    return h.kwargs_for(SPEC, case["inputs"])


def test_the_result_carries_generated_rows() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    result = shell_and_tube_design(**kwargs(case))
    assert isinstance(result, ShellAndTubeDesignResult)
    assert isinstance(result.applied_standards[0], AppliedStandardRow)
    assert result.tube_count == 276.0
    assert result.shell_inside_diameter.magnitude == pytest.approx(508.0)
    assert result.mawp_shell_side.magnitude == pytest.approx(43.019559902200484)


def test_an_unknown_grade_falls_back_to_the_class_defaults() -> None:
    case = next(c for c in ACTIVE if c["id"] == "an_unknown_grade_falls_back_to_the_class_defaults")
    result = shell_and_tube_design(**kwargs(case))
    assert result.shell_material == "Carbon Steel SA-516-70"
    assert result.tube_material == "Carbon Steel SA-179"
    assert result.tube_thermal_conductivity.magnitude == pytest.approx(52.0)


def test_the_nace_text_is_never_compliant() -> None:
    case = next(
        c for c in ACTIVE if c["id"] == "the_nace_assessment_reads_every_grade_as_not_compliant"
    )
    result = shell_and_tube_design(**kwargs(case))
    assert result.sour_service_required is True
    assert result.shell_nace_compliant is False
    assert result.tube_nace_compliant is False
    assert len(result.nace_issues) == 2


def test_a_stainless_tube_prices_as_carbon() -> None:
    base = next(c for c in ACTIVE if c["id"] == "the_defaults_worked_example")
    stainless = next(
        c
        for c in ACTIVE
        if c["id"] == "a_stainless_tube_moves_the_pressure_design_and_not_the_cost"
    )
    base_result = shell_and_tube_design(**kwargs(base))
    stainless_result = shell_and_tube_design(**kwargs(stainless))
    assert base_result.material_cost == stainless_result.material_cost
    assert base_result.total_dry_weight == stainless_result.total_dry_weight
    assert stainless_result.mawp_tube_side.magnitude > base_result.mawp_tube_side.magnitude


def test_a_partial_thermal_side_is_refused() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["tube_side_density"] = Q(998.0, "kg/m**3")
    with pytest.raises(InvalidInputError) as excinfo:
        shell_and_tube_design(**arguments)
    assert "tube_side" in str(excinfo.value)


def test_a_bad_designation_is_refused() -> None:
    case = next(c for c in ACTIVE if c["type"] == "worked_example")
    arguments = kwargs(case)
    arguments["tema_designation"] = "AE"
    with pytest.raises(InvalidInputError) as excinfo:
        shell_and_tube_design(**arguments)
    assert excinfo.value.field() == "tema_designation"
