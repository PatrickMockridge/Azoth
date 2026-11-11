"""The agreement and expected-value comparisons, including their two exceptions."""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any, cast

from azoth._registry_gen import spec
from azoth.rgen import compare, inputs, runner

ANTOINE = {
    "A": 5.23243,
    "B": 891.0098,
    "C": 332.0975,
    "D": 0.0,
    "E": 0.0,
    "form": "pow10",
    "Tc": 190.56,
    "Pc": 4_599_000.0,
    "T": 150.0,
}


def _run(kwargs: Mapping[str, Any]) -> dict[str, Any]:
    calc_id = "eos.antoine_vapor_pressure"
    decoded = inputs.decode(calc_id, inputs.encode(calc_id, kwargs))
    return cast("dict[str, Any]", runner.run(calc_id, decoded, "python")["result"])


def test_identical_documents_agree() -> None:
    result = _run(ANTOINE)
    report = compare.between_runs(result, result, tolerance=1e-12, diagnostic_bound=None)
    assert report["disagreements"] == []
    assert report["fields_compared"] >= 2


def test_a_tampered_value_is_a_disagreement() -> None:
    left = _run(ANTOINE)
    right: dict[str, Any] = dict(left)
    right["p_sat"] = {"magnitude": 1.0, "unit": "Pa"}
    report = compare.between_runs(left, right, tolerance=1e-12, diagnostic_bound=None)
    assert [item["path"] for item in report["disagreements"]] == ["p_sat"]


def test_quantities_meet_in_si() -> None:
    left = {"dp": {"magnitude": 1000.0, "unit": "mm"}}
    right = {"dp": {"magnitude": 1.0, "unit": "m"}}
    report = compare.between_runs(left, right, tolerance=1e-12, diagnostic_bound=None)
    assert report["disagreements"] == []


def test_absent_on_both_sides_is_not_a_disagreement() -> None:
    left = {"re": None, "regime": None}
    right = {"re": None, "regime": None}
    report = compare.between_runs(left, right, tolerance=1e-12, diagnostic_bound=None)
    assert report["disagreements"] == []


def test_absent_on_one_side_is_a_disagreement() -> None:
    report = compare.between_runs({"re": None}, {"re": 1.0}, tolerance=1e-12, diagnostic_bound=None)
    assert len(report["disagreements"]) == 1


def test_an_uncompared_field_is_recorded_not_judged() -> None:
    report = compare.between_runs(
        {"iterations": 3}, {"iterations": 7}, tolerance=1e-12, diagnostic_bound=None
    )
    assert report["disagreements"] == []
    assert report["uncompared"][0]["path"] == "iterations"


def test_a_diagnostic_is_held_to_the_solver_bound() -> None:
    agree = compare.between_runs(
        {"residual": 1e-13}, {"residual": 2e-13}, tolerance=1e-9, diagnostic_bound=1e-12
    )
    assert agree["disagreements"] == []
    disagree = compare.between_runs(
        {"residual": 1e-6}, {"residual": 2e-6}, tolerance=1e-9, diagnostic_bound=1e-12
    )
    assert len(disagree["disagreements"]) == 1


def test_warnings_compare_by_code_and_field_not_message() -> None:
    left = {"warnings": [{"code": "TRANSITIONAL_FLOW", "field": "re", "message": "0.0"}]}
    right = {"warnings": [{"code": "TRANSITIONAL_FLOW", "field": "re", "message": "1e-06"}]}
    report = compare.between_runs(left, right, tolerance=1e-12, diagnostic_bound=None)
    assert report["disagreements"] == []


def test_a_warning_code_that_differs_is_a_disagreement() -> None:
    left = {"warnings": [{"code": "TRANSITIONAL_FLOW", "field": "re", "message": "x"}]}
    right: dict[str, Any] = {"warnings": []}
    report = compare.between_runs(left, right, tolerance=1e-12, diagnostic_bound=None)
    assert len(report["disagreements"]) == 1


def test_expected_values_pass_at_the_case_tolerance() -> None:
    result = _run(ANTOINE)
    report = compare.against_expected(
        result,
        {"p_sat": 929295.9523003073},
        spec("eos.antoine_vapor_pressure")["outputs"],
        tolerance=1e-12,
        diagnostic_bound=None,
    )
    assert report["disagreements"] == []


def test_the_fifteenth_digit_is_within_the_case_tolerance() -> None:
    """`223_tm_c4_wagner` expects 10274.724598532119 and both kernels answer ...53211."""
    kwargs = {
        "A": -7.76451,
        "B": 1.45838,
        "C": -2.7758,
        "D": -1.23303,
        "E": 0.0,
        "form": "wagner",
        "Tc": 531.1,
        "Pc": 2_950_000.0,
        "T": 300.0,
    }
    calc_id = "eos.antoine_vapor_pressure"
    decoded = inputs.decode(calc_id, inputs.encode(calc_id, kwargs))
    result = runner.run(calc_id, decoded, "python")["result"]
    report = compare.against_expected(
        result,
        {"p_sat": 10274.724598532119},
        spec("eos.antoine_vapor_pressure")["outputs"],
        tolerance=1e-12,
        diagnostic_bound=None,
    )
    assert report["disagreements"] == []


def test_an_expected_value_outside_the_tolerance_fails() -> None:
    result = _run(ANTOINE)
    report = compare.against_expected(
        result,
        {"p_sat": 1.0},
        spec("eos.antoine_vapor_pressure")["outputs"],
        tolerance=1e-12,
        diagnostic_bound=None,
    )
    assert len(report["disagreements"]) == 1
