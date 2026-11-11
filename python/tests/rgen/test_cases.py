"""The declared case set: the union, the deduplication and the non-value records."""

from __future__ import annotations

from azoth.rgen import cases


def test_the_antoine_union_is_five_distinct_cases() -> None:
    declared = cases.declared_cases("eos.antoine_vapor_pressure")
    assert len(declared) == 5
    identifiers = {case["id"] for case in declared.cases}
    assert "i_pentane_antoine_dippr101_against_neqsim" in identifiers


def test_a_spec_reference_case_is_folded_into_the_external_copy() -> None:
    declared = cases.declared_cases("eos.antoine_vapor_pressure")
    ethylene = next(case for case in declared.cases if "ethylene" in case["id"])
    assert ethylene["source"] == "validation"
    assert "ethylene_pow10kpa" in ethylene.get("spec_ids", [])


def test_a_property_test_is_not_a_value_case() -> None:
    declared = cases.declared_cases("eos.antoine_vapor_pressure")
    not_value = {record["id"]: record for record in declared.not_value}
    assert not_value["round_trip_units"]["status"] == "not_a_value_case"


def test_a_calc_with_no_external_case_uses_its_spec_cases() -> None:
    declared = cases.declared_cases("thermal.conduction_plane_wall")
    assert len(declared) == 3
    assert {case["source"] for case in declared.cases} == {"spec"}
    not_value = {record["id"] for record in declared.not_value}
    assert {"monotonic", "round_trip_units"} <= not_value


def test_the_strictest_declared_tolerance_is_the_smallest() -> None:
    declared = cases.declared_cases("eos.antoine_vapor_pressure")
    assert cases.strictest_tolerance(declared) == 1e-12
