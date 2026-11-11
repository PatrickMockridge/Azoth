//! `thermal.incremental_zone_analysis`, against
//! `validation/neqsim/captures/incremental_zone_analysis_probe.tsv`.
//!
//! The capture drives `IncrementalZoneAnalysis` on one liquid zone, the phase cycle, a
//! `TWO_PHASE` tube regime, the LMTD branches, a single-point quality range and the
//! zone-count baffle scaling. This file also holds the empty list and the zero-flow refusal.

use azoth_core::units::{
    joules_per_kilogram_kelvin, kilograms_per_cubic_meter, kilograms_per_second, pascal_seconds,
    watts, watts_per_meter_kelvin,
};
use azoth_test_support as common;
use azoth_thermal::incremental_zone_analysis::{IncrementalZoneAnalysisRequest, rate};
use azoth_thermal::results::IncrementalZoneInput;
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.incremental_zone_analysis";

fn zones(case: &azoth_core::spec::TestCase) -> Vec<IncrementalZoneInput> {
    common::record_input(case, "zones")
        .iter()
        .map(|row| IncrementalZoneInput {
            zone_name: common::record_name(row, "zone_name")
                .expect("a stated name")
                .to_string(),
            duty: watts(common::record_number(row, "duty").expect("a stated duty")),
            hot_inlet_temp: azoth_core::units::kelvins(
                common::record_number(row, "hot_inlet_temp").expect("a stated temperature"),
            ),
            hot_outlet_temp: azoth_core::units::kelvins(
                common::record_number(row, "hot_outlet_temp").expect("a stated temperature"),
            ),
            cold_inlet_temp: azoth_core::units::kelvins(
                common::record_number(row, "cold_inlet_temp").expect("a stated temperature"),
            ),
            cold_outlet_temp: azoth_core::units::kelvins(
                common::record_number(row, "cold_outlet_temp").expect("a stated temperature"),
            ),
            tube_phase_regime: common::record_name(row, "tube_phase_regime")
                .expect("a stated regime")
                .to_string(),
            tube_density: kilograms_per_cubic_meter(
                common::record_number(row, "tube_density").expect("a stated density"),
            ),
            tube_viscosity: pascal_seconds(
                common::record_number(row, "tube_viscosity").expect("a stated viscosity"),
            ),
            tube_cp: joules_per_kilogram_kelvin(
                common::record_number(row, "tube_cp").expect("a stated heat capacity"),
            ),
            tube_conductivity: watts_per_meter_kelvin(
                common::record_number(row, "tube_conductivity").expect("a stated conductivity"),
            ),
            tube_mass_flow: kilograms_per_second(
                common::record_number(row, "tube_mass_flow").expect("a stated flow"),
            ),
            shell_phase_regime: common::record_name(row, "shell_phase_regime")
                .expect("a stated regime")
                .to_string(),
            shell_density: kilograms_per_cubic_meter(
                common::record_number(row, "shell_density").expect("a stated density"),
            ),
            shell_viscosity: pascal_seconds(
                common::record_number(row, "shell_viscosity").expect("a stated viscosity"),
            ),
            shell_cp: joules_per_kilogram_kelvin(
                common::record_number(row, "shell_cp").expect("a stated heat capacity"),
            ),
            shell_conductivity: watts_per_meter_kelvin(
                common::record_number(row, "shell_conductivity").expect("a stated conductivity"),
            ),
            shell_mass_flow: kilograms_per_second(
                common::record_number(row, "shell_mass_flow").expect("a stated flow"),
            ),
            tube_liquid_density: common::record_number(row, "tube_liquid_density")
                .map(kilograms_per_cubic_meter),
            tube_vapor_density: common::record_number(row, "tube_vapor_density")
                .map(kilograms_per_cubic_meter),
            tube_liquid_viscosity: common::record_number(row, "tube_liquid_viscosity")
                .map(pascal_seconds),
            tube_vapor_viscosity: common::record_number(row, "tube_vapor_viscosity")
                .map(pascal_seconds),
            tube_liquid_cp: common::record_number(row, "tube_liquid_cp")
                .map(joules_per_kilogram_kelvin),
            tube_liquid_conductivity: common::record_number(row, "tube_liquid_conductivity")
                .map(watts_per_meter_kelvin),
            tube_quality_in: common::record_number(row, "tube_quality_in"),
            tube_quality_out: common::record_number(row, "tube_quality_out"),
            tube_reduced_pressure: common::record_number(row, "tube_reduced_pressure"),
            tube_surface_tension: common::record_number(row, "tube_surface_tension")
                .map(azoth_core::units::newtons_per_meter),
            tube_heat_of_vaporization: common::record_number(row, "tube_heat_of_vaporization"),
        })
        .collect()
}

fn request_for(case: &azoth_core::spec::TestCase) -> IncrementalZoneAnalysisRequest {
    IncrementalZoneAnalysisRequest {
        tube_inside_diameter: common::input(case, "tube_inside_diameter"),
        tube_outside_diameter: common::input(case, "tube_outside_diameter"),
        tube_length: common::input(case, "tube_length"),
        tube_count: common::input(case, "tube_count"),
        tube_passes: common::input(case, "tube_passes"),
        shell_inside_diameter: common::input(case, "shell_inside_diameter"),
        baffle_spacing: common::input(case, "baffle_spacing"),
        tube_pitch: common::input(case, "tube_pitch"),
        triangular_pitch: case.flag("triangular_pitch").unwrap_or(true),
        tube_wall_conductivity: common::input(case, "tube_wall_conductivity"),
        fouling_tube: common::input(case, "fouling_tube"),
        fouling_shell: common::input(case, "fouling_shell"),
        zones: zones(case),
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        if case.kind != "worked_example" && case.kind != "reference" {
            continue;
        }
        let out = rate(&request_for(case))
            .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id));
        let context = format!("{}::{}", spec.id, case.id);
        for (field, value) in [
            ("number_of_zones", out.number_of_zones),
            ("total_required_area", out.total_required_area.value),
            ("total_duty", out.total_duty.value),
            ("weighted_average_u", out.weighted_average_u.value),
            (
                "total_tube_side_pressure_drop",
                out.total_tube_side_pressure_drop.value,
            ),
            (
                "total_shell_side_pressure_drop",
                out.total_shell_side_pressure_drop.value,
            ),
            (
                "minimum_approach_temperature",
                out.minimum_approach_temperature.value,
            ),
        ] {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(
                    value,
                    expected,
                    case.tolerance,
                    &format!("{context}::{field}"),
                );
            }
        }
        let Some(expected) = case.expected_record("zones") else {
            executed += 1;
            continue;
        };
        assert_eq!(out.zones.len(), expected.len(), "{context}: the zone count");
        for (index, (row, claims)) in out.zones.iter().zip(expected).enumerate() {
            let at = format!("{context}::zones[{index}]");
            if let Some(expected) = common::record_number(claims, "zone_index") {
                common::assert_close(
                    row.zone_index,
                    expected,
                    case.tolerance,
                    &format!("{at}::zone_index"),
                );
            }
            if let Some(name) = common::record_name(claims, "zone_name") {
                assert_eq!(row.zone_name, name, "{at}::zone_name");
            }
            if let Some(name) = common::record_name(claims, "tube_phase_regime") {
                assert_eq!(row.tube_phase_regime, name, "{at}::tube_phase_regime");
            }
            if let Some(name) = common::record_name(claims, "shell_phase_regime") {
                assert_eq!(row.shell_phase_regime, name, "{at}::shell_phase_regime");
            }
            for (field, value) in [
                (
                    "tube_side_heat_transfer_coefficient",
                    row.tube_side_heat_transfer_coefficient.value,
                ),
                (
                    "shell_side_heat_transfer_coefficient",
                    row.shell_side_heat_transfer_coefficient.value,
                ),
                (
                    "overall_heat_transfer_coefficient",
                    row.overall_heat_transfer_coefficient.value,
                ),
                ("lmtd", row.lmtd.value),
                ("duty", row.duty.value),
                ("required_area", row.required_area.value),
                ("tube_side_pressure_drop", row.tube_side_pressure_drop.value),
                (
                    "shell_side_pressure_drop",
                    row.shell_side_pressure_drop.value,
                ),
                ("hot_inlet_temp", row.hot_inlet_temp.value),
                ("hot_outlet_temp", row.hot_outlet_temp.value),
                ("cold_inlet_temp", row.cold_inlet_temp.value),
                ("cold_outlet_temp", row.cold_outlet_temp.value),
            ] {
                if let Some(expected) = common::record_number(claims, field) {
                    common::assert_close(
                        value,
                        expected,
                        case.tolerance,
                        &format!("{at}::{field}"),
                    );
                }
            }
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn an_empty_zone_list_answers_the_sentinels() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");
    let mut request = request_for(case);
    request.zones.clear();
    let out = rate(&request).expect("an empty list is not an error");
    assert_eq!(out.number_of_zones, 0.0);
    assert_eq!(out.total_required_area.value, 0.0);
    assert_eq!(out.weighted_average_u.value, 0.0);
    assert_eq!(out.minimum_approach_temperature.value, f64::MAX);
    assert!(out.zones.is_empty());
}

#[test]
fn a_zero_flow_is_refused() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");
    let mut request = request_for(case);
    request.zones[0].tube_mass_flow = kilograms_per_second(0.0);
    let error = rate(&request).expect_err("a zero flow is not a coefficient");
    assert!(error.to_string().contains("tube_mass_flow"), "{error}");

    let mut request = request_for(case);
    request.zones[0].shell_mass_flow = kilograms_per_second(0.0);
    assert!(rate(&request).is_err(), "a zero shell flow");
}

#[test]
fn a_non_integral_pass_count_is_refused() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");
    let mut request = request_for(case);
    request.tube_passes = 2.5;
    assert!(rate(&request).is_err());
}
