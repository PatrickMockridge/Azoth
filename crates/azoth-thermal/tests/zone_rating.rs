//! `thermal.zone_rating`, against `validation/neqsim/captures/zone_analysis_probe.tsv`.
//!
//! The capture drives `ThermalDesignCalculator.calculateZones` on its documented example, a
//! three-zone case, Bell-Delaware with sealing, an inner deposit, degenerate zone properties,
//! zero and negative LMTD, and the restore row. The cases hold the numbers; this file also
//! holds the refusals the class's permissive path does not make.

use azoth_core::units::{
    joules_per_kilogram_kelvin, kelvin_intervals, kilograms_per_cubic_meter, pascal_seconds, watts,
    watts_per_meter_kelvin,
};
use azoth_test_support as common;
use azoth_thermal::results::ZoneDefinitionRow;
use azoth_thermal::shell_and_tube_rating::{RatingRequest, ShellSideMethod};
use azoth_thermal::spec_gen;
use azoth_thermal::zone_rating::{ZoneRatingRequest, rate};

const CALC_ID: &str = "thermal.zone_rating";

fn zones(case: &azoth_core::spec::TestCase) -> Vec<ZoneDefinitionRow> {
    common::record_input(case, "zones")
        .iter()
        .map(|row| ZoneDefinitionRow {
            zone_name: common::record_name(row, "zone_name")
                .expect("a stated zone name")
                .to_string(),
            duty_fraction: common::record_number(row, "duty_fraction").expect("a stated fraction"),
            total_duty: watts(common::record_number(row, "total_duty").expect("a stated duty")),
            lmtd: kelvin_intervals(common::record_number(row, "lmtd").expect("a stated LMTD")),
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
        })
        .collect()
}

fn request_for(case: &azoth_core::spec::TestCase) -> ZoneRatingRequest {
    ZoneRatingRequest {
        base: RatingRequest {
            tube_outside_diameter: common::input(case, "tube_outside_diameter"),
            tube_inside_diameter: common::input(case, "tube_inside_diameter"),
            tube_length: common::input(case, "tube_length"),
            tube_count: common::input(case, "tube_count"),
            tube_passes: common::input(case, "tube_passes"),
            tube_pitch: common::input(case, "tube_pitch"),
            triangular_pitch: case.flag("triangular_pitch").unwrap_or(true),
            shell_inside_diameter: common::input(case, "shell_inside_diameter"),
            baffle_spacing: common::input(case, "baffle_spacing"),
            baffle_count: common::input(case, "baffle_count"),
            baffle_cut: common::input(case, "baffle_cut"),
            tube_wall_conductivity: common::input(case, "tube_wall_conductivity"),
            tube_mass_flow: common::input(case, "tube_mass_flow"),
            shell_mass_flow: common::input(case, "shell_mass_flow"),
            fouling_tube: common::input(case, "fouling_tube"),
            fouling_shell: common::input(case, "fouling_shell"),
            shell_side_method: common::input_str(case, "shell_side_method")
                .parse()
                .expect("a declared shell-side method"),
            tube_fouling_layer_thickness: case.input("tube_fouling_layer_thickness"),
            tube_fouling_layer_conductivity: case.input("tube_fouling_layer_conductivity"),
            shell_viscosity_wall: case.input("shell_viscosity_wall"),
            has_sealing: case.flag("has_sealing"),
            sealing_pairs: case.input("sealing_pairs"),
            tube_to_baffle_clearance: case.input("tube_to_baffle_clearance"),
            shell_to_baffle_clearance: case.input("shell_to_baffle_clearance"),
            bypass_area: case.input("bypass_area"),
            // No zone reads the base transport properties; the class defaults keep the
            // request well formed.
            ..RatingRequest::default()
        },
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
        common::assert_close(
            out.zone_count,
            case.expected_value("zone_count")
                .expect("every case states its zone count"),
            case.tolerance,
            &format!("{context}::zone_count"),
        );
        let expected = common::expected_record(case, "zones");
        assert_eq!(out.zones.len(), expected.len(), "{context}: the zone count");
        for (index, (row, claims)) in out.zones.iter().zip(expected).enumerate() {
            let at = format!("{context}::zones[{index}]");
            assert_eq!(
                row.zone_name,
                common::record_name(claims, "zone_name").expect("a stated name")
            );
            common::assert_close(
                row.duty_fraction,
                common::record_number(claims, "duty_fraction").expect("a stated fraction"),
                case.tolerance,
                &format!("{at}::duty_fraction"),
            );
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
                ("required_area", row.required_area.value),
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
fn a_zone_with_a_degenerate_property_is_refused() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");
    let mut request = request_for(case);
    request.zones[0].tube_viscosity = pascal_seconds(0.0);
    let error = rate(&request).expect_err("a zero viscosity is not an answer");
    assert!(
        error.to_string().contains("tube_viscosity"),
        "the refusal should name the field: {error}"
    );

    let mut request = request_for(case);
    request.zones[0].shell_conductivity = watts_per_meter_kelvin(0.0);
    let error = rate(&request).expect_err("a zero conductivity is not an answer");
    assert!(error.to_string().contains("shell_conductivity"), "{error}");
}

#[test]
fn a_base_the_strict_rating_refuses_is_refused_here() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");

    let mut request = request_for(case);
    request.base.tube_passes = 3.0;
    assert!(rate(&request).is_err(), "a non-dividing pass count");

    let mut request = request_for(case);
    request.base.tube_pitch = 0.019;
    assert!(
        rate(&request).is_err(),
        "a pitch not above the outside diameter"
    );

    let mut request = request_for(case);
    request.base.tube_fouling_layer_thickness = Some(0.008);
    request.base.tube_fouling_layer_conductivity = Some(0.2);
    assert!(rate(&request).is_err(), "a deposit that blocks the bore");

    let mut request = request_for(case);
    request.base.shell_side_method = ShellSideMethod::BellDelaware;
    request.base.has_sealing = Some(true);
    request.base.sealing_pairs = Some(0.0);
    assert!(rate(&request).is_err(), "a stated seal with no pair");
}

#[test]
fn an_empty_zone_list_answers_no_zones() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    let case = spec
        .all_tests()
        .find(|case| case.is_active() && case.kind == "worked_example")
        .expect("the worked example");
    let mut request = request_for(case);
    request.zones.clear();
    let out = rate(&request).expect("an empty list is not an error");
    assert_eq!(out.zone_count, 0.0);
    assert!(out.zones.is_empty());
}
