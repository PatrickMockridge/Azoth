//! `thermal.shell_and_tube_rating`, against `validation/neqsim/captures/thermal_rating_probe.tsv`.
//!
//! The capture drives `ThermalDesignCalculator` directly on its own defaults, the fouling test's
//! geometry, both pitches, both shell-side methods, all three tube-side regimes, the fouling pair,
//! the deposit layer, a wall viscosity and a stated bypass. The cases hold the numbers; this file
//! also holds the strict refusals the class's own tests assert.

use azoth_core::units::{
    joules_per_kilogram_kelvin, kilograms_per_cubic_meter, kilograms_per_second, meters,
    pascal_seconds, watts_per_meter_kelvin,
};
use azoth_test_support as common;
use azoth_thermal::spec_gen;
use azoth_thermal::{RatingRequest, ShellSideMethod, rate, shell_and_tube_rating};

const CALC_ID: &str = "thermal.shell_and_tube_rating";

fn request(case: &azoth_core::spec::TestCase) -> RatingRequest {
    RatingRequest {
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
        tube_density: common::input(case, "tube_density"),
        tube_viscosity: common::input(case, "tube_viscosity"),
        tube_cp: common::input(case, "tube_cp"),
        tube_conductivity: common::input(case, "tube_conductivity"),
        tube_mass_flow: common::input(case, "tube_mass_flow"),
        shell_density: common::input(case, "shell_density"),
        shell_viscosity: common::input(case, "shell_viscosity"),
        shell_cp: common::input(case, "shell_cp"),
        shell_conductivity: common::input(case, "shell_conductivity"),
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
        let out = rate(&request(case))
            .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id));
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, value, tolerance) in [
            (
                "tube_side_heat_transfer_coefficient",
                out.tube_side_heat_transfer_coefficient.value,
                case.tolerance,
            ),
            (
                "shell_side_heat_transfer_coefficient",
                out.shell_side_heat_transfer_coefficient.value,
                case.tolerance,
            ),
            (
                "overall_heat_transfer_coefficient",
                out.overall_heat_transfer_coefficient.value,
                case.tolerance,
            ),
            (
                "tube_side_pressure_drop",
                out.tube_side_pressure_drop.value,
                case.tolerance,
            ),
            (
                "shell_side_pressure_drop",
                out.shell_side_pressure_drop.value,
                case.tolerance,
            ),
            (
                "tube_side_velocity",
                out.tube_side_velocity.value,
                case.tolerance,
            ),
            (
                "shell_side_velocity",
                out.shell_side_velocity.value,
                case.tolerance,
            ),
            ("tube_side_reynolds", out.tube_side_reynolds, case.tolerance),
            (
                "shell_side_reynolds",
                out.shell_side_reynolds,
                case.tolerance,
            ),
            (
                "effective_tube_inside_diameter",
                out.effective_tube_inside_diameter.value,
                case.tolerance,
            ),
            (
                "outside_heat_transfer_area",
                out.outside_heat_transfer_area.value,
                case.tolerance,
            ),
            (
                "tube_fouling_layer_resistance_outside",
                out.tube_fouling_layer_resistance_outside,
                case.tolerance,
            ),
            ("jc", out.jc, case.tolerance),
            ("jl", out.jl, case.tolerance),
            ("jb", out.jb, case.tolerance),
            ("js", out.js, case.tolerance),
            ("jr", out.jr, case.tolerance),
        ] {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(value, expected, tolerance, &format!("{context}::{field}"));
            }
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_defaults_are_the_class_constructor() {
    let out = rate(&RatingRequest::default()).expect("the defaults compute");
    assert!((out.tube_side_heat_transfer_coefficient.value - 2353.2086219276475).abs() < 1e-9);
    assert!((out.shell_side_heat_transfer_coefficient.value - 505.55308270672384).abs() < 1e-9);
    assert!((out.overall_heat_transfer_coefficient.value - 336.4897827229193).abs() < 1e-9);
}

#[test]
fn the_flat_entry_point_is_the_request() {
    let request = RatingRequest::default();
    let from_request = rate(&request).expect("the request computes");
    let from_flat = shell_and_tube_rating(
        meters(request.tube_outside_diameter),
        meters(request.tube_inside_diameter),
        meters(request.tube_length),
        request.tube_count,
        request.tube_passes,
        meters(request.tube_pitch),
        request.triangular_pitch,
        meters(request.shell_inside_diameter),
        meters(request.baffle_spacing),
        request.baffle_count,
        request.baffle_cut,
        watts_per_meter_kelvin(request.tube_wall_conductivity),
        kilograms_per_cubic_meter(request.tube_density),
        pascal_seconds(request.tube_viscosity),
        joules_per_kilogram_kelvin(request.tube_cp),
        watts_per_meter_kelvin(request.tube_conductivity),
        kilograms_per_second(request.tube_mass_flow),
        kilograms_per_cubic_meter(request.shell_density),
        pascal_seconds(request.shell_viscosity),
        joules_per_kilogram_kelvin(request.shell_cp),
        watts_per_meter_kelvin(request.shell_conductivity),
        kilograms_per_second(request.shell_mass_flow),
        request.fouling_tube,
        request.fouling_shell,
        ShellSideMethod::Kern,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
        None,
    )
    .expect("the flat call computes");
    assert!(
        (from_flat.overall_heat_transfer_coefficient.value
            - from_request.overall_heat_transfer_coefficient.value)
            .abs()
            < 1e-12
    );
}

#[test]
fn strict_refusals_are_named() {
    let zero_tubes = RatingRequest {
        tube_count: 0.0,
        ..Default::default()
    };
    assert!(rate(&zero_tubes).is_err());

    let bad_pitch = RatingRequest {
        tube_pitch: RatingRequest::default().tube_outside_diameter,
        ..Default::default()
    };
    assert!(rate(&bad_pitch).is_err());

    let bad_cut = RatingRequest {
        baffle_cut: 0.5,
        ..Default::default()
    };
    assert!(rate(&bad_cut).is_err());

    let blocked_bore = RatingRequest {
        tube_fouling_layer_thickness: Some(0.02),
        ..Default::default()
    };
    assert!(rate(&blocked_bore).is_err());

    let zero_wall = RatingRequest {
        shell_viscosity_wall: Some(0.0),
        ..Default::default()
    };
    assert!(rate(&zero_wall).is_err());
}
