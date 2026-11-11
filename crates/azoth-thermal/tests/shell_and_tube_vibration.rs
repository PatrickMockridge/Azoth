//! `thermal.shell_and_tube_vibration`, against `validation/neqsim/captures/tube_vibration_probe.tsv`.
//!
//! The capture drives `VibrationAnalysis.performScreening` on a passing state and one state
//! for each of the three critical mechanisms, with the four helper intermediates printed.
//! The cases hold the numbers and the four flags; this file also holds the strict refusals.

use azoth_test_support as common;
use azoth_thermal::shell_and_tube_vibration::{VibrationRequest, rate};
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.shell_and_tube_vibration";

fn request(case: &azoth_core::spec::TestCase) -> VibrationRequest {
    VibrationRequest {
        tube_outside_diameter: common::input(case, "tube_outside_diameter"),
        tube_inside_diameter: common::input(case, "tube_inside_diameter"),
        unsupported_span: common::input(case, "unsupported_span"),
        tube_pitch: common::input(case, "tube_pitch"),
        tube_material_e: common::input(case, "tube_material_e"),
        tube_density: common::input(case, "tube_density"),
        crossflow_velocity: common::input(case, "crossflow_velocity"),
        shell_fluid_density: common::input(case, "shell_fluid_density"),
        tube_fluid_density: common::input(case, "tube_fluid_density"),
        shell_inside_diameter: common::input(case, "shell_inside_diameter"),
        sonic_velocity: common::input(case, "sonic_velocity"),
        damping_ratio: common::input(case, "damping_ratio"),
        triangular_pitch: case.flag("triangular_pitch").unwrap_or(false),
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
                "natural_frequency",
                out.natural_frequency.value,
                case.tolerance,
            ),
            (
                "vortex_shedding_frequency",
                out.vortex_shedding_frequency.value,
                case.tolerance,
            ),
            (
                "vortex_shedding_ratio",
                out.vortex_shedding_ratio,
                case.tolerance,
            ),
            (
                "critical_velocity",
                out.critical_velocity.value,
                case.tolerance,
            ),
            ("velocity_ratio", out.velocity_ratio, case.tolerance),
            (
                "acoustic_frequency",
                out.acoustic_frequency.value,
                case.tolerance,
            ),
            ("acoustic_ratio", out.acoustic_ratio, case.tolerance),
        ] {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(value, expected, tolerance, &format!("{context}::{field}"));
            }
        }
        for (field, value) in [
            ("vortex_shedding_critical", out.vortex_shedding_critical),
            ("fluid_elastic_critical", out.fluid_elastic_critical),
            ("acoustic_critical", out.acoustic_critical),
            ("passed", out.passed),
        ] {
            if let Some(expected) = case.flag(field) {
                assert_eq!(value, expected, "{context}::{field}");
            }
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_pass_state_is_the_one_in_the_capture() {
    let out = rate(&VibrationRequest {
        tube_outside_diameter: 0.01905,
        tube_inside_diameter: 0.01483,
        unsupported_span: 0.4,
        tube_pitch: 0.02381,
        tube_material_e: 2.0e11,
        tube_density: 7800.0,
        crossflow_velocity: 0.5,
        shell_fluid_density: 50.0,
        tube_fluid_density: 800.0,
        shell_inside_diameter: 0.5,
        sonic_velocity: 340.0,
        damping_ratio: 0.03,
        triangular_pitch: true,
    })
    .expect("the pass state computes");
    assert!((out.natural_frequency.value - 276.920_465_257_460_3).abs() < 1e-9);
    assert!((out.vortex_shedding_frequency.value - 8.341_372_730_975_722).abs() < 1e-9);
    assert!((out.critical_velocity.value - 56.897_887_868_483_934).abs() < 1e-9);
    assert!(out.passed);
    assert!(!out.vortex_shedding_critical);
    assert!(!out.fluid_elastic_critical);
    assert!(!out.acoustic_critical);
}

#[test]
fn the_layout_and_damping_guards_are_named() {
    let base = VibrationRequest {
        tube_outside_diameter: 0.01905,
        tube_inside_diameter: 0.01483,
        unsupported_span: 0.4,
        tube_pitch: 0.02381,
        tube_material_e: 2.0e11,
        tube_density: 7800.0,
        crossflow_velocity: 0.5,
        shell_fluid_density: 50.0,
        tube_fluid_density: 800.0,
        shell_inside_diameter: 0.5,
        sonic_velocity: 340.0,
        damping_ratio: 0.03,
        triangular_pitch: true,
    };
    let inside_above_outside = VibrationRequest {
        tube_inside_diameter: 0.02,
        ..base.clone()
    };
    assert!(rate(&inside_above_outside).is_err());

    let negative_damping = VibrationRequest {
        damping_ratio: -0.01,
        ..base.clone()
    };
    assert!(rate(&negative_damping).is_err());

    let negative_velocity = VibrationRequest {
        crossflow_velocity: -1.0,
        ..base.clone()
    };
    assert!(rate(&negative_velocity).is_err());

    let zero_outside = VibrationRequest {
        tube_outside_diameter: 0.0,
        ..base
    };
    assert!(rate(&zero_outside).is_err());
}
