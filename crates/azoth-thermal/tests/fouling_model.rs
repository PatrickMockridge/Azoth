//! `thermal.fouling_model`, against `validation/neqsim/captures/fouling_model_probe.tsv`.
//!
//! The capture drives `FoulingModel` on its three model types, each row stating the operating
//! conditions the port always applies: the fixed resistance, the Kern-Seaton asymptote, and
//! the Ebert-Panchal threshold model across flow, wall temperature, removal rate and the
//! stateful Euler path.

use azoth_test_support as common;
use azoth_thermal::fouling_model::{FoulingRequest, conditions, rate};
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.fouling_model";

fn request(case: &azoth_core::spec::TestCase) -> FoulingRequest {
    FoulingRequest {
        model: common::input_str(case, "model")
            .parse()
            .expect("a declared model"),
        velocity: common::input(case, "velocity"),
        fluid_density: common::input(case, "fluid_density"),
        fluid_viscosity: common::input(case, "fluid_viscosity"),
        wall_temperature: common::input(case, "wall_temperature"),
        tube_inner_diameter: common::input(case, "tube_inner_diameter"),
        alpha: common::input(case, "alpha"),
        beta: common::input(case, "beta"),
        activation_energy: common::input(case, "activation_energy"),
        gamma: common::input(case, "gamma"),
        fixed_fouling_resistance: common::input(case, "fixed_fouling_resistance"),
        rf_max: common::input(case, "rf_max"),
        time_constant: common::input(case, "time_constant") * 3600.0,
        time_hours: common::input(case, "time_hours") * 3600.0,
        target_fouling_resistance: common::input(case, "target_fouling_resistance"),
        threshold_rate: common::input(case, "threshold_rate"),
        euler_steps: case.input("euler_steps"),
        euler_step_hours: case.input("euler_step_hours").map(|value| value * 3600.0),
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
        for (field, value) in [
            ("fouling_resistance_at_time", out.fouling_resistance_at_time),
            ("current_fouling_resistance", out.current_fouling_resistance),
            (
                "asymptotic_fouling_resistance",
                out.asymptotic_fouling_resistance,
            ),
            (
                "threshold_wall_temperature",
                out.threshold_wall_temperature.value,
            ),
            ("threshold_velocity", out.threshold_velocity.value),
            ("time_to_target", out.time_to_target.value / 3600.0),
            ("operating_time", out.operating_time.value / 3600.0),
            ("fouling_rate", out.fouling_rate),
            ("reynolds_number", out.reynolds_number),
            ("wall_shear_stress", out.wall_shear_stress.value),
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
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_conditions_are_the_stated_ones() {
    let (reynolds, shear) = conditions(1.5, 800.0, 5e-4, 0.01905);
    common::assert_close(reynolds, 45720.0, 1e-12, "reynolds");
    common::assert_close(shear, 4.76675238185679, 1e-9, "shear");
    assert_eq!(conditions(1.5, 800.0, 0.0, 0.01905), (0.0, 0.0));
}

#[test]
fn a_non_integral_or_negative_step_count_is_refused() {
    let base = FoulingRequest {
        model: "ebert_panchal".parse().expect("a declared model"),
        velocity: 1.5,
        fluid_density: 800.0,
        fluid_viscosity: 5e-4,
        wall_temperature: 473.15,
        tube_inner_diameter: 0.01905,
        alpha: 0.5,
        beta: -0.66,
        activation_energy: 48000.0,
        gamma: 5e-9,
        fixed_fouling_resistance: 0.000176,
        rf_max: 0.0005,
        time_constant: 2000.0,
        time_hours: 1000.0,
        target_fouling_resistance: 1e-4,
        threshold_rate: 1e-10,
        euler_steps: Some(2.5),
        euler_step_hours: Some(1.0),
    };
    assert!(rate(&base).is_err());

    let negative = FoulingRequest {
        euler_steps: Some(-1.0),
        ..base.clone()
    };
    assert!(rate(&negative).is_err());

    let zero_velocity = FoulingRequest {
        velocity: 0.0,
        euler_steps: None,
        ..base
    };
    assert!(rate(&zero_velocity).is_err());
}
