//! Spec-driven tests for `eos.erosion_prediction`.

use azoth_core::units::{
    kilograms_per_cubic_meter, kilograms_per_square_meter_second, meters, meters_per_second,
    pascal_seconds, seconds,
};
use azoth_eos::erosion_prediction::erosion_prediction;
use azoth_eos::results::ErosionPredictionResult;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.erosion_prediction";

fn call(case: &azoth_core::spec::TestCase) -> ErosionPredictionResult {
    let mass = |name: &str| kilograms_per_cubic_meter(common::input(case, name));
    erosion_prediction(
        mass("mixture_density"),
        meters_per_second(common::input(case, "mixture_velocity")),
        meters(common::input(case, "pipe_diameter")),
        meters(common::input(case, "wall_thickness")),
        kilograms_per_square_meter_second(common::input(case, "sand_rate")),
        meters(common::input(case, "sand_particle_diameter")),
        mass("sand_particle_density"),
        pascal_seconds(common::input(case, "mixture_viscosity")),
        common::input_str(case, "geometry"),
        common::input(case, "bend_radius_ratio"),
        common::input_str(case, "pipe_material"),
        common::input(case, "api_c_factor"),
        meters(common::input(case, "corrosion_allowance")),
        seconds(common::input(case, "design_life")),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
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
        match case.kind {
            "worked_example" | "reference" => {
                executed += 1;
                let result = call(case);
                let context = &format!("{}::{}", spec.id, case.id);
                if let Some(expected) = case.expected_string("risk_level") {
                    assert_eq!(result.risk_level.as_str(), expected, "{context}: the risk");
                }
                for (field, value) in [
                    ("erosional_velocity", result.erosional_velocity.value),
                    ("velocity_ratio", result.velocity_ratio),
                    ("erosion_rate", result.erosion_rate.value),
                    ("cumulative_erosion", result.cumulative_erosion.value),
                    (
                        "remaining_wall_thickness",
                        result.remaining_wall_thickness.value,
                    ),
                ] {
                    if let Some(expected) = case.expected_value(field) {
                        common::assert_close(
                            value,
                            expected,
                            case.tolerance,
                            &format!("{context} ({field})"),
                        );
                    }
                }
                // The two flags are self-asserted: `TestCase` carries no boolean expectations.
                assert_eq!(
                    result.within_api_limits,
                    case.id != "api_v20",
                    "{context}: API limit"
                );
                assert_eq!(
                    result.within_erosion_limits,
                    !matches!(case.id, "eroded" | "rate_high" | "rate_critical"),
                    "{context}: erosion limit"
                );
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "mixture_density" => Some(common::input(case, "mixture_density")),
                        "mixture_velocity" => Some(common::input(case, "mixture_velocity")),
                        "pipe_diameter" => Some(common::input(case, "pipe_diameter")),
                        "wall_thickness" => Some(common::input(case, "wall_thickness")),
                        "sand_rate" => Some(common::input(case, "sand_rate")),
                        "sand_particle_diameter" => {
                            Some(common::input(case, "sand_particle_diameter"))
                        }
                        "sand_particle_density" => {
                            Some(common::input(case, "sand_particle_density"))
                        }
                        "mixture_viscosity" => Some(common::input(case, "mixture_viscosity")),
                        "bend_radius_ratio" => Some(common::input(case, "bend_radius_ratio")),
                        "api_c_factor" => Some(common::input(case, "api_c_factor")),
                        "corrosion_allowance" => Some(common::input(case, "corrosion_allowance")),
                        "design_life" => Some(common::input(case, "design_life")),
                        _ => None,
                    },
                    context,
                );
                common::assert_consistent(&result, context);
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{} has an unhandled kind `{other}`", spec.id, case.id),
        }
    }
    assert!(executed >= 2, "the spec should exercise more than one row");
}

/// **An unknown geometry is refused** rather than defaulted to the class's fall-through.
#[test]
fn an_unknown_geometry_is_refused() {
    let mass = kilograms_per_cubic_meter;
    let error = erosion_prediction(
        mass(100.0),
        meters_per_second(10.0),
        meters(0.1524),
        meters(0.01),
        kilograms_per_square_meter_second(0.0),
        meters(0.00025),
        mass(2650.0),
        pascal_seconds(0.001),
        "manifold",
        1.5,
        "carbon_steel",
        100.0,
        meters(0.003),
        seconds(788_940_000.0),
    )
    .expect_err("the geometry is not declared");
    assert!(error.to_string().contains("geometry"));
}
