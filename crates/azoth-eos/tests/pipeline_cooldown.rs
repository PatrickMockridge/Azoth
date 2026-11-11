//! Spec-driven tests for `eos.pipeline_cooldown`.

use azoth_core::units::{
    joules_per_kilogram_kelvin, kelvins, kilograms_per_cubic_meter, meters, seconds,
    watts_per_meter_kelvin, watts_per_square_meter_kelvin,
};
use azoth_eos::model_gen;
use azoth_eos::pipeline_cooldown::pipeline_cooldown;
use azoth_eos::results::PipelineCooldownResult;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pipeline_cooldown";

fn call(case: &azoth_core::spec::TestCase) -> PipelineCooldownResult {
    pipeline_cooldown(
        meters(common::input(case, "internal_diameter")),
        meters(common::input(case, "wall_thickness")),
        meters(common::input(case, "insulation_thickness")),
        meters(common::input(case, "coating_thickness")),
        kilograms_per_cubic_meter(common::input(case, "steel_density")),
        joules_per_kilogram_kelvin(common::input(case, "steel_specific_heat")),
        watts_per_meter_kelvin(common::input(case, "insulation_conductivity")),
        kilograms_per_cubic_meter(common::input(case, "fluid_density")),
        joules_per_kilogram_kelvin(common::input(case, "fluid_specific_heat")),
        kelvins(common::input(case, "initial_fluid_temperature")),
        kelvins(common::input(case, "ambient_temperature")),
        seconds(common::input(case, "time_step")),
        seconds(common::input(case, "total_time")),
        case.input("overall_u_value")
            .map(watts_per_square_meter_kelvin),
        case.flag("use_layer_calculation"),
        case.input("external_htc")
            .map(watts_per_square_meter_kelvin),
        case.input("target_temperature").map(kelvins),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn assert_vector(case: &azoth_core::spec::TestCase, context: &str, name: &str, got: &[f64]) {
    if let Some(expected) = case.expected_vector(name) {
        assert_eq!(got.len(), expected.len(), "{context}: {name} length");
        for (index, (value, want)) in got.iter().zip(expected).enumerate() {
            common::assert_close(
                *value,
                *want,
                case.tolerance,
                &format!("{context}: {name}[{index}]"),
            );
        }
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        assert_vector(case, context, "time", &result.time);
        assert_vector(
            case,
            context,
            "fluid_temperature",
            &result.fluid_temperature,
        );
        for (field, value) in [
            ("time_constant", result.time_constant.value),
            ("calculated_u_value", result.calculated_u_value.value),
            ("effective_u_value", result.effective_u_value.value),
            (
                "time_to_reach_target",
                result
                    .time_to_reach_target
                    .map_or(f64::NAN, |time| time.value),
            ),
        ] {
            if let Some(want) = case.expected_value(field) {
                common::assert_close(value, want, case.tolerance, &format!("{context} ({field})"));
            }
        }
        // A target the march never reaches is an absent output, and `TestCase` has no
        // expectation for absence.
        if case.id == "target_unreached" {
            assert!(
                result.time_to_reach_target.is_none(),
                "{context}: the 270 K target is not reached"
            );
        }
        common::assert_consistent(&result, context);
    }
}
