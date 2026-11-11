//! Spec-driven tests for `eos.emulsion_viscosity`.

use azoth_core::units::{kilograms_per_cubic_meter, pascal_seconds};
use azoth_eos::emulsion_viscosity::emulsion_viscosity;
use azoth_eos::results::EmulsionViscosityResult;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.emulsion_viscosity";

fn call(case: &azoth_core::spec::TestCase) -> EmulsionViscosityResult {
    emulsion_viscosity(
        pascal_seconds(common::input(case, "oil_viscosity")),
        pascal_seconds(common::input(case, "water_viscosity")),
        kilograms_per_cubic_meter(common::input(case, "oil_density")),
        kilograms_per_cubic_meter(common::input(case, "water_density")),
        common::input(case, "water_cut"),
        common::input(case, "max_packing_fraction"),
        common::input_str(case, "model"),
        common::input(case, "tightness_factor"),
        case.flag("demulsifier_present").unwrap_or(false),
        common::input(case, "demulsifier_efficiency"),
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
                if let Some(expected) = case.expected_string("emulsion_type") {
                    assert_eq!(
                        result.emulsion_type.as_str(),
                        expected,
                        "{context}: the phase type"
                    );
                }
                for (field, value) in [
                    ("effective_viscosity", result.effective_viscosity.value),
                    ("relative_viscosity", result.relative_viscosity),
                    ("viscosity_ratio", result.viscosity_ratio),
                    ("inversion_water_cut", result.inversion_water_cut),
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
                // The flag is self-asserted: `TestCase` carries no boolean expectations.
                let inverted = matches!(
                    case.id,
                    "woelflin_loose" | "richardson_low" | "inverted" | "low_clamp" | "wc_one"
                );
                assert_eq!(result.is_inverted, inverted, "{context}: inverted");
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "oil_viscosity" => Some(common::input(case, "oil_viscosity")),
                        "water_viscosity" => Some(common::input(case, "water_viscosity")),
                        "oil_density" => Some(common::input(case, "oil_density")),
                        "water_density" => Some(common::input(case, "water_density")),
                        "water_cut" => Some(common::input(case, "water_cut")),
                        "max_packing_fraction" => Some(common::input(case, "max_packing_fraction")),
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

/// **An unknown model is refused** rather than defaulted to Pal-Rhodes.
#[test]
fn an_unknown_model_is_refused() {
    let error = emulsion_viscosity(
        pascal_seconds(0.01),
        pascal_seconds(0.0005),
        kilograms_per_cubic_meter(850.0),
        kilograms_per_cubic_meter(1020.0),
        0.30,
        0.74,
        "bingham",
        0.5,
        false,
        0.0,
    )
    .expect_err("the model is not declared");
    assert!(error.to_string().contains("model"));
}
