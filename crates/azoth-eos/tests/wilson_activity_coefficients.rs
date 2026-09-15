//! Spec-driven tests for the `eos.wilson_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::{model_gen, wilson_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.wilson_activity_coefficients";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::WilsonActivityCoefficientsResult {
    wilson_activity_coefficients(
        common::input(case, "T"),
        case.vector("x").expect("the case declares x"),
        case.vector("M").expect("the case declares M"),
        case.vector("Tc").expect("the case declares Tc"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for (name, actual) in [("ln_gamma", &result.ln_gamma), ("gamma", &result.gamma)] {
            let expected = case
                .expected_vector(name)
                .unwrap_or_else(|| panic!("the case declares {name}"));
            for (i, (&got, &want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    got,
                    want,
                    case.tolerance,
                    &format!("{context} ({name}[{i}])"),
                );
            }
        }
        common::assert_consistent(&result, context);
    }
}

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let r =
        wilson_activity_coefficients(298.15, &[0.5, 0.5], &[0.058123, 0.1703], &[425.12, 658.0])
            .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_length_mismatch_is_refused() {
    let err = wilson_activity_coefficients(298.15, &[0.5, 0.5], &[0.058123], &[425.12, 658.0])
        .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("M"));
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err =
        wilson_activity_coefficients(298.15, &[0.6, 0.6], &[0.058123, 0.1703], &[425.12, 658.0])
            .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
