//! Spec-driven tests for the `eos.nrtl_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::{model_gen, nrtl_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.nrtl_activity_coefficients";

/// A flattened matrix case value, reshaped to the `N x N` nested form the model takes.
fn reshape(flat: &[f64], n: usize) -> Vec<Vec<f64>> {
    flat.chunks(n).map(<[f64]>::to_vec).collect()
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::NrtlActivityCoefficientsResult {
    let x = case.vector("x").expect("the case declares x");
    let n = x.len();
    nrtl_activity_coefficients(
        common::input(case, "T"),
        x,
        &reshape(case.matrix("Dij").expect("the case declares Dij"), n),
        &reshape(case.matrix("alpha").expect("the case declares alpha"), n),
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
    let r = nrtl_activity_coefficients(
        350.0,
        &[0.5, 0.5],
        &[vec![0.0, -48.68], vec![610.6, 0.0]],
        &[vec![0.0, 0.303], vec![0.303, 0.0]],
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = nrtl_activity_coefficients(
        350.0,
        &[0.6, 0.6],
        &[vec![0.0, -48.68], vec![610.6, 0.0]],
        &[vec![0.0, 0.303], vec![0.303, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

#[test]
fn a_non_symmetric_alpha_is_refused() {
    let err = nrtl_activity_coefficients(
        350.0,
        &[0.5, 0.5],
        &[vec![0.0, -48.68], vec![610.6, 0.0]],
        &[vec![0.0, 0.303], vec![0.4, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("alpha"));
}

#[test]
fn a_nonzero_diagonal_is_refused() {
    let err = nrtl_activity_coefficients(
        350.0,
        &[0.5, 0.5],
        &[vec![1.0, -48.68], vec![610.6, 0.0]],
        &[vec![0.0, 0.303], vec![0.303, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("Dij"));
}
