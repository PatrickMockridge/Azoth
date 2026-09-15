//! Spec-driven tests for the `eos.uniquac_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::{model_gen, uniquac_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.uniquac_activity_coefficients";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::UniquacActivityCoefficientsResult {
    let x = case.vector("x").expect("the case declares x");
    let n = x.len();
    uniquac_activity_coefficients(
        common::input(case, "T"),
        x,
        case.vector("r").expect("the case declares r"),
        case.vector("q").expect("the case declares q"),
        &reshape(case.matrix("aij").expect("the case declares aij"), n),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

fn reshape(flat: &[f64], n: usize) -> Vec<Vec<f64>> {
    flat.chunks(n).map(<[f64]>::to_vec).collect()
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
    let r = uniquac_activity_coefficients(
        298.15,
        &[0.5, 0.5],
        &[1.4311, 0.92],
        &[1.432, 1.4],
        &[vec![0.0, -71.0], vec![209.0, 0.0]],
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_nonzero_diagonal_is_refused() {
    let err = uniquac_activity_coefficients(
        298.15,
        &[0.5, 0.5],
        &[1.4311, 0.92],
        &[1.432, 1.4],
        &[vec![1.0, -71.0], vec![209.0, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("aij"));
}

#[test]
fn a_length_mismatch_is_refused() {
    let err = uniquac_activity_coefficients(
        298.15,
        &[0.5, 0.5],
        &[1.4311],
        &[1.432, 1.4],
        &[vec![0.0, -71.0], vec![209.0, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("r"));
}
