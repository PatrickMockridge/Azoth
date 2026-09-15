//! Spec-driven tests for the `eos.unifac_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::{model_gen, unifac_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.unifac_activity_coefficients";

/// A flattened matrix case value, reshaped to `rows x cols`.
fn reshape(flat: &[f64], rows: usize, cols: usize) -> Vec<Vec<f64>> {
    flat.chunks(cols).take(rows).map(<[f64]>::to_vec).collect()
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::UnifacActivityCoefficientsResult {
    let x = case.vector("x").expect("the case declares x");
    let group_r = case.vector("group_r").expect("the case declares group_r");
    let n = x.len();
    let g = group_r.len();
    unifac_activity_coefficients(
        common::input(case, "T"),
        x,
        &reshape(
            case.matrix("groups").expect("the case declares groups"),
            n,
            g,
        ),
        group_r,
        case.vector("group_q").expect("the case declares group_q"),
        &reshape(case.matrix("aij").expect("the case declares aij"), g, g),
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
    let r = unifac_activity_coefficients(
        298.15,
        &[0.5, 0.5],
        &[vec![1.0, 0.0], vec![0.0, 1.0]],
        &[1.4311, 0.92],
        &[1.432, 1.4],
        &[vec![0.0, -181.0], vec![289.6, 0.0]],
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_shape_mismatch_is_refused() {
    let err = unifac_activity_coefficients(
        298.15,
        &[0.5, 0.5],
        &[vec![1.0, 0.0, 0.0], vec![0.0, 1.0, 0.0]],
        &[1.4311, 0.92],
        &[1.432, 1.4],
        &[vec![0.0, -181.0], vec![289.6, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("groups"));
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = unifac_activity_coefficients(
        298.15,
        &[0.6, 0.6],
        &[vec![1.0, 0.0], vec![0.0, 1.0]],
        &[1.4311, 0.92],
        &[1.432, 1.4],
        &[vec![0.0, -181.0], vec![289.6, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
