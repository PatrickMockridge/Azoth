//! Spec-driven tests for the `eos.uniquac_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{self, UniquacParameters};
use azoth_eos::{model_gen, uniquac_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.uniquac_activity_coefficients";

/// Reshape a flattened `N x N` case value to the nested form the model takes.
fn reshape(flat: &[f64], n: usize) -> Vec<Vec<f64>> {
    flat.chunks(n).map(<[f64]>::to_vec).collect()
}

/// The methanol/water parameters the databank resolves, which are also the doc-test's
/// and the refusal tests' stand-ins.
fn methanol_water() -> UniquacParameters {
    UniquacParameters {
        r: vec![1.4311, 0.92],
        q: vec![1.432, 1.4],
    }
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::UniquacActivityCoefficientsResult {
    let x = case.vector("x").expect("the case declares x");
    let n = x.len();
    let params = databank::uniquac_parameters(case.list("components").expect("components"))
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    uniquac_activity_coefficients(
        &params,
        common::input(case, "T"),
        x,
        &reshape(case.matrix("aij").expect("the case declares aij"), n),
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
    let r = uniquac_activity_coefficients(
        &methanol_water(),
        298.15,
        &[0.5, 0.5],
        &[vec![0.0, -71.0], vec![209.0, 0.0]],
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The resolution itself: `r`/`q` are the UNIFAC group sums, and the case's recorded
/// values are exactly the ones the databank produces - which is the migration's whole
/// claim about this model.
#[test]
fn the_databank_resolves_the_methanol_water_group_sums() {
    assert_eq!(
        databank::uniquac_parameters(&["methanol", "water"]).unwrap(),
        methanol_water()
    );
}

/// Methanol is the CH3OH group alone (`R = 1.4311`, `Q = 1.432`) and water the H2O group
/// alone (`R = 0.92`, `Q = 1.4`), so each component's `r`/`q` is its single group's.
/// A component of several groups sums them, which is the part a one-group test misses.
#[test]
fn a_component_of_several_groups_sums_them() {
    let acetone = databank::uniquac_parameters(&["acetone"]).unwrap();
    let groups = databank::unifac_parameters(&["acetone"]).unwrap();
    let expected_r: f64 = groups
        .groups
        .iter()
        .zip(&groups.group_r)
        .map(|(&count, &r)| count * r)
        .sum();
    let expected_q: f64 = groups
        .groups
        .iter()
        .zip(&groups.group_q)
        .map(|(&count, &q)| count * q)
        .sum();
    assert!(groups.groups.iter().filter(|&&count| count > 0.0).count() > 1);
    assert_eq!(acetone.r, vec![expected_r]);
    assert_eq!(acetone.q, vec![expected_q]);
}

/// The source is the group decomposition, not NeqSim's `rUNIQUAQ`/`qUNIQUAQ` columns.
/// Those are zero for this component, and a resolver that read them would give `r = 0`
/// and divide by zero downstream.
#[test]
fn the_source_is_the_group_sums_and_not_the_zero_uniquac_columns() {
    let acetone = databank::uniquac_parameters(&["acetone"]).unwrap();
    assert_ne!(acetone.r, vec![0.0]);
    assert_ne!(acetone.q, vec![0.0]);
}

/// A name with no group decomposition is refused rather than given zero parameters.
#[test]
fn a_name_without_a_group_decomposition_is_refused() {
    let err = databank::uniquac_parameters(&["methanol", "unobtainium"]).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

#[test]
fn a_nonzero_diagonal_is_refused() {
    let err = uniquac_activity_coefficients(
        &methanol_water(),
        298.15,
        &[0.5, 0.5],
        &[vec![1.0, -71.0], vec![209.0, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("aij"));
}

#[test]
fn a_length_mismatch_is_refused() {
    let err = uniquac_activity_coefficients(
        &UniquacParameters {
            r: vec![1.4311],
            q: vec![1.432, 1.4],
        },
        298.15,
        &[0.5, 0.5],
        &[vec![0.0, -71.0], vec![209.0, 0.0]],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}
