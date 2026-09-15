//! Spec-driven tests for the `eos.wilke_viscosity` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::{model_gen, wilke_viscosity};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.wilke_viscosity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::WilkeViscosityResult {
    wilke_viscosity(
        case.vector("Tc").expect("the case declares Tc"),
        case.vector("Vc").expect("the case declares Vc"),
        case.vector("M").expect("the case declares M"),
        case.vector("omega").expect("the case declares omega"),
        case.vector("dipole").expect("the case declares dipole"),
        case.vector("kappa").expect("the case declares kappa"),
        common::input(case, "T"),
        common::input(case, "V"),
        case.vector("z").expect("the case declares z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        common::assert_close(
            result.mu.value,
            common::expected(case, "mu"),
            case.tolerance,
            &format!("{}::{} (mu)", spec.id, case.id),
        );
        common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
    }
}

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let r = wilke_viscosity(
        &[190.56, 369.83],
        &[9.9e-5, 0.000203],
        &[0.016043, 0.044097],
        &[0.0115, 0.1523],
        &[0.0, 0.0],
        &[0.0, 0.0],
        300.0,
        2.2987856717900145e-3,
        &[0.5, 0.5],
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = wilke_viscosity(
        &[190.56, 369.83],
        &[9.9e-5, 0.000203],
        &[0.016043, 0.044097],
        &[0.0115, 0.1523],
        &[0.0, 0.0],
        &[0.0, 0.0],
        300.0,
        2.2987856717900145e-3,
        &[0.6, 0.6],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("z"));
}

#[test]
fn a_length_mismatch_is_refused() {
    let err = wilke_viscosity(
        &[190.56, 369.83],
        &[9.9e-5],
        &[0.016043, 0.044097],
        &[0.0115, 0.1523],
        &[0.0, 0.0],
        &[0.0, 0.0],
        300.0,
        2.2987856717900145e-3,
        &[0.5, 0.5],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("Tc"));
}

#[test]
fn a_negative_molar_volume_is_refused() {
    let err = wilke_viscosity(
        &[190.56, 369.83],
        &[9.9e-5, 0.000203],
        &[0.016043, 0.044097],
        &[0.0115, 0.1523],
        &[0.0, 0.0],
        &[0.0, 0.0],
        300.0,
        -1.0,
        &[0.5, 0.5],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::OutOfRange { .. }), "{err:?}");
    assert_eq!(err.field(), Some("V"));
}
