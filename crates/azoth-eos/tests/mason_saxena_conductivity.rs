//! Spec-driven tests for the `eos.mason_saxena_conductivity` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::{mason_saxena_conductivity, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.mason_saxena_conductivity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::MasonSaxenaConductivityResult {
    mason_saxena_conductivity(
        case.vector("Cv0").expect("the case declares Cv0"),
        case.vector("M").expect("the case declares M"),
        case.vector("omega").expect("the case declares omega"),
        case.vector("Tc").expect("the case declares Tc"),
        case.vector("Vc").expect("the case declares Vc"),
        case.vector("dipole").expect("the case declares dipole"),
        case.vector("kappa").expect("the case declares kappa"),
        common::input(case, "T"),
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
            result.k.value,
            common::expected(case, "k"),
            case.tolerance,
            &format!("{}::{} (k)", spec.id, case.id),
        );
        common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
    }
}

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let r = mason_saxena_conductivity(
        &[27.544151394, 65.809299386],
        &[0.016043, 0.044097],
        &[0.0115, 0.1523],
        &[190.56, 369.83],
        &[9.9e-5, 0.000203],
        &[0.0, 0.0],
        &[0.0, 0.0],
        300.0,
        &[0.5, 0.5],
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = mason_saxena_conductivity(
        &[27.544151394, 65.809299386],
        &[0.016043, 0.044097],
        &[0.0115, 0.1523],
        &[190.56, 369.83],
        &[9.9e-5, 0.000203],
        &[0.0, 0.0],
        &[0.0, 0.0],
        300.0,
        &[0.6, 0.6],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("z"));
}

#[test]
fn a_length_mismatch_is_refused() {
    let err = mason_saxena_conductivity(
        &[27.544151394],
        &[0.016043, 0.044097],
        &[0.0115, 0.1523],
        &[190.56, 369.83],
        &[9.9e-5, 0.000203],
        &[0.0, 0.0],
        &[0.0, 0.0],
        300.0,
        &[0.5, 0.5],
    )
    .unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("Cv0"));
}
