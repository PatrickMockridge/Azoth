//! Spec-driven tests for the `eos.nrtl_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{self, NrtlParameters};
use azoth_eos::{model_gen, nrtl_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.nrtl_activity_coefficients";

/// The methanol/water matrices the databank carries, which are also the 350 K
/// doc-test's and the refusal tests' stand-ins.
fn methanol_water() -> NrtlParameters {
    NrtlParameters {
        alpha: vec![0.0, 0.303, 0.303, 0.0],
        dij: vec![0.0, -48.68, 610.6, 0.0],
    }
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::NrtlActivityCoefficientsResult {
    let params = databank::nrtl_parameters(case.list("components").expect("components"), None)
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    nrtl_activity_coefficients(
        &params,
        common::input(case, "T"),
        case.vector("x").expect("x"),
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
    let r = nrtl_activity_coefficients(&methanol_water(), 350.0, &[0.5, 0.5]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The resolution itself, which is what the migration moved into the databank: if this
/// returns the wrong numbers every case above still passes against a case file that was
/// regenerated from the same mistake.
#[test]
fn the_databank_resolves_the_methanol_water_matrices() {
    assert_eq!(
        databank::nrtl_parameters(&["methanol", "water"], None).unwrap(),
        methanol_water()
    );
}

/// `alpha` is symmetric and `dij` is not, and both come from one `INTER.csv` row. A
/// resolver that read only the first ordering would make NRTL a symmetric model, which
/// it is not, and this is the test that says so.
#[test]
fn alpha_is_symmetric_and_dij_is_directional() {
    // Row-major, so [1] is (0, 1) and [2] is (1, 0).
    let p = databank::nrtl_parameters(&["methanol", "water"], None).unwrap();
    assert_eq!(p.alpha[1], p.alpha[2]);
    assert_ne!(p.dij[1], p.dij[2]);
    // The diagonal is zero by construction, not by a check: nothing writes it.
    assert_eq!(p.alpha[0], 0.0);
    assert_eq!(p.dij[0], 0.0);
    assert_eq!(p.alpha[3], 0.0);
    assert_eq!(p.dij[3], 0.0);
}

/// A pair the interaction table does not carry is an ideal interaction, which is
/// NeqSim's behaviour and a quiet answer the spec states. A name in no table at all is
/// a different thing, and is refused.
#[test]
fn a_name_outside_the_databank_is_refused() {
    let err = databank::nrtl_parameters(&["methanol", "unobtainium"], None).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = nrtl_activity_coefficients(&methanol_water(), 350.0, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

/// The matrices are sized for the mixture the caller resolved them for. Handing them a
/// composition of a different length is the one shape mistake left reachable now that
/// the databank builds them.
#[test]
fn a_matrix_that_does_not_match_the_composition_is_refused() {
    let err = nrtl_activity_coefficients(&methanol_water(), 350.0, &[0.5, 0.25, 0.25]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}
