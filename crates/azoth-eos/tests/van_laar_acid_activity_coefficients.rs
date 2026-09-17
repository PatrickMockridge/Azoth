//! Spec-driven tests for the `eos.van_laar_acid_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{self, VanLaarAcidParameters};
use azoth_eos::{model_gen, van_laar_acid_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.van_laar_acid_activity_coefficients";

/// The ternary the databank resolves for the first case, which is also the doc-test's
/// stand-in.
fn ternary() -> VanLaarAcidParameters {
    VanLaarAcidParameters {
        acid_index: vec![1, 2, 3],
    }
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::VanLaarAcidActivityCoefficientsResult {
    let params =
        databank::van_laar_acid_parameters(case.list("components").expect("components"), None)
            .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    van_laar_acid_activity_coefficients(
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
    let r = van_laar_acid_activity_coefficients(&ternary(), 250.0, &[0.5, 0.3, 0.2]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The acid identities, both spellings of each - the formulae are what a caller writing
/// `hno3` means, and NeqSim accepts them.
#[test]
fn the_databank_resolves_each_spelling_to_the_same_acid() {
    let named =
        databank::van_laar_acid_parameters(&["water", "nitric acid", "sulfuric acid"], None)
            .unwrap();
    let formulae = databank::van_laar_acid_parameters(&["water", "hno3", "h2so4"], None).unwrap();
    assert_eq!(named.acid_index, vec![1, 2, 3]);
    assert_eq!(formulae.acid_index, vec![1, 2, 3]);
}

/// A component the model does not cover is a species it does not model, not an error:
/// the databank resolves it and the model gives it the penalty.
#[test]
fn a_component_outside_the_three_gets_index_zero() {
    let params = databank::van_laar_acid_parameters(&["water", "hno3", "nitrogen"], None).unwrap();
    assert_eq!(params.acid_index, vec![1, 2, 0]);
}

/// A name the *databank* does not carry is a different thing from one it does and the
/// model does not, and it is refused.
#[test]
fn a_name_outside_the_databank_is_refused() {
    let err = databank::van_laar_acid_parameters(&["water", "unobtainium"], None).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

/// A species outside the model does not enter the acids' answer, and pays the penalty
/// itself: the two acids of a water/nitric/nitrogen mixture at `0.5/0.3/0.2` read as a
/// water/nitric mixture at `0.625/0.375`, which is their basis over the acids alone.
///
/// **This is not a test of the renormalisation.** Measured, dropping it changes nothing:
/// each of the three expressions is homogeneous of degree zero, so it returns the same
/// value whether the basis is `0.5/0.3` or `0.625/0.375` or any scaling of either. What
/// this pins is that the carrier gas stays out of the ternary - which is a property of
/// the basis being built only from indices 1, 2 and 3, not of the division.
#[test]
fn a_carrier_gas_stays_out_of_the_acids_answer() {
    let with_gas =
        databank::van_laar_acid_parameters(&["water", "hno3", "nitrogen"], None).unwrap();
    let without = databank::van_laar_acid_parameters(&["water", "hno3"], None).unwrap();

    let diluted = van_laar_acid_activity_coefficients(&with_gas, 250.0, &[0.5, 0.3, 0.2]).unwrap();
    let neat = van_laar_acid_activity_coefficients(&without, 250.0, &[0.625, 0.375]).unwrap();

    // Not `assert_eq!`: 0.5/0.8 and 0.625 are the same number to a person and not to a
    // double, and the two bases reach the same gamma by different roundings.
    common::assert_close(diluted.gamma[0], neat.gamma[0], 1e-12, "diluted vs neat");
    common::assert_close(diluted.gamma[1], neat.gamma[1], 1e-12, "diluted vs neat");
    // The third component is the carrier gas, and pays the penalty.
    assert_eq!(diluted.gamma[2], 1.0e12);
}

/// With no acid at all the ternary basis is undefined, and NeqSim reads that as pure
/// water rather than refusing.
#[test]
fn a_mixture_with_no_acid_reads_as_pure_water() {
    let params = databank::van_laar_acid_parameters(&["nitrogen", "methane"], None).unwrap();
    let r = van_laar_acid_activity_coefficients(&params, 250.0, &[0.5, 0.5]).unwrap();
    assert_eq!(r.gamma, vec![1.0e12, 1.0e12]);
}

/// The Taleb fit covers 190-298 K; outside it the model extrapolates, and that is
/// reported rather than refused.
#[test]
fn a_temperature_outside_the_fitted_range_warns() {
    let r = van_laar_acid_activity_coefficients(&ternary(), 320.0, &[0.5, 0.3, 0.2]).unwrap();
    assert!(!r.is_clean(), "an extrapolation should carry a warning");
    let inside = van_laar_acid_activity_coefficients(&ternary(), 250.0, &[0.5, 0.3, 0.2]).unwrap();
    assert!(inside.is_clean(), "{:?}", inside.warnings);
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let err = van_laar_acid_activity_coefficients(&ternary(), 250.0, &[0.6, 0.6, 0.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

#[test]
fn a_length_mismatch_is_refused() {
    let err = van_laar_acid_activity_coefficients(&ternary(), 250.0, &[1.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}
