//! Spec-driven tests for the `eos.unifac_psrk_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{UnifacParameters, unifac_parameters, unifac_psrk_parameters};
use azoth_eos::unifac_activity_coefficients::unifac_activity_coefficients;
use azoth_eos::{model_gen, unifac_psrk_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.unifac_psrk_activity_coefficients";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::UnifacPsrkActivityCoefficientsResult {
    let params = unifac_psrk_parameters(case.list("components").expect("components"))
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    unifac_psrk_activity_coefficients(
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
    let params = unifac_psrk_parameters(&["methanol", "water"]).unwrap();
    let r = unifac_psrk_activity_coefficients(&params, 298.15, &[0.5, 0.5]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The whole of the difference between this model and plain UNIFAC is where `aij`
/// comes from, so the two must agree exactly wherever the PSRK tables carry no
/// temperature dependence - and the model is then checked by UNIFAC's own cases, which
/// are externally validated.
#[test]
fn a_pair_with_no_temperature_dependence_is_plain_unifac() {
    let psrk = unifac_psrk_parameters(&["methanol", "water"]).unwrap();
    let plain = unifac_parameters(&["methanol", "water"]).unwrap();
    assert!(
        psrk.bij.iter().all(|&b| b == 0.0) && psrk.cij.iter().all(|&c| c == 0.0),
        "methanol/water is chosen because its b and c are zero"
    );
    // Three temperatures, because the claim is that this holds at all of them.
    for t in [200.0, 298.15, 450.0] {
        let a = unifac_psrk_activity_coefficients(&psrk, t, &[0.5, 0.5]).unwrap();
        let b = unifac_activity_coefficients(&plain, t, &[0.5, 0.5]).unwrap();
        assert_eq!(a.gamma, b.gamma, "at {t} K");
    }
}

/// And where the tables do carry a temperature dependence, the model must be UNIFAC
/// evaluated on `a + b T + c T^2` - which is what the resolver's split is for. This is
/// the derivation the cases' own values came from, made executable.
#[test]
fn the_effective_interaction_is_unifac_on_the_temperature_adjusted_matrix() {
    let psrk = unifac_psrk_parameters(&["water", "methane"]).unwrap();
    for t in [250.0, 400.0] {
        let adjusted = UnifacParameters {
            groups: psrk.groups.clone(),
            group_r: psrk.group_r.clone(),
            group_q: psrk.group_q.clone(),
            aij: (0..psrk.aij.len())
                .map(|i| psrk.aij[i] + psrk.bij[i] * t + psrk.cij[i] * t * t)
                .collect(),
        };
        let a = unifac_psrk_activity_coefficients(&psrk, t, &[0.5, 0.5]).unwrap();
        let b = unifac_activity_coefficients(&adjusted, t, &[0.5, 0.5]).unwrap();
        assert_eq!(a.gamma, b.gamma, "at {t} K");
        // The adjustment is not a no-op for this pair, or the test above proves nothing.
        assert!(psrk.bij.iter().any(|&v| v != 0.0));
        assert!(psrk.cij.iter().any(|&v| v != 0.0));
    }
}

/// `water`/`methane` is a pair whose `b` and `c` are non-zero, so the answer moves with
/// temperature. If the model ignored the temperature - the one way it differs from
/// plain UNIFAC - these would be equal.
#[test]
fn the_temperature_actually_moves_the_answer() {
    let params = unifac_psrk_parameters(&["water", "methane"]).unwrap();
    let cold = unifac_psrk_activity_coefficients(&params, 250.0, &[0.5, 0.5]).unwrap();
    let hot = unifac_psrk_activity_coefficients(&params, 400.0, &[0.5, 0.5]).unwrap();
    assert_ne!(cold.gamma, hot.gamma);
}

/// The resolver reads the `B` and `C` tables and not just `A`: a resolver that left the
/// two zero would make this model plain UNIFAC, and the two tables are the only source.
#[test]
fn the_resolver_reads_the_second_and_third_tables() {
    let params = unifac_psrk_parameters(&["water", "methane"]).unwrap();
    assert!(params.bij.iter().any(|&v| v != 0.0), "b is all zero");
    assert!(params.cij.iter().any(|&v| v != 0.0), "c is all zero");
    // And `a` still comes from the first table.
    let plain = unifac_parameters(&["water", "methane"]).unwrap();
    assert_eq!(params.aij, plain.aij);
    assert_eq!(params.groups, plain.groups);
}

#[test]
fn a_name_without_a_group_decomposition_is_refused() {
    let err = unifac_psrk_parameters(&["methanol", "unobtainium"]).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

#[test]
fn a_shape_mismatch_is_refused() {
    let mut params = unifac_psrk_parameters(&["methanol", "water"]).unwrap();
    params.cij = vec![0.0; 3];
    let err = unifac_psrk_activity_coefficients(&params, 298.15, &[0.5, 0.5]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let params = unifac_psrk_parameters(&["methanol", "water"]).unwrap();
    let err = unifac_psrk_activity_coefficients(&params, 298.15, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
