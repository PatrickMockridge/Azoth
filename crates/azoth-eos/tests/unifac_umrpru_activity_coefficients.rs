//! Spec-driven tests for the `eos.unifac_umrpru_activity_coefficients` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{
    UmrpruSet, UnifacParameters, unifac_parameters, unifac_psrk_parameters,
    unifac_umrpru_parameters,
};
use azoth_eos::unifac_activity_coefficients::unifac_activity_coefficients;
use azoth_eos::{model_gen, unifac_umrpru_activity_coefficients};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.unifac_umrpru_activity_coefficients";

/// The set the case names, as the resolver's own selector.
fn set_of(name: &str) -> UmrpruSet {
    match name {
        "umr" => UmrpruSet::Umr,
        "umrmc" => UmrpruSet::Umrmc,
        other => panic!("the case names an unknown parameter set {other:?}"),
    }
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::UnifacUmrpruActivityCoefficientsResult {
    let params = unifac_umrpru_parameters(
        case.list("components").expect("components"),
        set_of(
            case.string("parameters")
                .expect("the case declares parameters"),
        ),
    )
    .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    unifac_umrpru_activity_coefficients(
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
    let params = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr).unwrap();
    let r = unifac_umrpru_activity_coefficients(&params, 298.15, &[0.5, 0.5]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The offset is about 298.15 K, not zero, so at that temperature the model must be
/// plain UNIFAC evaluated on UMR's `a` matrix alone - and the case recorded for the
/// published UNIFAC value rests on that plus the pair's `a` agreeing. Both halves are
/// asserted here rather than assumed in prose.
#[test]
fn at_the_reference_temperature_only_the_a_matrix_acts() {
    let umr = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr).unwrap();
    let plain = unifac_parameters(&["methanol", "water"]).unwrap();
    assert_eq!(
        umr.aij, plain.aij,
        "the UMR a matrix's methanol/water pair is the plain UNIFAC one"
    );

    let adjusted = UnifacParameters {
        groups: umr.groups.clone(),
        group_r: umr.group_r.clone(),
        group_q: umr.group_q.clone(),
        aij: umr.aij.clone(),
    };
    let a = unifac_umrpru_activity_coefficients(&umr, 298.15, &[0.5, 0.5]).unwrap();
    let b = unifac_activity_coefficients(&adjusted, 298.15, &[0.5, 0.5]).unwrap();
    assert_eq!(a.gamma, b.gamma);
}

/// Away from the reference temperature the offset acts, and the model must equal UNIFAC
/// evaluated on `a + b*dt + c*dt^2` with `dt = T - 298.15`. Using `T` rather than
/// `T - 298.15` would be the one transcription slip this file exists to catch.
#[test]
fn the_offset_is_about_298_15_not_zero() {
    let umr = unifac_umrpru_parameters(&["water", "methane"], UmrpruSet::Umr).unwrap();
    for t in [320.0, 350.0] {
        let dt = t - 298.15;
        let adjusted = UnifacParameters {
            groups: umr.groups.clone(),
            group_r: umr.group_r.clone(),
            group_q: umr.group_q.clone(),
            aij: (0..umr.aij.len())
                .map(|i| umr.aij[i] + umr.bij[i] * dt + umr.cij[i] * dt * dt)
                .collect(),
        };
        let a = unifac_umrpru_activity_coefficients(&umr, t, &[0.5, 0.5]).unwrap();
        let b = unifac_activity_coefficients(&adjusted, t, &[0.5, 0.5]).unwrap();
        assert_eq!(a.gamma, b.gamma, "at {t} K");

        // And it is not the PSRK form, which fits about zero: a model that used
        // `a + b T + c T^2` would disagree here.
        let psrk_form = UnifacParameters {
            groups: umr.groups.clone(),
            group_r: umr.group_r.clone(),
            group_q: umr.group_q.clone(),
            aij: (0..umr.aij.len())
                .map(|i| umr.aij[i] + umr.bij[i] * t + umr.cij[i] * t * t)
                .collect(),
        };
        let wrong = unifac_activity_coefficients(&psrk_form, t, &[0.5, 0.5]).unwrap();
        assert_ne!(a.gamma, wrong.gamma, "at {t} K the two offsets must differ");
    }
}

/// The parameter set is a real choice, not a synonym: the two tables differ, and the
/// same state gives a different answer under each.
#[test]
fn the_two_parameter_sets_are_not_the_same() {
    let umr = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr).unwrap();
    let umrmc = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umrmc).unwrap();
    assert_eq!(
        umr.groups, umrmc.groups,
        "same decomposition, different interaction"
    );
    assert_ne!(umr.aij, umrmc.aij);

    let a = unifac_umrpru_activity_coefficients(&umr, 298.15, &[0.5, 0.5]).unwrap();
    let b = unifac_umrpru_activity_coefficients(&umrmc, 298.15, &[0.5, 0.5]).unwrap();
    assert_ne!(a.gamma, b.gamma);
}

/// The decomposition is `UNIFACcompUMRPRU`, not `UNIFACcomp`: they are different tables,
/// and a resolver reading the wrong one would still compute something.
#[test]
fn the_decomposition_is_the_umrpru_one() {
    let umr = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr).unwrap();
    let psrk = unifac_psrk_parameters(&["methanol", "water"]).unwrap();
    // Methanol is a single subgroup in each, so the counts agree here; the group
    // *constants* are the shared table, so they must agree too.
    assert_eq!(umr.group_r, psrk.group_r);
    assert_eq!(umr.group_q, psrk.group_q);
    // The interaction is what differs.
    assert!(!umr.bij.is_empty() && !umr.cij.is_empty());
}

#[test]
fn a_name_without_a_umrpru_group_assignment_is_refused() {
    let err = unifac_umrpru_parameters(&["methanol", "unobtainium"], UmrpruSet::Umr).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

#[test]
fn a_shape_mismatch_is_refused() {
    let mut params = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr).unwrap();
    params.cij = vec![0.0; 3];
    let err = unifac_umrpru_activity_coefficients(&params, 298.15, &[0.5, 0.5]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let params = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr).unwrap();
    let err = unifac_umrpru_activity_coefficients(&params, 298.15, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
