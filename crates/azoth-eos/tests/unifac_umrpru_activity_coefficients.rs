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

/// The offset is about 298.15 K, not zero, so at that temperature the resolved matrix
/// must be the table's `a` itself - and this pair's `a` is the plain UNIFAC one, which
/// is what makes the two models comparable at that state.
///
/// **Their answers are still not equal, and this test says so.** The combinatorial term
/// differs, and on this pair - where the matrices agree and the offset vanishes - that
/// difference is the whole of the difference between the two models' cases.
#[test]
fn at_the_reference_temperature_only_the_a_matrix_acts() {
    let umr = unifac_umrpru_parameters(&["methanol", "water"], UmrpruSet::Umr).unwrap();
    let plain = unifac_parameters(&["methanol", "water"]).unwrap();
    assert_eq!(
        umr.aij, plain.aij,
        "the UMR a matrix's methanol/water pair is the plain UNIFAC one"
    );
    assert_eq!(
        azoth_eos::unifac_umrpru_activity_coefficients::umrpru_aij(&umr, 298.15),
        umr.aij,
        "at 298.15 K the offset is zero, so the resolved matrix is the table's a"
    );

    let adjusted = UnifacParameters {
        groups: umr.groups.clone(),
        group_r: umr.group_r.clone(),
        group_q: umr.group_q.clone(),
        aij: umr.aij.clone(),
    };
    let a = unifac_umrpru_activity_coefficients(&umr, 298.15, &[0.5, 0.5]).unwrap();
    let b = unifac_activity_coefficients(&adjusted, 298.15, &[0.5, 0.5]).unwrap();
    assert_ne!(
        a.gamma, b.gamma,
        "Flory-Huggins and Staverman-Guggenheim differ on the same matrix"
    );
}

/// The resolved matrix at a temperature is `a + b*dt + c*dt^2` with `dt = T - 298.15`.
///
/// Using `T` rather than `T - 298.15` is the one transcription slip this pins: the PSRK
/// tables fit about zero, and this model's do not, so the two forms of the offset are
/// what `parameters` selects between at a state away from the reference.
#[test]
fn the_offset_is_about_298_15_not_zero() {
    let umr = unifac_umrpru_parameters(&["water", "methane"], UmrpruSet::Umr).unwrap();
    for t in [320.0, 350.0] {
        let dt = t - 298.15;
        let expected: Vec<f64> = (0..umr.aij.len())
            .map(|i| umr.aij[i] + umr.bij[i] * dt + umr.cij[i] * dt * dt)
            .collect();
        assert_eq!(
            azoth_eos::unifac_umrpru_activity_coefficients::umrpru_aij(&umr, t),
            expected,
            "at {t} K"
        );

        // And the form that fits about zero gives a different matrix at this state, so
        // the two offsets are distinguishable here rather than coinciding.
        let psrk_form: Vec<f64> = (0..umr.aij.len())
            .map(|i| umr.aij[i] + umr.bij[i] * t + umr.cij[i] * t * t)
            .collect();
        assert_ne!(expected, psrk_form, "at {t} K the two offsets must differ");
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

/// The differential oracle, from `validation/neqsim/UmrCpaProbe.java`.
///
/// NeqSim's `SystemUMRCPAEoS` at 298.15 K and 70 bara, methane/water 0.98/0.02, on the
/// `_umrmc` tables: the GE phase the UMR mixing rule reads. `TPflashUMRCPADehydration`
/// `LifecycleTest` pins its water-in-gas envelope at that state, so this is the one
/// state the model's own component can be driven to.
///
/// **What this catches that no case can.** The registered cases are this library's own
/// numbers; this is NeqSim's, on the pair whose `ln gamma` the UMR-CPA model actually
/// reads. Before the combinatorial term was separated, the two differed by
/// `6.9e-06` and `1.9e-02` in `ln gamma` - small enough on the first component to read
/// as rounding, and it is not rounding.
#[test]
fn neqsims_umr_cpa_is_a_differential_oracle() {
    let params = unifac_umrpru_parameters(&["methane", "water"], UmrpruSet::Umrmc).unwrap();
    let result = unifac_umrpru_activity_coefficients(&params, 298.15, &[0.98, 0.02]).unwrap();

    // The probe's `lnGamma[0]`/`lnGamma[1]`.
    let expected = [0.000_245_756_047_888_7_f64, 0.887_875_427_408_146];
    for (i, &want) in expected.iter().enumerate() {
        let got = result.ln_gamma[i];
        assert!(
            (got - want).abs() <= 1e-12 * want.abs().max(1.0),
            "ln_gamma[{i}] = {got}, NeqSim's {want}"
        );
    }
}
