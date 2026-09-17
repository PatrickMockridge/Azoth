//! Spec-driven tests for the `eos.ge_uniquac_phase` model.

use azoth_core::AzothError;
use azoth_eos::databank::{ge_nrtl_phase_parameters, ge_uniquac_phase_parameters};
use azoth_eos::ge_uniquac_phase::ge_uniquac_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_uniquac_phase";

/// The interaction matrix the spec's first case states.
const AIJ: [[f64; 2]; 2] = [[0.0, -71.0], [209.0, 0.0]];

fn aij_of(case: &azoth_core::spec::TestCase) -> Vec<Vec<f64>> {
    let flat = case.matrix("aij").expect("the case declares aij");
    let n = case.vector("x").expect("x").len();
    flat.chunks(n).map(<[f64]>::to_vec).collect()
}

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GeUniquacPhaseResult {
    let params = ge_uniquac_phase_parameters(case.list("components").expect("components"), None)
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    ge_uniquac_phase(
        &params,
        common::input(case, "T"),
        common::input(case, "P"),
        case.vector("x").expect("x"),
        &aij_of(case),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        for name in ["gamma", "ln_gamma", "ln_phi"] {
            let actual = match name {
                "gamma" => &result.gamma,
                "ln_gamma" => &result.ln_gamma,
                _ => &result.ln_phi,
            };
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
        let expected = case
            .expected_vector("p_sat")
            .expect("the case declares p_sat");
        for (i, (&got, &want)) in result.p_sat.iter().zip(expected).enumerate() {
            common::assert_close(
                got.value,
                want,
                case.tolerance,
                &format!("{context} (p_sat[{i}])"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// The phase is `gamma_i P0_i / P` and nothing else - the identity every GE phase
/// satisfies, because `ComponentGE.fugcoef` is one method and none of them overrides it.
#[test]
fn the_fugacity_coefficient_is_gamma_times_p0_over_p() {
    let params = ge_uniquac_phase_parameters(&["methanol", "water"], None).unwrap();
    let (t, p) = (298.15, 100_000.0);
    let aij: Vec<Vec<f64>> = AIJ.iter().map(|row| row.to_vec()).collect();
    let r = ge_uniquac_phase(&params, t, p, &[0.5, 0.5], &aij).unwrap();
    for i in 0..2 {
        let composed = r.gamma[i].ln() + (r.p_sat[i].value / p).ln();
        assert!(
            (r.ln_phi[i] - composed).abs() < 1e-15,
            "component {i}: ln_phi is {} but ln(gamma P0 / P) is {composed}",
            r.ln_phi[i]
        );
    }
}

/// The activity coefficients are `eos.uniquac_activity_coefficients`', not a second copy.
#[test]
fn the_activity_coefficients_are_the_uniquac_models() {
    let names = ["methanol", "water"];
    let params = ge_uniquac_phase_parameters(&names, None).unwrap();
    let aij: Vec<Vec<f64>> = AIJ.iter().map(|row| row.to_vec()).collect();
    let r = ge_uniquac_phase(&params, 298.15, 100_000.0, &[0.5, 0.5], &aij).unwrap();

    let activity = azoth_eos::uniquac_activity_coefficients::uniquac_activity_coefficients(
        &azoth_eos::databank::uniquac_parameters(&names).unwrap(),
        298.15,
        &[0.5, 0.5],
        &aij,
    )
    .unwrap();
    assert_eq!(r.gamma, activity.gamma);
}

/// `aij` is directional, and the phase reads it that way.
///
/// Transposing the matrix must change the answer. A phase that built `tau_ij` from
/// `aij[j][i]`, or that symmetrised the matrix on the way in, would pass every other
/// test here and be wrong for every real pair - no UNIQUAC matrix is symmetric.
#[test]
fn the_interaction_matrix_is_directional() {
    let names = ["methanol", "water"];
    let params = ge_uniquac_phase_parameters(&names, None).unwrap();
    let aij: Vec<Vec<f64>> = AIJ.iter().map(|row| row.to_vec()).collect();
    let transposed: Vec<Vec<f64>> = (0..2)
        .map(|i| (0..2).map(|j| AIJ[j][i]).collect())
        .collect();

    let forward = ge_uniquac_phase(&params, 298.15, 100_000.0, &[0.5, 0.5], &aij).unwrap();
    let reversed = ge_uniquac_phase(&params, 298.15, 100_000.0, &[0.5, 0.5], &transposed).unwrap();
    assert_ne!(
        forward.gamma, reversed.gamma,
        "transposing aij must change the activity coefficients"
    );
}

/// The area and volume parameters are the UNIFAC group sums, not NeqSim's own columns.
///
/// `eos.uniquac_activity_coefficients` carries the argument; this asserts the *phase* did
/// not resolve a second set. `rUNIQUAQ`/`qUNIQUAQ` are 0.0 for 109 of `UNIFACcomp.csv`'s
/// 112 rows, so a phase reading them would divide by zero for almost every mixture.
#[test]
fn the_area_and_volume_parameters_are_the_group_sums() {
    let names = ["methanol", "water"];
    let phase = ge_uniquac_phase_parameters(&names, None).unwrap();
    let model = azoth_eos::databank::uniquac_parameters(&names).unwrap();
    assert_eq!(phase.r, model.r);
    assert_eq!(phase.q, model.q);
    assert_eq!(phase.r, vec![1.4311, 0.92]);
    assert_eq!(phase.q, vec![1.432, 1.4]);
}

#[test]
fn a_matrix_that_is_not_n_by_n_is_refused() {
    let params = ge_uniquac_phase_parameters(&["methanol", "water"], None).unwrap();
    for bad in [
        vec![vec![0.0, 1.0]],            // one row for two components
        vec![vec![0.0, 1.0], vec![1.0]], // a short row
        vec![vec![0.0], vec![1.0, 0.0]], // a short first row
    ] {
        let err = ge_uniquac_phase(&params, 298.15, 100_000.0, &[0.5, 0.5], &bad).unwrap_err();
        assert!(
            matches!(err, AzothError::InvalidInput { .. }),
            "{bad:?}: {err:?}"
        );
    }
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let params = ge_uniquac_phase_parameters(&["methanol", "water"], None).unwrap();
    let aij: Vec<Vec<f64>> = AIJ.iter().map(|row| row.to_vec()).collect();
    let err = ge_uniquac_phase(&params, 298.15, 100_000.0, &[0.6, 0.6], &aij).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

/// The four phases that share `ComponentGE.fugcoef` agree on `P0` and disagree on `gamma`.
///
/// NRTL, UNIFAC and Wilson were the first three; UNIQUAC is the fourth. One identity at
/// one state, four different activity models, which is what says the arithmetic is shared
/// and the models are not.
#[test]
fn four_ge_phases_share_the_vapour_pressure_and_not_the_activity() {
    let names = ["methanol", "water"];
    let (t, p) = (298.15, 100_000.0);
    let x = [0.5, 0.5];
    let aij: Vec<Vec<f64>> = AIJ.iter().map(|row| row.to_vec()).collect();

    let nrtl = azoth_eos::ge_nrtl_phase::ge_nrtl_phase(
        &ge_nrtl_phase_parameters(&names, None).unwrap(),
        t,
        p,
        &x,
    )
    .unwrap();
    let unifac = azoth_eos::ge_unifac_phase::ge_unifac_phase(
        &azoth_eos::databank::ge_unifac_phase_parameters(&names, None).unwrap(),
        t,
        p,
        &x,
    )
    .unwrap();
    let (mixture, _) = azoth_eos::databank::mixture_of(&names, None).unwrap();
    let wilson = azoth_eos::ge_wilson_phase::ge_wilson_phase(
        &azoth_eos::databank::ge_wilson_phase_parameters(&names, None).unwrap(),
        &mixture,
        t,
        p,
        &x,
    )
    .unwrap();
    let uniquac = ge_uniquac_phase(
        &ge_uniquac_phase_parameters(&names, None).unwrap(),
        t,
        p,
        &x,
        &aij,
    )
    .unwrap();

    for other in [&nrtl.p_sat, &unifac.p_sat, &wilson.p_sat] {
        assert_eq!(&uniquac.p_sat, other, "every GE phase reads the same P0");
    }
    assert_ne!(uniquac.gamma, nrtl.gamma);
    assert_ne!(uniquac.gamma, unifac.gamma);
    assert_ne!(uniquac.gamma, wilson.gamma);
}
