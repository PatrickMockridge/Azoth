//! Spec-driven tests for the `eos.ge_unifac_phase` model.

use azoth_core::AzothError;
use azoth_eos::databank::{ge_nrtl_phase_parameters, ge_unifac_phase_parameters};
use azoth_eos::ge_unifac_phase::ge_unifac_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_unifac_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GeUnifacPhaseResult {
    let params = ge_unifac_phase_parameters(case.list("components").expect("components"), None)
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    ge_unifac_phase(
        &params,
        common::input(case, "T"),
        common::input(case, "P"),
        case.vector("x").expect("x"),
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

/// The phase is `gamma_i P0_i / P` and nothing else.
///
/// The identity is the whole of `ComponentGE.fugcoef`, and it is the same method every GE
/// phase inherits - which is why this test and its NRTL counterpart look alike rather
/// than merely similar.
#[test]
fn the_fugacity_coefficient_is_gamma_times_p0_over_p() {
    let params = ge_unifac_phase_parameters(&["methanol", "water"], None).unwrap();
    let (t, p) = (298.15, 100_000.0);
    let r = ge_unifac_phase(&params, t, p, &[0.5, 0.5]).unwrap();
    for i in 0..2 {
        let composed = r.gamma[i].ln() + (r.p_sat[i].value / p).ln();
        assert!(
            (r.ln_phi[i] - composed).abs() < 1e-15,
            "component {i}: ln_phi is {} but ln(gamma P0 / P) is {composed}",
            r.ln_phi[i]
        );
    }
}

/// The activity coefficients are `eos.unifac_activity_coefficients`', not a second copy.
#[test]
fn the_activity_coefficients_are_the_unifac_models() {
    let names = ["methanol", "water"];
    let phase = ge_unifac_phase_parameters(&names, None).unwrap();
    let r = ge_unifac_phase(&phase, 298.15, 100_000.0, &[0.5, 0.5]).unwrap();

    let activity = azoth_eos::unifac_activity_coefficients::unifac_activity_coefficients(
        &azoth_eos::databank::unifac_parameters(&names).unwrap(),
        298.15,
        &[0.5, 0.5],
    )
    .unwrap();
    assert_eq!(r.gamma, activity.gamma);
}

/// The phase and the NRTL phase differ only in where `gamma` comes from.
///
/// Both are `ComponentGE.fugcoef`, so at one state the two must agree on `P0` exactly and
/// disagree on `gamma` - which is what says the shared arithmetic is shared and the
/// activity models are not.
#[test]
fn two_ge_phases_share_the_vapour_pressure_and_not_the_activity() {
    let names = ["methanol", "water"];
    let t = 298.15;
    let p = 100_000.0;
    let x = [0.5, 0.5];
    let unifac =
        ge_unifac_phase(&ge_unifac_phase_parameters(&names, None).unwrap(), t, p, &x).unwrap();
    let nrtl = azoth_eos::ge_nrtl_phase::ge_nrtl_phase(
        &ge_nrtl_phase_parameters(&names, None).unwrap(),
        t,
        p,
        &x,
    )
    .unwrap();

    assert_eq!(unifac.p_sat, nrtl.p_sat);
    assert_ne!(unifac.gamma, nrtl.gamma);
}

/// A component NeqSim's database tags a Henry's-law solute is refused, as for NRTL.
#[test]
fn a_henrys_law_component_is_refused_rather_than_computed() {
    for name in ["CO2", "methane", "n-hexane"] {
        let err = ge_unifac_phase_parameters(&[name, "water"], None).unwrap_err();
        assert!(
            matches!(err, AzothError::InvalidInput { .. }),
            "{name}: {err:?}"
        );
    }
}

/// A component with no UNIFAC group assignment is refused rather than given `r = q = 0`.
///
/// `ComponentGEUnifac` refuses it too, and for the reason it states: without groups `R`
/// and `Q` are zero and the activity coefficient is `NaN`.
#[test]
fn a_component_with_no_group_assignment_is_refused() {
    // Ethanol has no row in NeqSim's `UNIFACcomp.csv`.
    let err = ge_unifac_phase_parameters(&["ethanol", "water"], None).unwrap_err();
    assert!(
        matches!(
            err,
            AzothError::PropertyUnavailable { .. } | AzothError::InvalidInput { .. }
        ),
        "{err:?}"
    );
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let params = ge_unifac_phase_parameters(&["methanol", "water"], None).unwrap();
    let err = ge_unifac_phase(&params, 298.15, 100_000.0, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
