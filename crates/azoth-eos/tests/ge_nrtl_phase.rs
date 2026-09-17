//! Spec-driven tests for the `eos.ge_nrtl_phase` model.

use azoth_core::{AzothError, CalcResult};
use azoth_eos::databank::{self, ge_nrtl_phase_parameters};
use azoth_eos::ge_nrtl_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.ge_nrtl_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::GeNrtlPhaseResult {
    let params = ge_nrtl_phase_parameters(case.list("components").expect("components"), None)
        .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    ge_nrtl_phase(
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

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let r = ge_nrtl_phase(&params, 298.15, 100_000.0, &[0.5, 0.5]).unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}

/// The phase is `gamma_i P0_i / P` and nothing else.
///
/// The identity is the whole of `ComponentGE.fugcoef`, so asserting it against the two
/// factors the result reports is what says the composition happened rather than some
/// other arrangement of the same three numbers.
#[test]
fn the_fugacity_coefficient_is_gamma_times_p0_over_p() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let (t, p) = (298.15, 100_000.0);
    let r = ge_nrtl_phase(&params, t, p, &[0.5, 0.5]).unwrap();
    for i in 0..2 {
        let composed = r.gamma[i].ln() + (r.p_sat[i].value / p).ln();
        assert!(
            (r.ln_phi[i] - composed).abs() < 1e-15,
            "component {i}: ln_phi is {} but ln(gamma P0 / P) is {composed}",
            r.ln_phi[i]
        );
    }
}

/// The activity coefficients are `eos.nrtl_activity_coefficients`', not a second copy.
///
/// The phase resolves its NRTL matrices through the same function that model does, and
/// this is what says so: a phase that built its own matrices could drift from the model
/// whose spec the matrices belong to.
#[test]
fn the_activity_coefficients_are_the_nrtl_models() {
    let names = ["methanol", "water"];
    let phase = ge_nrtl_phase_parameters(&names, None).unwrap();
    let model = databank::nrtl_parameters(&names, None).unwrap();
    assert_eq!(phase.alpha, model.alpha);
    assert_eq!(phase.dij, model.dij);

    let r = ge_nrtl_phase(&phase, 298.15, 100_000.0, &[0.5, 0.5]).unwrap();
    let pure = azoth_eos::nrtl_activity_coefficients(&model, 298.15, &[0.5, 0.5]).unwrap();
    assert_eq!(r.gamma, pure.gamma);
}

/// A component the databank carries no Antoine correlation for is refused.
///
/// The phase's fugacity coefficient is `gamma_i P0_i / P`, and a component without a
/// correlation has no `P0` - a card-added substance is exactly that case, since a card
/// states the parameters a cubic reads.
#[test]
fn a_component_without_a_vapour_pressure_correlation_is_refused() {
    let overlay = {
        let mut card = azoth_eos::databank::Overlay::new();
        card.set_component(
            "unobtainium",
            azoth_eos::databank::ComponentOverride {
                tc: Some(500.0),
                pc: Some(4_000_000.0),
                omega: Some(0.2),
                cp: None,
            },
        );
        card
    };
    let err = ge_nrtl_phase_parameters(&["methanol", "unobtainium"], Some(&overlay)).unwrap_err();
    assert!(
        matches!(err, AzothError::PropertyUnavailable { .. }),
        "{err:?}"
    );
}

/// A component NeqSim's database tags as a Henry's-law solute is refused.
///
/// `ComponentGE.fugcoef` takes one of two branches on the component's
/// `REFERENCESTATETYPE`: `gamma_i P0_i / P` for `solvent`, and a Henry's-law coefficient
/// for anything else. Only the first is ported, so computing the second's expression for
/// a solute would be a wrong number with nothing to show that it is wrong. 49 of the
/// databank's 173 substances are tagged this way, and the two acids carry the literal
/// `0.0`.
#[test]
fn a_henrys_law_component_is_refused_rather_than_computed() {
    for name in ["CO2", "methane", "n-hexane", "nitric acid"] {
        let err = ge_nrtl_phase_parameters(&[name, "water"], None).unwrap_err();
        assert!(
            matches!(err, AzothError::InvalidInput { .. }),
            "{name}: {err:?}"
        );
        assert_eq!(err.field(), Some("components"), "{name}");
    }
}

/// The substances the phase does describe are the ones the database calls `solvent`.
///
/// Stated as a test rather than left to the case files, because it is what makes those
/// cases legitimate: a case whose components were not solvent-tagged would be exercising
/// a branch the model does not implement.
#[test]
fn the_solvent_tagged_substances_are_the_ones_that_resolve() {
    for name in ["water", "methanol", "ethanol", "MEG", "acetone"] {
        let params = ge_nrtl_phase_parameters(&[name, "water"], None).unwrap_or_else(|e| {
            panic!("{name} is tagged solvent, so the phase should resolve it: {e}")
        });
        assert_eq!(params.antoine.len(), 2, "{name}");
    }
}

#[test]
fn a_composition_that_does_not_sum_to_one_is_refused() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let err = ge_nrtl_phase(&params, 298.15, 100_000.0, &[0.6, 0.6]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}

#[test]
fn a_length_mismatch_is_refused() {
    let params = ge_nrtl_phase_parameters(&["methanol", "water"], None).unwrap();
    let err = ge_nrtl_phase(&params, 298.15, 100_000.0, &[1.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("components"));
}
