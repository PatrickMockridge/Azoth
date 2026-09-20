//! Spec-driven tests for the `eos.furst_electrolyte_phase` model.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::furst_electrolyte_phase::furst_electrolyte_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.furst_electrolyte_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::FurstElectrolytePhaseResult {
    let components: Vec<String> = case
        .list("components")
        .expect("components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    furst_electrolyte_phase(
        &components,
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        case.vector("x").expect("x"),
        common::input_str(case, "compressed_phase"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.z_factor,
            case.expected_value("z_factor").expect("z_factor"),
            case.tolerance,
            &format!("{context} (z_factor)"),
        );
        let expected = case.expected_vector("ln_phi").expect("ln_phi");
        for (i, (&got, &want)) in result.ln_phi.iter().zip(expected).enumerate() {
            common::assert_close(
                got,
                want,
                case.tolerance,
                &format!("{context} (ln_phi[{i}])"),
            );
        }
        common::assert_consistent(&result, context);
    }
}

/// **The root is not the cubic's.** The electrolyte terms carry a pressure, and this is
/// the claim that the residual's third term is doing something: the phase's `Z` must differ
/// from the root of the *same* cubic at the *same* reduced parameters.
///
/// Asserted as a difference and not as a value, because the case set already pins the
/// values. A port that had solved the cubic and added the electrolyte to `ln phi` alone
/// would pass every layer of the term and fail this.
#[test]
fn the_electrolyte_terms_move_the_root() {
    let case = azoth_eos::model_gen::model(MODEL_ID)
        .expect("in the table")
        .cases
        .iter()
        .find(|c| c.id == "the_shipped_tests_aqueous_phase")
        .expect("the case exists");
    let components: Vec<String> = case
        .list("components")
        .expect("components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    let brine = call(case);

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = azoth_eos::furst_electrolyte::furst_mixture_of(
        &names,
        azoth_eos::furst_dielectric::MixingRule::default_for_the_model(),
        None,
    )
    .expect("the brine builds");
    let reduced = mixture
        .reduced_parameters(kelvins(298.15), pascals(1001325.0))
        .expect("reduces");
    let x = case.vector("x").expect("x");
    let (a_mix, b_mix) = mixture.mixture_parameters(&reduced, x);
    let cubic = azoth_eos::srk_z_factor(a_mix, b_mix).expect("the cubic has a root");
    assert!(
        (brine.z_factor - cubic.z_min).abs() > 1.0e-7,
        "the brine sits at {} and its own cubic at {}, so the electrolyte's pressure is \
         what separates them",
        brine.z_factor,
        cubic.z_min
    );
}

/// An ion is a component here and not a scalar, which is the whole difference from the
/// Soreide-Whitson model - so a composition naming one has to be a composition.
#[test]
fn a_composition_that_is_not_one_is_refused() {
    let err = furst_electrolyte_phase(
        &["water".to_string(), "Na+".to_string()],
        kelvins(313.15),
        pascals(500_000.0),
        &[0.6, 0.6],
        "liquid",
    )
    .expect_err("the fractions sum to 1.2");
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("x"));
}
