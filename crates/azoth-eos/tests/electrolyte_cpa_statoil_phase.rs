//! Spec-driven tests for the `eos.electrolyte_cpa_statoil_phase` model.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::electrolyte_cpa_statoil_phase::electrolyte_cpa_statoil_phase;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.electrolyte_cpa_statoil_phase";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::ElectrolyteCpaStatoilPhaseResult {
    let components: Vec<String> = case
        .list("components")
        .expect("components")
        .iter()
        .map(|name| (*name).to_string())
        .collect();
    electrolyte_cpa_statoil_phase(
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

/// **The distribution function is the statoil one, and it is not the interface default.**
///
/// `PhaseElectrolyteCPAstatoil.calc_g` overrides `PhaseCPAInterface`'s simplified
/// Carnahan-Starling `(2 - eta)/(2(1 - eta)^3)` with `1/(1 - 1.9 eta)`. The two agree only at
/// `eta -> 0`, so a port that read the default would match at every dilute state and diverge
/// exactly where the model is used.
///
/// Asserted against the value the jar prints rather than against the formula alone: the
/// aqueous phase at 40 bara and 10 % salt reports `gcpa = 1.65969673735988`, and the
/// Carnahan-Starling form at the same packing fraction is 1.77% away.
#[test]
fn the_statoil_distribution_function_is_not_the_interface_default() {
    let names = ["CO2", "water", "Na+", "Cl-"];
    let x = [
        0.019_983_819_489_548_1,
        0.917_183_084_665_173,
        0.031_416_547_922_639_7,
        0.031_416_547_922_639_7,
    ];
    let (mixture, _) =
        azoth_eos::furst_electrolyte::electrolyte_cpa_statoil_mixture_of(&names, None)
            .expect("the brine builds");
    let association = mixture.association().expect("it associates");
    assert_eq!(
        association.rdf_form(),
        azoth_eos::association::RdfForm::Chapman,
        "the statoil phase is the classical form"
    );

    let reduced = mixture
        .reduced_parameters(kelvins(283.15), pascals(4.0e6))
        .expect("reduces");
    let state = mixture
        .phase_state(&reduced, &x, azoth_eos::mixture::RootSide::Liquid)
        .expect("solves");
    assert!(
        (state.z - 0.031_404_654_632_074_5).abs() < 1.0e-8,
        "Z = {}, NeqSim gives 0.0314046546320745",
        state.z
    );

    // The same state under the interface default, which must not agree.
    let (default, _) = azoth_eos::furst_electrolyte::furst_mixture_of(
        &names,
        azoth_eos::furst_dielectric::MixingRule::default_for_the_model(),
        None,
    )
    .and_then(|(m, g)| m.with_association().map(|m| (m, g)))
    .expect("the fluid builds");
    let reduced_default = default
        .reduced_parameters(kelvins(283.15), pascals(4.0e6))
        .expect("reduces");
    let state_default = default
        .phase_state(&reduced_default, &x, azoth_eos::mixture::RootSide::Liquid)
        .expect("solves");
    assert!(
        (state.ln_phi[1] - state_default.ln_phi[1]).abs() > 1.0e-3,
        "water's ln phi is {} under the statoil form and {} under Carnahan-Starling; the \
         distribution function has to separate them",
        state.ln_phi[1],
        state_default.ln_phi[1]
    );
}

/// An ion is a component here and not a scalar, so a composition naming one has to be one.
#[test]
fn a_composition_that_is_not_one_is_refused() {
    let err = electrolyte_cpa_statoil_phase(
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
