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

/// **The 2004 revision, against its own probe.**
///
/// `SystemFurstElectrolyteEosMod2004` is the base with five quantities zeroed - the solvent
/// dielectric constant's two temperature derivatives, the shielding parameter's, `XLR`'s, and
/// the solvent's composition derivative. Measured: the base's
/// `getSolventDiElectricConstantdT` is `-0.359218709298880` and this variant's is
/// `-0.00000000000000`, which is the `0 *` the source multiplies it by.
///
/// The two models are close - `Z` differs in the ninth digit - so what separates them is the
/// fugacity coefficients, and that is what this asserts.
#[test]
fn the_2004_revision_matches_its_own_oracle() {
    let names = ["methane", "water", "Na+", "Cl-"];
    let x = [
        0.000_226_524_776_743_935,
        0.997_777_272_904_135,
        0.000_998_101_159_560_386,
        0.000_998_101_159_560_386,
    ];
    let (mod2004, _) = azoth_eos::furst_electrolyte::furst_mod2004_mixture_of(
        &names,
        azoth_eos::furst_dielectric::MixingRule::default_for_the_model(),
        None,
    )
    .expect("the 2004 mixture builds");
    let reduced = mod2004
        .reduced_parameters(kelvins(298.15), pascals(1001325.0))
        .expect("reduces");
    let state = mod2004
        .phase_state(&reduced, &x, azoth_eos::mixture::RootSide::Liquid)
        .expect("solves");

    assert!(
        (state.z - 0.009_635_856_372_047_19).abs() < 1.0e-8,
        "Z = {}, NeqSim gives 0.00963585637204719",
        state.z
    );
    // Master's numbers, at the `1e-7` this settles to. They were `1e-6`-apart from its
    // 3.20.0 values before the 2004 revision stopped adding the extensive `FBornD` term.
    let want: [f64; 4] = [
        8.372_595_908_242_34,
        -5.748_783_296_884_28,
        -275.908_793_598_670,
        -166.578_354_735_104,
    ];
    for (i, &expected) in want.iter().enumerate() {
        let scale = expected.abs().max(1.0);
        assert!(
            (state.ln_phi[i] - expected).abs() < 1.0e-7 * scale,
            "ln phi[{i}] = {}, NeqSim gives {expected}",
            state.ln_phi[i]
        );
    }
}
