//! Spec-driven tests for the `eos.parachor_mixture_surface_tension` model.
//!
//! The oracle is `validation/neqsim/ParachorProbe.java` and its capture
//! `captures/parachor_probe.tsv`, which prints every input the sum reads - both phases' densities,
//! molar masses and compositions, and each component's parachor - so the case pins the arithmetic
//! and not a band.

use azoth_core::AzothError;
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_eos::parachor_mixture_surface_tension::parachor_mixture_surface_tension;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.parachor_mixture_surface_tension";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::ParachorMixtureSurfaceTensionResult {
    parachor_mixture_surface_tension(
        case.vector("parachors")
            .expect("the case states the parachors"),
        kilograms_per_cubic_meter(case.input("rho_gas").expect("the case states rho_gas")),
        kilograms_per_mole(case.input("M_gas").expect("the case states M_gas")),
        case.vector("x_gas")
            .expect("the case states the gas fractions"),
        kilograms_per_cubic_meter(
            case.input("rho_liquid")
                .expect("the case states rho_liquid"),
        ),
        kilograms_per_mole(case.input("M_liquid").expect("the case states M_liquid")),
        case.vector("x_liquid")
            .expect("the case states the liquid fractions"),
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
            result.sigma.value,
            case.expected_value("sigma").expect("the case states sigma"),
            case.tolerance,
            &format!("{context} (sigma)"),
        );
        common::assert_consistent(&result, context);
    }
}

/// **Both phases are in the sum, and each with its own molar density** - which is what makes this
/// the mixture form rather than the pure-component one.
///
/// The two are told apart by moving a density alone: the liquid's is the sum's positive term and
/// the gas's the negative one, so raising the gas's density *lowers* the tension and raising the
/// liquid's raises it. A port that used one phase's molar density for both would move both the
/// same way.
#[test]
fn the_two_phases_enter_with_opposite_signs() {
    let parachors = [77.3, 191.7];
    let x_gas = [0.8356249351028912, 0.1643750648971088];
    let x_liquid = [0.09542082642782022, 0.9045791735721797];
    let run = |rho_gas: f64, rho_liquid: f64| {
        parachor_mixture_surface_tension(
            &parachors,
            kilograms_per_cubic_meter(rho_gas),
            kilograms_per_mole(0.022959902730870334),
            &x_gas,
            kilograms_per_cubic_meter(rho_liquid),
            kilograms_per_mole(0.054107691623917306),
            &x_liquid,
        )
        .expect("the sum runs")
        .sigma
        .value
    };
    let base = run(19.938328330315997, 549.3851668858096);
    assert!(
        run(39.876656660631994, 549.3851668858096) < base,
        "a denser gas"
    );
    assert!(
        run(19.938328330315997, 1098.7703337716) > base,
        "a denser liquid"
    );
}

/// The class answers `0.0` where its arithmetic throws, and this refuses instead: a density or a
/// molar mass that is not positive, and three vectors that are not the same length.
#[test]
fn the_refusals_are_named_rather_than_a_silent_zero() {
    let parachors = [77.3, 191.7];
    let x_gas = [0.8356249351028912, 0.1643750648971088];
    let x_liquid = [0.09542082642782022, 0.9045791735721797];

    let with = |rho_gas: f64, m_gas: f64| {
        parachor_mixture_surface_tension(
            &parachors,
            kilograms_per_cubic_meter(rho_gas),
            kilograms_per_mole(m_gas),
            &x_gas,
            kilograms_per_cubic_meter(549.3851668858096),
            kilograms_per_mole(0.054107691623917306),
            &x_liquid,
        )
    };
    assert!(matches!(
        with(0.0, 0.022959902730870334),
        Err(AzothError::OutOfRange { .. })
    ));
    assert!(matches!(
        with(19.9, 0.0),
        Err(AzothError::OutOfRange { .. })
    ));
    assert!(matches!(
        parachor_mixture_surface_tension(
            &parachors[..1],
            kilograms_per_cubic_meter(19.938328330315997),
            kilograms_per_mole(0.022959902730870334),
            &x_gas,
            kilograms_per_cubic_meter(549.3851668858096),
            kilograms_per_mole(0.054107691623917306),
            &x_liquid,
        ),
        Err(AzothError::InvalidInput { .. })
    ));
}
