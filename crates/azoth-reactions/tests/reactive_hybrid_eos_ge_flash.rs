//! The reactive hybrid flash against its oracle, `HybridEosGeReactiveProbe`'s capture.
//!
//! The fluid is `SystemHybridEosGeFlashTest.createReactiveScaleSystem`'s - methane 5 /
//! CO2 0.05 / [n-heptane 2] / water 55.5 / Ca++ 6e-4 / Cl- 2e-4 / HCO3- 1e-3 at 313.15 K and
//! 50 bar - and the numbers are `captures/hybrid_eos_ge_reactive_probe.tsv`'s, beside
//! `captures/hybrid_eos_ge_reactive_loop_probe.tsv`'s for the loop's own interior.

use azoth_core::spec::ModelAlgorithm;
use azoth_eos::Cubic;
use azoth_reactions::reactive_hybrid_eos_ge_flash::{
    CoupledAlgorithm, log_activities, reactive_components, reactive_hybrid_eos_ge_flash,
};

/// The two-phase fluid, in the capture's own component order.
const TWO_PHASE: [&str; 9] = [
    "methane", "CO2", "water", "Ca++", "Cl-", "HCO3-", "CO3--", "OH-", "H3O+",
];

/// The three-phase fluid: the same feed with the oil former.
const THREE_PHASE: [&str; 10] = [
    "methane",
    "n-heptane",
    "CO2",
    "water",
    "Ca++",
    "Cl-",
    "HCO3-",
    "CO3--",
    "OH-",
    "H3O+",
];

const TWO_PHASE_MOLES: [f64; 9] = [
    5.0, 0.05, 55.5, 6.0e-4, 2.0e-4, 1.0e-3, 1.0e-10, 1.0e-10, 1.0e-10,
];
const THREE_PHASE_MOLES: [f64; 10] = [
    5.0, 2.0, 0.05, 55.5, 6.0e-4, 2.0e-4, 1.0e-3, 1.0e-10, 1.0e-10, 1.0e-10,
];

/// `ChemicalEquilibrium`'s own defaults, which `solveChemEq`'s caller passes straight through.
const CHEMISTRY_MAX_ITERATIONS: u32 = 100;
const CHEMISTRY_TOLERANCE: f64 = 1.0e-8;

/// The band the two answers agree in, over `beta`, `x` and the brine's species amounts.
///
/// Measured: `beta` lands `2e-11` relative from the capture, the brine's compositions
/// `2e-9` to `4e-9`, and its species amounts `4e-10`. The band is stated once here rather than
/// per number because it is one measurement - the two codes solve the same coupled fixed point
/// by their own routes and stop at their own points inside it.
const BAND: f64 = 1.0e-8;

fn algorithm() -> &'static ModelAlgorithm {
    azoth_eos::algorithm_of(&azoth_eos::model_gen::HYBRID_EOS_GE_FLASH_SPEC)
        .expect("the hybrid flash's own algorithm")
}

fn coupled_algorithm() -> CoupledAlgorithm<'static> {
    CoupledAlgorithm {
        algorithm: algorithm(),
        chemistry_max_iterations: CHEMISTRY_MAX_ITERATIONS,
        chemistry_tolerance: CHEMISTRY_TOLERANCE,
    }
}

fn relative(value: f64, expected: f64) -> f64 {
    if expected == 0.0 {
        (value - expected).abs()
    } else {
        ((value - expected) / expected).abs()
    }
}

/// **The two-phase reactive flash, against the endpoint capture.**
///
/// Three passes and a certified state, which is the class's own minimum and the whole of what
/// the loop's convergence test is for. The fractions, the brine's composition and its species
/// amounts are the capture's: this is the answer NeqSim reaches on
/// `reactiveGasAqueousFlashSupportsCalciteScalePotential`'s fluid.
#[test]
fn the_two_phase_reactive_flash_reproduces_the_capture() {
    let result = reactive_hybrid_eos_ge_flash(
        &TWO_PHASE,
        Cubic::Srk,
        313.15,
        50.0e5,
        &TWO_PHASE_MOLES,
        coupled_algorithm(),
    )
    .expect("the two-phase reactive flash");

    assert!(result.converged);
    assert_eq!(result.passes, 3);

    // The gas role, and the brine that takes the rest. The oil role is the vanished one, and
    // its fraction is the solver's floor rather than zero.
    assert!(relative(result.beta[0], 0.08346174755576197).abs() < BAND);
    assert!(relative(result.beta[2], 0.916538252444238).abs() < BAND);
    assert!(result.beta[1] < 1.0e-9);
    assert!((result.beta.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);

    // The brine, component by component, in the capture's own order.
    let captured = [
        (0, 8.337190524932367e-13),
        (1, 9.194988724623066e-5),
        (2, 0.9998755149894019),
        (3, 1.0811192924231683e-5),
        (4, 3.603730974743895e-6),
        (5, 1.8068514254397906e-5),
        (6, 4.720668890247729e-10),
        (7, 2.0619739366820892e-10),
        (8, 5.100610059954797e-8),
    ];
    for (index, expected) in captured {
        let value = result.x[2][index];
        assert!(
            relative(value, expected).abs() < BAND,
            "{}: {value} against the capture's {expected}",
            TWO_PHASE[index]
        );
    }

    // The species amounts, which are what `getNumberOfMolesInPhase` reports for the brine and
    // what the calcite scale potential is computed from.
    let species = [
        ("CO2", 0.005103038372766726),
        ("water", 55.491129720661775),
        ("HCO3-", 0.0010027670978232206),
        ("CO3--", 2.6198786331897105e-8),
        ("OH-", 1.1443550870656359e-8),
        ("H3O+", 2.8307385294304776e-6),
    ];
    assert_eq!(
        result.reactive_components,
        species
            .iter()
            .map(|(name, _)| (*name).to_string())
            .collect::<Vec<_>>()
    );
    for ((name, expected), value) in species.iter().zip(&result.aqueous_moles) {
        assert!(
            relative(*value, *expected).abs() < BAND,
            "{name}: {value} against the capture's {expected}"
        );
    }

    // The inventory the passes adjusted, which is the feed where nothing reacted and these
    // amounts where something did.
    assert!(relative(result.coupled_moles[1], 0.049997206913475076).abs() < BAND);
    assert!(relative(result.coupled_moles[2], 55.499994364714574).abs() < BAND);
    assert!((result.coupled_moles[0] - 5.0).abs() < 1.0e-12);
    assert!((result.coupled_moles[3] - 6.0e-4).abs() < 1.0e-15);

    // The class's own gates on the state: every element conserved and the brine neutral. The
    // element residual is the split's against the capture's `2.19e-10`, and the charge is
    // exact to rounding.
    assert!(
        result.element_residual <= 1.0e-8,
        "{}",
        result.element_residual
    );
    assert!(
        result.charge_residual <= 1.0e-8,
        "{}",
        result.charge_residual
    );
    assert!(result.max_material_balance_residual <= 1.0e-7);
    assert!(result.max_log_fugacity_residual <= 1.0e-5);

    // And the two quantities the loop itself stops on.
    assert!(result.chemical_deviation <= 1.0e-10);
    assert!(result.residual <= 1.0e-10);
}

/// **The three-phase reactive flash**: the same feed with an oil former, so all three roles are
/// active and the answer is the capture's `reactive-gas-oil-aqueous` block.
#[test]
fn the_three_phase_reactive_flash_reproduces_the_capture() {
    let result = reactive_hybrid_eos_ge_flash(
        &THREE_PHASE,
        Cubic::Srk,
        313.15,
        50.0e5,
        &THREE_PHASE_MOLES,
        coupled_algorithm(),
    )
    .expect("the three-phase reactive flash");

    assert!(result.converged);
    assert_eq!(result.passes, 3);

    let captured = [
        (0, 0.07289352057955786),
        (1, 0.039906308796565454),
        (2, 0.8872001706238768),
    ];
    for (role, expected) in captured {
        let value = result.beta[role];
        assert!(
            relative(value, expected).abs() < BAND,
            "role {role}: {value} against the capture's {expected}"
        );
    }
    assert!((result.beta.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);

    // The oil role exists here, and the brine's species are the three-phase block's.
    assert!(result.beta[1] > 0.01);
    let species = [
        ("CO2", 0.004506916293757877),
        ("HCO3-", 0.0010024287512929388),
        ("CO3--", 2.9644120982970915e-8),
        ("H3O+", 2.5007922770332012e-6),
        ("OH-", 1.2952594065994403e-8),
    ];
    for (name, expected) in species {
        let index = result
            .reactive_components
            .iter()
            .position(|component| component == name)
            .expect("a reactive component");
        let value = result.aqueous_moles[index];
        assert!(
            relative(value, expected).abs() < BAND,
            "{name}: {value} against the capture's {expected}"
        );
    }
    assert!(
        result.element_residual <= 1.0e-8,
        "{}",
        result.element_residual
    );
    assert!(
        result.charge_residual <= 1.0e-8,
        "{}",
        result.charge_residual
    );
}

/// The components the fluid drives, which is `getAllComponents`' set in this library's order.
///
/// `Ca++` and `Cl-` are both in it and neither is reactive: no surviving reaction of the pitzer
/// source names them for a carbonate brine, so they stay out of its chemistry. A rule that took
/// every component the element table covers would put them in, and the calcite the fluid holds
/// no reaction for would move.
#[test]
fn the_reactive_set_is_the_captured_one() {
    assert_eq!(
        reactive_components(&TWO_PHASE).expect("the set"),
        vec!["CO2", "water", "HCO3-", "CO3--", "OH-", "H3O+"]
    );
    assert_eq!(
        reactive_components(&THREE_PHASE).expect("the set"),
        vec!["CO2", "water", "HCO3-", "CO3--", "OH-", "H3O+"]
    );
}

/// **The activity vector the chemistry is handed, against the capture's.**
///
/// At the brine the first coupled pass is solved from. The rule is
/// [`log_activities`](azoth_reactions::reactive_hybrid_eos_ge_flash::log_activities)'s, and this
/// is the measurement it rests on: water as `ln gamma`, an ion as `ln gamma - ln gamma_inf`, and
/// a neutral that is not water with the reference phase's own `ln(m/x)` shift on top.
#[test]
fn the_activity_vector_is_the_captured_one() {
    // The brine at the first coupled pass, from the loop capture's `seed_pass_aqueous_moles`.
    let moles = [
        4.626975367632675e-11,
        0.005103321335336981,
        55.491135348959475,
        6.0e-4,
        2.0e-4,
        1.0e-3,
        1.0e-10,
        1.0e-10,
        1.0e-10,
    ];
    let total: f64 = moles.iter().sum();
    let x: Vec<f64> = moles.iter().map(|value| value / total).collect();
    let activity = log_activities(&TWO_PHASE, 313.15, 50.0e5, &x).expect("the vector");

    // The capture's `log_activity_before` for the reactive set, in the fluid's order.
    let captured = [
        (1, 1.243962971440027e-4),
        (2, 9.293683004418796e-5),
        (5, -0.04569936620212545),
        (6, -0.20257719812796537),
        (7, -0.04953084620998993),
        (8, -0.051797969636428576),
    ];
    for (index, expected) in captured {
        let value = activity[index];
        // Absolute and not relative: two of these are differences of logarithms near zero, so a
        // relative comparison measures the printout's own rounding.
        assert!(
            (value - expected).abs() < 1.0e-10,
            "{}: {value} against the capture's {expected}",
            TWO_PHASE[index]
        );
    }
}
