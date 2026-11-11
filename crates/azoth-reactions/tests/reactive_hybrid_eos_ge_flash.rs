//! Spec-driven tests for `reactions.reactive_hybrid_eos_ge_flash`.
//!
//! The oracle is `validation/neqsim/HybridEosGeReactiveProbe.java` and its capture
//! `captures/hybrid_eos_ge_reactive_probe.tsv`, beside
//! `HybridEosGeReactiveLoopProbe.java`'s and `captures/hybrid_eos_ge_reactive_loop_probe.tsv`
//! for the loop's own interior - the per-pass deviation, residual, inventory and activity
//! vector.
//!
//! **The fluid is `SystemHybridEosGeFlashTest.createReactiveScaleSystem`'s**, in both forms, at
//! 313.15 K and 50 bar. Its three roles are `[gas, oil, aqueous]` in the solver's order, and
//! with no oil former the oil role is the vanished one.
//!
//! The case's matrix expectation - the per-role compositions - is a `SELF_ASSERTED_EXPECTATION`
//! in `tools/gen_registry.py`: the generator carries matrices as *inputs* and not as
//! expectations, so what pins `x` here is the capture transcription below and, on the Python
//! side, the case's own numbers.

use azoth_core::spec::TestCase;
use azoth_eos::Cubic;
use azoth_reactions::model_gen;
use azoth_reactions::reactive_hybrid_eos_ge_flash::{
    HYBRID_SOLVER_TOLERANCE, MAXIMUM_REACTIVE_PASSES, MINIMUM_REACTIVE_PASSES,
    REACTIVE_COMPOSITION_TOLERANCE, ReactiveHybridEosGeFlashResult, log_activities,
    reactive_components, reactive_hybrid_eos_ge_flash,
};
use azoth_test_support as common;

const MODEL_ID: &str = "reactions.reactive_hybrid_eos_ge_flash";

/// The two-phase fluid, in the capture's own component order.
const TWO_PHASE: [&str; 9] = [
    "methane", "CO2", "water", "Ca++", "Cl-", "HCO3-", "CO3--", "OH-", "H3O+",
];

/// The three-phase fluid: the same feed with the oil former.
const THREE_PHASE: [&str; 10] = [
    "methane",
    "CO2",
    "n-heptane",
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
    5.0, 0.05, 2.0, 55.5, 6.0e-4, 2.0e-4, 1.0e-3, 1.0e-10, 1.0e-10, 1.0e-10,
];

/// The band the two answers agree in, over `beta`, `x` and the brine's species amounts.
///
/// Measured: `beta` lands `2e-11` relative from the capture, the brine's compositions `2e-9` to
/// `4e-9`, its species amounts `4e-10`, and the coupled inventory `2e-9`. The band is stated
/// once because it is one measurement: the two codes solve the same coupled fixed point by their
/// own routes and stop at their own points inside it.
const BAND: f64 = 1.0e-8;

fn names(components: &[&str]) -> Vec<String> {
    components.iter().map(|name| (*name).to_string()).collect()
}

fn relative(value: f64, expected: f64) -> f64 {
    if expected == 0.0 {
        (value - expected).abs()
    } else {
        ((value - expected) / expected).abs()
    }
}

fn run(components: &[&str], moles: &[f64]) -> ReactiveHybridEosGeFlashResult {
    reactive_hybrid_eos_ge_flash(&names(components), Cubic::Srk, 313.15, 50.0e5, moles)
        .expect("the reactive hybrid flash")
}

fn call(case: &TestCase) -> ReactiveHybridEosGeFlashResult {
    let components = case.list("components").expect("components");
    let moles = case.vector("moles").expect("moles");
    reactive_hybrid_eos_ge_flash(
        &names(components),
        common::input_str(case, "cubic").parse().expect("cubic"),
        common::input(case, "T"),
        common::input(case, "P"),
        moles,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

/// **Every case the spec declares**, which is the same two fluids the capture holds.
#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model has no cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);

        // Three passes, because a two-pass answer may not certify - and every captured fluid
        // stops at the floor rather than below it.
        assert!(
            result.passes >= MINIMUM_REACTIVE_PASSES,
            "{context}: {} passes",
            result.passes
        );
        assert!(
            result.chemical_deviation <= REACTIVE_COMPOSITION_TOLERANCE,
            "{context}: deviation {}",
            result.chemical_deviation
        );
        assert!(
            result.residual <= HYBRID_SOLVER_TOLERANCE,
            "{context}: residual {}",
            result.residual
        );
        common::assert_close(
            result.phase_fractions.iter().sum::<f64>(),
            1.0,
            1.0e-12,
            &format!("{context} (beta sum)"),
        );

        if let Some(expected) = case.expected_value("passes") {
            common::assert_close(
                result.passes as f64,
                expected,
                case.tolerance,
                &format!("{context} (passes)"),
            );
        }
        for name in ["coupled_moles", "aqueous_moles", "phase_fractions"] {
            let Some(expected) = case.expected_vector(name) else {
                continue;
            };
            let got = match name {
                "coupled_moles" => &result.coupled_moles,
                "aqueous_moles" => &result.aqueous_moles,
                _ => &result.phase_fractions,
            };
            assert_eq!(got.len(), expected.len(), "{context}: {name} length");
            for (index, (value, want)) in got.iter().zip(expected).enumerate() {
                common::assert_close(
                    *value,
                    *want,
                    case.tolerance,
                    &format!("{context} ({name}[{index}])"),
                );
            }
        }
    }
}

/// **The two-phase reactive flash, against the endpoint capture.**
///
/// Three passes and a brine that retains its carbonate: this is the answer NeqSim reaches on
/// `reactiveGasAqueousFlashSupportsCalciteScalePotential`'s fluid, and it is the state the
/// calcite scale potential is computed from.
#[test]
fn the_two_phase_reactive_flash_reproduces_the_capture() {
    let result = run(&TWO_PHASE, &TWO_PHASE_MOLES);

    // **The pass count is a path quantity and is not pinned**, as the three-phase case's
    // source says: the class takes three and this port four on this state, and the state is
    // what is asserted below.
    assert!(result.passes >= 3);
    assert!(relative(result.phase_fractions[0], 0.08340628923382257).abs() < BAND);
    assert!(relative(result.phase_fractions[2], 0.9165937107661775).abs() < BAND);
    // The oil role is the one the solver drove to nothing; it carries the fraction floor.
    assert!(result.phase_fractions[1] < 1.0e-9);

    // The brine, component by component, in the capture's own order.
    let captured = [
        (0, 8.342220336303274e-13),
        (1, 1.5227074776605413e-4),
        (2, 0.9998151292793369),
        (3, 1.0810539127903152e-5),
        (4, 3.603513042634384e-6),
        (5, 1.81011904868991e-5),
        (6, 2.8609138495807395e-10),
        (7, 1.247240121639693e-10),
        (8, 8.431859003812431e-8),
    ];
    for (index, expected) in captured {
        let value = result.x[2][index];
        assert!(
            relative(value, expected).abs() < BAND,
            "{}: {value} against the capture's {expected}",
            TWO_PHASE[index]
        );
    }

    // The brine's species amounts, which are what `getNumberOfMolesInPhase` reports.
    let species = [
        ("CO2", 0.008451238886302745),
        ("water", 55.49113420432702),
        ("HCO3-", 0.0010046413193313182),
        ("CO3--", 1.5878470901768904e-8),
        ("OH-", 6.922356638553393e-9),
        ("H3O+", 4.679799353604247e-6),
    ];
    assert_eq!(
        reactive_components(&TWO_PHASE).expect("the set"),
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

    // The inventory the passes adjusted: the feed where nothing reacted, and the chemistry's
    // own amounts where something did.
    assert!(relative(result.coupled_moles[1], 0.049995343002476114).abs() < BAND);
    assert!(relative(result.coupled_moles[2], 55.49999065627943).abs() < BAND);
    assert!(relative(result.coupled_moles[5], 0.0010046413193313182).abs() < BAND);
    assert!((result.coupled_moles[0] - 5.0).abs() < 1.0e-12);
    assert!((result.coupled_moles[3] - 6.0e-4).abs() < 1.0e-15);

    // The class's own gates on the state, against the capture's `2.19e-10` and `4.17e-13`.
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
}

/// **The three-phase reactive flash**: the same feed with an oil former, so all three roles are
/// active and every row of the capture's three phase blocks is pinned.
#[test]
fn the_three_phase_reactive_flash_reproduces_the_capture() {
    let result = run(&THREE_PHASE, &THREE_PHASE_MOLES);

    assert!(result.passes >= 3);
    let captured_beta = [0.07285239347314516, 0.03989899384231154, 0.8872486126845432];
    for (role, expected) in captured_beta.iter().enumerate() {
        let value = result.phase_fractions[role];
        assert!(
            relative(value, *expected).abs() < BAND,
            "role {role}: {value} against the capture's {expected}"
        );
    }

    let captured_x = [
        [
            0.9841252651179889,
            0.007347302356230273,
            0.00677100239483458,
            0.0017564301309462372,
            1.0000000000010392e-50,
            1.0000000000010392e-50,
            1.0000000000010392e-50,
            1.0000000000010392e-50,
            1.0000000000010392e-50,
            1.0000000000010392e-50,
        ],
        [
            0.20646824314902654,
            0.0036003076983455277,
            0.7889978097022985,
            9.336394503294806e-4,
            1.0000000000004426e-50,
            1.0000000000004426e-50,
            1.0000000000004426e-50,
            1.0000000000004426e-50,
            1.0000000000004426e-50,
            1.0000000000004426e-50,
        ],
        [
            8.29390510640891e-13,
            1.3564865226317645e-4,
            2.3888644848789504e-15,
            0.9998317683467991,
            1.0811006754295511e-5,
            3.6036689180985036e-6,
            1.809271509414047e-5,
            3.2084883009135925e-10,
            1.3994654782381029e-10,
            7.514854401968292e-8,
        ],
    ];
    for (role, row) in captured_x.iter().enumerate() {
        for (index, expected) in row.iter().enumerate() {
            let value = result.x[role][index];
            assert!(
                relative(value, *expected).abs() < BAND,
                "x[{role}][{}]: {value} against the capture's {expected}",
                THREE_PHASE[index]
            );
        }
    }

    let species = [
        ("CO2", 0.007528363750340491),
        ("HCO3-", 0.0010041274881426791),
        ("CO3--", 1.7806787326104148e-8),
        ("H3O+", 4.1706686006735345e-6),
        ("OH-", 7.766892631060786e-9),
    ];
    let set = reactive_components(&THREE_PHASE).expect("the set");
    for (name, expected) in species {
        let index = set
            .iter()
            .position(|component| component == name)
            .expect("a reactive component");
        let value = result.aqueous_moles[index];
        assert!(
            relative(value, expected).abs() < BAND,
            "{name}: {value} against the capture's {expected}"
        );
    }

    // The coupled inventory, which the loop capture prints per pass.
    let coupled = [
        (0, 5.0),
        (1, 0.049995854905294075),
        (2, 2.0),
        (3, 55.499991676467666),
        (6, 0.0010041274881426785),
        (7, 1.7806787326104138e-8),
    ];
    for (index, expected) in coupled {
        let value = result.coupled_moles[index];
        assert!(
            relative(value, expected).abs() < BAND,
            "coupled[{index}]: {value} against the capture's {expected}"
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

/// **The activity vector the chemistry is handed, against the capture's.**
///
/// At the brine the first coupled pass is solved from. The rule is
/// [`log_activities`](azoth_reactions::reactive_hybrid_eos_ge_flash::log_activities)'s, and this
/// is the measurement it rests on: water as `ln gamma`, an ion as `ln gamma - ln gamma_inf`, and
/// a neutral that is not water with the reference phase's own `ln(m/x)` shift on top.
#[test]
fn the_activity_vector_is_the_captured_one() {
    // The captured composition is the phase's own; rebuild it from the loop capture's
    // `seed_pass_aqueous_moles` rather than from the rounded numbers above.
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
        // relative comparison would measure the printout's own rounding.
        assert!(
            (value - expected).abs() < 1.0e-10,
            "{}: {value} against the capture's {expected}",
            TWO_PHASE[index]
        );
    }
}

/// The components the fluid drives, which is `getAllComponents`' set in this library's order.
///
/// `Ca++` and `Cl-` are both in the element table and neither is reactive: no surviving reaction
/// of the pitzer source names them for a carbonate brine. A rule that took every component the
/// element table covers would put them into the chemistry - and would then keep the reactions
/// that name them.
#[test]
fn the_reactive_set_is_the_captured_one() {
    let expected = vec!["CO2", "water", "HCO3-", "CO3--", "OH-", "H3O+"];
    assert_eq!(reactive_components(&TWO_PHASE).expect("the set"), expected);
    assert_eq!(
        reactive_components(&THREE_PHASE).expect("the set"),
        expected
    );
}

/// **The spec's algorithm and the kernel's constants are the same numbers.**
///
/// The loop's four constants are the class's, and the spec states the two of them that are a
/// stopping rule. Nothing reads the spec at run time, so this is what keeps the pair honest.
#[test]
fn the_spec_states_the_loops_own_numbers() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    let algorithm = spec.algorithm.expect("the model is a procedure");
    assert_eq!(algorithm.max_iterations, MAXIMUM_REACTIVE_PASSES);
    assert_eq!(algorithm.tolerance, HYBRID_SOLVER_TOLERANCE);
}
