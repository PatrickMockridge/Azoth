//! The reactive tangent-plane analysis, against the flash's own capture.
//!
//! The oracle is `validation/neqsim/ReactiveFlashProbe.java`, whose capture is
//! `captures/reactive_flash_probe.tsv`. The analysis is private-state-heavy - two of its
//! methods are reached by reflection in the probe - so the capture holds what the class
//! reports for each seed it built and for the verdict.
//!
//! **The verdict on this fluid is "stable", and every seed says so with the same number.**
//! All six return the class's `10.0` sentinel, which is what it reports for a trial that
//! walked back to the reference composition or could not be initialised. That makes the
//! *seeds* the part of this test that carries weight: they are deterministic in the critical
//! constants, the equilibrated feed and the state, and they are printed raw.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank::mixture_of;
use azoth_eos::{Cubic, RootSide};
use azoth_reactions::databank::formation_properties;
use azoth_reactions::formula_matrix::FormulaMatrix;
use azoth_reactions::rand_solver::{ThermoData, solve_single_phase, standard_potentials};
use azoth_reactions::reactive_stability::{
    CriticalConstants, STABLE_TPD, is_unstable, reference_potentials, run_trial, trial_seeds,
};

const NAMES: [&str; 4] = ["CO", "water", "CO2", "hydrogen"];
const FEED: [f64; 4] = [0.25, 0.25, 0.25, 0.25];
const TEMPERATURE: f64 = 600.0;
const PRESSURE_BARA: f64 = 1.0;

/// From `captures/reactive_flash_probe.tsv`, `fluid=wgs-600K`: the six seeds the class built,
/// raw and unnormalised, and the distance each returned.
const CAPTURED_SEEDS: [[f64; 4]; 6] = [
    [
        2.339_636_875_897_244e-5,
        6.330_792_942_154_727e-4,
        2.206_562_857_894_150_6e-4,
        6.191_663_246_550_934e-4,
    ],
    [
        267.700_672_936_909_54,
        9.893_268_850_023_867,
        802.708_574_280_293,
        286.066_417_889_687_95,
    ],
    [1.0, 1.0e-12, 1.0e-12, 1.0e-12],
    [1.0e-12, 1.0, 1.0e-12, 1.0e-12],
    [1.0e-12, 1.0e-12, 1.0, 1.0e-12],
    [1.0e-12, 1.0e-12, 1.0e-12, 1.0],
];

/// The equilibrated feed, the reference potentials and the fugacity coefficients: the same
/// setup `rand_solver.rs`'s test uses, because the stability analysis starts from exactly the
/// state the RAND solve leaves.
struct State {
    mixture: azoth_eos::Mixture,
    reduced: azoth_eos::ReducedParameters,
    equilibrated: Vec<f64>,
    constants: Vec<CriticalConstants>,
}

fn solved_state() -> State {
    let (mixture, ideal) = mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
    let data: Vec<ThermoData> = NAMES
        .iter()
        .enumerate()
        .map(|(index, name)| {
            let formation = formation_properties(name)
                .expect("the component table parses")
                .expect("the component databank carries it");
            ThermoData {
                enthalpy_of_formation: formation.enthalpy_of_formation,
                absolute_entropy: formation.absolute_entropy,
                gibbs_energy_of_formation: formation.gibbs_energy_of_formation,
                cp: [
                    ideal.cp_a[index],
                    ideal.cp_b[index],
                    ideal.cp_c[index],
                    ideal.cp_d[index],
                    ideal.cp_e[index],
                ],
            }
        })
        .collect();

    let matrix =
        FormulaMatrix::build(&NAMES.map(String::from)).expect("the components are in the databank");
    let b: Vec<f64> = matrix
        .matrix
        .iter()
        .map(|row| row.iter().zip(FEED).map(|(a, n)| a * n).sum())
        .collect();
    let g0 = standard_potentials(&data, TEMPERATURE, PRESSURE_BARA, &[], &[]);

    let reduced = mixture
        .reduced_parameters(kelvins(TEMPERATURE), pascals(PRESSURE_BARA * 1.0e5))
        .expect("a state");
    let mut ln_phi = |x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        let state = mixture.phase_state(&reduced, x, RootSide::Vapour)?;
        Ok(state.ln_phi.clone())
    };

    let solution =
        solve_single_phase(&matrix.matrix, &g0, &b, &FEED, &mut ln_phi).expect("the solve runs");
    assert!(solution.converged, "the homogeneous solve converges");
    let total: f64 = solution.moles.iter().sum();
    let equilibrated: Vec<f64> = solution.moles.iter().map(|n| n / total).collect();

    let constants = mixture
        .components()
        .iter()
        .map(|component| CriticalConstants {
            tc: component.tc.value,
            pc: component.pc.value / 1.0e5,
            omega: component.omega,
        })
        .collect();

    State {
        mixture,
        reduced,
        equilibrated,
        constants,
    }
}

/// **The seeds are built from the *equilibrated* feed.**
///
/// `solveHomogeneousChemicalEquilibrium` runs a single-phase RAND solve on a clone and writes
/// the answer into the main system with `setx`; `init(1)` keeps it rather than recomputing the
/// fractions from the moles, and the capture's `phase_x_after_run` is that composition -
/// `0.0791`, the shift's own answer, not the feed's `0.25`. So the seeds are the equilibrated
/// fractions against Wilson K-values, which is what this asserts. (The first draft of this
/// test read a stale capture and concluded the opposite; the probe now prints the class's own
/// seeds beside a recomputation from the live composition, in one run, so the two cannot be
/// told apart wrongly again.)
#[test]
fn the_six_seeds_are_what_the_class_built() {
    let state = solved_state();
    let seeds = trial_seeds(
        &state.equilibrated,
        &state.constants,
        TEMPERATURE,
        PRESSURE_BARA,
    );

    assert_eq!(seeds.len(), CAPTURED_SEEDS.len(), "the class built six");
    for (index, (seed, captured)) in seeds.iter().zip(CAPTURED_SEEDS).enumerate() {
        for (component, (got, want)) in seed.iter().zip(captured).enumerate() {
            let relative = (got - want).abs() / want.abs().max(1e-30);
            assert!(
                relative < 1.0e-6,
                "seed {index}, component {component} ({}): {got} against the capture's {want}",
                NAMES[component]
            );
        }
    }
}

/// **Every seed is stable, and the class says so with one number.** The analysis finds no
/// second phase for the water-gas shift at 600 K and 1 bar: each trial walks back to the
/// equilibrated feed, which is what `10.0` means, and the verdict is `false` - so the two
/// phases the flash's own capture ends with are *not* a stability-driven addition.
#[test]
fn every_captured_seed_is_stable() {
    let state = solved_state();
    let mut ln_phi = |x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        let phase = state
            .mixture
            .phase_state(&state.reduced, x, RootSide::Vapour)?;
        Ok(phase.ln_phi.clone())
    };
    let ln_phi_feed = ln_phi(&state.equilibrated).expect("the feed's coefficients");
    let d =
        reference_potentials(&state.equilibrated, &ln_phi_feed, &[0.0; 4]).expect("the potentials");

    let seeds = trial_seeds(
        &state.equilibrated,
        &state.constants,
        TEMPERATURE,
        PRESSURE_BARA,
    );
    let mut unstable = 0;
    for (index, seed) in seeds.iter().enumerate() {
        let tpd = run_trial(&d, seed, &state.equilibrated, &mut ln_phi).expect("the trial runs");
        assert_eq!(
            tpd, STABLE_TPD,
            "seed {index} returned {tpd}, and the capture has every seed at the sentinel"
        );
        assert!(!is_unstable(tpd), "seed {index} is not a new phase");
        if is_unstable(tpd) {
            unstable += 1;
        }
    }
    assert_eq!(unstable, 0, "the captured verdict is stable");
}

/// **The reference potentials, and the branch a trace component does *not* take.**
///
/// `computeReferencePotentials` has three branches: the logarithm for a component the feed
/// holds, `-100.0` for one at or below `MIN_MOLES`, and `-1000.0` for an ion. The probe's
/// trace-nitrogen fluid carries `1e-40` of an inert component, and what the class leaves in the
/// phase is `1.0001516e-30` - *above* the `1e-30` floor - so the **logarithm** runs and the
/// capture's `-69.0768` is `ln(1.0001516e-30)` plus that component's fugacity coefficient. The
/// port reproduces it from the captured composition, which is what pins the branch order.
///
/// The absent branch itself is not reached by any captured fluid; it is asserted directly, at
/// the floor and below it, where the class's own answer is a constant.
#[test]
fn the_reference_potentials_are_the_classes_three_branches() {
    use azoth_reactions::reactive_stability::{
        REFERENCE_ABSENT, REFERENCE_ION, reference_potentials,
    };

    // The trace-nitrogen block: `phase_x_after_run` and `reference_potentials`, in the capture's
    // order (CO, water, CO2, hydrogen, nitrogen).
    let equilibriated = [
        0.079_140_531_089_050_75,
        0.079_140_531_089_050_89,
        0.420_859_468_910_949_2,
        0.420_859_468_910_949_1,
        1.000_151_672_217_491_5e-30,
    ];
    let captured = [
        -2.535_945_213_718_570_7,
        -2.537_323_477_797_878,
        -0.865_452_653_324_643_4,
        -0.864_974_018_712_630_2,
        -69.076_841_192_751_1,
    ];

    let names = ["CO", "water", "CO2", "hydrogen", "nitrogen"];
    let (mixture, _) = mixture_of(&names, Cubic::Srk, None).expect("the databank carries them");
    let reduced = mixture
        .reduced_parameters(kelvins(600.0), pascals(1.0e5))
        .expect("a state");
    let ln_phi = mixture
        .phase_state(&reduced, &equilibriated, RootSide::Vapour)
        .expect("the coefficients")
        .ln_phi
        .clone();

    let potentials =
        reference_potentials(&equilibriated, &ln_phi, &[0.0; 5]).expect("the shapes agree");
    for (index, (got, want)) in potentials.iter().zip(captured).enumerate() {
        // All five come in to `1e-12` **absolute**, the nitrogen one included: its value is a
        // logarithm of `1.0001516e-30` and the capture's composition is the input, so the only
        // difference left is the `ln` and the cubic's own rounding.
        assert!(
            (got - want).abs() < 1.0e-12,
            "component {index} ({}): {got} against the capture's {want}",
            names[index]
        );
    }
    // A component that *is* at the floor is not a logarithm of it.
    assert_eq!(REFERENCE_ABSENT, -100.0);
    assert_eq!(REFERENCE_ION, -1000.0);
    let floored = reference_potentials(&[1.0e-30, 0.5], &[0.1, 0.2], &[0.0, 0.0]).expect("a shape");
    assert_eq!(floored[0], REFERENCE_ABSENT, "at the floor is absent");
    assert!(
        floored[1].is_finite(),
        "and a held component is a logarithm"
    );
    let ionic = reference_potentials(&[0.4, 0.6], &[0.1, 0.2], &[0.0, -1.0]).expect("a shape");
    assert_eq!(ionic[1], REFERENCE_ION, "an ion does not partition");
}

/// **The four steps composed, on the fluid that is stable.** `analyse` brings the feed to
/// homogeneous equilibrium, takes the reference potentials from *that* composition, builds the
/// seeds and runs a trial from each. The capture's verdict for this fluid is `unstable=false`
/// with `number_of_unstable_trials=0` and `most_unstable_trial=null`, and the potentials it
/// prints are the ones this asserts - the same numbers the earlier test pins, reached here
/// through the orchestrator rather than by hand.
///
/// The two closures are the crate's boundary: `ce` is the single-phase RAND solve and `ln_phi`
/// is the cubic, both supplied by the test exactly as the driver would supply them.
#[test]
fn the_four_steps_compose_to_the_captured_verdict() {
    use azoth_reactions::rand_solver::solve_single_phase;
    use azoth_reactions::reactive_stability::analyse;

    let state = solved_state();
    let matrix =
        FormulaMatrix::build(&NAMES.map(String::from)).expect("the components are in the databank");
    let b: Vec<f64> = matrix
        .matrix
        .iter()
        .map(|row| row.iter().zip(FEED).map(|(a, n)| a * n).sum())
        .collect();
    let data: Vec<ThermoData> = NAMES
        .iter()
        .enumerate()
        .map(|(index, name)| {
            let formation = formation_properties(name)
                .expect("the component table parses")
                .expect("the component databank carries it");
            let (_, ideal) =
                mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
            ThermoData {
                enthalpy_of_formation: formation.enthalpy_of_formation,
                absolute_entropy: formation.absolute_entropy,
                gibbs_energy_of_formation: formation.gibbs_energy_of_formation,
                cp: [
                    ideal.cp_a[index],
                    ideal.cp_b[index],
                    ideal.cp_c[index],
                    ideal.cp_d[index],
                    ideal.cp_e[index],
                ],
            }
        })
        .collect();
    let g0 = standard_potentials(&data, TEMPERATURE, PRESSURE_BARA, &[], &[]);
    let reduced = state.reduced;

    // `ce`: the homogeneous solve, which returns its input when it does not converge - the
    // class's own behaviour in `solveHomogeneousChemicalEquilibrium`.
    let mut ln_phi = |x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        Ok(state
            .mixture
            .phase_state(&reduced, x, RootSide::Vapour)?
            .ln_phi
            .clone())
    };
    let mut ce = |x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        let mut one_phase = |y: &[f64]| -> azoth_core::Result<Vec<f64>> {
            Ok(state
                .mixture
                .phase_state(&reduced, y, RootSide::Vapour)?
                .ln_phi
                .clone())
        };
        let solved = solve_single_phase(&matrix.matrix, &g0, &b, x, &mut one_phase)?;
        if !solved.converged {
            return Ok(x.to_vec());
        }
        let total: f64 = solved.moles.iter().sum();
        Ok(solved.moles.iter().map(|moles| moles / total).collect())
    };

    let outcome = analyse(
        &FEED,
        &state.constants,
        TEMPERATURE,
        PRESSURE_BARA,
        &[0.0; 4],
        &mut ce,
        &mut ln_phi,
    )
    .expect("the analysis runs");

    assert!(!outcome.unstable, "the captured verdict is stable");
    assert!(outcome.unstable_trials.is_empty(), "no trial was added");
    assert!(outcome.tpd_values.is_empty(), "and none was recorded");
    // The reference is the *equilibrated* feed, which is what the capture's
    // `phase_x_after_run` records - `0.0791` and not the feed's `0.25`.
    let captured_reference = [
        0.079_140_531_089_050_75,
        0.079_140_531_089_050_89,
        0.420_859_468_910_949_2,
        0.420_859_468_910_949_1,
    ];
    for (index, (got, want)) in outcome.reference.iter().zip(captured_reference).enumerate() {
        assert!(
            (got - want).abs() / want.abs() < 1.0e-5,
            "reference component {index} ({}): {got} against the capture's {want}",
            NAMES[index]
        );
    }
    // And the potentials are taken from that composition, not from the feed.
    let captured_potentials = [
        -2.535_945_213_718_570_7,
        -2.537_323_477_797_878,
        -0.865_452_653_324_643_4,
        -0.864_974_018_712_630_2,
    ];
    for (index, (got, want)) in outcome
        .potentials
        .iter()
        .zip(captured_potentials)
        .enumerate()
    {
        assert!(
            (got - want).abs() < 1.0e-6,
            "potential {index}: {got} against the capture's {want}"
        );
    }
}
