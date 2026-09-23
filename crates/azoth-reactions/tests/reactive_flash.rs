//! The reactive flash driver's phase list, against its own capture.
//!
//! The oracle is `validation/neqsim/ReactiveFlashProbe.java`, whose capture is
//! `captures/reactive_flash_probe.tsv`. The blocks read here are the `-forced-one-phase` ones,
//! which are the only states that reach the driver's VLE initialisation: `SystemSrkEos`
//! constructs with two phases, so a system that has only been initialised hands the driver
//! `np = 2`, `maxPhases = 2` and its `skipStability` rule - and the phase list has to be cut
//! back with `setNumberOfPhases(1)` *after* the last `init` for the other branch to be
//! reachable at all.
//!
//! **What the two branches report.** Forced to one phase, the water-gas shift at 600 K takes
//! `solveSinglePhaseChemicalEquilibrium` and reports `final_gibbs_energy = 0.0`, because that
//! branch returns before the driver computes one. At 300 K the Rachford-Rice root comes out
//! interior, the VLE initialisation adds a phase, and the driver reports the beta-weighted sum
//! over both - `-0.7162112`, which is one phase's worth because the two betas sum to one.
//! The class's own default states, two phases at `beta = 1.0` each, report **twice** one
//! phase's worth; the mechanism is measured in the capture and stated in the module docs.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank::mixture_of;
use azoth_eos::{Cubic, RootSide};
use azoth_reactions::databank::formation_properties;
use azoth_reactions::formula_matrix::FormulaMatrix;
use azoth_reactions::rand_solver::{ThermoData, standard_potentials};
use azoth_reactions::reactive_flash::{PhaseGibbs, total_gibbs_energy, vle_initialization};
use azoth_reactions::reactive_stability::CriticalConstants;

const NAMES: [&str; 4] = ["CO", "water", "CO2", "hydrogen"];
const FEED: [f64; 4] = [0.25, 0.25, 0.25, 0.25];

/// From `fluid=wgs-300K-forced-one-phase`, `rachford_rice_v`, and `vle_x[0]`/`vle_x[1]`.
const CAPTURED_VAPOUR_FRACTION: f64 = 0.789_403_682_026_609;
const CAPTURED_LIQUID: [f64; 4] = [
    3.438_771_947_499_286_7e-4,
    0.994_383_058_112_622_2,
    0.004_685_834_312_447_280_6,
    5.872_303_801_806_66e-4,
];
const CAPTURED_VAPOUR: [f64; 4] = [
    0.316_603_008_599_469_05,
    0.051_414_238_658_794_65,
    0.315_444_665_659_592_06,
    0.316_538_087_082_144_26,
];

/// From `fluid=wgs-300K-forced-one-phase`: the betas the VLE initialisation set, which the
/// driver's `computeGibbsEnergy` weighs the phases by, and the energy it reports.
const CAPTURED_BETA_LIQUID: f64 = 0.210_596_317_973_390_96;
const CAPTURED_BETA_VAPOUR: f64 = 0.789_403_682_026_609;
const CAPTURED_GIBBS_300K: f64 = -0.716_211_225_797_898_8;

/// From `fluid=wgs-600K`, `equilibrium_moles[0]`: the moles the multiphase solve leaves, which
/// for a split of two identical phases is the whole fluid's composition.
const CAPTURED_MOLES: [f64; 4] = [
    0.079_140_502_548_023_4,
    0.079_140_502_560_693_7,
    0.420_859_061_906_977_14,
    0.420_859_061_880_687_67,
];

/// From `fluid=wgs-600K`, `phase[0]` and `phase[1]`: the driver's default two-phase state,
/// both at `beta = 1.0`, and the energy their sum reports.
const CAPTURED_PHASE_600K: [[f64; 4]; 2] = [
    [
        0.079_140_571_487_661_56,
        0.079_140_571_500_331_87,
        0.420_859_428_519_148,
        0.420_859_428_492_858_5,
    ],
    [
        0.079_140_575_963_793_04,
        0.079_140_576_002_380_47,
        0.420_859_424_057_090_15,
        0.420_859_423_976_736_2,
    ],
];
const CAPTURED_GIBBS_600K: f64 = -2.259_535_542_715_054_7;

/// From `fluid=wgs-300K-forced-one-phase`, `equilibrium_moles[0]` and `[1]`: the moles the
/// multiphase solve leaves in each phase, whose sum is the fluid's composition.
const CAPTURED_MOLES_300K: [[f64; 4]; 2] = [
    [
        2.301_530_064_017_200_4e-5,
        2.301_531_734_980_803_5e-5,
        0.007_045_898_122_781_015,
        0.007_045_897_177_273_817,
    ],
    [
        0.001_604_927_761_908_962_3,
        0.001_604_927_764_738_454_3,
        0.491_332_069_020_232_46,
        0.491_332_082_246_425_8,
    ],
];

/// The 300 K converged phase compositions, `phase_x[0]` and `phase_x[1]`.
const CONVERGED_300K: [[f64; 4]; 2] = [
    [
        0.001_627_923_612_413_166_7,
        0.001_627_924_794_323_039_4,
        0.498_372_109_235_549_7,
        0.498_372_042_357_714_1,
    ],
    [
        0.001_627_923_802_483_865_3,
        0.001_627_923_805_353_899_3,
        0.498_372_069_488_229_4,
        0.498_372_082_903_932_87,
    ],
];

/// The component databank's three formation columns and its heat-capacity polynomial, in the
/// fluid's order - the same data `computeG0` reads.
fn formation_data() -> Vec<ThermoData> {
    let (_, ideal) = mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
    NAMES
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
        .collect()
}

fn mixture_and_constants() -> (azoth_eos::Mixture, Vec<CriticalConstants>) {
    let (mixture, _) = mixture_of(&NAMES, Cubic::Srk, None).expect("the databank carries them");
    let constants = mixture
        .components()
        .iter()
        .map(|component| CriticalConstants {
            tc: component.tc.value,
            pc: component.pc.value / 1.0e5,
            omega: component.omega,
        })
        .collect();
    (mixture, constants)
}

fn ln_phi_at(
    mixture: &azoth_eos::Mixture,
    temperature: f64,
    pressure_bara: f64,
    x: &[f64],
) -> Vec<f64> {
    let reduced = mixture
        .reduced_parameters(kelvins(temperature), pascals(pressure_bara * 1.0e5))
        .expect("a state");
    mixture
        .phase_state(&reduced, x, RootSide::Vapour)
        .expect("the coefficients")
        .ln_phi
        .clone()
}

/// **The driver's VLE initialisation, on the one fluid that reaches it.** `initializeWithVLEFlash`
/// is the second entry point into the driver's multiphase machinery: at 300 K the water-gas
/// shift's Rachford-Rice root is `0.7894`, the class adds a phase, and it sets the water-rich
/// liquid and the light vapour that the reactive solve then starts from. Unlike the Wilson
/// seeds of the stability analysis, these are not a recomputation - the probe calls the
/// class's own private method on a fresh system and prints what it left.
#[test]
fn the_vle_initialisation_is_the_classes_own() {
    let (_, constants) = mixture_and_constants();
    let initialisation =
        vle_initialization(&FEED, &constants, 300.0, 1.0).expect("the 300 K root is interior");

    assert!(
        (initialisation.vapour_fraction - CAPTURED_VAPOUR_FRACTION).abs()
            / CAPTURED_VAPOUR_FRACTION
            < 1.0e-12,
        "V is {} against the capture's {CAPTURED_VAPOUR_FRACTION}",
        initialisation.vapour_fraction
    );
    for (name, got, want) in [("liquid", &initialisation.liquid, CAPTURED_LIQUID)]
        .into_iter()
        .chain([("vapour", &initialisation.vapour, CAPTURED_VAPOUR)])
    {
        for (component, (got, want)) in got.iter().zip(want).enumerate() {
            assert!(
                (got - want).abs() / want.abs() < 1.0e-9,
                "{name} component {component} ({}): {got} against the capture's {want}",
                NAMES[component]
            );
        }
    }
}

/// **A root at the bound is no initialisation at all, and the driver adds no phase.** At 600 K
/// every component's Wilson K exceeds one, so the Rachford-Rice iteration is driven to `V = 1`
/// and `initializeWithVLEFlash` returns before it touches the system - which the capture's
/// `after_vle_init_phases=1` and its untouched `vle_x[0] = 0.25 0.25 0.25 0.25` both show.
#[test]
fn the_vapour_bound_is_no_initialisation() {
    let (_, constants) = mixture_and_constants();
    assert!(
        vle_initialization(&FEED, &constants, 600.0, 1.0).is_none(),
        "the 600 K root is V = 1, which the class returns from"
    );
}

/// **The captured Gibbs measure is the phase list's own sum, weights and all.** `computeGibbsEnergy`
/// weighs each phase by `phase.getBeta()`, and the capture holds both regimes: the class's
/// default states carry two phases at `beta = 1.0` and report twice one phase's worth, while
/// the VLE-initialised split carries `(0.2106, 0.7894)` and reports one phase's worth. Both
/// numbers are reproduced here from the captured compositions, so the doubling is a measured
/// consequence of the weights rather than a ratio that happens to come out at two.
#[test]
fn the_captured_gibbs_measure_is_the_phase_list_sum() {
    let (mixture, _) = mixture_and_constants();

    let double = total_gibbs_energy(&CAPTURED_PHASE_600K.map(|x| PhaseGibbs {
        beta: 1.0,
        ln_phi: ln_phi_at(&mixture, 600.0, 1.0, &x),
        fractions: x.to_vec(),
    }));
    assert!(
        (double - CAPTURED_GIBBS_600K).abs() / CAPTURED_GIBBS_600K.abs() < 1.0e-5,
        "the default state's two phases at beta one sum to {double} against the capture's \
         {CAPTURED_GIBBS_600K}"
    );

    let split = total_gibbs_energy(&[
        PhaseGibbs {
            beta: CAPTURED_BETA_LIQUID,
            fractions: CONVERGED_300K[0].to_vec(),
            ln_phi: ln_phi_at(&mixture, 300.0, 1.0, &CONVERGED_300K[0]),
        },
        PhaseGibbs {
            beta: CAPTURED_BETA_VAPOUR,
            fractions: CONVERGED_300K[1].to_vec(),
            ln_phi: ln_phi_at(&mixture, 300.0, 1.0, &CONVERGED_300K[1]),
        },
    ]);
    assert!(
        (split - CAPTURED_GIBBS_300K).abs() / CAPTURED_GIBBS_300K.abs() < 1.0e-5,
        "the initialised split's phases sum to {split} against the capture's {CAPTURED_GIBBS_300K}"
    );

    // The weights are what makes the difference: the same two 600 K phases weighted by the
    // 300 K split's betas report half as much, which is the claim the doubling rests on.
    let reweighted = total_gibbs_energy(&CAPTURED_PHASE_600K.map(|x| PhaseGibbs {
        beta: 0.5,
        ln_phi: ln_phi_at(&mixture, 600.0, 1.0, &x),
        fractions: x.to_vec(),
    }));
    assert!(
        (reweighted - CAPTURED_GIBBS_600K / 2.0).abs() / (CAPTURED_GIBBS_600K / 2.0).abs() < 1.0e-5,
        "halving the weights halves the sum: {reweighted}"
    );
}

/// **What the multiphase solve reproduces, and what it does not.**
///
/// The class's `solve()` is reached with two phases on both captured states - the pair
/// `SystemSrkEos` constructs, which is what the driver's `skipStability` rule leaves in place -
/// and it drives them with the feed's frozen element balance. Its answer has two parts and this
/// port reproduces one of them:
///
/// * the **composition** - the moles summed over the phases - to `1e-6` on the class's own
///   fluid, which is the part the element balance and the equilibrium determine;
/// * not the **split**. Both phases converge to the *same* composition, so the Gibbs energy is
///   flat along the direction that trades moles between them: any split satisfies the
///   equilibrium conditions, and which one a run reports is decided by its path. NeqSim's
///   600 K run ends at `(1.0, ~0)` and this port at `(0.5, 0.5)`; at 300 K NeqSim reports
///   `(0.014138, 0.985862)` and this port `(0.016, 0.984)` - all four on the same flat line,
///   and neither code converges the split tighter than its own `1e-4` relaxed tolerance.
///
/// The two phases are both of the class's `gas` type, so both take the cubic's vapour root;
/// rooting the water-rich phase liquid instead is what collapses the split to a single phase
/// here, which is why the closure is handed the phase's index.
#[test]
fn the_multiphase_solve_reproduces_the_composition() {
    use azoth_reactions::rand_solver::PhaseFeed;

    let (mixture, constants) = mixture_and_constants();
    let matrix =
        FormulaMatrix::build(&NAMES.map(String::from)).expect("the components are in the databank");
    let b: Vec<f64> = matrix
        .matrix
        .iter()
        .map(|row| row.iter().zip(FEED).map(|(a, n)| a * n).sum())
        .collect();
    let formation = formation_data();
    let initialisation =
        vle_initialization(&FEED, &constants, 300.0, 1.0).expect("the 300 K root is interior");
    let vapour_fraction = initialisation.vapour_fraction;

    // The class's own two-phase state: `SystemSrkEos`' pair, both holding the whole feed at
    // `beta = 1.0`, which is what the captured `wgs-600K` block is.
    let constructor_pair = [
        PhaseFeed {
            fractions: FEED.to_vec(),
            beta: 1.0,
        },
        PhaseFeed {
            fractions: FEED.to_vec(),
            beta: 1.0,
        },
    ];
    // And the VLE-initialised pair the forced-one-phase 300 K block reaches.
    let initialised = [
        PhaseFeed {
            fractions: initialisation.liquid.clone(),
            beta: 1.0 - vapour_fraction,
        },
        PhaseFeed {
            fractions: initialisation.vapour.clone(),
            beta: vapour_fraction,
        },
    ];

    // The captured `wgs-600K` overall: `equilibrium_moles[0]` is the first phase's moles, and
    // the second phase's are the same number for a split of the same composition.
    let captured_600k = CAPTURED_MOLES;
    let solution = solved(
        &mixture,
        &matrix.matrix,
        &b,
        &formation,
        600.0,
        &constructor_pair,
    );
    assert!(solution.converged, "the captured state converges");
    assert_eq!(solution.phase_amounts.len(), 2, "the multiphase branch ran");
    let overall: Vec<f64> = (0..NAMES.len())
        .map(|i| solution.phase_moles.iter().map(|phase| phase[i]).sum())
        .collect();
    for (index, (got, want)) in overall.iter().zip(captured_600k).enumerate() {
        let relative = (got - want).abs() / want.abs();
        assert!(
            relative < 1.0e-5,
            "600 K component {index} ({}): {got} against the capture's {want}, a relative \
             {relative:.3e}",
            NAMES[index]
        );
    }
    // Both phases hold the same composition, which is what makes the split a flat direction.
    for (index, (first, second)) in solution.phase_moles[0]
        .iter()
        .zip(&solution.phase_moles[1])
        .enumerate()
    {
        assert!(
            (first - second).abs() < 1.0e-9,
            "600 K component {index} ({}) is {first} in one phase and {second} in the other",
            NAMES[index]
        );
    }

    // The 300 K state, where the VLE initialisation set the two phases apart. The capture
    // prints both phases' moles, so the overall is their sum.
    let mut captured_300k = [0.0_f64; 4];
    for phase in CAPTURED_MOLES_300K {
        for (total, moles) in captured_300k.iter_mut().zip(phase) {
            *total += moles;
        }
    }
    // **The 300 K band is the measurement, and it is wider than the 600 K one**: the worst
    // component lands `1.2e-5` from the capture here against `1e-6` there, because this state's
    // own residual is `2e-5` - NeqSim stops it at the relaxed multiphase tolerance - and the
    // split it stops on is a point on the flat direction rather than a converged one.
    let solution = solved(
        &mixture,
        &matrix.matrix,
        &b,
        &formation,
        300.0,
        &initialised,
    );
    let overall: Vec<f64> = (0..NAMES.len())
        .map(|i| solution.phase_moles.iter().map(|phase| phase[i]).sum())
        .collect();
    for (index, (got, want)) in overall.iter().zip(captured_300k).enumerate() {
        let relative = (got - want).abs() / want.abs();
        assert!(
            relative < 5.0e-5,
            "300 K component {index} ({}): {got} against the capture's {want}, a relative \
             {relative:.3e}",
            NAMES[index]
        );
    }
    // The element balance, which is the constraint the whole solve rests on.
    assert!(
        solution.element_residual < 1.0e-4,
        "the element residual is {}",
        solution.element_residual
    );
}

/// One RAND solve at a state, with both phases on the vapour root - which is what the class's
/// two `gas` phases do.
fn solved(
    mixture: &azoth_eos::Mixture,
    a_matrix: &[Vec<f64>],
    b: &[f64],
    formation: &[ThermoData],
    temperature: f64,
    phases: &[azoth_reactions::rand_solver::PhaseFeed],
) -> azoth_reactions::rand_solver::RandSolution {
    use azoth_reactions::rand_solver::solve;

    let g0 = standard_potentials(formation, temperature, 1.0);
    let reduced = mixture
        .reduced_parameters(kelvins(temperature), pascals(1.0e5))
        .expect("a state");
    let mut ln_phi = |_phase: usize, x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        Ok(mixture
            .phase_state(&reduced, x, RootSide::Vapour)?
            .ln_phi
            .clone())
    };
    solve(
        a_matrix,
        &g0,
        b,
        FEED.iter().sum::<f64>(),
        phases,
        &mut ln_phi,
    )
    .expect("the solve runs")
}

/// **The captured states take the outer loop's "already at the ceiling" branch**, and the
/// stability analysis is never consulted on them.
///
/// `SystemSrkEos` hands the driver two phases and `maxPhases = 2`, so the loop solves once,
/// finds nothing to remove (both fractions are `1.0`), and accepts the phase list at the
/// ceiling - which is `skipStability` showing through: the class never reaches
/// `ReactiveStabilityAnalysis` on these fluids, and the capture's `unstable=false` came from
/// the separate probe that runs the analysis on a forced one-phase state.
#[test]
fn the_captured_state_converges_at_the_phase_ceiling() {
    use azoth_reactions::rand_solver::PhaseFeed;
    use azoth_reactions::reactive_flash::outer_loop;

    let (mixture, _) = mixture_and_constants();
    let matrix =
        FormulaMatrix::build(&NAMES.map(String::from)).expect("the components are in the databank");
    let b: Vec<f64> = matrix
        .matrix
        .iter()
        .map(|row| row.iter().zip(FEED).map(|(a, n)| a * n).sum())
        .collect();
    let g0 = standard_potentials(&formation_data(), 600.0, 1.0);
    let reduced = mixture
        .reduced_parameters(kelvins(600.0), pascals(1.0e5))
        .expect("a state");

    let pair = [
        PhaseFeed {
            fractions: FEED.to_vec(),
            beta: 1.0,
        },
        PhaseFeed {
            fractions: FEED.to_vec(),
            beta: 1.0,
        },
    ];
    let mut ln_phi = |_phase: usize, x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        Ok(mixture
            .phase_state(&reduced, x, RootSide::Vapour)?
            .ln_phi
            .clone())
    };
    let mut asked = 0_u32;
    let mut stability = |_phases: &[PhaseFeed]| -> azoth_core::Result<Vec<Vec<f64>>> {
        asked += 1;
        Ok(Vec::new())
    };

    let outcome = outer_loop(
        &matrix.matrix,
        &g0,
        &b,
        FEED.iter().sum::<f64>(),
        pair.to_vec(),
        2,
        &mut ln_phi,
        &mut stability,
    )
    .expect("the loop runs");

    assert!(outcome.converged, "the captured state converges");
    assert_eq!(asked, 0, "the stability analysis is never consulted");
    assert_eq!(outcome.phases.len(), 2, "both phases are kept");
    assert_eq!(
        outcome.total_iterations, outcome.solution.iterations,
        "one solve ran, and the loop counted its passes"
    );
    assert_eq!(
        outcome.phases[0].beta, 1.0,
        "the stale fractions are untouched"
    );

    let overall: Vec<f64> = (0..NAMES.len())
        .map(|i| {
            outcome
                .solution
                .phase_moles
                .iter()
                .map(|phase| phase[i])
                .sum()
        })
        .collect();
    for (index, (got, want)) in overall.iter().zip(CAPTURED_MOLES).enumerate() {
        let relative = (got - want).abs() / want.abs();
        assert!(
            relative < 1.0e-5,
            "component {index} ({}): {got} against the capture's {want}, a relative {relative:.3e}",
            NAMES[index]
        );
    }
}

/// **The removal step's floor and its renormalisation.** A phase under `1e-12` goes, and the
/// fractions that remain are divided by their sum. The class calls `normalizeBeta` only when
/// something was removed, which is why nothing moves on the captured states.
///
/// **No capture reaches this with a removal** - the fraction it tests is the *system's* stale
/// array, which on the captured states is `(1.0, 1.0)` or `(0.2106, 0.7894)`, both far above
/// the floor. The test is of the port's own construction and says so.
#[test]
fn the_removal_floor_and_the_renormalisation() {
    use azoth_reactions::rand_solver::PhaseFeed;
    use azoth_reactions::reactive_flash::{
        MIN_PHASE_FRACTION, normalise_betas, remove_negligible_phases,
    };

    let phase = |beta: f64| PhaseFeed {
        fractions: FEED.to_vec(),
        beta,
    };
    assert_eq!(MIN_PHASE_FRACTION, 1.0e-12);

    let mut phases = vec![phase(0.5), phase(1.0e-13), phase(0.25)];
    assert!(remove_negligible_phases(&mut phases), "one was removed");
    assert_eq!(phases.len(), 2);
    assert!(
        (phases[0].beta - 0.5 / 0.75).abs() < 1.0e-15
            && (phases[1].beta - 0.25 / 0.75).abs() < 1.0e-15,
        "the remaining fractions renormalise: {:?}",
        phases.iter().map(|phase| phase.beta).collect::<Vec<f64>>()
    );

    // Nothing to remove leaves the list alone, and the last phase is never dropped.
    let mut phases = vec![phase(0.5), phase(0.5)];
    assert!(!remove_negligible_phases(&mut phases));
    assert_eq!(phases.len(), 2);
    let mut only = vec![phase(1.0e-13)];
    assert!(!remove_negligible_phases(&mut only), "the last phase stays");
    assert_eq!(only.len(), 1);

    let mut zero = vec![phase(0.0), phase(0.0)];
    assert!(normalise_betas(&mut zero).is_err(), "a zero sum is refused");
}

/// **The trial-phase addition**: one phase, the first of the trials, at `0.01`, with every
/// fraction renormalised - and a refusal at the ceiling that leaves the list untouched.
///
/// **No capture reaches this.** The driver's own two-phase start skips the stability analysis,
/// so no trial is ever produced on the captured fluids; the class's own `addTrialPhases` is
/// guarded on `system.getMaxNumberOfPhases()`, which is what the port takes as the ceiling.
#[test]
fn the_trial_phase_is_added_at_the_guard() {
    use azoth_reactions::rand_solver::PhaseFeed;
    use azoth_reactions::reactive_flash::{TRIAL_PHASE_BETA, add_trial_phase};

    let pair = vec![
        PhaseFeed {
            fractions: FEED.to_vec(),
            beta: 1.0,
        },
        PhaseFeed {
            fractions: FEED.to_vec(),
            beta: 1.0,
        },
    ];
    let trial = vec![0.1, 0.2, 0.3, 0.4];

    let mut phases = pair.clone();
    assert!(
        add_trial_phase(&mut phases, std::slice::from_ref(&trial), 3).expect("the shapes agree"),
        "there is room for a third"
    );
    assert_eq!(phases.len(), 3);
    assert_eq!(TRIAL_PHASE_BETA, 0.01);
    let total: f64 = phases.iter().map(|phase| phase.beta).sum();
    assert!((total - 1.0).abs() < 1.0e-15, "the fractions renormalise");
    assert!(
        (phases[2].fractions[3] - 0.4).abs() < 1.0e-15,
        "the trial is kept"
    );

    // At the ceiling, and with no trials at all, the list is returned unchanged.
    let mut phases = pair.clone();
    assert!(!add_trial_phase(&mut phases, &[trial], 2).expect("the shapes agree"));
    assert_eq!(phases.len(), 2);
    let mut phases = pair;
    assert!(!add_trial_phase(&mut phases, &[], 3).expect("nothing to add"));
    assert_eq!(phases.len(), 2);
}

/// **The non-reactive fallback, on a fluid that has no reaction to run.** Methane and water
/// have rank 2 over three elements, so `NR = 0`; at 300 K and 50 bar they do not mix, and the
/// driver's `runNonReactiveFlash` is what the capture's block records. It is a conventional VLE
/// flash - Wilson K-values, then successive substitution with `K_i = phi_liq / phi_vap`.
///
/// Three things the capture pins beyond the split: the path reports `converged = true`, its
/// **`total_iterations` is `0`** because the driver never counts these passes, and the gas
/// phase is index 0 (`phase0_type=GAS`), so the *liquid* composition is written at index 1.
///
/// **This is the tightest match in the tranche**: both compositions and the vapour fraction
/// agree to `1e-12`, which is round-off, because this path is a Wilson seed and 6 successive
/// substitutions over the same cubic - there is no long Newton trajectory for two
/// implementations to drift apart on.
#[test]
fn the_non_reactive_fallback_is_the_captured_split() {
    use azoth_reactions::reactive_flash::non_reactive_flash;

    const FEED: [f64; 2] = [0.5, 0.5];
    let (mixture, constants) = mixture_and_constants_for(&["methane", "water"]);

    let reduced = mixture
        .reduced_parameters(kelvins(300.0), pascals(50.0e5))
        .expect("a state");
    let mut ln_phi = |phase: usize, x: &[f64]| -> azoth_core::Result<Vec<f64>> {
        // The class reads the phase *types*; index 0 is `GAS` in the capture, index 1 is the
        // liquid that takes the cubic's small root.
        let side = if phase == 0 {
            RootSide::Vapour
        } else {
            RootSide::Liquid
        };
        Ok(mixture.phase_state(&reduced, x, side)?.ln_phi.clone())
    };

    let outcome =
        non_reactive_flash(&FEED, &constants, 300.0, 50.0, 1, &mut ln_phi).expect("the flash runs");

    assert!(outcome.converged, "the capture has this one converged");
    // The passes are reported here and **not** by the driver: the capture's
    // `total_iterations` is `0` because `runNonReactiveFlash` never adds its own passes to the
    // counter. The port's number is the successive substitution's, and it is not zero.
    assert!(
        outcome.iterations >= 1,
        "the substitution ran {} passes",
        outcome.iterations
    );
    assert!(
        !outcome.all_supercritical,
        "water is below its critical point"
    );
    assert_eq!(outcome.phases.len(), 2);

    // `phase_x[0]` is the gas (methane-rich) and `phase_x[1]` the liquid (water-rich).
    let captured_vapour = [0.999_330_366_296_782_f64, 6.696_337_032_181_276e-4];
    let captured_liquid = [2.920_908_649_624_794_4e-7_f64, 0.999_999_707_909_135_1];
    for (index, (got, want)) in outcome.phases[0]
        .fractions
        .iter()
        .zip(captured_vapour)
        .enumerate()
    {
        assert!(
            (got - want).abs() / want.abs() < 1.0e-12,
            "vapour component {index}: {got} against the capture's {want}"
        );
    }
    for (index, (got, want)) in outcome.phases[1]
        .fractions
        .iter()
        .zip(captured_liquid)
        .enumerate()
    {
        assert!(
            (got - want).abs() / want.abs() < 1.0e-12,
            "liquid component {index}: {got} against the capture's {want}"
        );
    }

    // `beta[0]` is the gas fraction and `beta[1]` the liquid's, and they sum to one.
    let captured_vapour_fraction = 0.500_334_895_161_083_2_f64;
    assert!(
        (outcome.phases[0].beta - captured_vapour_fraction).abs() / captured_vapour_fraction
            < 1.0e-12,
        "the vapour fraction is {} against the capture's {captured_vapour_fraction}",
        outcome.phases[0].beta
    );
    assert!((outcome.phases[0].beta + outcome.phases[1].beta - 1.0).abs() < 1.0e-12);
}

/// A mixture built from an arbitrary component list, with its critical constants.
fn mixture_and_constants_for(names: &[&str]) -> (azoth_eos::Mixture, Vec<CriticalConstants>) {
    let (mixture, _) = mixture_of(names, Cubic::Srk, None).expect("the databank carries them");
    let constants = mixture
        .components()
        .iter()
        .map(|component| CriticalConstants {
            tc: component.tc.value,
            pc: component.pc.value / 1.0e5,
            omega: component.omega,
        })
        .collect();
    (mixture, constants)
}
