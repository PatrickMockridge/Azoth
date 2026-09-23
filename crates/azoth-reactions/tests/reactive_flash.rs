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
