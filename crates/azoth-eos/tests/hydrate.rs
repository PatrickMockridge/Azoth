//! The hydrate fugacity, checked against the capture's own layers.
//!
//! `validation/neqsim/captures/hydrate_probe.tsv` prints the gas phase's per-component
//! fugacities at the converged state, which are exactly what NeqSim's `setFug` copies into
//! the hydrate's `reffug`. So the arithmetic can be driven from those four numbers alone and
//! compared with the occupancies and the water fugacity the same capture reports - no fluid
//! flash, and therefore no way for a fluid-side difference to hide a hydrate-side one.
//!
//! **The structure is II at every one of these states, and the capture said I for a while.**
//! The stable structure lives on the *water* component: `ComponentHydratePVTsim.fugcoef` runs
//! its two-structure selection inside the water branch and writes the winner into water's own
//! field, so a probe that read whichever component came first read methane's stale zero - and
//! every occupancy printed from it was the other structure's. `cavprwat` is the tell: `2/17`
//! and `1/17`, not `1/23` and `3/23`.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::mixture::RootSide;
use azoth_eos::{Cubic, databank, hydrate};

/// `(T in K, P in bara, [f_methane, f_ethane, f_propane, f_water] in bara, the stable
/// structure, and six occupancies as `[small: methane, ethane, propane; large: ...]`)`.
type State = (f64, f64, [f64; 4], usize, [f64; 6]);

const STATES: [State; 3] = [
    (
        293.227_691_103_840,
        100.0,
        [
            73.913_006_100_490_8,
            5.623_707_032_765_60,
            0.757_593_084_435_071,
            0.018_831_198_125_747_6,
        ],
        1,
        [
            0.883_388_125_313_154,
            0.0,
            0.0,
            0.098_687_644_680_376_3,
            0.473_812_737_652_016,
            0.423_505_120_488_215,
        ],
    ),
    (
        288.530_460_404_555,
        50.0,
        [
            39.572_555_564_715_4,
            3.877_475_364_990_52,
            0.635_251_541_845_998,
            0.013_077_400_912_688_7,
        ],
        1,
        [
            0.827_853_003_809_279,
            0.0,
            0.0,
            0.061_464_211_812_117_9,
            0.440_926_974_790_581,
            0.493_360_317_509_043,
        ],
    ),
    (
        297.061_002_920_514,
        200.0,
        [
            134.872_480_485_926,
            7.373_163_724_510_20,
            0.788_096_602_127_644,
            0.026_587_509_996_694_3,
        ],
        1,
        [
            0.923_519_664_857_266,
            0.0,
            0.0,
            0.160_891_597_042_569,
            0.493_336_214_283_545,
            0.341_939_918_179_603,
        ],
    ),
];

/// The reference water phase's fugacity at the capture's first state, in bara.
///
/// NeqSim builds that phase as a pure-water phase of the *host's* class, so it is an SRK water
/// state here - and `reference_water_fugacity` below computes it from this crate's own SRK
/// water rather than taking the number, with this constant asserting the two agree.
const REFERENCE_WATER_FUGACITY_BARA: f64 = 0.018_831_203_388_759_6;

/// The fugacity of pure water at a state, on the cubic the hydrate's fluid runs, in pascals:
/// what NeqSim's `refPhase` contributes to the coefficient.
fn reference_water_fugacity(t: f64, p: f64) -> f64 {
    let (mixture, _) = databank::mixture_of(&["water"], Cubic::Srk, None).expect("water resolves");
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("the reduction");
    let state = mixture
        .phase_state(&reduced, &[1.0], RootSide::Liquid)
        .expect("a liquid water root");
    state.ln_phi[0].exp() * p
}

/// The four substances' hydrate records, as the resolver builds them.
fn guests() -> Vec<hydrate::HydrateGuest> {
    ["methane", "ethane", "propane", "water"]
        .iter()
        .map(|name| {
            let entry = databank::entry(name, None).expect("the substance is in the databank");
            hydrate::HydrateGuest {
                name: entry.name.clone(),
                langmuir_a: entry.hydrate_langmuir_a,
                langmuir_b: entry.hydrate_langmuir_b,
                former: entry.hydrate_former,
            }
        })
        .collect()
}

#[test]
fn the_occupancies_reproduce_the_capture() {
    for (t, _p_bara, fugacities, structure, expected) in STATES {
        let records = guests();
        // The capture's fugacities are in bara, as NeqSim's internal ones are; the kernel
        // states its own in pascals and converts where the Langmuir product needs it.
        let refs: Vec<f64> = fugacities.iter().map(|value| value * 1.0e5).collect();

        for cavity in 0..2 {
            let occupied = hydrate::occupancy(&records, &refs, structure, cavity, t);
            for (index, guest) in ["methane", "ethane", "propane"].iter().enumerate() {
                let wanted = expected[cavity * 3 + index];
                assert!(
                    (occupied[index] - wanted).abs() < 1.0e-12,
                    "T = {t}: cavity {cavity} {guest} is {}, the capture says {wanted}",
                    occupied[index]
                );
            }
            assert_eq!(occupied[3], 0.0, "water occupies no cage");
        }
    }
}

#[test]
fn the_water_fugacity_coefficient_reproduces_the_capture() {
    // At the answer the coefficient is `alpha_water / P`, because the hydrate's water fugacity
    // *is* the fluid's there - which is the definition of the temperature being solved for.
    for (t, p_bara, fugacities, structure, _) in STATES {
        let records = guests();
        let refs: Vec<f64> = fugacities.iter().map(|value| value * 1.0e5).collect();

        let (found, coefficient) = hydrate::stable_structure(
            &records,
            &refs,
            t,
            p_bara * 1.0e5,
            reference_water_fugacity(t, p_bara * 1.0e5),
        )
        .expect("the sum");
        assert_eq!(found, structure, "T = {t}: the stable structure");
        let expected = fugacities[3] / p_bara;
        assert!(
            (coefficient / expected - 1.0).abs() < 1.0e-6,
            "T = {t}: the coefficient is {coefficient}, and the capture's water fugacity over \
             the pressure is {expected}"
        );
    }
}

/// The cavity sum, asserted on its own rather than only through the coefficient.
///
/// It is `-0.577685053424543` at the first state, against a chemical-potential change of
/// `+0.577684723785035` - they cancel to a part in `1e7`, which is *why* the coefficient comes
/// out as `alpha_w/P` at an answer. A port that dropped a sign in either term would still land
/// near an answer if it dropped the other one the same way, so the term is checked directly.
#[test]
fn the_cavity_sum_reproduces_the_probe() {
    let (t, _p_bara, fugacities, structure, _) = STATES[0];
    let records = guests();
    let refs: Vec<f64> = fugacities.iter().map(|value| value * 1.0e5).collect();

    let mut cavity_sum = 0.0;
    for cavity in 0..2 {
        let occupied: f64 = hydrate::occupancy(&records, &refs, structure, cavity, t)
            .iter()
            .sum();
        cavity_sum += hydrate::CAVITIES_PER_WATER[structure][cavity] * (1.0 - occupied).ln();
    }
    assert!(
        (cavity_sum + 0.577_685_053_424_543).abs() < 1.0e-12,
        "the cavity sum is {cavity_sum}, and the probe prints -0.577685053424543"
    );
}

#[test]
fn a_full_cavity_is_refused_rather_than_clamped() {
    let records = guests();
    // Fugacities large enough that the cavity sum passes one: at 100 bar the guests' are of
    // order `1e7` Pa and the Langmuir constants of order one, so `1e30` saturates.
    let saturated = vec![1.0e30, 1.0e30, 1.0e30, 1.0e30];
    let error = hydrate::stable_structure(
        &records,
        &saturated,
        293.0,
        1.0e7,
        reference_water_fugacity(293.0, 1.0e7),
    )
    .expect_err("a fully occupied cavity has no logarithm");
    assert!(matches!(error, AzothError::OutOfRange { .. }), "{error:?}");
}

/// This crate's pure-water SRK fugacity reproduces NeqSim's reference phase.
///
/// The coefficient's third term is `ln(f_w^ref/f_w^fluid)`, and `f_w^ref` is a *computed*
/// quantity rather than a constant: NeqSim builds it from the host's own equation. So the two
/// implementations have to agree on it, and this is where the capture's printed
/// `ref_water_fugacity` is held to this crate's water.
#[test]
fn the_reference_water_fugacity_reproduces_the_probe() {
    let (t, p_bara, _, _, _) = STATES[0];
    let p = p_bara * 1.0e5;
    let computed = reference_water_fugacity(t, p) / 1.0e5;
    assert!(
        (computed / REFERENCE_WATER_FUGACITY_BARA - 1.0).abs() < 1.0e-6,
        "the reference water fugacity is {computed} bar, and the probe prints \
         {REFERENCE_WATER_FUGACITY_BARA}"
    );
}
