//! The multiphase split, checked by reduction to the two-phase flash.
//!
//! A general *n*-phase solver has no registered calculation to reduce to at *n* > 2, so
//! the discipline is the one the rest of this crate uses: reduce it to the case that does.
//! At two phases the fractions and compositions must be the ones `eos.pt_flash` finds -
//! a different iteration, a different formulation and a different unknown, converging on
//! the same state.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::mixture::{RootSide, wilson_k};
use azoth_eos::multiphase::{MultiphasePhase, PhaseKind, solve_phase_fractions};

fn phases_from_wilson(
    mixture: &azoth_eos::Mixture,
    t: f64,
    p: f64,
    feed: &[f64],
) -> Vec<MultiphasePhase> {
    let k = wilson_k(mixture, kelvins(t), pascals(p));
    let vapour_raw: Vec<f64> = k.iter().zip(feed).map(|(&ki, &zi)| ki * zi).collect();
    let liquid_raw: Vec<f64> = k.iter().zip(feed).map(|(&ki, &zi)| zi / ki).collect();
    let vapour_total: f64 = vapour_raw.iter().sum();
    let liquid_total: f64 = liquid_raw.iter().sum();
    vec![
        MultiphasePhase {
            fraction: 0.5,
            composition: vapour_raw.iter().map(|v| v / vapour_total).collect(),
            kind: PhaseKind::Cubic(RootSide::Vapour),
        },
        MultiphasePhase {
            fraction: 0.5,
            composition: liquid_raw.iter().map(|v| v / liquid_total).collect(),
            kind: PhaseKind::Cubic(RootSide::Liquid),
        },
    ]
}

/// Two phases must be the two-phase flash's answer.
#[test]
fn two_phases_reduce_to_the_two_phase_flash() {
    let (mixture, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("the pair");
    let algorithm = test_algorithm();
    for (t, p_pa, feed) in [
        (330.0, 2_500_000.0, vec![0.6, 0.4]),
        (300.0, 3_000_000.0, vec![0.1, 0.9]),
        (350.0, 5_000_000.0, vec![0.5, 0.5]),
    ] {
        let reference = azoth_eos::pt_flash::pt_flash(&mixture, kelvins(t), pascals(p_pa), &feed)
            .expect("the flash computes");
        let reduced = mixture
            .reduced_parameters(kelvins(t), pascals(p_pa))
            .expect("a state");
        let mut phases = phases_from_wilson(&mixture, t, p_pa, &feed);
        let split = solve_phase_fractions(&mixture, &reduced, &feed, &mut phases, &algorithm)
            .unwrap_or_else(|e| panic!("T={t}, P={p_pa}: {e}"));

        let beta = split.fractions[0];
        let reference_beta = reference.beta.expect("a split has a vapour fraction");
        // Only where the flash's own answer is a *split*. A subcooled or superheated feed
        // has a negative or above-one vapour fraction, which is the negative flash and
        // not a phase amount; the multiphase solver has no negative phase to report and
        // drives that phase to the floor instead - a difference in what the two can say,
        // not in where they converge.
        if !(0.0..=1.0).contains(&reference_beta) {
            assert!(
                beta <= 1.0e-9 || beta >= 1.0 - 1.0e-9,
                "T={t}, P={p_pa}: the flash reports the negative flash {reference_beta}, so \
                 the absent phase should be at the floor and it is {beta}"
            );
            continue;
        }
        assert!(
            (beta - reference_beta).abs() < 1e-6,
            "T={t}, P={p_pa}: the multiphase solve gives {beta} and the flash {reference_beta}"
        );
        // The vapour phase is the first, and its composition is the flash's.
        for i in 0..feed.len() {
            assert!(
                (split.compositions[0][i] - reference.y[i]).abs() < 1e-5,
                "T={t}, P={p_pa}: y[{i}] is {} and the flash gives {}",
                split.compositions[0][i],
                reference.y[i]
            );
        }
    }
}

/// The phases account for the feed, to the normalisation's cost.
///
/// `sum_k beta_k x_ik = z_i` is an identity of the closed form `x_ik = z_i/(E_i phi_ik)`,
/// and the iteration then rescales each phase to sum to one - which the identity does not
/// survive exactly. Eleven digits is what that costs, measured, and the tolerance is that
/// rather than the arithmetic's last bit.
#[test]
fn the_phases_account_for_the_feed() {
    let (mixture, _) = databank::mixture_of(&["methane", "propane", "n-butane"], Cubic::Pr, None)
        .expect("the ternary");
    let algorithm = test_algorithm();
    let feed = vec![0.5, 0.3, 0.2];
    let (t, p_pa) = (330.0, 4_000_000.0);
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p_pa))
        .expect("a state");
    let mut phases = phases_from_wilson(&mixture, t, p_pa, &feed);
    let split = solve_phase_fractions(&mixture, &reduced, &feed, &mut phases, &algorithm)
        .expect("the solve converges");
    for i in 0..feed.len() {
        let accounted: f64 = split
            .fractions
            .iter()
            .zip(&split.compositions)
            .map(|(&beta, x)| beta * x[i])
            .sum();
        // The normalisation's cost and not the form's: the closed `x_ik = z_i/(E_i phi_ik)`
        // balances exactly, and rescaling each phase to sum to one is what these eleven
        // digits pay for.
        assert!(
            (accounted - feed[i]).abs() < 1e-10,
            "component {i}: the phases hold {accounted} and the feed has {}",
            feed[i]
        );
    }
    let total: f64 = split.fractions.iter().sum();
    assert!((total - 1.0).abs() < 1e-12, "the fractions sum to {total}");
}

/// The algorithm a test runs, since the solver reads its stopping rule from one.
fn test_algorithm() -> azoth_core::ModelAlgorithm {
    azoth_core::ModelAlgorithm {
        scheme: "multiphase_fraction_newton",
        convergence: "absolute",
        tolerance: 1e-10,
        max_iterations: 100,
        bracket: None,
        initialisation: None,
        initial_temperature: None,
        inner: None,
        fallback: None,
    }
}
