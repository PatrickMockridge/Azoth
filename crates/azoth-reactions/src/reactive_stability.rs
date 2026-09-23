//! The reactive tangent-plane stability analysis, from
//! `flashops/reactiveflash/ReactiveStabilityAnalysis.java`.
//!
//! It decides whether a second phase forms, and it does so in four steps: bring the feed to
//! **homogeneous** chemical equilibrium (a single-phase [`crate::rand_solver`] solve), take
//! the reference potentials `d_i = ln x_i + ln phi_i` from *that* composition, seed trial
//! phases with Wilson K-values and pure components, and run a tangent-plane trial from each.
//! A trial whose distance comes in under [`TPD_THRESHOLD`] is a phase the fluid is unstable
//! with respect to.
//!
//! # The trial is a successive substitution on the log-composition
//!
//! ```text
//! logW_i = d_i - ln phi_i(W)          unnormalised, iterated to SS_TOL
//! TPD    = 1 - sum_i exp(logW_i)      the tangent-plane distance at the stationary point
//! ```
//!
//! **Two answers are not distances and the class reports both as the same number.** A trial
//! that walks back to the feed's own composition is the *trivial* solution, and a trial whose
//! phase cannot be initialised is not a distance either; both come back as
//! [`STABLE_TPD`] = `10.0`, which is how the caller reads "this seed found nothing". The
//! captured states are all of that kind - every one of their six seeds returns `10.0` - so
//! the sentinel is not an edge case here but the whole of the answer.

use azoth_core::{AzothError, Result};

/// The successive-substitution tolerance, from `SS_TOL`.
pub const SS_TOL: f64 = 1.0e-9;

/// The pass cap, from `MAX_SS_ITER`.
pub const MAX_SS_ITERATIONS: usize = 100;

/// The distance below which a trial counts as a new phase, from `TPD_THRESHOLD`.
pub const TPD_THRESHOLD: f64 = -1.0e-8;

/// The floor a trial mole fraction is kept above, from `MIN_MOLES`.
pub const MIN_MOLES: f64 = 1.0e-30;

/// What a trial reports when it is stable, trivial, or cannot be initialised.
pub const STABLE_TPD: f64 = 10.0;

/// The composition difference below which a trial is the reference's own solution, from
/// `trivialCheck < 1e-4`.
pub const TRIVIAL_TOLERANCE: f64 = 1.0e-4;

/// One component's critical constants, which the Wilson seed reads.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CriticalConstants {
    /// The critical temperature, in K.
    pub tc: f64,
    /// The critical pressure, in bara - the unit the Wilson correlation is written in.
    pub pc: f64,
    /// The acentric factor.
    pub omega: f64,
}

/// The Wilson K-value `(pc / P) exp(5.373 (1 + omega) (1 - tc / T))`, floored at `1e-20`.
///
/// **A component with no critical constants gets a K of 1**, which is the class's own
/// fallback and the reason a fluid of pseudo-components still produces trials - they are just
/// not separation-shaped ones.
#[must_use]
pub fn wilson_k(constants: &CriticalConstants, temperature: f64, pressure: f64) -> f64 {
    if constants.pc > 0.0 && constants.tc > 0.0 {
        ((constants.pc / pressure)
            * (5.373 * (1.0 + constants.omega) * (1.0 - constants.tc / temperature)).exp())
        .max(1.0e-20)
    } else {
        1.0
    }
}

/// The trial seeds `generateTrialPhases` builds, unnormalised and in its own order:
///
/// 1. a liquid-like trial `z_i / K_i`, and
/// 2. a vapour-like trial `K_i z_i`,
///
/// both **only when some `|ln K|` exceeds `0.01`**, and then
/// 3. one pure-component trial per component - `1.0` on that component and `1e-12` on the
///    rest - which is built for every component the feed holds.
///
/// The returned vectors are what the trial *starts* from, not fractions: the class normalises
/// them as it sets composition, and the capture prints them raw.
#[must_use]
#[allow(clippy::needless_range_loop)] // the class's own loop over the components
pub fn trial_seeds(
    fractions: &[f64],
    constants: &[CriticalConstants],
    temperature: f64,
    pressure: f64,
) -> Vec<Vec<f64>> {
    let nc = fractions.len();
    let k: Vec<f64> = constants
        .iter()
        .map(|entry| wilson_k(entry, temperature, pressure))
        .collect();
    let all_near_one = k.iter().all(|value| value.ln().abs() <= 0.01);

    let mut seeds: Vec<Vec<f64>> = Vec::new();
    if !all_near_one {
        seeds.push(
            (0..nc)
                .map(|i| {
                    if fractions[i] > MIN_MOLES {
                        fractions[i] / k[i]
                    } else {
                        MIN_MOLES
                    }
                })
                .collect(),
        );
        seeds.push(
            (0..nc)
                .map(|i| {
                    if fractions[i] > MIN_MOLES {
                        k[i] * fractions[i]
                    } else {
                        MIN_MOLES
                    }
                })
                .collect(),
        );
    }
    for j in 0..nc {
        // The class skips a component the feed does not hold and pins the ions out; this
        // crate's caller has already refused an ionic fluid, so the presence test is the one
        // that bites.
        if fractions[j] > 1e-100 {
            seeds.push(
                (0..nc)
                    .map(|i| if i == j { 1.0 } else { 1.0e-12 })
                    .collect(),
            );
        }
    }
    seeds
}

/// The value `computeReferencePotentials` gives a component whose mole fraction is at or below
/// [`MIN_MOLES`]: "effectively absent", and **not** a logarithm of the floor.
pub const REFERENCE_ABSENT: f64 = -100.0;

/// The value it gives an ion, which does not participate in phase equilibrium.
pub const REFERENCE_ION: f64 = -1000.0;

/// The reference potentials `d_i = ln x_i + ln phi_i` of the equilibrated feed, from
/// `computeReferencePotentials`.
///
/// Three branches, and the middle one is easy to get wrong: a component the feed holds takes
/// the logarithm, one at or below [`MIN_MOLES`] takes [`REFERENCE_ABSENT`], and an ion takes
/// [`REFERENCE_ION`]. Writing the absent case as `ln(max(x, MIN_MOLES))` is a *different
/// number* - `-69.08` and not `-100` - and the difference does not show on any fluid whose
/// components are all present.
///
/// **The absent branch is not reached by any captured fluid.** The probe's trace-nitrogen case
/// carries `1e-40` and the class's own floor leaves the phase holding `1.0001516e-30`, which is
/// *above* `MIN_MOLES`, so the logarithm runs there and the capture pins that instead. `charges`
/// is the per-component ionic charge, zero for the neutral fluids this crate's callers admit.
///
/// # Errors
/// [`AzothError::InvalidInput`] on a shape disagreement.
pub fn reference_potentials(
    fractions: &[f64],
    ln_phi: &[f64],
    charges: &[f64],
) -> Result<Vec<f64>> {
    if fractions.len() != ln_phi.len() || fractions.len() != charges.len() {
        return Err(AzothError::InvalidInput {
            field: "fractions".to_string(),
            reason: format!(
                "{} composition entr(ies) against {} fugacity coefficient(s) and {} charge(s)",
                fractions.len(),
                ln_phi.len(),
                charges.len()
            ),
        });
    }
    Ok((0..fractions.len())
        .map(|i| {
            if charges[i] != 0.0 {
                // The ion override comes last in the class, so it wins over the logarithm.
                return REFERENCE_ION;
            }
            let x = fractions[i];
            if x > MIN_MOLES {
                x.ln() + ln_phi[i]
            } else {
                REFERENCE_ABSENT
            }
        })
        .collect())
}

/// One tangent-plane trial, returning its distance or [`STABLE_TPD`].
///
/// `ln_phi` is the same closure [`crate::rand_solver::solve`] takes: a mole-fraction
/// composition in, each component's `ln(phi_i)` out.
#[allow(clippy::needless_range_loop)] // the same index loops the class writes
pub fn run_trial(
    d: &[f64],
    seed: &[f64],
    fractions: &[f64],
    ln_phi: &mut dyn FnMut(&[f64]) -> Result<Vec<f64>>,
) -> Result<f64> {
    let nc = d.len();
    let mut log_w = vec![0.0_f64; nc];
    let mut previous = vec![0.0_f64; nc];
    let mut sum_w = 0.0_f64;
    for i in 0..nc {
        let w = seed[i].max(MIN_MOLES);
        log_w[i] = w.ln();
        sum_w += w;
    }

    let mut fractions_w: Vec<f64> = (0..nc).map(|i| log_w[i].exp() / sum_w).collect();

    for _ in 0..MAX_SS_ITERATIONS {
        previous.copy_from_slice(&log_w);
        let ln_phi_here = ln_phi(&fractions_w)?;

        let mut error = 0.0_f64;
        sum_w = 0.0;
        for i in 0..nc {
            log_w[i] = d[i] - ln_phi_here[i];
            error += (log_w[i] - previous[i]).abs();
            sum_w += log_w[i].exp();
        }
        if !sum_w.is_finite() || sum_w <= 0.0 {
            return Ok(STABLE_TPD);
        }
        for i in 0..nc {
            fractions_w[i] = log_w[i].exp() / sum_w;
        }
        if error < SS_TOL {
            break;
        }
    }

    // `TPD = 1 - sum_i exp(logW_i)`.
    let mut tpd = 1.0_f64;
    for i in 0..nc {
        tpd -= log_w[i].exp();
    }

    // A trial that walked back to the reference is the reference's own solution and not a
    // second phase.
    let trivial: f64 = fractions_w
        .iter()
        .zip(fractions)
        .map(|(w, x)| (w - x).abs())
        .sum();
    if trivial < TRIVIAL_TOLERANCE {
        return Ok(STABLE_TPD);
    }
    Ok(tpd)
}

/// Whether a distance counts as a new phase, and the composition to start from when it does.
#[must_use]
pub fn is_unstable(tpd: f64) -> bool {
    tpd < TPD_THRESHOLD
}
