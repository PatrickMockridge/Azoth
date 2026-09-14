//! `eos.stability_test` - the tangent-plane stability test.
//!
//! Whether a feed at a fixed temperature and pressure is stable as a single phase. The
//! question [`crate::pt_flash`] cannot ask: successive substitution finds *a* stationary
//! point, and a flash converging to `x = y = z` has shown that its own starting point was
//! not a split, not that the feed is single phase.
//!
//! Spec: `specs/models/eos/stability_test.yaml`, which carries the criterion, the two
//! Wilson trials and what they cost, the root the feed is placed on, and the iteration
//! cap.
//!
//! Michelsen's tangent-plane distance, at each of two trial compositions: the feed is
//! unstable if either trial reaches a stationary point below the tangent plane. An
//! **unconverged trial raises** rather than being discarded - a tangent-plane distance
//! bounds stability only at a stationary point, and discarding a partial one would turn
//! "could not tell" into "stable".

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::mixture::{Mixture, PhaseState, ReducedParameters, RootSide, normalise, wilson_k};
use crate::model_gen;
use crate::pr_z_factor;

use crate::results::{StabilityTestResult, StabilityVerdict};

/// Below this the feed is unstable.
///
/// Negative rather than zero because a trivial trial - one that converges to the
/// feed - reaches `tm` of order `1e-16`, and a strict `tm < 0` would call a
/// single-phase feed unstable on rounding alone.
const TM_LIMIT: f64 = -1.0e-08;

/// A trial whose mole numbers leave this range is diverging rather than converging.
///
/// `exp` of a large `ln W` overflows to infinity and the next normalisation is a
/// division by it, so the guard is what keeps "diverged" from being reported as a
/// number.
const LOG_W_CEILING: f64 = 700.0;

/// `ln`, with a non-positive argument giving negative infinity rather than raising.
///
/// A component absent from the feed has `z_i = 0` and a reference potential of
/// `-inf`, which is the correct value: its trial mole number is zero and stays zero.
/// Spelled out rather than left to `f64::ln` because the Python reference has to be
/// (`math.log` raises, and `log(-0.0)` differs between the two), and the two sides
/// are meant to be readable as the same function.
fn ln(value: f64) -> f64 {
    if value > 0.0 {
        value.ln()
    } else {
        f64::NEG_INFINITY
    }
}

/// The feed's phase state, on whichever admissible root has the lower Gibbs energy.
///
/// The comparison is `A^R / RT - ln Z + Z` and **not** `A^R / RT + Z`; the spec's notes
/// carry the reduction and the measurement behind it. A single admissible root is the
/// common case and is taken directly.
fn feed_state(mixture: &Mixture, reduced: &ReducedParameters, z: &[f64]) -> Result<PhaseState> {
    let (a_mix, b_mix) = mixture.mixture_parameters(reduced, z);
    let roots = pr_z_factor(a_mix, b_mix)?;
    let mut candidates = vec![roots.z_min];
    if roots.z_min != roots.z_max {
        candidates.push(roots.z_max);
    }

    let mut best: Option<f64> = None;
    let mut best_gibbs = f64::INFINITY;
    for compressibility in candidates {
        let residual = mixture.helmholtz_energy(reduced, z, compressibility)?;
        let gibbs = residual - compressibility.ln() + compressibility;
        if gibbs < best_gibbs {
            best_gibbs = gibbs;
            best = Some(compressibility);
        }
    }
    // `candidates` is never empty - `z_min` is always pushed - so this is a proof
    // obligation rather than a case, and it is returned rather than unwrapped because
    // this crate does not panic on a path a caller can reach.
    let compressibility = best.ok_or_else(|| AzothError::InvalidInput {
        field: "z".to_string(),
        reason: "the cubic returned no admissible root for the feed, so there is no \
                 reference state for the tangent plane to be measured from"
            .to_string(),
    })?;
    mixture.phase_state_at(reduced, z, compressibility)
}

/// What one trial found: its composition, its distance, and the steps it took.
struct Trial {
    w: Vec<f64>,
    tm: f64,
    iterations: u32,
}

/// Iterate one trial to a stationary point of the tangent-plane distance.
///
/// `liquid` names which root the trial's phase claims: the smallest for the
/// liquid-like seed, the largest for the vapour-like one. Selecting by *ordering* is
/// the rule [`crate::pr_z_factor`] fixes, and both implementations follow it rather
/// than re-deriving a root per iteration.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap, or if a mole
///   number leaves the representable range. **Not** a trial discarded and the feed
///   called stable: see the module documentation.
/// * Propagates the mixture's and the cubic's range checks.
fn trial(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    d: &[f64],
    seed: &[f64],
    liquid: bool,
    tolerance: f64,
    cap: u32,
) -> Result<Trial> {
    let mut w = seed.to_vec();
    normalise(&mut w);
    let mut ln_w: Vec<f64> = w.iter().map(|&value| ln(value)).collect();
    let mut residual = f64::INFINITY;
    let side = if liquid {
        RootSide::Liquid
    } else {
        RootSide::Vapour
    };

    for step in 1..=cap {
        let state = mixture.phase_state(reduced, &w, side)?;
        let ln_w_new: Vec<f64> = d
            .iter()
            .zip(&state.ln_phi)
            .map(|(&di, &lp)| di - lp)
            .collect();

        if ln_w_new.iter().any(|&value| value > LOG_W_CEILING) {
            return Err(AzothError::SolverNotConverged {
                iterations: step,
                residual,
                tolerance,
            });
        }

        // The unnormalised mole numbers, which is what `tm` is written against.
        let mole_numbers: Vec<f64> = ln_w_new.iter().map(|&value| value.exp()).collect();
        let totals: f64 = mole_numbers.iter().sum();

        residual = (ln_w_new
            .iter()
            .zip(&ln_w)
            .map(|(&new, &old)| (new - old) * (new - old))
            .sum::<f64>()
            / ln_w_new.len() as f64)
            .sqrt();
        w = mole_numbers;
        normalise(&mut w);
        ln_w = ln_w_new;

        if residual <= tolerance {
            return Ok(Trial {
                w,
                tm: 1.0 - totals,
                iterations: step,
            });
        }
    }

    Err(AzothError::SolverNotConverged {
        iterations: cap,
        residual,
        tolerance,
    })
}

/// Whether a mixture at a temperature and pressure is stable as a single phase.
///
/// `z` is the overall composition, in mole fractions, and is **checked rather than
/// renormalised** - silently rescaling a caller's composition would make their error
/// invisible in a way that changes every number downstream.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `z` is the wrong length, has a negative entry,
///   or does not sum to one.
/// * [`AzothError::SolverNotConverged`] if a trial hits its cap, or if a trial's mole
///   numbers leave the representable range - which is the iteration diverging rather
///   than a state to diagnose.
/// * Propagates the kernels' range checks.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::{StabilityVerdict, databank, stability_test};
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"])
///     .expect("the pair resolves")
///     .0;
/// let r = stability_test(&mixture, kelvins(330.0), pascals(2_500_000.0), &[0.6, 0.4])?;
/// assert_eq!(r.verdict, StabilityVerdict::Unstable);
/// assert!((r.tm[1] - -0.21403386436968574).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn stability_test(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<StabilityTestResult> {
    let spec = &model_gen::STABILITY_TEST_SPEC;
    let mut warnings = Vec::new();

    let n = mixture.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!("a feed for {n} components has {} entries", z.len()),
        ));
    }
    if let Some(bad) = z.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "z[{bad}] is {} but a mole fraction cannot be negative",
                z[bad]
            ),
        ));
    }
    let sum: f64 = z.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the feed's mole fractions sum to {sum}, not to one. Renormalising it \
                 here would make a composition error invisible in every number \
                 downstream, so it is refused instead"
            ),
        ));
    }

    let min_t_over_tc = mixture
        .components()
        .iter()
        .map(|c| t.value / c.tc.value)
        .fold(f64::INFINITY, f64::min);
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "min_t_over_tc").then_some(min_t_over_tc),
        &mut warnings,
    )?;

    let reduced = mixture.reduced_parameters(t, p)?;
    warnings.extend(reduced.warnings.iter().cloned());

    let algorithm = algorithm_of(spec)?;

    // The reference potentials, from the feed on its lower-Gibbs root.
    let feed = feed_state(mixture, &reduced, z)?;
    let d: Vec<f64> = z
        .iter()
        .zip(&feed.ln_phi)
        .map(|(&zi, &lp)| ln(zi) + lp)
        .collect();

    let k = wilson_k(mixture, t, p);
    let vapour_seed: Vec<f64> = z.iter().zip(&k).map(|(&zi, &ki)| zi * ki).collect();
    let liquid_seed: Vec<f64> = z.iter().zip(&k).map(|(&zi, &ki)| zi / ki).collect();

    let mut tm = Vec::with_capacity(2);
    let mut w = Vec::with_capacity(2);
    let mut iterations = Vec::with_capacity(2);

    // Ordered vapour-like first, then liquid-like, and the order is the contract:
    // `tm` and `w` are positional and a caller reads `tm[0]` as the vapour-like
    // trial. Swapping them would silently relabel the two.
    for (seed, liquid) in [(&vapour_seed, false), (&liquid_seed, true)] {
        let outcome = trial(
            mixture,
            &reduced,
            &d,
            seed,
            liquid,
            algorithm.tolerance,
            algorithm.max_iterations,
        )?;
        w.push(outcome.w);
        tm.push(outcome.tm);
        iterations.push(outcome.iterations);
    }

    let unstable = tm.iter().any(|&distance| distance < TM_LIMIT);

    Ok(StabilityTestResult {
        verdict: if unstable {
            StabilityVerdict::Unstable
        } else {
            StabilityVerdict::Stable
        },
        tm,
        w,
        iterations,
        min_t_over_tc,
        warnings,
    })
}
