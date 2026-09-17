//! `eos.pvf_flash` - the temperature at which a feed's vapour fraction is a given value.
//!
//! Spec: `specs/models/eos/pvf_flash.toml`. NeqSim's `PVFflash`.
//!
//! At a fixed pressure the vapour fraction rises monotonically through the two-phase
//! region, from zero at the bubble point to one at the dew point, so a specified
//! fraction is a temperature. The search is Illinois' - a bracketing regula falsi - and
//! it needs no derivative of the vapour fraction, which is what makes it robust where a
//! Newton on a *flashed* quantity is not: the fraction is a maximum of zero and one at
//! the ends of the region, so its slope vanishes exactly where the iteration has to
//! cross out of one phase into the other.

use azoth_core::units::{Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::mixture::Mixture;
use crate::model_gen;
use crate::pt_flash::pt_flash;
use crate::results::{Phase, PtFlashResult, PvfFlashResult};

/// Where the bracket may be widened to, in kelvin, and by how much each step.
///
/// Upstream's `bracketAttempts` loop: ten kelvin a side, twenty times, floored at 50 K
/// and capped at 2000 K. It is what lets a feed whose two-phase region lies well away
/// from its own temperature be reached rather than reported as unanswerable.
const WIDEN_STEP: f64 = 10.0;
const WIDEN_LIMIT: u32 = 20;
const COLDEST: f64 = 50.0;
const HOTTEST: f64 = 2000.0;

/// The vapour fraction at a temperature, which is what the search brackets on.
///
/// A single-phase feed answers one or zero rather than an absence: the search is *for* a
/// fraction, and a trial outside the two-phase region has to return something ordered
/// against the specification for the bracket to move.
fn beta_at(mixture: &Mixture, p: Pressure, feed: &[f64], t: f64) -> Result<(f64, PtFlashResult)> {
    let flash = pt_flash(mixture, kelvins(t), p, feed)?;
    let beta = match flash.phase {
        Phase::TwoPhase => flash.beta.unwrap_or(0.0),
        Phase::AllVapour | Phase::Trivial => 1.0,
        Phase::AllLiquid => 0.0,
    };
    Ok((beta, flash))
}

/// The temperature at which a feed's vapour fraction at a pressure is `beta`.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `P` is not positive, or if `beta` is not strictly
///   inside `(0, 1)`.
/// * [`AzothError::SolverNotConverged`] if the search brackets nothing, or reaches its
///   cap with the fraction still short of the one asked for.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::databank;
/// use azoth_eos::pvf_flash::pvf_flash;
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves")
///     .0;
/// let r = pvf_flash(
///     &mixture,
///     pascals(2_500_000.0),
///     0.8422055475803881,
///     kelvins(330.0),
///     &[0.6, 0.4],
/// )?;
/// assert!((r.t.value - 330.0).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pvf_flash(
    mixture: &Mixture,
    p: Pressure,
    beta: f64,
    feed_temperature: ThermodynamicTemperature,
    feed: &[f64],
) -> Result<PvfFlashResult> {
    let spec = &model_gen::PVF_FLASH_SPEC;
    let mut warnings = Vec::new();

    // Checked before the range checks, because the spec's own bound on `beta` is the
    // same interval and would answer first with a message that names a bound rather
    // than the two models whose job the endpoints are.
    if beta <= 0.0 || beta >= 1.0 {
        return Err(AzothError::InvalidInput {
            field: "beta".to_string(),
            reason: format!(
                "a vapour fraction of {beta} is the {} point, and those are \
                 `eos.bubble_temperature` and `eos.dew_temperature` - calculations with \
                 their own procedures. This model solves the interior of the two-phase \
                 region, which is the part they do not",
                if beta <= 0.0 { "bubble" } else { "dew" }
            ),
        });
    }
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "beta" => Some(beta),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;

    let bracket = algorithm.bracket.ok_or_else(|| AzothError::InvalidInput {
        field: "algorithm.bracket".to_string(),
        reason: format!(
            "scheme `{}` brackets the vapour fraction by temperature but the spec \
             declares no bracket",
            algorithm.scheme
        ),
    })?;
    let centre = feed_temperature.value;
    let (mut cold, mut hot) = (bracket.lower * centre, bracket.upper * centre);

    let mut cold_beta = beta_at(mixture, p, feed, cold)?.0;
    let mut hot_beta = beta_at(mixture, p, feed, hot)?.0;
    // Widen until the specification is inside, which is what makes the search a
    // bracketing one rather than a descent from a guess.
    let mut attempts = 0;
    while cold_beta > beta && attempts < WIDEN_LIMIT {
        cold = (cold - WIDEN_STEP).max(COLDEST);
        cold_beta = beta_at(mixture, p, feed, cold)?.0;
        attempts += 1;
    }
    attempts = 0;
    while hot_beta < beta && attempts < WIDEN_LIMIT {
        hot = (hot + WIDEN_STEP).min(HOTTEST);
        hot_beta = beta_at(mixture, p, feed, hot)?.0;
        attempts += 1;
    }
    if !(cold_beta <= beta && beta <= hot_beta) {
        return Err(AzothError::SolverNotConverged {
            iterations: 0,
            residual: beta - cold_beta.min(hot_beta),
            tolerance: algorithm.tolerance,
        });
    }

    // Illinois' method: regula falsi, with the end that has not moved halved before the
    // next step. Plain regula falsi converges from one side only and crawls when the
    // function is curved; halving the stale end restores the bisection's guarantee
    // without giving up the secant's speed.
    let (mut t_a, mut t_b) = (cold, hot);
    let (mut f_a, mut f_b) = (cold_beta - beta, hot_beta - beta);
    let mut iterations = 0;
    let mut residual = f64::NAN;
    let mut answer = f64::NAN;

    for step in 1..=algorithm.max_iterations {
        iterations = step;
        let t_c = (t_a - f_a * (t_b - t_a) / (f_b - f_a)).clamp(COLDEST, HOTTEST);
        let f_c = beta_at(mixture, p, feed, t_c)?.0 - beta;
        residual = f_c.abs();
        if residual < algorithm.tolerance {
            answer = t_c;
            break;
        }
        if f_c * f_b < 0.0 {
            t_a = t_c;
            f_a = f_c;
        } else {
            t_b = t_c;
            f_b = f_c;
            f_a *= 0.5;
        }
        if (t_b - t_a).abs() < 1.0e-10 {
            answer = t_c;
            break;
        }
    }
    if answer.is_nan() {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual,
            tolerance: algorithm.tolerance,
        });
    }

    let (found, flash) = beta_at(mixture, p, feed, answer)?;
    warnings.extend(flash.warnings.iter().cloned());
    Ok(PvfFlashResult {
        t: kelvins(answer),
        beta: found,
        phase: flash.phase,
        x: flash.x,
        y: flash.y,
        k: flash.k,
        z_liquid: flash.z_liquid,
        z_vapour: flash.z_vapour,
        iterations,
        residual: (found - beta).abs(),
        warnings,
    })
}
