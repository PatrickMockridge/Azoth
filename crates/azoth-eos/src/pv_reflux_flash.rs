//! `eos.pv_reflux_flash` - the temperature at which a phase ratio is a given value.
//!
//! Spec: `specs/models/eos/pv_reflux_flash.toml`. NeqSim's `PVrefluxflash`.
//!
//! A distillation column's condenser fixes how much liquid it returns against how much
//! it takes off, and a reboiler fixes the same ratio the other way round. Both are a
//! ratio of the two phase amounts, and at a fixed pressure that ratio moves with the
//! temperature, so a specified ratio is a temperature.
//!
//! The derivative is a secant over the last two iterates rather than a slope of the
//! flashed ratio, which is what makes the first step a *probe*: there is no derivative
//! until there are two iterates, and inventing one from the first residual alone would
//! take a step sized by the temperature itself.

use azoth_core::units::{Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::mixture::Mixture;
use crate::model_gen;
use crate::pt_flash::pt_flash;
use crate::results::{PtFlashResult, PvRefluxFlashResult};

/// The largest temperature step one iteration may take, in kelvin.
///
/// Upstream's two-kelvin cap. It is what keeps a secant taken from two nearly-equal
/// residuals from throwing the iterate across the two-phase region.
const MAX_STEP: f64 = 2.0;

/// How far the first iteration probes, in kelvin, before any derivative exists.
const PROBE: f64 = 0.1;

/// Which phase the reflux ratio is of.
///
/// A condenser asks for the vapour's - the liquid it returns over the vapour it takes -
/// and a reboiler for the liquid's. The two ratios are reciprocals, so the choice is not
/// cosmetic: it decides which end of the column's operating line the answer sits on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RefluxPhase {
    /// The vapour phase: the ratio is `(1 - beta) / beta`.
    Vapour,
    /// The liquid phase: the ratio is `beta / (1 - beta)`.
    Liquid,
}

impl RefluxPhase {
    /// The ratio of the phases at a vapour fraction.
    ///
    /// `1 / beta_phase - 1`, which is upstream's expression and is the other phase's
    /// amount over this one's.
    fn ratio(self, beta: Option<f64>) -> f64 {
        // A single-phase feed has no ratio: one of the amounts is zero, and upstream
        // reads the phase fraction as exactly one or zero there rather than refusing.
        // `1/1 - 1 = 0` for an absent second phase, and the same at the other end.
        let fraction = match (self, beta) {
            (Self::Vapour, Some(b)) => b,
            (Self::Vapour, None) => 1.0,
            (Self::Liquid, Some(b)) => 1.0 - b,
            (Self::Liquid, None) => 1.0,
        };
        1.0 / fraction - 1.0
    }
}

/// The temperature at which a phase ratio at a pressure is `reflux`.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `P` or the starting temperature is not positive.
/// * [`AzothError::SolverNotConverged`] if the iteration reaches its cap without the
///   ratio meeting the one asked for - a ratio the two-phase region does not contain at
///   this pressure.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::databank;
/// use azoth_eos::pv_reflux_flash::RefluxPhase;
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves")
///     .0;
/// let r = azoth_eos::pv_reflux_flash::pv_reflux_flash(
///     &mixture,
///     pascals(2_500_000.0),
///     0.1873586001338361,
///     RefluxPhase::Vapour,
///     kelvins(330.0),
///     &[0.6, 0.4],
/// )?;
/// assert!((r.t.value - 330.0).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pv_reflux_flash(
    mixture: &Mixture,
    p: Pressure,
    reflux: f64,
    phase_of: RefluxPhase,
    feed_temperature: ThermodynamicTemperature,
    feed: &[f64],
) -> Result<PvRefluxFlashResult> {
    let spec = &model_gen::PV_REFLUX_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "temperature" => Some(feed_temperature.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let at = |t: f64| -> Result<PtFlashResult> { pt_flash(mixture, kelvins(t), p, feed) };

    // Upstream's own arrangement of the three iterates: `t_old` holds the temperature
    // the flash on hand was taken at, `t_older` the one before it, and the secant is
    // over the two. The shift happens before `t_old` is read, which is what makes the
    // first iteration's difference `T_start - 0` - and why the first step is a probe
    // rather than a secant.
    let mut current = feed_temperature.value;
    // Declared uninitialised: `t_old`'s and `f`'s first values are upstream's zeros, and
    // the loop reaches `t_older`/`f_old` through them before writing either.
    let mut t_old = 0.0;
    let mut t_older;
    let mut f_old;
    let mut f = 0.0;
    let mut iterations = 0;
    let mut flash = at(current)?;

    for step in 1..=algorithm.max_iterations {
        iterations = step;
        f_old = f;
        t_older = t_old;
        t_old = current;
        f = reflux - phase_of.ratio(flash.beta);

        let measured = t_old - t_older;
        let slope = if measured.abs() > 1.0e-12 {
            (f - f_old) / measured
        } else {
            0.0
        };

        if f.abs() <= algorithm.tolerance {
            break;
        }
        if step >= algorithm.max_iterations {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual: f.abs(),
                tolerance: algorithm.tolerance,
            });
        }

        current = if step < 2 || slope == 0.0 || !slope.is_finite() {
            // No derivative yet, or a degenerate one: step a fixed tenth of a kelvin
            // *away* from the residual's sign, so the second iterate brackets it.
            let probe = if f > 0.0 {
                PROBE
            } else if f < 0.0 {
                -PROBE
            } else {
                0.0
            };
            t_old + probe
        } else {
            let full = (f / slope).clamp(-MAX_STEP, MAX_STEP);
            let damping = (0.4 + 0.06 * f64::from(step)).min(1.0);
            t_old - full * damping
        };
        flash = at(current)?;
    }

    warnings.extend(flash.warnings.iter().cloned());
    let residual = (reflux - phase_of.ratio(flash.beta)).abs();
    Ok(PvRefluxFlashResult {
        t: kelvins(current),
        beta: flash.beta,
        phase: flash.phase,
        x: flash.x,
        y: flash.y,
        k: flash.k,
        z_liquid: flash.z_liquid,
        z_vapour: flash.z_vapour,
        iterations,
        residual,
        warnings,
    })
}
