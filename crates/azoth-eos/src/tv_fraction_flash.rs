//! `eos.tv_fraction_flash` - the pressure at which a gas volume fraction is a given value.
//!
//! Spec: `specs/models/eos/tv_fraction_flash.toml`. NeqSim's `TVfractionFlash`.
//!
//! A *volume* fraction, not a mole fraction: the gas phase's share of the mixture's
//! volume. It is what ASTM D6377's vapour pressure is defined on, where the four-to-one
//! ratio is volumes and not moles, and at a fixed temperature it rises monotonically as
//! the pressure falls - so a specified fraction is a pressure.
//!
//! The residual's derivative is a central difference of the residual itself rather than
//! the quotient rule upstream writes. That rule mixes a volume-*corrected* residual with
//! uncorrected derivatives of the mixture and the gas phase separately; differencing the
//! quantity the iteration actually drives is self-consistent, and it costs two more
//! flashes a step.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::algorithm_of;
use crate::mixture::Mixture;
use crate::model_gen;
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::pt_flash::pt_flash;
use crate::results::{Phase, PtFlashResult, TvFractionFlashResult};

/// How far the pressure is walked down when the feed is single phase, and how many times.
///
/// Upstream's preamble: a feed that is all liquid at the starting pressure has no gas to
/// take a fraction of, so the pressure is lowered by a tenth, up to twenty times, until
/// one appears. Without it the residual is a constant and the Newton has nothing to
/// follow.
const WALK_FACTOR: f64 = 0.9;
const WALK_LIMIT: u32 = 20;

/// The damping's floor and its start, upstream's adaptive parameters.
///
/// The first step is `1/(1 + 100)` of the Newton step - one per cent - and the damping
/// halves towards twenty as the error improves. It is what keeps a Newton taken from a
/// pressure far from the answer from crossing the two-phase region in one step.
const DAMPING_START: f64 = 100.0;
const DAMPING_FLOOR: f64 = 20.0;

/// The largest pressure step one iteration may take, and the fraction of the current
/// pressure beyond which it is capped anyway.
///
/// **Upstream's ten is ten *bar*** - NeqSim holds pressures in bar and this crate holds
/// them in pascals - so the constant is written in the unit the cap is expressed in and
/// converted once, here. Written as `10.0` in pascals it would cap every step at ten
/// pascals and the iteration would crawl towards an answer it could never reach.
const MAX_STEP_BAR: f64 = 10.0;
const MAX_STEP_PA: f64 = MAX_STEP_BAR * 1.0e5;
const MAX_STEP_FRACTION: f64 = 0.5;

/// The fewest steps the iteration takes whatever the residual does.
const MINIMUM_STEPS: u32 = 6;

/// The gas phase's share of the mixture's volume, at a pressure.
///
/// The phase volumes are the cubic's own, `z R T / P`. Upstream's default applies the
/// Peneloux translation as well, `v - sum_i x_i c_i`; this does not, because nothing
/// wires a translation onto a databank-built mixture on either side - `eos.pr_peneloux_shift`
/// computes one, and no caller sets it. The spec's assumption records the difference.
fn volume_fraction(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    feed: &[f64],
) -> Result<(f64, PtFlashResult)> {
    let flash = pt_flash(mixture, t, p, feed)?;
    let (beta, x, y) = match flash.beta {
        Some(beta) => (beta, flash.x.clone(), flash.y.clone()),
        // A single-phase feed has no gas to take a share of. Upstream reads the phase
        // fraction as one or zero there, and so does this: the residual is then a
        // constant, which is what tells the iteration the outlet is single phase.
        None => {
            let one = if flash.phase == Phase::AllLiquid {
                0.0
            } else {
                1.0
            };
            return Ok((one, flash));
        }
    };
    let reduced = mixture.reduced_parameters(t, p)?;
    let v_liquid = mixture
        .phase_state(&reduced, &x, crate::mixture::RootSide::Liquid)?
        .z
        * MOLAR_GAS_CONSTANT
        * t.value
        / p.value;
    let v_vapour = mixture
        .phase_state(&reduced, &y, crate::mixture::RootSide::Vapour)?
        .z
        * MOLAR_GAS_CONSTANT
        * t.value
        / p.value;
    let total = (1.0 - beta) * v_liquid + beta * v_vapour;
    Ok((beta * v_vapour / total, flash))
}

/// The pressure at which a feed's gas volume fraction at a temperature is `fraction`.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::InvalidInput`] if `fraction` is not strictly inside `(0, 1)`, or the
///   feed is the wrong length.
/// * [`AzothError::SolverNotConverged`] if the iteration reaches its cap without the
///   fraction meeting the one asked for.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::databank;
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves")
///     .0;
/// let r = azoth_eos::tv_fraction_flash::tv_fraction_flash(
///     &mixture,
///     kelvins(330.0),
///     0.9,
///     pascals(2_500_000.0),
///     &[0.6, 0.4],
/// )?;
/// assert!(r.pressure.value > 0.0);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tv_fraction_flash(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    fraction: f64,
    start_pressure: Pressure,
    feed: &[f64],
) -> Result<TvFractionFlashResult> {
    let spec = &model_gen::TV_FRACTION_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(start_pressure.value),
            _ => None,
        },
        &mut warnings,
    )?;
    if fraction <= 0.0 || fraction >= 1.0 {
        return Err(AzothError::invalid_input(
            "fraction",
            format!(
                "a volume fraction of {fraction} is not inside (0, 1). At zero the whole \
                 mixture is liquid and at one it is all gas, and neither has a \
                 two-phase pressure to find"
            ),
        ));
    }

    let algorithm = algorithm_of(spec)?;
    let mut pressure = start_pressure.value;
    let (mut found, mut flash) = volume_fraction(mixture, t, pascals(pressure), feed)?;

    // The preamble: a feed that is single phase at the starting pressure has no gas to
    // take a share of, so the pressure is walked down until one appears. Upstream walks
    // twenty times and gives up; giving up here is an error naming the state rather than
    // an answer at whatever the last pressure was.
    let mut attempts = 0;
    while flash.beta.is_none() && attempts < WALK_LIMIT {
        pressure *= WALK_FACTOR;
        attempts += 1;
        (found, flash) = volume_fraction(mixture, t, pascals(pressure), feed)?;
    }
    if flash.beta.is_none() {
        return Err(AzothError::OutOfRange {
            field: "fraction".to_string(),
            value: fraction,
            detail: format!(
                "the feed is single phase at every pressure from the {} asked for down to \
                 {pressure:.6e} Pa, so there is no volume fraction to solve for. A state \
                 with no gas phase at all has no gas volume share",
                start_pressure.value
            ),
        });
    }

    let mut error = 100.0_f64;
    let mut error_old = error;
    let mut damping = DAMPING_START;
    let mut iterations = 0;

    for step in 1..=algorithm.max_iterations {
        iterations = step;
        // The derivative is a central difference of the residual the iteration drives,
        // taken at a step relative to the pressure so it is the same fraction of the
        // state in either implementation.
        let h = 1.0e-4 * pressure;
        let above = volume_fraction(mixture, t, pascals(pressure + h), feed)?.0;
        let below = volume_fraction(mixture, t, pascals((pressure - h).max(1.0)), feed)?.0;
        let slope = (above - below) / (2.0 * h);
        if !slope.is_finite() || slope.abs() < 1.0e-30 {
            return Err(AzothError::SolverNotConverged {
                iterations,
                residual: error,
                tolerance: algorithm.tolerance,
            });
        }

        if step > 3 && error < error_old * 0.9 {
            damping = (damping * 0.9).max(DAMPING_FLOOR);
        }
        let step_fraction = f64::from(step) / (f64::from(step) + damping);
        let mut candidate = pressure - step_fraction * (found - fraction) / slope;
        if candidate <= 0.0 {
            candidate = pressure * 0.9;
        }
        let cap = MAX_STEP_PA.min(MAX_STEP_FRACTION * pressure.abs());
        if (candidate - pressure).abs() > cap {
            candidate = pressure + (candidate - pressure).signum() * cap;
        }

        let moved = (candidate - pressure).abs();
        pressure = candidate;
        (found, flash) = volume_fraction(mixture, t, pascals(pressure), feed)?;
        error_old = error;
        error = (found - fraction).abs();

        if error < 1.0e-8 && step > 3 {
            break;
        }
        let done = error <= algorithm.tolerance && moved <= 1.0e-6;
        if (done && step >= MINIMUM_STEPS) || step >= algorithm.max_iterations {
            break;
        }
    }

    if error > 1.0e-4 {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: error,
            tolerance: algorithm.tolerance,
        });
    }
    if error > algorithm.tolerance {
        warnings.push(azoth_core::Warning::new(
            azoth_core::WarningCode::SolverNotConverged,
            format!(
                "the volume fraction at the answer is {found}, against the {fraction} \
                 asked for. The iteration left by its cap rather than by its tolerance, \
                 and upstream accepts this to 1e-4 where its own loop tests 1e-6"
            ),
        ));
    }

    warnings.extend(flash.warnings.iter().cloned());
    Ok(TvFractionFlashResult {
        pressure: pascals(pressure),
        temperature: t,
        beta: flash.beta,
        volume_fraction: found,
        phase: flash.phase,
        x: flash.x,
        y: flash.y,
        k: flash.k,
        z_liquid: flash.z_liquid,
        z_vapour: flash.z_vapour,
        iterations,
        residual: error,
        warnings,
    })
}
