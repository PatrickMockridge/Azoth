//! `eos.ph_flash` - the temperature a mixture reaches at a given pressure and enthalpy.
//!
//! Spec: `specs/models/eos/ph_flash.yaml`
//!
//! # The procedure is an outer solve over two things that already exist
//!
//! The enthalpy of a mixture at a pressure is a **strictly increasing** function of
//! temperature, which is what makes a bisection well posed, and it is assembled from
//! parts this crate already has: the phase split at a trial temperature, from
//! [`pt_flash`], and each phase's enthalpy, from [`molar_enthalpy_entropy`].
//!
//! ```text
//! H(T) = (1 - beta) * H_liquid(x, Z_l) + beta * H_vapour(y, Z_v)
//! ```
//!
//! A single-phase state is the same expression with the whole of it on one side, and
//! the root that describes that phase. It is deliberately *not* expressed through
//! `beta` there: the flash reports a `beta` outside `[0, 1]` for a single-phase feed -
//! it is the extrapolated split, not a physical one - and multiplying an enthalpy by
//! 1.888 would be a wrong answer shaped exactly like a right one.
//!
//! # The mirror
//!
//! `python/src/azoth/eos/reference/ph_flash.py` runs the same procedure, and
//! `python/tests/models/test_ph_flash.py` compares them case by case with the
//! **iteration counts required to match**. That is why the scan is a fixed one over
//! the spec's bracket rather than an adaptive expansion: two expansions that stop at
//! different points take different numbers of steps, and this project treats that as a
//! disagreement rather than a detail.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::algorithm_of;
use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pt_flash::pt_flash;
use crate::results::{PhFlashResult, Phase, PtFlashResult};

/// The molar enthalpy of a mixture at a temperature and pressure, and its split.
///
/// The composition the flash settles on is the equilibrium one, so this is the
/// enthalpy of the *feed* at that state - which is what makes it comparable with a
/// duty a caller supplied.
///
/// # Errors
/// Propagates whatever [`pt_flash`] and [`molar_enthalpy_entropy`] raise. Neither is
/// swallowed: a trial state that cannot be evaluated is a bracket that does not
/// contain the answer, and reporting it as a number would hide that.
pub fn enthalpy_at(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<(f64, PtFlashResult)> {
    let flash = pt_flash(mixture, t, p, z)?;

    let h = match flash.phase {
        Phase::TwoPhase => {
            let liquid =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.x, flash.z_liquid)?;
            let vapour =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.y, flash.z_vapour)?;
            let beta = flash.beta.unwrap_or(0.0);
            (1.0 - beta) * liquid.h.value + beta * vapour.h.value
        }
        // One phase, so the whole feed is in it and its own root describes it. `beta`
        // is deliberately unused - see the module docstring.
        _ => {
            let root = if flash.phase == Phase::AllLiquid {
                flash.z_liquid
            } else {
                flash.z_vapour
            };
            molar_enthalpy_entropy(mixture, ideal_gas, t, p, z, root)?
                .h
                .value
        }
    };

    Ok((h, flash))
}

/// The narrowest interval of the spec's scan that contains the requested enthalpy.
///
/// Scanned from the bottom up, returning the *first* sign change, so the interval is
/// determined by the bracket alone rather than by where a search happened to start -
/// which is what lets the two implementations agree on the iteration count.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if no temperature on the scan produces an
///   enthalpy on the other side of the target, which means the requested state is
///   outside the range this model covers.
#[allow(clippy::too_many_arguments)] // The three bracket values are the spec's, passed
// separately so the signature reads like the block it mirrors.
pub fn bracket_by_scan(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    lower: f64,
    upper: f64,
    steps: u32,
) -> Result<(f64, f64)> {
    let mut previous_t = lower;
    let (mut previous_h, _) = enthalpy_at(mixture, ideal_gas, kelvins(previous_t), p, z)?;

    for index in 1..=steps {
        let t = lower + (upper - lower) * f64::from(index) / f64::from(steps);
        let (h, _) = enthalpy_at(mixture, ideal_gas, kelvins(t), p, z)?;
        if (h - target) * (previous_h - target) <= 0.0 {
            return Ok((previous_t, t));
        }
        previous_t = t;
        previous_h = h;
    }

    Err(AzothError::SolverNotConverged {
        iterations: steps,
        residual: (previous_h - target).abs() / target.abs().max(1.0),
        tolerance: 0.0,
    })
}

/// One of each distinct warning, in first-seen order.
///
/// The search evaluates the flash thousands of times, so the same caveat arrives
/// thousands of times. Returning them all would make the result depend on how many
/// iterations the search took, which is a property of the algorithm rather than of the
/// state the caller asked about.
fn distinct(warnings: &[Warning]) -> Vec<Warning> {
    let mut out: Vec<Warning> = Vec::new();
    for warning in warnings {
        let seen = out.iter().any(|kept| {
            kept.code == warning.code
                && kept.field == warning.field
                && kept.message == warning.message
        });
        if !seen {
            out.push(warning.clone());
        }
    }
    out
}

/// The temperature at which a mixture has a given molar enthalpy at a pressure.
///
/// `h` is a *difference* from the datum `ideal_gas` carries, not an absolute quantity:
/// two calls with different reference values are not comparable, and their difference
/// is a plausible number rather than an error.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `p` is not positive, or a range check on the answer
///   fails.
/// * [`AzothError::InvalidInput`] if the spec declares no bracket, which would be a
///   generator bug rather than a caller's.
/// * [`AzothError::SolverNotConverged`] if no temperature on the bracket covers the
///   requested enthalpy, or if the bisection reaches its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{joules_per_mole, pascals, kelvins};
/// use azoth_eos::mixture::{Component, Mixture};
/// use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
/// use azoth_eos::ph_flash;
///
/// let mixture = Mixture::new(
///     vec![
///         Component::new(kelvins(190.56), pascals(4_599_000.0), 0.0115)?,
///         Component::new(kelvins(425.12), pascals(3_796_000.0), 0.2002)?,
///     ],
///     vec![0.0, 0.01289789, 0.01289789, 0.0],
/// )?;
/// let ideal_gas = IdealGasModel {
///     cp_a: vec![3.0, 5.0],
///     cp_b: vec![0.0, 0.0],
///     cp_c: vec![0.0, 0.0],
///     cp_d: vec![0.0, 0.0],
///     h_ref: vec![0.0, 0.0],
///     s_ref: vec![0.0, 0.0],
///     t_ref: kelvins(300.0),
///     p_ref: pascals(100_000.0),
/// };
/// let r = ph_flash::ph_flash(
///     &mixture,
///     &ideal_gas,
///     pascals(2_000_000.0),
///     joules_per_mole(-6723.003844102389),
///     &[0.6, 0.4],
/// )?;
/// assert!((r.temperature.value - 300.0).abs() < 1.0e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn ph_flash(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    h: MolarEnergy,
    z: &[f64],
) -> Result<PhFlashResult> {
    let spec = &model_gen::PH_FLASH_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            "H" => Some(h.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let algorithm = algorithm_of(spec)?;
    let Some(bracket) = algorithm.bracket else {
        return Err(AzothError::InvalidInput {
            field: "algorithm.bracket".to_string(),
            reason: format!(
                "scheme `{}` scans for a bracket but the spec declares none",
                algorithm.scheme
            ),
        });
    };

    let (lower, upper) = bracket_by_scan(
        mixture,
        ideal_gas,
        p,
        h.value,
        z,
        bracket.lower,
        bracket.upper,
        bracket.steps,
    )?;

    let (mut lo, mut hi) = (lower, upper);
    let (mut lo_h, _) = enthalpy_at(mixture, ideal_gas, kelvins(lo), p, z)?;

    let mut iterations: u32 = 0;
    let mut mid = lo;
    let mut mid_h = lo_h;
    let mut state = None;

    while iterations < algorithm.max_iterations {
        iterations += 1;
        mid = 0.5 * (lo + hi);
        let (evaluated_h, flash) = enthalpy_at(mixture, ideal_gas, kelvins(mid), p, z)?;
        mid_h = evaluated_h;
        warnings.extend(flash.warnings.iter().cloned());
        state = Some(flash);

        // The bracket is narrowed on the *temperature*, not on the enthalpy. Relative
        // convergence on an enthalpy that crosses zero is ill-conditioned - it is
        // genuinely zero at some temperature for a datum that puts it there - while
        // the temperature interval is well behaved and is what the answer is.
        if (hi - lo) <= algorithm.tolerance * mid {
            break;
        }

        if (mid_h - h.value) * (lo_h - h.value) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            lo_h = mid_h;
        }
    }

    if iterations >= algorithm.max_iterations && (hi - lo) > algorithm.tolerance * mid {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: (mid_h - h.value).abs() / h.value.abs().max(1.0),
            tolerance: algorithm.tolerance,
        });
    }

    let Some(flash) = state else {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: f64::INFINITY,
            tolerance: algorithm.tolerance,
        });
    };

    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "T").then_some(mid),
        &mut warnings,
    )?;

    Ok(PhFlashResult {
        temperature: kelvins(mid),
        beta: flash.beta,
        x: flash.x,
        y: flash.y,
        k: flash.k,
        phase: flash.phase,
        z_liquid: flash.z_liquid,
        z_vapour: flash.z_vapour,
        iterations,
        residual: (mid_h - h.value).abs() / h.value.abs().max(1.0),
        warnings: distinct(&warnings),
    })
}
