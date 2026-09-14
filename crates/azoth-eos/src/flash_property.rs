//! Inverting a molar property for temperature, at a fixed pressure.
//!
//! The machinery behind [`crate::ph_flash`] and [`crate::ps_flash`], which differ in
//! exactly one thing: whether the property being inverted is the enthalpy or the
//! entropy. Everything else - how the state at a trial temperature is assembled, why the
//! bracket is a scan, why the answer is narrowed on the temperature rather than on the
//! property, why the warnings are deduplicated - is the same for both, and it is subtle
//! enough that writing it twice would invite the two copies to disagree.
//!
//! The precedent is [`crate::phase_boundary`], shared by `eos.bubble_pressure` and
//! `eos.dew_pressure` because the guard against the trivial solution has to be written
//! once. The same argument applies here, and it is stronger for the single-phase branch
//! below: a feed that is entirely one phase has no vapour fraction, and the flash's
//! value for it is an *extrapolation* rather than a number anybody should use.

use azoth_core::units::{Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, Warning};

use crate::mixture::{Mixture, RootSide};
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pt_flash::pt_flash;
use crate::results::{Phase, PtFlashResult};

/// Which molar property is being inverted.
///
/// An enum rather than a bool because the call sites read as the physics: a reader
/// seeing `Property::Entropy` knows what is being held constant without looking up what
/// `true` meant.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Property {
    /// The molar enthalpy, in J/mol. Inverted by `eos.ph_flash`.
    Enthalpy,
    /// The molar entropy, in J/(mol*K). Inverted by `eos.ps_flash`.
    Entropy,
}

impl Property {
    /// The property of a phase, in its base unit.
    fn of(self, phase: &crate::results::MolarEnthalpyEntropyResult) -> f64 {
        match self {
            Self::Enthalpy => phase.h.value,
            Self::Entropy => phase.s.value,
        }
    }
}

/// The molar property of a mixture at a temperature and pressure, and its split.
///
/// The composition the flash settles on is the equilibrium one, so this is the property
/// of the *feed* at that state - which is what makes it comparable with a duty or a
/// change a caller supplied.
///
/// # Errors
/// Propagates whatever [`pt_flash`] and [`molar_enthalpy_entropy`] raise. Neither is
/// swallowed: a trial state that cannot be evaluated is a bracket that does not contain
/// the answer, and reporting it as a number would hide that.
pub fn property_at(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    which: Property,
) -> Result<(f64, PtFlashResult)> {
    let flash = pt_flash(mixture, t, p, z)?;

    let value = match flash.phase {
        Phase::TwoPhase => {
            let liquid =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.x, flash.z_liquid)?;
            let vapour =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.y, flash.z_vapour)?;
            let beta = flash.beta.unwrap_or(0.0);
            (1.0 - beta) * which.of(&liquid) + beta * which.of(&vapour)
        }
        // One phase, so the whole feed is in it and *its* root describes it. `beta` is
        // deliberately unused - the flash's value there is an extrapolation - and so are
        // `z_liquid`/`z_vapour`, which belong to the extrapolated phase compositions.
        // Reading the root from one of those is what made `H(T)` discontinuous: where the
        // flash reports a negative flash its phase composition is not a state, its cubic
        // root is a different number from the feed's, and the step between the two showed
        // up as a hundred kelvin's worth of enthalpy - so `eos.ph_flash` inverted a
        // function with a jump in it and `process.heater` returned a wrong temperature.
        _ => {
            let side = if flash.phase == Phase::AllLiquid {
                RootSide::Liquid
            } else {
                RootSide::Vapour
            };
            let reduced = mixture.reduced_parameters(t, p)?;
            let root = mixture.phase_state(&reduced, z, side)?.z;
            which.of(&molar_enthalpy_entropy(mixture, ideal_gas, t, p, z, root)?)
        }
    };

    Ok((value, flash))
}

/// The narrowest interval of the spec's scan that contains the requested property.
///
/// Scanned from the bottom up, returning the *first* sign change, so the interval is
/// determined by the bracket alone rather than by where a search happened to start -
/// which is what lets the two implementations agree on the iteration count.
///
/// # A temperature with no state is skipped, not fatal
///
/// Some `(T, P)` pairs on the scan have no admissible liquid root: the cubic's smallest
/// root falls below the mixture's `B`, so `ln(Z - B)` is the logarithm of a negative
/// number and the state does not exist. That is not rare and it is not confined to the
/// ends of the range - on methane/n-butane at 15 bar it happens at 160 K and nowhere
/// else between 100 K and 400 K.
///
/// Such a point has no property, so it cannot bracket anything and cannot be compared
/// against the target. Aborting the whole search on one - which this did until the
/// process layer needed a valve at 15 bar - made `eos.ph_flash` unusable at ordinary
/// states, and the failure looked like a caller's error rather than a gap in the scan.
///
/// **Only an out-of-range state is skipped.** An [`AzothError::InvalidInput`] - a
/// composition that is not a composition, a vector of the wrong length - does not depend
/// on the temperature, so it is the caller's error at every point and it propagates
/// immediately. Skipping those too would turn "your `z` is wrong" into "the solver did
/// not converge", which is a worse answer to a question nobody asked: it reports a
/// failure of the search where the search was never the problem.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if no temperature on the scan produces a
///   property on the other side of the target, which means the requested state is
///   outside the range this model covers. When nothing on the scan was evaluable at all
///   the residual is infinite, which is the honest reading rather than a number.
/// * Anything [`property_at`] raises that is not an out-of-range state.
#[allow(clippy::too_many_arguments)] // The three bracket values are the spec's, passed
// separately so the signature reads like the block it mirrors.
pub fn bracket_by_scan(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    which: Property,
    lower: f64,
    upper: f64,
    steps: u32,
) -> Result<(f64, f64)> {
    let mut previous: Option<(f64, f64)> = None;

    for index in 0..=steps {
        let t = lower + (upper - lower) * f64::from(index) / f64::from(steps);
        let value = match property_at(mixture, ideal_gas, kelvins(t), p, z, which) {
            Ok((value, _)) => value,
            Err(AzothError::OutOfRange { .. }) => continue,
            Err(other) => return Err(other),
        };
        if let Some((previous_t, previous_value)) = previous
            && (value - target) * (previous_value - target) <= 0.0
        {
            return Ok((previous_t, t));
        }
        previous = Some((t, value));
    }

    Err(AzothError::SolverNotConverged {
        iterations: steps,
        residual: previous.map_or(f64::INFINITY, |(_, value)| {
            (value - target).abs() / target.abs().max(1.0)
        }),
        tolerance: 0.0,
    })
}

/// One of each distinct warning, in first-seen order.
///
/// The search evaluates the flash thousands of times, so the same caveat arrives
/// thousands of times. Returning them all would make the result depend on how many
/// iterations the search took, which is a property of the algorithm rather than of the
/// state the caller asked about.
pub fn distinct(warnings: &[Warning]) -> Vec<Warning> {
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

/// What a property inversion found.
pub struct Solved {
    /// The temperature that satisfies the property. This is the answer.
    pub temperature: f64,
    /// `|X(T) - X_target| / max(|X_target|, 1)` at the answer.
    pub residual: f64,
    /// The flash evaluated at the answer.
    pub flash: PtFlashResult,
    /// Bisection steps taken.
    pub iterations: u32,
    /// Caveats from every trial, deduplicated.
    pub warnings: Vec<Warning>,
}

/// Invert `X(T, P) = target` for `T` by bracketing, then bisecting.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the spec declares no bracket, which would be a
///   generator bug rather than a caller's.
/// * [`AzothError::SolverNotConverged`] if no temperature on the bracket covers the
///   requested property, or if the bisection reaches its cap.
pub fn solve_temperature(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    which: Property,
    algorithm: &azoth_core::spec::ModelAlgorithm,
) -> Result<Solved> {
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
        target,
        z,
        which,
        bracket.lower,
        bracket.upper,
        bracket.steps,
    )?;

    let (mut lo, mut hi) = (lower, upper);
    let (mut lo_value, _) = property_at(mixture, ideal_gas, kelvins(lo), p, z, which)?;

    let mut iterations: u32 = 0;
    let mut mid = lo;
    let mut residual = f64::INFINITY;
    let mut warnings: Vec<Warning> = Vec::new();
    let mut state = None;

    while iterations < algorithm.max_iterations {
        iterations += 1;
        mid = 0.5 * (lo + hi);
        let (mid_value, flash) = property_at(mixture, ideal_gas, kelvins(mid), p, z, which)?;
        warnings.extend(flash.warnings.iter().cloned());
        state = Some(flash);

        // Convergence is declared on the *residual*, not on the width of the temperature
        // bracket. The bracket bounds the residual through `|H'| * (hi - lo) / 2` only
        // while `H(T)` is continuous; where it is not, the bisection collapses onto the
        // jump and returns a temperature whose enthalpy is not the one asked for, with no
        // error. It did exactly that for 105 of 111 enthalpy targets before
        // `eos.pt_flash` stopped handing back an extrapolated vapour fraction.
        residual = (mid_value - target).abs() / target.abs().max(1.0);
        if residual <= algorithm.tolerance {
            break;
        }

        if (mid_value - target) * (lo_value - target) <= 0.0 {
            hi = mid;
        } else {
            lo = mid;
            lo_value = mid_value;
        }
    }

    if residual > algorithm.tolerance {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual,
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

    Ok(Solved {
        temperature: mid,
        residual,
        flash,
        iterations,
        warnings: distinct(&warnings),
    })
}
