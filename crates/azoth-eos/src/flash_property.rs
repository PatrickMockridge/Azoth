//! Inverting a molar property for temperature, at a fixed pressure.
//!
//! The shared machinery behind [`crate::ph_flash`] and [`crate::ps_flash`], which differ
//! in two things: which property is being inverted, and the variable the iteration runs
//! in. Everything else - how the state at a trial temperature is assembled, the damping,
//! the step clamp, what happens when a trial temperature cannot be evaluated, why the
//! warnings are deduplicated - is the same for both, and subtle enough that writing it
//! twice would invite the two copies to disagree.
//!
//! The iteration is upstream's: `thermodynamicoperations/flashops/PSFlash.java` and
//! `PHflash.java`, NeqSim 3.20.0. Both are quasi-Newton in the temperature - entropy in
//! `T`, enthalpy in `1/T` - damped by a factor that halves whenever the residual grows,
//! and neither is fatal when a trial temperature cannot be evaluated.

use azoth_core::spec::ModelAlgorithm;
use azoth_core::units::{Pressure, ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, Warning};

use crate::mixture::{Mixture, RootSide};
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pt_flash::pt_flash;
use crate::results::{Phase, PtFlashResult};

/// The largest temperature step one iteration may take, in kelvin.
///
/// Upstream clamps to ten in both solvers. It is what keeps a Newton step taken far from
/// the root - where the derivative is a poor local model - from throwing the iterate
/// into a region the cubic cannot describe.
const MAX_STEP: f64 = 10.0;

/// How much of a Newton step the first iteration takes, before any damping.
const INITIAL_FACTOR: f64 = 0.8;

/// Below this the damping stops halving: a factor smaller than this cannot make progress
/// and only drives the iteration into its cap.
const MIN_FACTOR: f64 = 0.1;

/// The relative term in the entropy solver's tolerance.
///
/// Upstream's `RELATIVE_ENTROPY_FLASH_TOLERANCE`. The tolerance is
/// `max(algorithm.tolerance, |target| * this)`, so a large target is not held to an
/// absolute residual that is vanishingly small relative to it.
const RELATIVE_ENTROPY_TOLERANCE: f64 = 1.0e-10;

/// A residual this small that has stopped improving is accepted rather than driven on.
///
/// Upstream's `STAGNANT_ENTROPY_TOLERANCE`, with its `STAGNANT_ITERATION_LIMIT` of five.
/// The rule is what lets a state near the edge of the cubic's validity return an answer
/// instead of exhausting the cap: below this the residual is moving by rounding, and
/// further iterations cannot improve it.
const STAGNANT_RESIDUAL: f64 = 1.0e-4;

/// Consecutive non-improving iterations at a low residual before the answer is accepted.
const STAGNANT_LIMIT: u32 = 5;

/// How many times the enthalpy solver may halve the gap back towards a temperature that
/// worked, before the failure is reported.
const RETRY_LIMIT: u32 = 15;

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

/// One evaluation of the property at a trial temperature.
///
/// Carries the heat capacity alongside the value because the Newton step needs it:
/// `dS/dT = cp/T` and, in reciprocal temperature, `dH/d(1/T) = -T**2 * cp`. Computing it
/// in the same call as the property is what keeps the step's derivative and the value it
/// corrects describing the same state.
struct Evaluation {
    value: f64,
    cp: f64,
    flash: PtFlashResult,
}

/// The molar property of a mixture at a temperature and pressure, and its split.
///
/// The composition the flash settles on is the equilibrium one, so this is the property
/// of the *feed* at that state - which is what makes it comparable with a duty or a
/// change a caller supplied.
///
/// # Errors
/// Propagates whatever [`pt_flash`] and [`molar_enthalpy_entropy`] raise. Neither is
/// swallowed here: the caller decides whether a trial temperature that cannot be
/// evaluated is fatal, and both solvers below decide it is not.
pub fn property_at(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    which: Property,
) -> Result<(f64, PtFlashResult)> {
    let evaluation = evaluate(mixture, ideal_gas, t, p, z, which)?;
    Ok((evaluation.value, evaluation.flash))
}

/// The property and the heat capacity at a state, on whichever branch the flash picks.
fn evaluate(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    which: Property,
) -> Result<Evaluation> {
    let flash = pt_flash(mixture, t, p, z)?;

    let (value, cp) = match flash.phase {
        Phase::TwoPhase => {
            let liquid =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.x, flash.z_liquid)?;
            let vapour =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.y, flash.z_vapour)?;
            let beta = flash.beta.unwrap_or(0.0);
            (
                (1.0 - beta) * which.of(&liquid) + beta * which.of(&vapour),
                (1.0 - beta) * liquid.cp.value + beta * vapour.cp.value,
            )
        }
        // One phase, so the whole feed is in it and *its* root describes it. `beta` is
        // deliberately unused, and so are `z_liquid`/`z_vapour`, which belong to the
        // extrapolated phase compositions: where the flash reports a negative flash its
        // phase composition is not a state, and its cubic root is a different number
        // from the feed's.
        _ => {
            let side = if flash.phase == Phase::AllLiquid {
                RootSide::Liquid
            } else {
                RootSide::Vapour
            };
            let reduced = mixture.reduced_parameters(t, p)?;
            let root = mixture.phase_state(&reduced, z, side)?.z;
            let state = molar_enthalpy_entropy(mixture, ideal_gas, t, p, z, root)?;
            (which.of(&state), state.cp.value)
        }
    };

    Ok(Evaluation { value, cp, flash })
}

/// The derivative of the property with respect to temperature, at constant pressure.
///
/// Taken from the property itself rather than from the heat capacity, because across a
/// phase boundary the two are not the same quantity. `cp` is the phase-fraction-weighted
/// heat capacity of the two phases; the equilibrium `dS/dT` at constant pressure also
/// carries the latent heat of the split changing with temperature, and that term is the
/// **larger** of the two. Measured on methane/n-butane at 5 bar, `dS/dT` reaches
/// 1.3 J/(mol*K**2) where `cp/T` is 0.5. Stepping on the smaller one makes the iteration
/// a fixed point with a gain of about three, which oscillates between two temperatures
/// instead of converging.
///
/// A central difference costs two evaluations and is exact to second order in `delta`,
/// which is chosen relative to the temperature so the step is the same fraction of the
/// state in either implementation. Where either side cannot be evaluated the heat
/// capacity's value is used instead, which is the frozen-composition slope and is close
/// to the right one wherever the state is single phase.
fn slope(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    z: &[f64],
    which: Property,
    temperature: f64,
    fallback: f64,
) -> f64 {
    let delta = (1.0e-4 * temperature).max(1.0e-6);
    let above = evaluate(
        mixture,
        ideal_gas,
        kelvins(temperature + delta),
        p,
        z,
        which,
    );
    let below = evaluate(
        mixture,
        ideal_gas,
        kelvins((temperature - delta).max(1.0)),
        p,
        z,
        which,
    );
    match (above, below) {
        (Ok(hi), Ok(lo)) => (hi.value - lo.value) / (2.0 * delta),
        _ => fallback,
    }
}

/// The outcome of inverting a property for temperature.
pub struct Solution {
    /// The temperature the property was found at.
    pub temperature: f64,
    /// How far the property at that temperature is from the one asked for.
    pub residual: f64,
    /// The flash at the answer.
    pub flash: PtFlashResult,
    /// Iterations taken.
    pub iterations: u32,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

/// Invert `X(T, P) = target` for `T`.
///
/// `Enthalpy` runs upstream's `PHflash.solveQ` and `Entropy` its `PSFlash.solveQ`. They
/// share the damping, the step clamp and the recovery, and differ in the variable the
/// step is taken in and in the residual that is tested.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the iteration reaches its cap without the
///   residual falling below the tolerance, or if the first evaluation at the starting
///   temperature fails outright - there is then no state to iterate from.
/// * [`AzothError::InvalidInput`] from [`property_at`], which does not depend on the
///   temperature and is therefore the caller's error wherever the iteration goes.
pub fn solve_temperature(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    which: Property,
    algorithm: &ModelAlgorithm,
) -> Result<Solution> {
    match which {
        Property::Entropy => solve_entropy(mixture, ideal_gas, p, target, z, algorithm),
        Property::Enthalpy => solve_enthalpy(mixture, ideal_gas, p, target, z, algorithm),
    }
}

/// The starting temperature, from the spec, floored where upstream floors it.
fn start_temperature(algorithm: &ModelAlgorithm) -> f64 {
    algorithm.initial_temperature.unwrap_or(300.0).max(50.0)
}

/// Whether a failed trial is the caller's error rather than a state this path cannot
/// describe.
///
/// An [`AzothError::InvalidInput`] - a composition that is not a composition, a vector of
/// the wrong length - does not depend on the temperature, so it is the caller's error at
/// every point and it must not be absorbed by the recovery. Everything else a trial
/// temperature can raise is a property of *that* temperature, and the solvers recover
/// from it.
fn is_temperature_dependent_failure(error: &AzothError) -> bool {
    !matches!(error, AzothError::InvalidInput { .. })
}

/// Invert the entropy, by upstream's `PSFlash.solveQ`.
fn solve_entropy(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    algorithm: &ModelAlgorithm,
) -> Result<Solution> {
    let tolerance = algorithm
        .tolerance
        .max(target.abs() * RELATIVE_ENTROPY_TOLERANCE);
    let stagnation = STAGNANT_RESIDUAL.min(tolerance);
    let cap = algorithm.max_iterations;

    let mut temperature = start_temperature(algorithm);
    let mut evaluation = evaluate(
        mixture,
        ideal_gas,
        kelvins(temperature),
        p,
        z,
        Property::Entropy,
    )?;
    let mut warnings = evaluation.flash.warnings.clone();

    // Upstream's initial values, kept because the damping rule reads them: a first
    // iteration counts as an improvement on `1.0e10` and so takes a half step, where
    // seeding the errors at infinity would make the relaxation below evaluate to zero
    // and the iteration would never move.
    let mut iterations = 1;
    let mut error = 1.0_f64;
    let mut error_old = 1.0e10_f64;
    let mut factor = INITIAL_FACTOR;
    let mut correct_factor = true;
    let mut stagnant = 0;

    loop {
        if error > error_old && factor > MIN_FACTOR && correct_factor {
            factor *= 0.5;
        } else if error < error_old && correct_factor {
            factor = 1.0;
        }
        iterations += 1;

        // The residual `target - S` falls as `S` rises, so its derivative is the
        // property's slope negated.
        let residual = target - evaluation.value;
        let derivative = -slope(
            mixture,
            ideal_gas,
            p,
            z,
            Property::Entropy,
            temperature,
            evaluation.cp / temperature,
        );
        if derivative.is_finite() && derivative != 0.0 {
            let mut candidate = temperature - factor * residual / derivative;

            if !candidate.is_finite() {
                candidate = temperature + 1.0;
                correct_factor = false;
            } else if candidate < 0.0 {
                candidate = (temperature - MAX_STEP).abs();
                correct_factor = false;
            } else if (temperature - candidate).abs() > MAX_STEP {
                candidate = temperature - (temperature - candidate).signum() * MAX_STEP;
                correct_factor = false;
            } else {
                correct_factor = true;
            }

            match evaluate(
                mixture,
                ideal_gas,
                kelvins(candidate),
                p,
                z,
                Property::Entropy,
            ) {
                Ok(next) => {
                    temperature = candidate;
                    evaluation = next;
                    warnings.extend(evaluation.flash.warnings.iter().cloned());
                }
                Err(ref e) if is_temperature_dependent_failure(e) => {
                    // The step is undone and the damping halved, which is upstream's
                    // response. The state that was good stays good, so nothing is lost
                    // but the progress this step would have made, and the residual is
                    // read from the state that remains.
                    factor *= 0.5;
                    correct_factor = false;
                }
                Err(e) => return Err(e),
            }
        }

        error_old = error;
        error = (target - evaluation.value).abs();
        if iterations > 3 && (error - error_old).abs() <= tolerance && error <= stagnation {
            stagnant += 1;
        } else {
            stagnant = 0;
        }

        if (error + error_old) <= tolerance && iterations >= 3 {
            break;
        }
        if stagnant >= STAGNANT_LIMIT || iterations >= cap {
            break;
        }
    }

    if error > tolerance {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: error,
            tolerance,
        });
    }

    Ok(Solution {
        temperature,
        residual: error,
        flash: evaluation.flash,
        iterations,
        warnings: distinct(&warnings),
    })
}

/// Invert the enthalpy, by upstream's `PHflash.solveQ`, in reciprocal temperature.
///
/// The variable is `1/T` rather than `T` because an enthalpy against temperature is close
/// to linear in the reciprocal, which makes the Newton step a good model over a much
/// wider range - the reason upstream solves it this way and the reason the step clamp is
/// applied to the temperature it converts back to.
fn solve_enthalpy(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    algorithm: &ModelAlgorithm,
) -> Result<Solution> {
    let tolerance = algorithm.tolerance;
    let cap = algorithm.max_iterations;
    // The residual is relative, so a target of zero has no scale to be relative to.
    // Upstream divides by `|Hspec|` unconditionally; refusing is the honest alternative
    // to returning whatever the division produces.
    let scale = target.abs();
    if scale == 0.0 {
        return Err(AzothError::invalid_input(
            "H",
            "an enthalpy of exactly zero has no scale for the solver's relative residual \
             to be measured against, so the inversion would divide by it"
                .to_string(),
        ));
    }

    let mut temperature = start_temperature(algorithm);
    let mut evaluation = evaluate(
        mixture,
        ideal_gas,
        kelvins(temperature),
        p,
        z,
        Property::Enthalpy,
    )?;
    let mut warnings = evaluation.flash.warnings.clone();

    // Upstream's initial values, kept because the damping rule reads them: a first
    // iteration counts as an improvement on `1.0e10` and so takes a half step.
    let mut iterations = 1;
    let mut error = 1.0_f64;
    let mut error_old = 1.0e10_f64;
    let mut factor = INITIAL_FACTOR;
    let mut correct_factor = true;
    let mut retries = 0;
    // The bracket upstream keeps from the sign of the residual. It starts open, because
    // the first iteration has seen one temperature and nothing bounds it from the other
    // side yet.
    let mut min_temperature = 0.0_f64;
    let mut max_temperature = 1.0e10_f64;

    loop {
        if error > error_old && factor > MIN_FACTOR && correct_factor {
            factor *= 0.5;
        } else if error < error_old && correct_factor {
            factor = f64::from(iterations) / (f64::from(iterations) + 1.0);
        }
        iterations += 1;

        // The step is taken in `1/T`, where the residual `(H - target)/scale` has
        // derivative `-T**2 * dH/dT / scale` - so the update is applied to the
        // reciprocal and the clamps below are applied to the temperature it converts
        // back to.
        let residual = (evaluation.value - target) / scale;
        let derivative = -temperature
            * temperature
            * slope(
                mixture,
                ideal_gas,
                p,
                z,
                Property::Enthalpy,
                temperature,
                evaluation.cp,
            )
            / scale;
        if derivative.is_finite() && derivative != 0.0 {
            let reciprocal = 1.0 / temperature - factor * residual / derivative;
            let mut candidate = if reciprocal == 0.0 {
                f64::INFINITY
            } else {
                1.0 / reciprocal
            };

            if !candidate.is_finite() {
                candidate = temperature + 1.0;
                correct_factor = false;
            } else if candidate < 0.0 {
                candidate = (temperature + MAX_STEP).abs();
                correct_factor = false;
            } else if (temperature - candidate).abs() > MAX_STEP {
                candidate = temperature - (temperature - candidate).signum() * MAX_STEP;
                correct_factor = false;
            } else {
                correct_factor = true;
            }
            candidate = candidate.clamp(min_temperature + 0.1, max_temperature - 0.1);

            // A trial temperature the inner flash cannot settle is backed off towards
            // the last one that worked rather than reported. Upstream's comment is the
            // reason: a trial temperature can land where the cubic has no valid root, and
            // aborting the whole inversion there would make the model fail wherever its
            // path crossed such a region.
            let next = loop {
                match evaluate(
                    mixture,
                    ideal_gas,
                    kelvins(candidate),
                    p,
                    z,
                    Property::Enthalpy,
                ) {
                    Ok(next) => break Some(next),
                    Err(ref e) if is_temperature_dependent_failure(e) => {
                        retries += 1;
                        if retries > RETRY_LIMIT {
                            return Err(AzothError::SolverNotConverged {
                                iterations,
                                residual: error,
                                tolerance,
                            });
                        }
                        candidate = 0.5 * (candidate + temperature);
                        if (candidate - temperature).abs() < 1.0e-09 {
                            break None;
                        }
                    }
                    Err(e) => return Err(e),
                }
            };

            if let Some(next) = next {
                temperature = candidate;
                evaluation = next;
                warnings.extend(evaluation.flash.warnings.iter().cloned());

                // The bracket tightens from the sign of the residual at the temperature
                // just evaluated, which is how a step that overshoots is detected.
                if residual > 0.0 && temperature > max_temperature {
                    max_temperature = temperature;
                } else if residual < 0.0 && temperature < min_temperature {
                    min_temperature = temperature;
                }
            }
        }

        error_old = error;
        error = ((evaluation.value - target) / scale).abs();

        if (error + error_old) <= tolerance && iterations >= 3 {
            break;
        }
        if iterations >= cap {
            break;
        }
    }

    if error > tolerance {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: error,
            tolerance,
        });
    }

    Ok(Solution {
        temperature,
        residual: error,
        flash: evaluation.flash,
        iterations,
        warnings: distinct(&warnings),
    })
}

/// One of each distinct warning, in first-seen order.
///
/// The search evaluates the flash many times, so the same caveat arrives many times.
/// Returning them all would make the result depend on how many iterations the search
/// took, which is a property of the algorithm rather than of the state the caller asked
/// about.
pub fn distinct(warnings: &[Warning]) -> Vec<Warning> {
    let mut out: Vec<Warning> = Vec::new();
    for warning in warnings {
        if !out.iter().any(|seen| {
            seen.code == warning.code
                && seen.field == warning.field
                && seen.message == warning.message
        }) {
            out.push(warning.clone());
        }
    }
    out
}
