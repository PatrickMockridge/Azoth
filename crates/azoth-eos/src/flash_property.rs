//! Inverting a molar property for temperature, at a fixed pressure.
//!
//! The shared machinery behind [`crate::ph_flash`] and [`crate::ps_flash`], which differ
//! in two things: which property is being inverted, and the variable the iteration runs
//! in. Everything else - how the state at a trial temperature is assembled, the damping,
//! the step clamp, what happens when a trial temperature cannot be evaluated, why the
//! warnings are deduplicated - is the same for both, and subtle enough that writing it
//! twice would invite the two copies to disagree.
//!
//! The iteration is upstream's `solveQ`: `thermodynamicoperations/flashops/PSFlash.java`
//! and `PHflash.java`, NeqSim 3.20.0, both running it for `type == 0` - which is every
//! default caller. It is a quasi-Newton in the temperature - entropy in `T`, enthalpy in
//! `1/T` - damped by a factor that halves whenever the residual grows, and it is not
//! fatal when a trial temperature cannot be evaluated.
//!
//! Two more of upstream's flash classes are **not ported because nothing reaches
//! them**. `dTPflash` extends `TPflash` and iterates only a named subset of components
//! toward iso-fugacity - a membrane model - and its sole caller in the whole checkout is
//! an example test. `TPgradientFlash` solves the composition and pressure profile of a
//! fluid column under a temperature gradient and gravity, and its only callers are two
//! tests. Neither is dead in the sense `QfuncFlash` is - both have a public
//! `ThermodynamicOperations` method - but neither is reached from production code, and
//! a model with no caller is a promise nothing keeps.
//!
//! `QfuncFlash` is where upstream declares the shape the subclasses share - the
//! `calcdQdP`/`calcdQdT`/`calcdQdPP`/`calcdQdTT` quartet and a `run()` that calls the
//! second-order solver. **Its own `run()` is never invoked**: every subclass
//! (`PSFlash`, `TSFlash`, `TUflash`, `PVflash`, `THflash`, and the reference-EOS
//! variants) overrides it and none calls `super.run()`, and nothing constructs the base
//! class. So the "2x2 Newton on the Q_P/Q_T pair" that class looks like it runs is code
//! no execution reaches, and the scheme that runs is each subclass's own.
//!
//! Upstream carries a second scheme behind `type != 0`, `SysNewtonRhapsonPHflash`, which
//! solves the isofugacity residuals and the energy residual together in `(u, ln T)`. It
//! is reachable - `ThrottlingValve` and `Compressor` pass `type = 1` - and it is **not**
//! ported. Its Jacobian entry for the energy row is `dH/dT = Cp`, and `cp` is the
//! phase-fraction-weighted heat capacity rather than the equilibrium `dH/dT` that
//! [`slope`] below documents the difference between; across a phase boundary the two
//! differ by the latent heat of the split moving with temperature, which is the larger
//! term. Measured on methane/n-butane, it converges in three steps where this one does
//! and diverges from a cold start where this one takes eighteen. Both specs record it.

use azoth_core::spec::ModelAlgorithm;
use azoth_core::units::{Pressure, ThermodynamicTemperature, kelvins, pascals};
use azoth_core::warning::WarningCode;
use azoth_core::{AzothError, Result, Warning};

use crate::mixture::{Mixture, RootSide};
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::pt_flash::pt_flash;
use crate::results::{MolarEnthalpyEntropyResult, Phase, PtFlashResult};

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
    /// The molar enthalpy, in J/mol. Inverted by `eos.ph_flash` and `eos.th_flash`.
    Enthalpy,
    /// The molar entropy, in J/(mol*K). Inverted by `eos.ps_flash` and `eos.ts_flash`.
    Entropy,
    /// The molar volume, in m**3/mol. Inverted by `eos.tv_flash` and `eos.pv_flash`.
    Volume,
    /// The molar internal energy `U = H - P*V`, in J/mol. Inverted by `eos.tu_flash`,
    /// `eos.pu_flash` and `eos.vu_flash`.
    InternalEnergy,
}

/// The molar properties of one phase at a state, the ones an inversion assembles.
///
/// `h`, `s` and `cp` come from [`molar_enthalpy_entropy`]; `v` is `z*R*T/P` for the
/// phase's compressibility factor, the same arithmetic `eos.pr_molar_volume` does.
struct PhaseProperties {
    h: f64,
    s: f64,
    cp: f64,
    v: f64,
}

impl PhaseProperties {
    fn of(state: &MolarEnthalpyEntropyResult, z: f64, t: f64, p: f64) -> Self {
        Self {
            h: state.h.value,
            s: state.s.value,
            cp: state.cp.value,
            v: z * MOLAR_GAS_CONSTANT * t / p,
        }
    }

    fn value(&self, which: Property, pressure: f64) -> f64 {
        match which {
            Property::Enthalpy => self.h,
            Property::Entropy => self.s,
            Property::Volume => self.v,
            Property::InternalEnergy => self.h - pressure * self.v,
        }
    }
}

/// One evaluation of the property at a trial state.
///
/// Carries the heat capacity and the molar volume alongside the value because the Newton
/// step needs a derivative: `cp` for a temperature step (`dS/dT = cp/T`,
/// `dH/d(1/T) = -T**2*cp`) and the volume for a pressure step (`dH/dP ~ V`). Computing
/// them in the same call as the property keeps the step's derivative and the value it
/// corrects describing the same state.
struct Evaluation {
    value: f64,
    cp: f64,
    volume: f64,
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
    let tk = t.value;
    let pk = p.value;

    let (value, cp, volume) = match flash.phase {
        Phase::TwoPhase => {
            let liquid =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.x, flash.z_liquid)?;
            let vapour =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.y, flash.z_vapour)?;
            let liquid_p = PhaseProperties::of(&liquid, flash.z_liquid, tk, pk);
            let vapour_p = PhaseProperties::of(&vapour, flash.z_vapour, tk, pk);
            let beta = flash.beta.unwrap_or(0.0);
            (
                (1.0 - beta) * liquid_p.value(which, pk) + beta * vapour_p.value(which, pk),
                (1.0 - beta) * liquid_p.cp + beta * vapour_p.cp,
                (1.0 - beta) * liquid_p.v + beta * vapour_p.v,
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
            let properties = PhaseProperties::of(&state, root, tk, pk);
            (properties.value(which, pk), properties.cp, properties.v)
        }
    };

    Ok(Evaluation {
        value,
        cp,
        volume,
        flash,
    })
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

/// The derivative of the property with respect to pressure, at constant temperature.
///
/// The pressure-side twin of [`slope`], a central difference over pressure. The fallback
/// is a crude single-phase estimate: `dV/dP ~ -V/P`, `dH/dP ~ V` and, by the Maxwell
/// relation `dS/dP = -dV/dT`, `dS/dP ~ -V/T`.
fn slope_pressure(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    z: &[f64],
    which: Property,
    pressure: f64,
    fallback: f64,
) -> f64 {
    let delta = (1.0e-4 * pressure).max(1.0);
    let above = evaluate(mixture, ideal_gas, t, pascals(pressure + delta), z, which);
    let below = evaluate(
        mixture,
        ideal_gas,
        t,
        pascals((pressure - delta).max(1.0)),
        z,
        which,
    );
    match (above, below) {
        (Ok(hi), Ok(lo)) => (hi.value - lo.value) / (2.0 * delta),
        _ => fallback,
    }
}

/// The outcome of inverting a property for temperature or pressure.
pub struct Solution {
    /// The temperature at the answer - the one solved for, or the fixed one.
    pub temperature: f64,
    /// The pressure at the answer - the one solved for, or the fixed one.
    pub pressure: f64,
    /// How far the property at that state is from the one asked for.
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
        Property::Enthalpy => solve_reciprocal(mixture, ideal_gas, p, target, z, which, algorithm),
        Property::Volume => solve_volume(mixture, ideal_gas, p, target, z, algorithm),
        Property::InternalEnergy => {
            solve_reciprocal(mixture, ideal_gas, p, target, z, which, algorithm)
        }
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
        pressure: p.value,
        residual: error,
        flash: evaluation.flash,
        iterations,
        warnings: distinct(&warnings),
    })
}

/// The relative term in the volume solver's tolerance, upstream's `PVflash` test
/// `|V - Vspec|/Vspec < 1e-9`.
const RELATIVE_VOLUME_TOLERANCE: f64 = 1.0e-9;

/// Invert the volume, by upstream's `PVflash.solveQ`, in temperature.
///
/// The same Newton-over-temperature structure as [`solve_entropy`]: volume rises with
/// temperature at fixed pressure, so the residual `target - V` has derivative `-dV/dT`.
fn solve_volume(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    algorithm: &ModelAlgorithm,
) -> Result<Solution> {
    let tolerance = algorithm
        .tolerance
        .max(target.abs() * RELATIVE_VOLUME_TOLERANCE);
    let stagnation = STAGNANT_RESIDUAL.min(tolerance);
    let cap = algorithm.max_iterations;

    let mut temperature = start_temperature(algorithm);
    let mut evaluation = evaluate(
        mixture,
        ideal_gas,
        kelvins(temperature),
        p,
        z,
        Property::Volume,
    )?;
    let mut warnings = evaluation.flash.warnings.clone();

    // Upstream's initial values, kept because the damping rule reads them.
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

        let residual = target - evaluation.value;
        let derivative = -slope(
            mixture,
            ideal_gas,
            p,
            z,
            Property::Volume,
            temperature,
            evaluation.volume / temperature,
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
                Property::Volume,
            ) {
                Ok(next) => {
                    temperature = candidate;
                    evaluation = next;
                    warnings.extend(evaluation.flash.warnings.iter().cloned());
                }
                Err(ref e) if is_temperature_dependent_failure(e) => {
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
        pressure: p.value,
        residual: error,
        flash: evaluation.flash,
        iterations,
        warnings: distinct(&warnings),
    })
}

/// Invert a property in reciprocal temperature, by upstream's `PHflash.solveQ` /
/// `PUflash.solveQ`.
///
/// The variable is `1/T` rather than `T` because a property against temperature is close
/// to linear in the reciprocal, which makes the Newton step a good model over a much
/// wider range - the reason upstream solves it this way and the reason the step clamp is
/// applied to the temperature it converts back to.
fn solve_reciprocal(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    p: Pressure,
    target: f64,
    z: &[f64],
    which: Property,
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
    let mut evaluation = evaluate(mixture, ideal_gas, kelvins(temperature), p, z, which)?;
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
            * slope(mixture, ideal_gas, p, z, which, temperature, evaluation.cp)
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
                match evaluate(mixture, ideal_gas, kelvins(candidate), p, z, which) {
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
        pressure: p.value,
        residual: error,
        flash: evaluation.flash,
        iterations,
        warnings: distinct(&warnings),
    })
}

/// Invert `X(T, P) = target` for `P`, at a fixed temperature.
///
/// Newton over pressure with a central-difference slope, the same damping, step clamp
/// and recovery as [`solve_temperature`]. Upstream's `TVflash`, `THflash`, `TSflash` and
/// `TUflash` all reduce to this; they differ only in which property is inverted.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the iteration reaches its cap without the
///   residual falling below the tolerance, or if the first evaluation fails outright.
/// * [`AzothError::InvalidInput`] from [`property_at`].
#[allow(clippy::too_many_arguments)] // The signature is the solver's contract.
pub fn solve_pressure(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    target: f64,
    z: &[f64],
    which: Property,
    algorithm: &ModelAlgorithm,
    start_pressure: f64,
) -> Result<Solution> {
    let tolerance = algorithm.tolerance;
    let cap = algorithm.max_iterations;
    // The residual is relative to the target's magnitude, so a value like a molar volume
    // - whose answer differs from the target only in the last bits - is not compared by a
    // subtraction that cancels to machine noise.
    let scale = target.abs().max(1.0);

    let mut pressure = start_pressure.max(1.0);
    let mut evaluation = evaluate(mixture, ideal_gas, t, pascals(pressure), z, which)?;
    let mut warnings = evaluation.flash.warnings.clone();

    // Upstream's initial values, kept because the damping rule reads them: a first
    // iteration counts as an improvement on `1.0e10` and so takes a half step.
    let mut iterations = 1;
    let mut error = 1.0_f64;
    let mut error_old = 1.0e10_f64;
    let mut factor = INITIAL_FACTOR;
    let mut correct_factor = true;
    let mut retries = 0;

    loop {
        if error > error_old && factor > MIN_FACTOR && correct_factor {
            factor *= 0.5;
        } else if error < error_old && correct_factor {
            factor = 1.0;
        }
        iterations += 1;

        let residual = (evaluation.value - target) / scale;
        let fallback = match which {
            Property::Volume => -evaluation.volume / pressure,
            Property::Enthalpy | Property::InternalEnergy => evaluation.volume,
            Property::Entropy => -evaluation.volume / t.value,
        };
        let derivative =
            slope_pressure(mixture, ideal_gas, t, z, which, pressure, fallback) / scale;
        if derivative.is_finite() && derivative != 0.0 {
            let mut candidate = pressure - factor * residual / derivative;

            if !candidate.is_finite() {
                candidate = pressure * 1.1;
                correct_factor = false;
            } else if candidate <= 0.0 {
                candidate = pressure / 2.0;
                correct_factor = false;
            } else if (pressure - candidate).abs() > 0.5 * pressure {
                candidate = pressure - (pressure - candidate).signum() * 0.5 * pressure;
                correct_factor = false;
            } else {
                correct_factor = true;
            }

            // A trial pressure the inner flash cannot settle is backed off towards the
            // last one that worked, as in [`solve_enthalpy`].
            let next = loop {
                match evaluate(mixture, ideal_gas, t, pascals(candidate), z, which) {
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
                        candidate = 0.5 * (candidate + pressure);
                        if (candidate - pressure).abs() < 1.0 {
                            break None;
                        }
                    }
                    Err(e) => return Err(e),
                }
            };

            if let Some(next) = next {
                pressure = candidate;
                evaluation = next;
                warnings.extend(evaluation.flash.warnings.iter().cloned());
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
        temperature: t.value,
        pressure,
        residual: error,
        flash: evaluation.flash,
        iterations,
        warnings: distinct(&warnings),
    })
}

/// What a `(V, X)` flash holds besides the volume.
///
/// `OptimizedVUflash` and `VHflashQfunc` are the same solver with a different energy
/// target: their `calcdQdP`, `calcdQdPP` and `calcdQdTT` are the same lines, and
/// `calcdQdT` reads `Hspec` where the VU form reads `Uspec + P Vspec`. That expression
/// *is* the enthalpy, so the two differ in nothing but whether the caller supplies it
/// or it is assembled from the volume.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EnergyTarget {
    /// An internal energy, at the volume the flash is already holding. The enthalpy
    /// target is `u + P v`, which moves with the pressure the iteration tries.
    InternalEnergy(f64),
    /// An enthalpy, supplied whole.
    Enthalpy(f64),
}

impl EnergyTarget {
    /// The enthalpy the state should have, at a trial pressure.
    fn enthalpy(self, pressure: f64, v_spec: f64) -> f64 {
        match self {
            Self::InternalEnergy(u) => u + pressure * v_spec,
            Self::Enthalpy(h) => h,
        }
    }
}

/// Invert `(V, X) = (v_spec, target)` for `(P, T)` together, by upstream's
/// `OptimizedVUflash`.
///
/// A decoupled 2x2 Newton in the Q-function form `Q_P = P(V - Vspec)/(R T)` and
/// `Q_T = (Uspec + P Vspec - H)/(T R)`, with the diagonal derivatives only. One inner
/// [`pt_flash`] per iteration. This is the one solver that moves both state variables, so
/// it does not share the single-variable `solve_temperature`/`solve_pressure` loop.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the iteration reaches its cap without both the
///   volume and the internal-energy specification being met to the upstream's `1e-3`
///   relative acceptance.
#[allow(clippy::too_many_arguments)] // The signature is the solver's contract.
pub fn solve_pressure_temperature(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    v_spec: f64,
    target: EnergyTarget,
    z: &[f64],
    algorithm: &ModelAlgorithm,
    start_pressure: f64,
    start_temperature: f64,
) -> Result<Solution> {
    let cap = algorithm.max_iterations;
    let mut pressure = start_pressure.max(1.0);
    let mut temperature = start_temperature.max(50.0);
    let mut damping = 0.8;
    let mut stagnation = 0;

    let mut evaluation = evaluate(
        mixture,
        ideal_gas,
        kelvins(temperature),
        pascals(pressure),
        z,
        Property::Enthalpy,
    )?;
    let mut warnings = evaluation.flash.warnings.clone();

    let mut iterations = 0;
    let mut last_error = f64::MAX;
    // Whether the loop left by its own convergence criterion rather than by the cap.
    // The two are different answers and only one of them is converged.
    let mut converged = false;

    loop {
        iterations += 1;
        let h = evaluation.value;
        let v = evaluation.volume;
        let cp = evaluation.cp;
        let tk = temperature;
        let pk = pressure;

        let q_p = pk * (v - v_spec) / (MOLAR_GAS_CONSTANT * tk);
        let q_t = (target.enthalpy(pk, v_spec) - h) / (tk * MOLAR_GAS_CONSTANT);
        let dvdp = slope_pressure(
            mixture,
            ideal_gas,
            kelvins(tk),
            z,
            Property::Volume,
            pk,
            -v / pk,
        );
        let mut dq_pp =
            (v - v_spec) / (MOLAR_GAS_CONSTANT * tk) + pk * dvdp / (MOLAR_GAS_CONSTANT * tk);
        let mut dq_tt = -cp / (tk * MOLAR_GAS_CONSTANT) - q_t / tk;
        if dq_pp.abs() < 1.0e-12 {
            dq_pp = dq_pp.signum() * 1.0e-12;
        }
        if dq_tt.abs() < 1.0e-12 {
            dq_tt = dq_tt.signum() * 1.0e-12;
        }

        let delta_p = (-damping * q_p / dq_pp).clamp(-0.3 * pk, 0.3 * pk);
        let delta_t = (-damping * q_t / dq_tt).clamp(-50.0, 50.0);
        let ny_p = (pk + delta_p).clamp(100.0, 2.0e8);
        let ny_t = (tk + delta_t).clamp(50.0, 5000.0);

        match evaluate(
            mixture,
            ideal_gas,
            kelvins(ny_t),
            pascals(ny_p),
            z,
            Property::Enthalpy,
        ) {
            Ok(next) => {
                pressure = ny_p;
                temperature = ny_t;
                evaluation = next;
                warnings.extend(evaluation.flash.warnings.iter().cloned());
            }
            Err(ref e) if is_temperature_dependent_failure(e) => {
                damping = (damping * 0.7).max(0.05);
            }
            Err(e) => return Err(e),
        }

        let pres_error = ((ny_p - pk) / ny_p.max(0.1)).abs();
        let temp_error = ((ny_t - tk) / ny_t.max(1.0)).abs();
        let total_error = pres_error + temp_error;
        let vol_err = ((evaluation.volume - v_spec) / v_spec).abs();
        let h_target = target.enthalpy(ny_p, v_spec);
        let h_err = ((evaluation.value - h_target) / h_target.abs().max(1.0)).abs();

        if total_error < 1.0e-6 && vol_err < 1.0e-6 && h_err < 1.0e-5 {
            converged = true;
            break;
        }
        if total_error < last_error {
            damping = (damping * 1.1).min(0.8);
            stagnation = 0;
        } else {
            damping = (damping * 0.7).max(0.05);
            stagnation += 1;
        }
        let _ = stagnation;
        last_error = total_error;

        if iterations >= cap {
            break;
        }
    }

    let vol_err = ((evaluation.volume - v_spec) / v_spec).abs();
    let h_target = target.enthalpy(pressure, v_spec);
    let h_err = ((evaluation.value - h_target) / h_target.abs().max(1.0)).abs();
    if vol_err >= 1.0e-3 || h_err >= 1.0e-3 {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: vol_err.max(h_err),
            tolerance: 1.0e-3,
        });
    }
    if !converged {
        // NeqSim's `lastRunConverged = false`. The acceptance below is the *specification*'s,
        // and it is loose enough to pass on an iterate the iteration never settled at:
        // a liquid's volume barely moves with pressure, so at pure propane 250 K the
        // 10 bar state comes back as 15.5 bar - 55% out - with a relative volume error
        // of 5.5e-4, under the 1e-3 the specification is judged by. NeqSim returns that
        // number too and sets a flag; returning it without a word is the one thing
        // neither does.
        warnings.push(Warning::new(
            WarningCode::SolverNotConverged,
            format!(
                "the iteration reached its cap of {cap} without both the volume and the \
                 energy residual falling below its own convergence criterion, so the \
                 pressure and temperature are the last iterate rather than a converged \
                 state. They satisfy the volume and internal energy asked for to \
                 {:.1e} and {:.1e} relative, which is the acceptance NeqSim 3.20.0 \
                 applies and is looser than the iteration's.",
                vol_err, h_err
            ),
        ));
    }

    Ok(Solution {
        temperature,
        pressure,
        residual: vol_err.max(h_err),
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
