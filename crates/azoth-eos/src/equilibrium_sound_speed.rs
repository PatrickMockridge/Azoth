//! `eos.equilibrium_sound_speed` - the homogeneous-equilibrium acoustic derivative.
//!
//! Spec: `specs/models/eos/equilibrium_sound_speed.toml`. The port source is NeqSim's
//! `thermo/util/EquilibriumSoundSpeed.java`.
//!
//! `c^2 = (dp/drho)` at constant specific entropy and total composition: each pressure
//! sample is re-equilibrated at fixed `s` (mixtures by a bracketed TP root, pure fluids by
//! the saturation-quality solve), accepted only when the entropy, inventory and fugacity
//! residuals close, and differentiated on a central or one-sided stencil that stays on the
//! centre state's phase branch.
//!
//! **This is not [`crate::sound_speed`].** That model is the frozen phase speed the class
//! reads off one phase; this one lets the phase split relax with pressure and is the value
//! the HEM release model uses.

use azoth_core::units::{
    Pressure, SpecificHeatCapacity, ThermodynamicTemperature, kelvins, pascals,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::cubic::Cubic;
use crate::databank;
use crate::entropy_state::{self, RootPolicy};
use crate::mixture::{Mixture, RootSide};
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pt_flash::pt_flash;
use crate::results::{EquilibriumSoundSpeedResult, EquilibriumSoundSpeedSample, Phase};

use crate::model_gen;

/// `EquilibriumSoundSpeed.ENTROPY_TOLERANCE`, J/(kg K).
const ENTROPY_TOLERANCE: f64 = 1.0e-7;
/// `BALANCE_TOLERANCE`.
const BALANCE_TOLERANCE: f64 = 1.0e-8;
/// `FUGACITY_TOLERANCE`.
const FUGACITY_TOLERANCE: f64 = 1.0e-6;
/// `DERIVATIVE_TOLERANCE`: the successive-derivative band that accepts a step.
const DERIVATIVE_TOLERANCE: f64 = 2.0e-3;
/// `MAX_REFINEMENTS`: the step is halved this many times after the first estimate.
const MAX_REFINEMENTS: usize = 10;
/// The class refuses a relative pressure step outside this band.
const MIN_RELATIVE_STEP: f64 = 1.0e-6;
/// See [`MIN_RELATIVE_STEP`].
const MAX_RELATIVE_STEP: f64 = 0.05;

/// `EquilibriumSoundSpeed.Status`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquilibriumStatus {
    /// Every closure and step-convergence check passed.
    Converged,
    /// A checked equilibrium state could not be obtained.
    FlashFailed,
    /// The requested physics or phase assemblage is unsupported.
    UnsupportedState,
    /// No positive, step-converged density derivative was found.
    DerivativeNotConverged,
}

impl EquilibriumStatus {
    /// The spec's value for this status.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Converged => "converged",
            Self::FlashFailed => "flash_failed",
            Self::UnsupportedState => "unsupported_state",
            Self::DerivativeNotConverged => "derivative_not_converged",
        }
    }
}

/// `EquilibriumSoundSpeed.Stencil`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Stencil {
    /// Symmetric three-point difference.
    Central,
    /// Second-order difference towards increasing pressure.
    Forward,
    /// Second-order difference towards decreasing pressure.
    Backward,
    /// No accepted stencil.
    Unavailable,
}

impl Stencil {
    /// The spec's value for this stencil.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Central => "central",
            Self::Forward => "forward",
            Self::Backward => "backward",
            Self::Unavailable => "unavailable",
        }
    }
}

/// One accepted pressure sample.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumSample {
    /// Absolute pressure, Pa.
    pub pressure: f64,
    /// Temperature, K.
    pub temperature: f64,
    /// Mass over EOS total volume, kg/m3.
    pub density: f64,
    /// Signed specific-entropy error, J/(kg K).
    pub entropy_residual: f64,
    /// The largest component-inventory error over the total moles.
    pub component_residual: f64,
    /// The largest absolute interphase log-fugacity ratio; zero for one phase.
    pub fugacity_residual: f64,
    /// What the flash found.
    pub phase: Phase,
    /// The vapour fraction where the state is two-phase.
    pub vapour_fraction: Option<f64>,
}

/// The class's `Result`.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumSoundSpeedOutcome {
    /// The convergence or failure category.
    pub status: EquilibriumStatus,
    /// The class's own message for that category.
    pub message: String,
    /// The acoustic speed, m/s; absent when no step converged.
    pub sound_speed: Option<f64>,
    /// The final pressure step, Pa.
    pub pressure_step: f64,
    /// The successive-derivative difference at that step, relative.
    pub relative_step_error: Option<f64>,
    /// The target specific entropy, J/(kg K).
    pub specific_entropy: f64,
    /// The accepted stencil.
    pub stencil: Stencil,
    /// Whether any sampled stencil changed the phase assemblage.
    pub phase_boundary: bool,
    /// `tp`/quality evaluations, including rejected trials.
    pub flash_evaluations: u32,
    /// The final stencil's states: centre first.
    pub samples: Vec<EquilibriumSample>,
}

/// One intermediate failure category: every refusal here is the class's `FLASH_FAILED`.
enum CalcError {
    Flash(String),
}

type Calc<T> = std::result::Result<T, CalcError>;

/// The equilibrium acoustic speed of a mixture at a state.
///
/// See [`equilibrium_speed_of`] for the arithmetic. `specific_entropy` is the fixed state
/// variable; when absent it is the equilibrium flash's own value at `(T, P, z)`. A pure fluid
/// in its two-phase region has no unique temperature at a pressure, so its entropy has to be
/// stated.
///
/// # Errors
/// `InvalidInput` for a step outside `[1e-6, 0.05]`, a non-positive state, a composition of
/// the wrong length, or a missing molar mass; and whatever the databank refuses.
pub fn equilibrium_sound_speed(
    components: &[String],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    relative_pressure_step: Option<f64>,
    specific_entropy: Option<SpecificHeatCapacity>,
) -> Result<EquilibriumSoundSpeedResult> {
    let spec = &model_gen::EQUILIBRIUM_SOUND_SPEED_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            "relative_pressure_step" => relative_pressure_step,
            _ => None,
        },
        &mut warnings,
    )?;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
    if z.len() != mixture.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {} components needs {} mole fractions, and {z:?} has {}",
                mixture.len(),
                mixture.len(),
                z.len()
            ),
        ));
    }

    let outcome = equilibrium_speed_of(
        &mixture,
        &ideal_gas,
        t,
        p,
        z,
        relative_pressure_step.unwrap_or(1.0e-3),
        specific_entropy,
    )?;

    Ok(EquilibriumSoundSpeedResult {
        sound_speed: outcome.sound_speed,
        status: outcome.status.as_str().to_string(),
        message: outcome.message,
        pressure_step: outcome.pressure_step,
        relative_step_error: outcome.relative_step_error,
        specific_entropy: outcome.specific_entropy,
        stencil: outcome.stencil.as_str().to_string(),
        phase_boundary_encountered: outcome.phase_boundary,
        flash_evaluations: outcome.flash_evaluations,
        samples: outcome
            .samples
            .into_iter()
            .map(|sample| EquilibriumSoundSpeedSample {
                pressure: pascals(sample.pressure),
                temperature: kelvins(sample.temperature),
                density: sample.density,
                entropy_residual: sample.entropy_residual,
                phase: sample.phase,
                vapour_fraction: sample.vapour_fraction,
            })
            .collect(),
        warnings,
    })
}

/// The equilibrium acoustic speed of a mixture at a state.
///
/// `specific_entropy` is the fixed state variable; when absent it is taken from the
/// equilibrium flash at `(t, p, z)`. A pure fluid in its two-phase region has no unique
/// temperature at a pressure, so its entropy has to be stated.
///
/// # Errors
/// `InvalidInput` for a step outside `[1e-6, 0.05]`, a non-positive state, a composition
/// of the wrong length, or a missing molar mass; and whatever the flash or the databank
/// refuses while an evaluation is being attempted - those becoming the outcome's
/// `flash_failed`, not an error, once inside the calculation.
pub fn equilibrium_speed_of(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    relative_step: f64,
    specific_entropy: Option<SpecificHeatCapacity>,
) -> Result<EquilibriumSoundSpeedOutcome> {
    if !relative_step.is_finite()
        || !(MIN_RELATIVE_STEP..=MAX_RELATIVE_STEP).contains(&relative_step)
    {
        return Err(AzothError::invalid_input(
            "relative_pressure_step",
            format!("the step must be in [1e-6, 0.05], and {relative_step} is not"),
        ));
    }
    if p.value <= 0.0 || !p.value.is_finite() || t.value <= 0.0 || !t.value.is_finite() {
        return Err(AzothError::invalid_input(
            "state",
            format!(
                "a finite positive pressure and temperature are required, and the state is \
                 {:.6} K and {:.6} Pa",
                t.value, p.value
            ),
        ));
    }
    if z.len() != mixture.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {} components needs that many mole fractions, and {z:?} has {}",
                mixture.len(),
                z.len()
            ),
        ));
    }
    crate::phase_boundary::check_composition(z, mixture.len(), "z")?;
    let molar_mass = molar_mass_of(mixture, z)?;

    // The entropy the class reads from the state it was handed: ours is the flash at the
    // stated `(T, P, z)`, which is the state a caller's fluid is in after its own flash.
    let target = match specific_entropy {
        Some(value) => value.value,
        None => {
            let flash = pt_flash(mixture, t, p, z)?;
            specific_entropy_at(mixture, ideal_gas, t, p, z, &flash, molar_mass)?
        }
    };

    let mut calculation = Calculation {
        mixture,
        ideal_gas,
        z,
        molar_mass,
        target,
        evaluations: 0,
        boundary: false,
        start_temperature: t.value,
    };
    Ok(calculation.run(p.value, relative_step))
}

impl Calculation<'_> {
    fn run(&mut self, pressure: f64, relative_step: f64) -> EquilibriumSoundSpeedOutcome {
        let mut samples: Vec<EquilibriumSample> = Vec::new();
        let mut stencil = Stencil::Unavailable;
        let mut step = f64::NAN;
        let mut error = None;

        let centre = match self.solve(pressure) {
            Ok(sample) => sample,
            Err(CalcError::Flash(message)) => {
                return self.outcome(
                    EquilibriumStatus::FlashFailed,
                    message,
                    None,
                    step,
                    error,
                    stencil,
                    samples,
                );
            }
        };
        samples = vec![centre.clone()];
        let mut previous: Option<(Stencil, f64)> = None;
        step = pressure * relative_step;

        for refinement in 0..=MAX_REFINEMENTS {
            if refinement > 0 {
                step *= 0.5;
            }
            let current = match self.estimate(&centre, step) {
                Ok(Some(current)) => current,
                Ok(None) => {
                    previous = None;
                    continue;
                }
                Err(CalcError::Flash(message)) => {
                    return self.outcome(
                        EquilibriumStatus::FlashFailed,
                        message,
                        None,
                        step,
                        error,
                        stencil,
                        samples,
                    );
                }
            };
            samples = current.states.clone();
            stencil = current.stencil;
            if let Some((last_stencil, last_derivative)) = previous {
                if last_stencil == current.stencil
                    && last_derivative > 0.0
                    && current.derivative > 0.0
                    && last_derivative.is_finite()
                    && current.derivative.is_finite()
                {
                    error = Some(
                        (current.derivative - last_derivative).abs() / current.derivative.abs(),
                    );
                    let speed = (1.0 / current.derivative).sqrt();
                    if error.is_some_and(|e| e <= DERIVATIVE_TOLERANCE)
                        && speed.is_finite()
                        && speed > 0.0
                    {
                        return self.outcome(
                            EquilibriumStatus::Converged,
                            "entropy, inventory, fugacity and step convergence checks passed"
                                .to_string(),
                            Some(speed),
                            step,
                            error,
                            stencil,
                            samples,
                        );
                    }
                }
            }
            previous = Some((current.stencil, current.derivative));
        }

        self.outcome(
            EquilibriumStatus::DerivativeNotConverged,
            "no positive step-converged equilibrium density derivative".to_string(),
            None,
            step * 2.0,
            error,
            stencil,
            samples,
        )
    }

    #[allow(clippy::too_many_arguments)] // one argument per field of the Result
    fn outcome(
        &self,
        status: EquilibriumStatus,
        message: String,
        sound_speed: Option<f64>,
        pressure_step: f64,
        relative_step_error: Option<f64>,
        stencil: Stencil,
        samples: Vec<EquilibriumSample>,
    ) -> EquilibriumSoundSpeedOutcome {
        EquilibriumSoundSpeedOutcome {
            status,
            message,
            sound_speed,
            pressure_step,
            relative_step_error,
            specific_entropy: self.target,
            stencil,
            phase_boundary: self.boundary,
            flash_evaluations: self.evaluations,
            samples,
        }
    }

    /// One stencil around `centre` at a pressure step, Pa.
    fn estimate(&mut self, centre: &EquilibriumSample, step: f64) -> Calc<Option<Estimate>> {
        let pressure = centre.pressure;
        let low = self.solve(pressure - step)?;
        let high = self.solve(pressure + step)?;
        let same_low = same_assemblage(centre.phase, low.phase);
        let same_high = same_assemblage(centre.phase, high.phase);
        self.boundary |= !same_low || !same_high;
        if same_low && same_high {
            return Ok(Some(Estimate {
                stencil: Stencil::Central,
                derivative: (high.density - low.density) / (2.0 * step),
                states: vec![centre.clone(), low, high],
            }));
        }
        if same_high {
            let far = self.solve(pressure + 2.0 * step)?;
            self.boundary |= !same_assemblage(centre.phase, far.phase);
            if !same_assemblage(centre.phase, far.phase) {
                return Ok(None);
            }
            return Ok(Some(Estimate {
                stencil: Stencil::Forward,
                derivative: (-3.0 * centre.density + 4.0 * high.density - far.density)
                    / (2.0 * step),
                states: vec![centre.clone(), high, far],
            }));
        }
        if same_low {
            let far = self.solve(pressure - 2.0 * step)?;
            self.boundary |= !same_assemblage(centre.phase, far.phase);
            if !same_assemblage(centre.phase, far.phase) {
                return Ok(None);
            }
            return Ok(Some(Estimate {
                stencil: Stencil::Backward,
                derivative: (3.0 * centre.density - 4.0 * low.density + far.density) / (2.0 * step),
                states: vec![centre.clone(), low, far],
            }));
        }
        Ok(None)
    }

    /// The checked equilibrium state at a pressure and the target entropy.
    fn solve(&mut self, pressure: f64) -> Calc<EquilibriumSample> {
        let state = entropy_state::isentropic_state(
            self.mixture,
            self.ideal_gas,
            self.z,
            self.molar_mass,
            self.target,
            pressure,
            self.start_temperature,
            &mut self.evaluations,
            RootPolicy::SOUND_SPEED,
        )
        .map_err(|error| CalcError::Flash(error.0))?;
        self.checked(EquilibriumSample {
            pressure: state.pressure,
            temperature: state.temperature,
            density: state.density,
            entropy_residual: state.entropy_residual,
            component_residual: state.component_residual,
            fugacity_residual: state.fugacity_residual,
            phase: state.phase,
            vapour_fraction: state.vapour_fraction,
        })
    }

    /// The class's closure check: a sample the result may carry.
    fn checked(&self, sample: EquilibriumSample) -> Calc<EquilibriumSample> {
        let closure = sample.entropy_residual.abs() <= ENTROPY_TOLERANCE
            && sample.component_residual <= BALANCE_TOLERANCE
            && sample.component_residual.is_finite()
            && sample.fugacity_residual <= FUGACITY_TOLERANCE
            && sample.fugacity_residual.is_finite()
            && sample.density > 0.0
            && sample.density.is_finite()
            && sample.temperature > 0.0
            && sample.temperature.is_finite()
            && sample.pressure > 0.0
            && sample.pressure.is_finite();
        if closure {
            Ok(sample)
        } else {
            Err(CalcError::Flash(format!(
                "the equilibrium closure failed at {:.6} Pa: entropy residual = {:.3e} \
                 J/(kg K), component residual = {:.3e}, log fugacity residual = {:.3e}",
                sample.pressure,
                sample.entropy_residual,
                sample.component_residual,
                sample.fugacity_residual
            )))
        }
    }
}

struct Calculation<'a> {
    mixture: &'a Mixture,
    ideal_gas: &'a IdealGasModel,
    z: &'a [f64],
    molar_mass: f64,
    target: f64,
    evaluations: u32,
    boundary: bool,
    start_temperature: f64,
}

struct Estimate {
    stencil: Stencil,
    states: Vec<EquilibriumSample>,
    derivative: f64,
}

/// The class compares the *sets* of phase type names; ours has one two-phase and three
/// single-phase variants.
fn same_assemblage(first: Phase, second: Phase) -> bool {
    matches!(
        (first, second),
        (Phase::TwoPhase, Phase::TwoPhase)
            | (Phase::AllVapour, Phase::AllVapour)
            | (Phase::AllLiquid, Phase::AllLiquid)
            | (Phase::Trivial, Phase::Trivial)
    )
}

/// The mixture's molar mass, kg/mol.
fn molar_mass_of(mixture: &Mixture, composition: &[f64]) -> Result<f64> {
    let mut total = 0.0;
    for (index, component) in mixture.components().iter().enumerate() {
        let mass = component.molar_mass.ok_or_else(|| {
            AzothError::invalid_input(
                "molar_mass",
                format!(
                    "component {index} carries no molar mass, and a specific entropy divides \
                     by it"
                ),
            )
        })?;
        total += composition[index] * mass;
    }
    if total <= 0.0 || !total.is_finite() {
        return Err(AzothError::invalid_input(
            "molar_mass",
            format!("the mixture's molar mass is {total} kg/mol"),
        ));
    }
    Ok(total)
}

/// The flashed state's specific entropy, J/(kg K).
fn specific_entropy_at(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
    flash: &crate::results::PtFlashResult,
    molar_mass: f64,
) -> Result<f64> {
    let molar = match flash.phase {
        Phase::TwoPhase => {
            let beta = flash.vapour_fraction.ok_or_else(|| {
                AzothError::invalid_input(
                    "state",
                    "the flash split the feed and reported no vapour fraction".to_string(),
                )
            })?;
            let liquid =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.x, flash.z_liquid)?;
            let vapour =
                molar_enthalpy_entropy(mixture, ideal_gas, t, p, &flash.y, flash.z_vapour)?;
            (1.0 - beta) * liquid.s.value + beta * vapour.s.value
        }
        Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => {
            let side = if matches!(flash.phase, Phase::AllVapour) {
                RootSide::Vapour
            } else {
                RootSide::Liquid
            };
            let reduced = mixture.reduced_parameters(t, p)?;
            let root = mixture.phase_state(&reduced, z, side)?.z;
            molar_enthalpy_entropy(mixture, ideal_gas, t, p, z, root)?
                .s
                .value
        }
    };
    Ok(molar / molar_mass)
}
