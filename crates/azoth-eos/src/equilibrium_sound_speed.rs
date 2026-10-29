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
use crate::mixture::{Mixture, RootSide};
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::pt_flash::pt_flash;
use crate::pure_saturation::pure_saturation;
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
/// The bracket's temperature tolerance, relative.
const ROOT_TOLERANCE: f64 = 1.0e-12;
/// The entropy bracket may widen this many times before it is refused.
const BRACKET_STEPS: usize = 30;
/// The bisection cap.
const BISECTION_STEPS: usize = 100;
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
        if pressure <= 0.0 || !pressure.is_finite() {
            return Err(CalcError::Flash(format!(
                "the pressure sample {pressure} Pa is not positive"
            )));
        }
        if self.mixture.len() == 1 {
            if let Some(state) = self.saturation_state(pressure)? {
                return Ok(state);
            }
        }
        self.tp(pressure)
    }

    /// The pure fluid's saturation-quality state, where the target entropy lies between the
    /// saturated phases. `None` sends the caller to the TP root.
    fn saturation_state(&mut self, pressure: f64) -> Calc<Option<EquilibriumSample>> {
        let component = &self.mixture.components()[0];
        let tc = component.tc.value;
        let omega = component.omega;
        // **A tolerant bracket.** The cold end is where the cubic's liquid branch is hardest
        // to find, so the first fraction that answers is the bound rather than 0.3 Tc by fiat.
        let mut low_t = 0.0;
        let mut p_low = f64::NAN;
        for fraction in [0.3, 0.4, 0.5, 0.6] {
            let candidate = (1.0_f64).max(fraction * tc);
            if let Ok(result) =
                pure_saturation(component.tc, component.pc, omega, kelvins(candidate))
            {
                low_t = candidate;
                p_low = result.p_sat.value;
                break;
            }
        }
        // **The warm end is bounded by `eos.pure_saturation` itself**, which declines above
        // about 0.94 Tc; the scan takes the first temperature whose saturation pressure is
        // above the sample, and a state above all of them falls back to the TP root.
        let mut high_t = 0.0;
        let mut p_high = f64::NAN;
        for fraction in [0.94, 0.92, 0.90, 0.85, 0.80] {
            if let Ok(result) =
                pure_saturation(component.tc, component.pc, omega, kelvins(fraction * tc))
            {
                if result.p_sat.value > pressure {
                    high_t = fraction * tc;
                    p_high = result.p_sat.value;
                    break;
                }
            }
        }
        if high_t <= 0.0 || low_t <= 0.0 || pressure <= p_low || pressure >= p_high {
            return Ok(None);
        }
        let mut lower = low_t;
        let mut upper = high_t;
        for _ in 0..BISECTION_STEPS {
            let middle = 0.5 * (lower + upper);
            let trial = match pure_saturation(component.tc, component.pc, omega, kelvins(middle)) {
                Ok(result) => result.p_sat.value,
                Err(_) => return Ok(None),
            };
            if trial < pressure {
                lower = middle;
            } else {
                upper = middle;
            }
            if upper - lower <= ROOT_TOLERANCE * 1.0_f64.max(middle) {
                break;
            }
        }
        let temperature = 0.5 * (lower + upper);
        let reduced = match self
            .mixture
            .reduced_parameters(kelvins(temperature), pascals(pressure))
        {
            Ok(reduced) => reduced,
            Err(_) => return Ok(None),
        };
        let liquid_root = match self.mixture.phase_state(&reduced, self.z, RootSide::Liquid) {
            Ok(state) => state.z,
            Err(_) => return Ok(None),
        };
        let vapour_root = match self.mixture.phase_state(&reduced, self.z, RootSide::Vapour) {
            Ok(state) => state.z,
            Err(_) => return Ok(None),
        };
        if (vapour_root - liquid_root).abs() < 1.0e-12 {
            return Ok(None);
        }
        let liquid = molar_enthalpy_entropy(
            self.mixture,
            self.ideal_gas,
            kelvins(temperature),
            pascals(pressure),
            self.z,
            liquid_root,
        );
        let vapour = molar_enthalpy_entropy(
            self.mixture,
            self.ideal_gas,
            kelvins(temperature),
            pascals(pressure),
            self.z,
            vapour_root,
        );
        let (Ok(liquid), Ok(vapour)) = (liquid, vapour) else {
            return Ok(None);
        };
        let s_liquid = liquid.s.value / self.molar_mass;
        let s_vapour = vapour.s.value / self.molar_mass;
        if s_vapour <= s_liquid || !(s_liquid..=s_vapour).contains(&self.target) {
            return Ok(None);
        }
        self.evaluations += 1;
        let quality = (self.target - s_liquid) / (s_vapour - s_liquid);
        // **The EOS volume, unshifted.** The class uses `getVolume` with volume correction
        // off in its own tests - a density correction is not part of the acoustic derivative.
        let v_liquid = liquid_root * MOLAR_GAS_CONSTANT * temperature / pressure;
        let v_vapour = vapour_root * MOLAR_GAS_CONSTANT * temperature / pressure;
        let volume = (1.0 - quality) * v_liquid + quality * v_vapour;
        if volume <= 0.0 || !volume.is_finite() {
            return Err(CalcError::Flash(format!(
                "the saturation-quality volume at {pressure} Pa is {volume} m3/mol"
            )));
        }
        self.checked(EquilibriumSample {
            pressure,
            temperature,
            density: self.molar_mass / volume,
            entropy_residual: (1.0 - quality) * s_liquid + quality * s_vapour - self.target,
            component_residual: 0.0,
            fugacity_residual: 0.0,
            phase: Phase::TwoPhase,
            vapour_fraction: Some(quality),
        })
        .map(Some)
    }

    /// The equilibrium state at a pressure by the bracketed TP entropy root.
    fn tp(&mut self, pressure: f64) -> Calc<EquilibriumSample> {
        let guess = self.start_temperature;
        let initial = self.tp_at(pressure, guess)?;
        if initial.entropy_residual.abs() <= ENTROPY_TOLERANCE {
            return self.checked(initial);
        }
        let mut lower = guess;
        let mut upper = guess;
        let mut low_error = initial.entropy_residual;
        let mut high_error = initial.entropy_residual;
        let mut span = 1.0_f64.max(guess * 0.005);
        for _ in 0..BRACKET_STEPS {
            if !(low_error > 0.0 || high_error < 0.0) {
                break;
            }
            if low_error > 0.0 {
                lower = (guess - span).max(1.0);
                low_error = self.tp_at(pressure, lower)?.entropy_residual;
            }
            if high_error < 0.0 {
                upper = guess + span;
                high_error = self.tp_at(pressure, upper)?.entropy_residual;
            }
            span *= 1.6;
        }
        if !(low_error <= 0.0 && high_error >= 0.0) {
            return Err(CalcError::Flash(format!(
                "could not bracket the specific entropy at {:.6} Pa",
                pressure
            )));
        }
        let mut last_error = f64::NAN;
        for _ in 0..BISECTION_STEPS {
            let temperature = 0.5 * (lower + upper);
            let trial = self.tp_at(pressure, temperature)?;
            last_error = trial.entropy_residual;
            if last_error.abs() <= ENTROPY_TOLERANCE {
                return self.checked(trial);
            }
            if last_error > 0.0 {
                upper = temperature;
            } else {
                lower = temperature;
            }
            if upper - lower <= ROOT_TOLERANCE * 1.0_f64.max(temperature) {
                break;
            }
        }
        Err(CalcError::Flash(format!(
            "the entropy root did not close at {:.6} Pa; residual = {last_error} J/(kg K)",
            pressure
        )))
    }

    /// One TP flash and its checked sample.
    fn tp_at(&mut self, pressure: f64, temperature: f64) -> Calc<EquilibriumSample> {
        self.evaluations += 1;
        let flash = pt_flash(
            self.mixture,
            kelvins(temperature),
            pascals(pressure),
            self.z,
        )
        .map_err(|error| CalcError::Flash(format!("{error}")))?;
        self.assemble(pressure, temperature, &flash)
    }

    /// Turn a flash into the checked sample, or refuse it.
    fn assemble(
        &self,
        pressure: f64,
        temperature: f64,
        flash: &crate::results::PtFlashResult,
    ) -> Calc<EquilibriumSample> {
        let (entropy, density, balance, fugacity) = match flash.phase {
            Phase::TwoPhase => {
                let beta = flash.vapour_fraction.ok_or_else(|| {
                    CalcError::Flash(
                        "the flash split the feed and reported no vapour fraction".to_string(),
                    )
                })?;
                if !(0.0..=1.0).contains(&beta) {
                    return Err(CalcError::Flash(format!(
                        "the vapour fraction {beta} is outside [0, 1]"
                    )));
                }
                let liquid = molar_enthalpy_entropy(
                    self.mixture,
                    self.ideal_gas,
                    kelvins(temperature),
                    pascals(pressure),
                    &flash.x,
                    flash.z_liquid,
                );
                let vapour = molar_enthalpy_entropy(
                    self.mixture,
                    self.ideal_gas,
                    kelvins(temperature),
                    pascals(pressure),
                    &flash.y,
                    flash.z_vapour,
                );
                let (Ok(liquid), Ok(vapour)) = (liquid, vapour) else {
                    return Err(CalcError::Flash(
                        "the phase entropy is not finite".to_string(),
                    ));
                };
                let entropy =
                    ((1.0 - beta) * liquid.s.value + beta * vapour.s.value) / self.molar_mass;
                let v_liquid = flash.z_liquid * MOLAR_GAS_CONSTANT * temperature / pressure;
                let v_vapour = flash.z_vapour * MOLAR_GAS_CONSTANT * temperature / pressure;
                let volume = (1.0 - beta) * v_liquid + beta * v_vapour;
                if volume <= 0.0 || !volume.is_finite() {
                    return Err(CalcError::Flash(format!(
                        "the two-phase volume at {pressure} Pa is {volume} m3/mol"
                    )));
                }
                let mut balance = 0.0_f64;
                let mut fugacity = 0.0_f64;
                for index in 0..self.z.len() {
                    let recovered = (1.0 - beta) * flash.x[index] + beta * flash.y[index];
                    balance = balance.max((recovered - self.z[index]).abs());
                    if self.z[index] > 1.0e-12 {
                        if flash.x[index] <= 0.0 || flash.y[index] <= 0.0 {
                            return Err(CalcError::Flash(
                                "a phase composition is not positive".to_string(),
                            ));
                        }
                        let ln_ratio = flash.x[index].ln() + flash.ln_phi_liquid[index]
                            - flash.y[index].ln()
                            - flash.ln_phi_vapour[index];
                        fugacity = fugacity.max(ln_ratio.abs());
                    }
                }
                (entropy, self.molar_mass / volume, balance, fugacity)
            }
            Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => {
                let side = if matches!(flash.phase, Phase::AllVapour) {
                    RootSide::Vapour
                } else {
                    RootSide::Liquid
                };
                let reduced = self
                    .mixture
                    .reduced_parameters(kelvins(temperature), pascals(pressure))
                    .map_err(|error| CalcError::Flash(format!("{error}")))?;
                let root = self
                    .mixture
                    .phase_state(&reduced, self.z, side)
                    .map_err(|error| CalcError::Flash(format!("{error}")))?
                    .z;
                let state = molar_enthalpy_entropy(
                    self.mixture,
                    self.ideal_gas,
                    kelvins(temperature),
                    pascals(pressure),
                    self.z,
                    root,
                )
                .map_err(|error| CalcError::Flash(format!("{error}")))?;
                let volume = root * MOLAR_GAS_CONSTANT * temperature / pressure;
                if volume <= 0.0 || !volume.is_finite() {
                    return Err(CalcError::Flash(format!(
                        "the single-phase volume at {pressure} Pa is {volume} m3/mol"
                    )));
                }
                (
                    state.s.value / self.molar_mass,
                    self.molar_mass / volume,
                    0.0,
                    0.0,
                )
            }
        };
        Ok(EquilibriumSample {
            pressure,
            temperature,
            density,
            entropy_residual: entropy - self.target,
            component_residual: balance,
            fugacity_residual: fugacity,
            phase: flash.phase,
            vapour_fraction: flash.vapour_fraction,
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
