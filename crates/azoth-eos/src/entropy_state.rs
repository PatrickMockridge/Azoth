//! One equilibrium state at a fixed specific entropy.
//!
//! `eos.equilibrium_sound_speed` and the HEM release model both re-equilibrate a fluid at a
//! fixed specific entropy: a pure fluid through its saturation quality, a mixture through a
//! bracketed TP root. This module is that one arithmetic, with the two roots' stopping
//! policies as parameters rather than two copies.
//!
//! **The policies are measured differences between the two classes.**
//! `EquilibriumSoundSpeed` returns early when its first trial already closes, stops a
//! bisection on a `1e-12 max(1, T)` bracket, and lets a non-finite residual fall through to
//! the "did not close" report; `HomogeneousEquilibriumReleaseModel` does none of the three
//! and throws on a non-finite residual. The acceptance thresholds stay with the callers:
//! the sound speed accepts an entropy residual of `1e-7` J/(kg K), the release model `1e-5`.

use azoth_core::units::{kelvins, pascals};

use crate::mixture::{Mixture, RootSide};
use crate::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use crate::pr_molar_volume::MOLAR_GAS_CONSTANT;
use crate::pt_flash::pt_flash;
use crate::pure_saturation::pure_saturation;
use crate::results::Phase;

/// The root's target residual, J/(kg K).
pub const ENTROPY_TOLERANCE: f64 = 1.0e-7;
/// The bisection's relative temperature tolerance.
pub const ROOT_TOLERANCE: f64 = 1.0e-12;
/// Entropy-bracket widenings before the bracket is refused.
pub const BRACKET_STEPS: usize = 30;
/// The bisection cap.
pub const BISECTION_STEPS: usize = 100;

/// One equilibrium state at the target specific entropy.
#[derive(Debug, Clone, PartialEq)]
pub struct EntropyState {
    /// Absolute pressure, Pa.
    pub pressure: f64,
    /// Temperature, K.
    pub temperature: f64,
    /// Mass over EOS total volume, kg/m3.
    pub density: f64,
    /// Molar enthalpy, J/mol.
    pub molar_enthalpy: f64,
    /// Molar entropy, J/(mol K).
    pub molar_entropy: f64,
    /// Signed specific-entropy error, J/(kg K).
    pub entropy_residual: f64,
    /// The largest component-inventory error over the total moles.
    pub component_residual: f64,
    /// The largest absolute interphase log-fugacity ratio; zero for one phase.
    pub fugacity_residual: f64,
    /// The largest phase-composition sum error; zero for one phase.
    pub phase_composition_residual: f64,
    /// What the flash found.
    pub phase: Phase,
    /// The vapour fraction where the state is two-phase.
    pub vapour_fraction: Option<f64>,
}

/// How one class stops its entropy root.
#[derive(Debug, Clone, Copy)]
pub struct RootPolicy {
    /// Return the first trial when it already meets [`ENTROPY_TOLERANCE`].
    pub accept_initial: bool,
    /// Stop a stalled bisection at this relative bracket width, or run to the cap.
    pub bracket_width_tolerance: Option<f64>,
    /// Refuse a non-finite residual as its own failure rather than as a missed root.
    pub nonfinite_is_error: bool,
}

impl RootPolicy {
    /// `EquilibriumSoundSpeed`'s root.
    pub const SOUND_SPEED: Self = Self {
        accept_initial: true,
        bracket_width_tolerance: Some(ROOT_TOLERANCE),
        nonfinite_is_error: false,
    };
    /// `HomogeneousEquilibriumReleaseModel`'s root.
    pub const HEM: Self = Self {
        accept_initial: false,
        bracket_width_tolerance: None,
        nonfinite_is_error: true,
    };
}

/// Why an entropy state could not be obtained.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct EntropyStateError(pub String);

impl std::fmt::Display for EntropyStateError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        formatter.write_str(&self.0)
    }
}

impl std::error::Error for EntropyStateError {}

/// The equilibrium state at `pressure` and the target specific entropy.
///
/// `evaluations` counts every TP flash and every accepted saturation-quality state, including
/// rejected trials, as both classes do.
///
/// # Errors
/// A non-positive pressure, a bracket that cannot be widened to contain the entropy, a
/// non-finite residual under a policy that refuses one, or a root that does not close.
#[allow(clippy::too_many_arguments)] // one argument per input the state reads
pub fn isentropic_state(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    z: &[f64],
    molar_mass: f64,
    target: f64,
    pressure: f64,
    start_temperature: f64,
    evaluations: &mut u32,
    policy: RootPolicy,
) -> Result<EntropyState, EntropyStateError> {
    if pressure <= 0.0 || !pressure.is_finite() {
        return Err(EntropyStateError(format!(
            "the pressure sample {pressure} Pa is not positive"
        )));
    }
    if mixture.len() == 1 {
        if let Some(state) = saturation_state(
            mixture,
            ideal_gas,
            z,
            molar_mass,
            target,
            pressure,
            evaluations,
        )? {
            return Ok(state);
        }
    }
    tp_root(
        mixture,
        ideal_gas,
        z,
        molar_mass,
        target,
        pressure,
        start_temperature,
        evaluations,
        policy,
    )
}

/// The pure fluid's saturation-quality state, where the target entropy lies between the
/// saturated phases. `None` sends the caller to the TP root.
#[allow(clippy::too_many_arguments)] // one argument per input the assembly reads
fn saturation_state(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    z: &[f64],
    molar_mass: f64,
    target: f64,
    pressure: f64,
    evaluations: &mut u32,
) -> Result<Option<EntropyState>, EntropyStateError> {
    let component = &mixture.components()[0];
    let tc = component.tc.value;
    let omega = component.omega;
    // **A tolerant bracket.** The cold end is where the cubic's liquid branch is hardest
    // to find, so the first fraction that answers is the bound rather than 0.3 Tc by fiat.
    let mut low_t = 0.0;
    let mut p_low = f64::NAN;
    for fraction in [0.3, 0.4, 0.5, 0.6] {
        let candidate = (1.0_f64).max(fraction * tc);
        if let Ok(result) = pure_saturation(component.tc, component.pc, omega, kelvins(candidate)) {
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
    let reduced = match mixture.reduced_parameters(kelvins(temperature), pascals(pressure)) {
        Ok(reduced) => reduced,
        Err(_) => return Ok(None),
    };
    let liquid_root = match mixture.phase_state(&reduced, z, RootSide::Liquid) {
        Ok(state) => state.z,
        Err(_) => return Ok(None),
    };
    let vapour_root = match mixture.phase_state(&reduced, z, RootSide::Vapour) {
        Ok(state) => state.z,
        Err(_) => return Ok(None),
    };
    if (vapour_root - liquid_root).abs() < 1.0e-12 {
        return Ok(None);
    }
    let liquid = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        kelvins(temperature),
        pascals(pressure),
        z,
        liquid_root,
    );
    let vapour = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        kelvins(temperature),
        pascals(pressure),
        z,
        vapour_root,
    );
    let (Ok(liquid), Ok(vapour)) = (liquid, vapour) else {
        return Ok(None);
    };
    let s_liquid = liquid.s.value / molar_mass;
    let s_vapour = vapour.s.value / molar_mass;
    if s_vapour <= s_liquid || !(s_liquid..=s_vapour).contains(&target) {
        return Ok(None);
    }
    *evaluations += 1;
    let quality = (target - s_liquid) / (s_vapour - s_liquid);
    // **The EOS volume, unshifted.** The class uses `getVolume` with volume correction
    // off in its own tests - a density correction is not part of the acoustic derivative.
    let v_liquid = liquid_root * MOLAR_GAS_CONSTANT * temperature / pressure;
    let v_vapour = vapour_root * MOLAR_GAS_CONSTANT * temperature / pressure;
    let volume = (1.0 - quality) * v_liquid + quality * v_vapour;
    if volume <= 0.0 || !volume.is_finite() {
        return Err(EntropyStateError(format!(
            "the saturation-quality volume at {pressure} Pa is {volume} m3/mol"
        )));
    }
    Ok(Some(EntropyState {
        pressure,
        temperature,
        density: molar_mass / volume,
        molar_enthalpy: (1.0 - quality) * liquid.h.value + quality * vapour.h.value,
        molar_entropy: (1.0 - quality) * liquid.s.value + quality * vapour.s.value,
        entropy_residual: (1.0 - quality) * s_liquid + quality * s_vapour - target,
        component_residual: 0.0,
        fugacity_residual: 0.0,
        phase_composition_residual: 0.0,
        phase: Phase::TwoPhase,
        vapour_fraction: Some(quality),
    }))
}

/// The equilibrium state at a pressure by the bracketed TP entropy root.
#[allow(clippy::too_many_arguments)] // one argument per input the root reads
fn tp_root(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    z: &[f64],
    molar_mass: f64,
    target: f64,
    pressure: f64,
    start_temperature: f64,
    evaluations: &mut u32,
    policy: RootPolicy,
) -> Result<EntropyState, EntropyStateError> {
    let guess = start_temperature;
    let initial = tp_at(
        mixture,
        ideal_gas,
        z,
        molar_mass,
        target,
        pressure,
        guess,
        evaluations,
    )?;
    if policy.accept_initial && initial.entropy_residual.abs() <= ENTROPY_TOLERANCE {
        return Ok(initial);
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
            low_error = tp_at(
                mixture,
                ideal_gas,
                z,
                molar_mass,
                target,
                pressure,
                lower,
                evaluations,
            )?
            .entropy_residual;
        }
        if high_error < 0.0 {
            upper = guess + span;
            high_error = tp_at(
                mixture,
                ideal_gas,
                z,
                molar_mass,
                target,
                pressure,
                upper,
                evaluations,
            )?
            .entropy_residual;
        }
        span *= 1.6;
    }
    if !(low_error <= 0.0 && high_error >= 0.0) {
        return Err(EntropyStateError(format!(
            "could not bracket the specific entropy at {pressure:.6} Pa"
        )));
    }
    let mut last_error = f64::NAN;
    for _ in 0..BISECTION_STEPS {
        let temperature = 0.5 * (lower + upper);
        let trial = tp_at(
            mixture,
            ideal_gas,
            z,
            molar_mass,
            target,
            pressure,
            temperature,
            evaluations,
        )?;
        last_error = trial.entropy_residual;
        if policy.nonfinite_is_error && !last_error.is_finite() {
            return Err(EntropyStateError("nonfinite entropy root".to_string()));
        }
        if last_error.abs() <= ENTROPY_TOLERANCE {
            return Ok(trial);
        }
        if last_error > 0.0 {
            upper = temperature;
        } else {
            lower = temperature;
        }
        if policy
            .bracket_width_tolerance
            .is_some_and(|tolerance| upper - lower <= tolerance * 1.0_f64.max(temperature))
        {
            break;
        }
    }
    Err(EntropyStateError(format!(
        "the entropy root did not close at {pressure:.6} Pa; residual = {last_error} J/(kg K)"
    )))
}

/// One TP flash and its assembled state.
#[allow(clippy::too_many_arguments)] // one argument per input the flash reads
fn tp_at(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    z: &[f64],
    molar_mass: f64,
    target: f64,
    pressure: f64,
    temperature: f64,
    evaluations: &mut u32,
) -> Result<EntropyState, EntropyStateError> {
    *evaluations += 1;
    let flash = pt_flash(mixture, kelvins(temperature), pascals(pressure), z)
        .map_err(|error| EntropyStateError(format!("{error}")))?;
    assemble(
        mixture,
        ideal_gas,
        z,
        molar_mass,
        target,
        pressure,
        temperature,
        &flash,
    )
}

/// Turn a flash into the assembled state, or refuse it.
#[allow(clippy::too_many_arguments)] // one argument per input the assembly reads
fn assemble(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    z: &[f64],
    molar_mass: f64,
    target: f64,
    pressure: f64,
    temperature: f64,
    flash: &crate::results::PtFlashResult,
) -> Result<EntropyState, EntropyStateError> {
    match flash.phase {
        Phase::TwoPhase => {
            let beta = flash.vapour_fraction.ok_or_else(|| {
                EntropyStateError(
                    "the flash split the feed and reported no vapour fraction".to_string(),
                )
            })?;
            if !(0.0..=1.0).contains(&beta) {
                return Err(EntropyStateError(format!(
                    "the vapour fraction {beta} is outside [0, 1]"
                )));
            }
            let liquid = molar_enthalpy_entropy(
                mixture,
                ideal_gas,
                kelvins(temperature),
                pascals(pressure),
                &flash.x,
                flash.z_liquid,
            );
            let vapour = molar_enthalpy_entropy(
                mixture,
                ideal_gas,
                kelvins(temperature),
                pascals(pressure),
                &flash.y,
                flash.z_vapour,
            );
            let (Ok(liquid), Ok(vapour)) = (liquid, vapour) else {
                return Err(EntropyStateError(
                    "the phase entropy is not finite".to_string(),
                ));
            };
            let molar_entropy = (1.0 - beta) * liquid.s.value + beta * vapour.s.value;
            let molar_enthalpy = (1.0 - beta) * liquid.h.value + beta * vapour.h.value;
            let v_liquid = flash.z_liquid * MOLAR_GAS_CONSTANT * temperature / pressure;
            let v_vapour = flash.z_vapour * MOLAR_GAS_CONSTANT * temperature / pressure;
            let volume = (1.0 - beta) * v_liquid + beta * v_vapour;
            if volume <= 0.0 || !volume.is_finite() {
                return Err(EntropyStateError(format!(
                    "the two-phase volume at {pressure} Pa is {volume} m3/mol"
                )));
            }
            let mut balance = 0.0_f64;
            let mut fugacity = 0.0_f64;
            for (index, fraction) in z.iter().enumerate() {
                let recovered = (1.0 - beta) * flash.x[index] + beta * flash.y[index];
                balance = balance.max((recovered - fraction).abs());
                if *fraction > 1.0e-12 {
                    if flash.x[index] <= 0.0 || flash.y[index] <= 0.0 {
                        return Err(EntropyStateError(
                            "a phase composition is not positive".to_string(),
                        ));
                    }
                    let ln_ratio = flash.x[index].ln() + flash.ln_phi_liquid[index]
                        - flash.y[index].ln()
                        - flash.ln_phi_vapour[index];
                    fugacity = fugacity.max(ln_ratio.abs());
                }
            }
            let phase_residual = (flash.x.iter().sum::<f64>() - 1.0)
                .abs()
                .max((flash.y.iter().sum::<f64>() - 1.0).abs());
            Ok(EntropyState {
                pressure,
                temperature,
                density: molar_mass / volume,
                molar_enthalpy,
                molar_entropy,
                entropy_residual: molar_entropy / molar_mass - target,
                component_residual: balance,
                fugacity_residual: fugacity,
                phase_composition_residual: phase_residual,
                phase: flash.phase,
                vapour_fraction: flash.vapour_fraction,
            })
        }
        Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => {
            let side = if matches!(flash.phase, Phase::AllVapour) {
                RootSide::Vapour
            } else {
                RootSide::Liquid
            };
            let reduced = mixture
                .reduced_parameters(kelvins(temperature), pascals(pressure))
                .map_err(|error| EntropyStateError(format!("{error}")))?;
            let root = mixture
                .phase_state(&reduced, z, side)
                .map_err(|error| EntropyStateError(format!("{error}")))?
                .z;
            let state = molar_enthalpy_entropy(
                mixture,
                ideal_gas,
                kelvins(temperature),
                pascals(pressure),
                z,
                root,
            )
            .map_err(|error| EntropyStateError(format!("{error}")))?;
            let volume = root * MOLAR_GAS_CONSTANT * temperature / pressure;
            if volume <= 0.0 || !volume.is_finite() {
                return Err(EntropyStateError(format!(
                    "the single-phase volume at {pressure} Pa is {volume} m3/mol"
                )));
            }
            Ok(EntropyState {
                pressure,
                temperature,
                density: molar_mass / volume,
                molar_enthalpy: state.h.value,
                molar_entropy: state.s.value,
                entropy_residual: state.s.value / molar_mass - target,
                component_residual: 0.0,
                fugacity_residual: 0.0,
                phase_composition_residual: 0.0,
                phase: flash.phase,
                vapour_fraction: None,
            })
        }
    }
}
