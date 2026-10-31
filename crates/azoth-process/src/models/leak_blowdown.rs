//! `process.leak_blowdown` - NeqSim's legacy vessel-blowdown time loop.
//!
//! Spec: `specs/models/process/leak_blowdown.toml`. The port source is
//! `process/safety/release/LeakModel.java`, through `calculateSourceTerm(duration, timeStep)`.
//!
//! The class records the simplified orifice rate, the capped jet velocity, the unclamped
//! momentum and the droplet SMD, then depletes the inventory and updates the pressure and
//! temperature isentropically. Its scalar methods are this arithmetic, not separate ids; the
//! modern `calculateReleaseFlow` path is `process.hem_release`.

use std::f64::consts::PI;

use azoth_core::units::{
    Force, Length, Mass, MassRate, Pressure, ThermodynamicTemperature, Time, Velocity, Volume,
    kelvins, kilograms, kilograms_per_second, meters, meters_per_second, newtons, pascals, seconds,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use azoth_eos::cubic::Cubic;
use azoth_eos::databank;
use azoth_eos::entropy_state::molar_mass_of;
use azoth_eos::heat_capacity::phase_capacity;
use azoth_eos::mixture::{Mixture, ReducedParameters, RootSide};
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::pr_molar_volume::MOLAR_GAS_CONSTANT;
use azoth_eos::pt_flash::pt_flash;
use azoth_eos::results::{Phase, PtFlashResult};

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;

/// The class's fixed gas constant, J/(kmol K).
const R: f64 = 8314.0;
/// The class's surrogate surface tension for the droplet correlation, N/m.
const SURFACE_TENSION: f64 = 0.02;
/// The class's Nukiyama-Tanasawa coefficient.
const SMD_COEFFICIENT: f64 = 585.0;
/// The class's liquid-fraction exponent.
const SMD_LIQUID_EXPONENT: f64 = 0.45;
/// The class's liquid-fraction gate.
const SMD_LIQUID_FLOOR: f64 = 0.001;
/// The inventory floor, relative to the initial mass.
const MASS_FLOOR: f64 = 0.01;
/// The declared-point cap; a finer request is refused rather than run for hours.
const MAX_POINTS: u32 = 100_000;
/// The velocity cap, as a multiple of the ideal-gas sound speed.
const VELOCITY_CAP: f64 = 1.2;

/// One recorded step of the blowdown.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LeakBlowdownPoint {
    /// The step's time.
    #[serde(serialize_with = "scalar")]
    pub time: Time,
    /// The simplified orifice rate at the state.
    #[serde(serialize_with = "scalar")]
    pub mass_flow_rate: MassRate,
    /// The state's temperature.
    #[serde(serialize_with = "scalar")]
    pub temperature: ThermodynamicTemperature,
    /// The state's pressure.
    #[serde(serialize_with = "scalar")]
    pub pressure: Pressure,
    /// The physical vapour mole fraction: 1 all-vapour, 0 all-liquid, beta two-phase.
    pub vapour_mole_fraction: f64,
    /// `mdot/(rho A Cd)`, capped at 1.2 times the ideal-gas sound speed.
    #[serde(serialize_with = "scalar")]
    pub jet_velocity: Velocity,
    /// `mdot v + (P - back) A`, unclamped.
    #[serde(serialize_with = "scalar")]
    pub jet_momentum: Force,
    /// The liquid Sauter mean diameter, zero for a single phase or a trace liquid.
    #[serde(serialize_with = "scalar")]
    pub droplet_smd: Length,
}

/// One diagnostic row.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LeakBlowdownDiagnostic {
    /// The stable machine code.
    pub code: String,
    /// The engineering detail.
    pub message: String,
}

/// Result of `process.leak_blowdown`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct LeakBlowdownResult {
    /// `valid`, `valid_with_warnings`, `invalid` or `unsupported`.
    pub status: String,
    /// The scenario label, echoed.
    pub scenario_name: String,
    /// The opening the calculation used.
    #[serde(serialize_with = "scalar")]
    pub hole_diameter: Length,
    /// The declared point count, `ceil(duration/time_step)+1`.
    pub number_of_points: u32,
    /// The visited rows.
    pub recorded_points: u32,
    /// The recorded series.
    pub points: Vec<LeakBlowdownPoint>,
    /// The summed `mdot dt`; absent on failure.
    #[serde(serialize_with = "optional_scalar")]
    pub total_mass_released: Option<Mass>,
    /// The largest recorded rate; absent on failure.
    #[serde(serialize_with = "optional_scalar")]
    pub peak_mass_flow_rate: Option<MassRate>,
    /// The class's `totalReleased/(0.5 peakRate)`; absent on failure.
    #[serde(serialize_with = "optional_scalar")]
    pub time_to_empty: Option<Time>,
    /// The model's records, in order.
    pub diagnostics: Vec<LeakBlowdownDiagnostic>,
    /// The diagnostics' machine codes alone.
    pub diagnostic_codes: Vec<String>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// A state failure the class would throw out of the loop; the port reports it.
#[derive(Debug, Clone)]
enum BlowdownFailure {
    Flash(String),
    Density(String),
    MolarMass(String),
    Gamma(String),
    Droplet(String),
}

impl BlowdownFailure {
    fn code(&self) -> &'static str {
        match self {
            Self::Flash(_) => "FLASH_FAILED",
            Self::Density(_) => "DENSITY_UNAVAILABLE",
            Self::MolarMass(_) => "MOLAR_MASS_UNAVAILABLE",
            Self::Gamma(_) => "GAMMA_OUT_OF_RANGE",
            Self::Droplet(_) => "DROPLET_SMD_NOT_FINITE",
        }
    }

    fn message(&self) -> &str {
        match self {
            Self::Flash(message)
            | Self::Density(message)
            | Self::MolarMass(message)
            | Self::Gamma(message)
            | Self::Droplet(message) => message,
        }
    }
}

/// One flashed state and the quantities the class reads from it.
struct State {
    flash: PtFlashResult,
    pressure: f64,
    temperature: f64,
    /// The corrected total density, kg/m3.
    density: f64,
    /// The class's `Cp/Cv`, summed over the phases.
    gamma: f64,
    /// The two-phase vapour fraction; `None` for a single phase.
    vapour_fraction: Option<f64>,
    /// The vapour-phase corrected density, when two phases exist.
    vapour_density: Option<f64>,
    /// The liquid-phase corrected density, when two phases exist.
    liquid_density: Option<f64>,
}

impl State {
    /// The physical vapour mole fraction the row reports.
    fn vapour_mole_fraction(&self) -> f64 {
        match self.flash.phase {
            Phase::TwoPhase => self.vapour_fraction.unwrap_or(0.0),
            Phase::AllVapour => 1.0,
            Phase::AllLiquid | Phase::Trivial => 0.0,
        }
    }
}

/// Blow a vessel down through a circular opening.
///
/// # Errors
/// `InvalidInput` for a state or opening the class refuses at construction, a composition
/// that is not one entry per component and sums to one, a non-positive duration or step, or a
/// declared point count over the port's cap.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn leak_blowdown(
    components: &[String],
    initial_z: &[f64],
    initial_t: ThermodynamicTemperature,
    initial_p: Pressure,
    duration: Time,
    vessel_volume: Option<Volume>,
    hole_diameter: Option<Length>,
    discharge_coefficient: Option<f64>,
    back_pressure: Option<Pressure>,
    time_step: Option<Time>,
    scenario_name: Option<&str>,
) -> Result<LeakBlowdownResult> {
    let spec = &model_gen::LEAK_BLOWDOWN_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "initial_t" => Some(initial_t.value),
            "initial_p" => Some(initial_p.value),
            "vessel_volume" => vessel_volume.map(|v| v.value),
            "hole_diameter" => hole_diameter.map(|d| d.value),
            "discharge_coefficient" => discharge_coefficient,
            "back_pressure" => back_pressure.map(|p| p.value),
            "duration" => Some(duration.value),
            "time_step" => time_step.map(|d| d.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let volume = vessel_volume.map(|v| v.value).unwrap_or(1.0);
    let diameter = hole_diameter.unwrap_or_else(|| meters(0.01));
    let coefficient = discharge_coefficient.unwrap_or(0.62);
    let back = back_pressure.unwrap_or_else(|| pascals(101_325.0)).value;
    let span = duration;
    let step = time_step.unwrap_or_else(|| seconds(1.0));
    let scenario = scenario_name.unwrap_or("Release").to_string();

    if !volume.is_finite() || volume <= 0.0 {
        return Err(AzothError::invalid_input(
            "vessel_volume",
            format!("a vessel has a positive volume, and {volume} m3 is not"),
        ));
    }
    if !step.value.is_finite() || step.value <= 0.0 || !span.value.is_finite() || span.value <= 0.0
    {
        return Err(AzothError::invalid_input(
            "duration",
            format!(
                "a positive finite duration and step are required, and they are {} s and {} s",
                span.value, step.value
            ),
        ));
    }
    if !coefficient.is_finite() || coefficient <= 0.0 || coefficient > 1.0 {
        return Err(AzothError::invalid_input(
            "discharge_coefficient",
            format!("the effective-area factor must be in (0, 1], and {coefficient} is not"),
        ));
    }
    let declared = (span.value / step.value).ceil();
    if !declared.is_finite() || declared < 0.0 || declared > f64::from(MAX_POINTS - 1) {
        return Err(AzothError::invalid_input(
            "duration",
            format!("the declared point count is {declared}+1, and the port's cap is {MAX_POINTS}"),
        ));
    }
    let number_of_points = declared as u32 + 1;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
    if initial_z.len() != mixture.len() {
        return Err(AzothError::invalid_input(
            "initial_z",
            format!(
                "a mixture of {} components needs that many mole fractions, and {initial_z:?} \
                 has {}",
                mixture.len(),
                initial_z.len()
            ),
        ));
    }
    azoth_eos::phase_boundary::check_composition(initial_z, mixture.len(), "initial_z")?;
    if initial_z.iter().any(|value| !value.is_finite()) {
        return Err(AzothError::invalid_input(
            "initial_z",
            "a mole fraction must be finite",
        ));
    }

    let mut diagnostics = vec![LeakBlowdownDiagnostic {
        code: "MODEL_ASSUMPTIONS".to_string(),
        message: "Rigid adiabatic vessel; equilibrium TP flashes; simplified orifice; no slip \
                  or friction in the jet"
            .to_string(),
    }];
    if initial_p.value <= back {
        diagnostics.push(LeakBlowdownDiagnostic {
            code: "NO_FORWARD_FLOW".to_string(),
            message: "The initial pressure is at or below the receiving pressure".to_string(),
        });
    }

    let mut state = match flash_state(
        &mixture,
        &ideal_gas,
        initial_t.value,
        initial_p.value,
        initial_z,
    ) {
        Ok(state) => state,
        Err(failure) => {
            return Ok(failure_result(
                &scenario,
                diameter,
                number_of_points,
                failure,
                warnings,
            ));
        }
    };
    let initial_mass = state.density * volume;
    if !initial_mass.is_finite() || initial_mass <= 0.0 {
        return Ok(failure_result(
            &scenario,
            diameter,
            number_of_points,
            BlowdownFailure::Density(format!(
                "the initial density is {} kg/m3 and the volume {volume} m3",
                state.density
            )),
            warnings,
        ));
    }

    let area = PI * diameter.value * diameter.value / 4.0;
    let mut mass = initial_mass;
    let mut total = 0.0;
    let mut peak: f64 = 0.0;
    let mut points: Vec<LeakBlowdownPoint> = Vec::new();
    let mut floor_reached = false;
    let mut mass_floor = false;

    for index in 0..number_of_points {
        if mass <= MASS_FLOOR * initial_mass || mass.is_nan() {
            mass_floor = true;
            break;
        }
        let time = f64::from(index) * step.value;
        let (row, mdot) = match record_point(
            &mixture,
            &state,
            initial_z,
            back,
            diameter.value,
            coefficient,
            area,
            time,
        ) {
            Ok(row) => row,
            Err(failure) => {
                return Ok(failure_result(
                    &scenario,
                    diameter,
                    number_of_points,
                    failure,
                    warnings,
                ));
            }
        };
        peak = peak.max(mdot);
        points.push(row);

        mass -= mdot * step.value;
        total += mdot * step.value;
        if mass <= MASS_FLOOR * initial_mass || mass.is_nan() {
            mass_floor = true;
            break;
        }

        let ratio = (mass / initial_mass).powf(state.gamma);
        let new_pressure = state.pressure * ratio;
        let floored = new_pressure.max(back);
        if floored > new_pressure {
            floor_reached = true;
        }
        let new_temperature =
            state.temperature * (floored / state.pressure).powf((state.gamma - 1.0) / state.gamma);
        state = match flash_state(&mixture, &ideal_gas, new_temperature, floored, initial_z) {
            Ok(state) => state,
            Err(failure) => {
                return Ok(failure_result(
                    &scenario,
                    diameter,
                    number_of_points,
                    failure,
                    warnings,
                ));
            }
        };
    }

    if floor_reached {
        diagnostics.push(LeakBlowdownDiagnostic {
            code: "PRESSURE_FLOOR_REACHED".to_string(),
            message: "The isentropic update was floored at the receiving pressure".to_string(),
        });
    }
    if mass_floor {
        diagnostics.push(LeakBlowdownDiagnostic {
            code: "MASS_FLOOR_REACHED".to_string(),
            message: "The loop ended below 1% of the initial inventory".to_string(),
        });
    }
    let recorded_points = points.len() as u32;
    let time_to_empty = if peak > 0.0 {
        total / (0.5 * peak)
    } else {
        0.0
    };
    Ok(LeakBlowdownResult {
        status: "valid".to_string(),
        scenario_name: scenario,
        hole_diameter: diameter,
        number_of_points,
        recorded_points,
        points,
        total_mass_released: Some(kilograms(total)),
        peak_mass_flow_rate: Some(kilograms_per_second(peak)),
        time_to_empty: Some(seconds(time_to_empty)),
        diagnostic_codes: diagnostics.iter().map(|row| row.code.clone()).collect(),
        diagnostics,
        warnings,
    })
}

/// One step's row and the rate, from the current flashed state.
#[allow(clippy::too_many_arguments)] // one argument per state the step reads
fn record_point(
    mixture: &Mixture,
    state: &State,
    z: &[f64],
    back: f64,
    diameter: f64,
    coefficient: f64,
    area: f64,
    time: f64,
) -> std::result::Result<(LeakBlowdownPoint, f64), BlowdownFailure> {
    let molar_mass = molar_mass_of(mixture, z)
        .map_err(|error| BlowdownFailure::MolarMass(format!("{error}")))?;
    let mw = 1000.0 * molar_mass;
    let mdot = if state.pressure <= back {
        0.0
    } else {
        let critical = (2.0 / (state.gamma + 1.0)).powf(state.gamma / (state.gamma - 1.0));
        let ratio = back / state.pressure;
        if ratio <= critical {
            let factor =
                (2.0 / (state.gamma + 1.0)).powf((state.gamma + 1.0) / (2.0 * (state.gamma - 1.0)));
            coefficient
                * area
                * state.pressure
                * (state.gamma * mw / (R * state.temperature)).sqrt()
                * factor
        } else {
            let drop = state.pressure - back;
            coefficient * area * (2.0 * state.density * drop).sqrt()
        }
    };
    if !mdot.is_finite() || mdot < 0.0 {
        return Err(BlowdownFailure::Flash(format!(
            "the rate at {:.6} Pa is {mdot} kg/s",
            state.pressure
        )));
    }
    let velocity = if mdot <= 0.0 {
        0.0
    } else {
        let raw = mdot / (state.density * area * coefficient);
        let sound = (state.gamma * R * state.temperature / mw).sqrt();
        raw.min(VELOCITY_CAP * sound)
    };
    let momentum = mdot * velocity + (state.pressure - back) * area;
    let smd = droplet_smd(state, diameter, velocity)?;
    Ok((
        LeakBlowdownPoint {
            time: seconds(time),
            mass_flow_rate: kilograms_per_second(mdot),
            temperature: kelvins(state.temperature),
            pressure: pascals(state.pressure),
            vapour_mole_fraction: state.vapour_mole_fraction(),
            jet_velocity: meters_per_second(velocity),
            jet_momentum: newtons(momentum),
            droplet_smd: meters(smd),
        },
        mdot,
    ))
}

/// The class's liquid droplet Sauter mean diameter, m.
fn droplet_smd(
    state: &State,
    diameter: f64,
    velocity: f64,
) -> std::result::Result<f64, BlowdownFailure> {
    let Some(beta) = state.vapour_fraction else {
        return Ok(0.0);
    };
    let liquid_fraction = 1.0 - beta;
    if liquid_fraction < SMD_LIQUID_FLOOR {
        return Ok(0.0);
    }
    let (Some(vapour), Some(liquid)) = (state.vapour_density, state.liquid_density) else {
        return Ok(0.0);
    };
    if vapour <= 0.0 || liquid <= 0.0 || vapour.is_nan() || liquid.is_nan() {
        return Ok(0.0);
    }
    let weber = vapour * velocity * velocity * diameter / SURFACE_TENSION;
    if weber.is_nan() {
        return Err(BlowdownFailure::Droplet(
            "the droplet Weber number is NaN".to_string(),
        ));
    }
    let smd = SMD_COEFFICIENT * diameter / weber.sqrt()
        * (liquid / vapour).sqrt()
        * liquid_fraction.powf(SMD_LIQUID_EXPONENT);
    if smd.is_nan() {
        return Err(BlowdownFailure::Droplet(
            "the droplet Sauter mean diameter is NaN".to_string(),
        ));
    }
    Ok(smd.min(diameter / 2.0))
}

/// Flash one state and read the density, gamma and phase densities the class uses.
fn flash_state(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    temperature: f64,
    pressure: f64,
    z: &[f64],
) -> std::result::Result<State, BlowdownFailure> {
    let t = kelvins(temperature);
    let p = pascals(pressure);
    let reduced = mixture
        .reduced_parameters(t, p)
        .map_err(|error| BlowdownFailure::Flash(format!("{error}")))?;
    let flash =
        pt_flash(mixture, t, p, z).map_err(|error| BlowdownFailure::Flash(format!("{error}")))?;
    let (density, gamma, vapour_fraction, vapour_density, liquid_density) = match flash.phase {
        Phase::TwoPhase => {
            let beta = flash.vapour_fraction.ok_or_else(|| {
                BlowdownFailure::Flash(
                    "the flash split the feed and reported no vapour fraction".to_string(),
                )
            })?;
            let vapour = phase_of(mixture, ideal_gas, &reduced, t, p, &flash.y, flash.z_vapour)?;
            let liquid = phase_of(mixture, ideal_gas, &reduced, t, p, &flash.x, flash.z_liquid)?;
            let mass = beta * vapour.molar_mass + (1.0 - beta) * liquid.molar_mass;
            let volume = beta * vapour.molar_volume + (1.0 - beta) * liquid.molar_volume;
            let density = mass / volume;
            let cp = beta * vapour.cp + (1.0 - beta) * liquid.cp;
            let cv = beta * vapour.cv + (1.0 - beta) * liquid.cv;
            (
                density,
                cp / cv,
                Some(beta),
                Some(vapour.density),
                Some(liquid.density),
            )
        }
        Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => {
            if matches!(flash.phase, Phase::Trivial) {
                return Err(BlowdownFailure::Flash(format!(
                    "the flash found a trivial solution at {pressure} Pa"
                )));
            }
            let side = if matches!(flash.phase, Phase::AllVapour) {
                RootSide::Vapour
            } else {
                RootSide::Liquid
            };
            let root = mixture
                .phase_state(&reduced, z, side)
                .map_err(|error| BlowdownFailure::Flash(format!("{error}")))?
                .z;
            let phase = phase_of(mixture, ideal_gas, &reduced, t, p, z, root)?;
            (phase.density, phase.cp / phase.cv, None, None, None)
        }
    };
    if !density.is_finite() || density <= 0.0 {
        return Err(BlowdownFailure::Density(format!(
            "the density at {pressure} Pa is {density} kg/m3"
        )));
    }
    if !gamma.is_finite() || gamma <= 1.0 || gamma > 2.0 {
        return Err(BlowdownFailure::Gamma(format!(
            "the capacity ratio at {pressure} Pa is {gamma}, and the port refuses one outside \
             (1, 2]"
        )));
    }
    Ok(State {
        flash,
        pressure,
        temperature,
        density,
        gamma,
        vapour_fraction,
        vapour_density,
        liquid_density,
    })
}

/// One phase's molar mass, corrected molar volume, density and capacities.
struct PhaseProperties {
    molar_mass: f64,
    molar_volume: f64,
    density: f64,
    cp: f64,
    cv: f64,
}

#[allow(clippy::too_many_arguments)] // one argument per phase input
fn phase_of(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    reduced: &ReducedParameters,
    t: ThermodynamicTemperature,
    p: Pressure,
    composition: &[f64],
    root: f64,
) -> std::result::Result<PhaseProperties, BlowdownFailure> {
    let molar_mass = molar_mass_of(mixture, composition)
        .map_err(|error| BlowdownFailure::MolarMass(format!("{error}")))?;
    let shift = mixture.volume_shift(composition);
    let molar_volume = root * MOLAR_GAS_CONSTANT * t.value / p.value - shift;
    if !molar_volume.is_finite() || molar_volume <= 0.0 {
        return Err(BlowdownFailure::Density(format!(
            "the corrected molar volume at {:.6} Pa is {molar_volume} m3/mol",
            p.value
        )));
    }
    let capacity = phase_capacity(mixture, ideal_gas, reduced, t, p, composition, root)
        .map_err(|error| BlowdownFailure::Gamma(format!("{error}")))?;
    Ok(PhaseProperties {
        molar_mass,
        molar_volume,
        density: molar_mass / molar_volume,
        cp: capacity.cp,
        cv: capacity.cv,
    })
}

/// The class's failure result: no series, no totals, one diagnostic.
fn failure_result(
    scenario: &str,
    diameter: Length,
    number_of_points: u32,
    failure: BlowdownFailure,
    warnings: Vec<Warning>,
) -> LeakBlowdownResult {
    LeakBlowdownResult {
        status: "invalid".to_string(),
        scenario_name: scenario.to_string(),
        hole_diameter: diameter,
        number_of_points,
        recorded_points: 0,
        points: Vec::new(),
        total_mass_released: None,
        peak_mass_flow_rate: None,
        time_to_empty: None,
        diagnostics: vec![LeakBlowdownDiagnostic {
            code: failure.code().to_string(),
            message: failure.message().to_string(),
        }],
        diagnostic_codes: vec![failure.code().to_string()],
        warnings,
    }
}
