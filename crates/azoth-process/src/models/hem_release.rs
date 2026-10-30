//! `process.hem_release` - NeqSim's homogeneous-equilibrium orifice release.
//!
//! Spec: `specs/models/process/hem_release.toml`. The port source is
//! `process/safety/release/HomogeneousEquilibriumReleaseModel.java`, through its `calculate`.
//!
//! The model samples an isentropic expansion from the upstream stagnation state down to the
//! receiving pressure, maximises `rho sqrt(2 (h0 - h))` over pressure, and reports the throat
//! as the flow-limiting station. It is the consumer the equilibrium sound speed was ported
//! for: the throat's Mach check reads `eos.equilibrium_sound_speed`.
//!
//! **The arithmetic is the shared [`azoth_eos::entropy_state`] helper**, with the release
//! model's own root policy and the class's own closure thresholds.

use azoth_core::units::{
    Length, MassDensity, MassFlux, MassRate, MolarEnergy, MolarHeatCapacity, Pressure,
    ThermodynamicTemperature, Velocity, joules_per_kilogram_kelvin, joules_per_mole,
    joules_per_mole_kelvin, kelvins, kilograms_per_cubic_meter, kilograms_per_second,
    kilograms_per_square_meter_second, meters_per_second, pascals,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use azoth_eos::cubic::Cubic;
use azoth_eos::databank;
use azoth_eos::entropy_state::{self, EntropyState, RootPolicy};
use azoth_eos::equilibrium_sound_speed::{EquilibriumStatus, equilibrium_speed_of};
use azoth_eos::mixture::Mixture;
use azoth_eos::molar_enthalpy_entropy::IdealGasModel;
use azoth_eos::phase_boundary::check_composition;
use azoth_eos::results::Phase;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::model_gen;

/// The class's pressure samples between the receiving pressure and the upstream state.
const SAMPLES: usize = 48;
/// `PRESSURE_TOLERANCE`: the golden section's bracket as a fraction of the upstream pressure.
const PRESSURE_TOLERANCE: f64 = 1.0e-7;
/// The golden section's iteration cap.
const REFINE_ITERATIONS: usize = 80;
/// `sample()`'s accepted specific-entropy residual, J/(kg K).
const ENTROPY_CLOSURE: f64 = 1.0e-5;
/// `checkEquilibrium`'s component-inventory residual, over the total moles.
const BALANCE_CLOSURE: f64 = 1.0e-8;
/// `checkEquilibrium`'s interphase log-fugacity residual.
const FUGACITY_CLOSURE: f64 = 1.0e-5;
/// The smallest enthalpy drop that is still an expansion, J/kg.
const ENTHALPY_FLOOR: f64 = -1.0e-5;
/// The choked-flow test's margin over the receiving pressure.
const CHOKED_MARGIN: f64 = 1.0e-5;
/// The class's Mach band for the throat sound speed.
const MACH_BAND: f64 = 0.05;
/// The CO2 solid-capable guard's temperature, K.
const CO2_SOLID_K: f64 = 216.592;
/// The relative pressure step the throat sound speed is evaluated at.
const SOUND_SPEED_RELATIVE_STEP: f64 = 1.0e-3;

/// The four stations the class reports.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ReleaseStation {
    /// The upstream stagnation state, at zero velocity.
    UpstreamStagnation,
    /// The maximum-flux state, or the receiving-pressure state when the flow is not choked.
    ThroatCritical,
    /// The physical opening; a zero-length orifice aliases it to the throat.
    OrificeExit,
    /// The isentropic receiving-pressure state.
    AmbientExpanded,
}

impl ReleaseStation {
    /// The spec's value for this station.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::UpstreamStagnation => "upstream_stagnation",
            Self::ThroatCritical => "throat_critical",
            Self::OrificeExit => "orifice_exit",
            Self::AmbientExpanded => "ambient_expanded",
        }
    }
}

/// The station's phase assemblage, in this port's three words.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ReleasePhase {
    /// One vapour phase.
    AllVapour,
    /// One liquid phase.
    AllLiquid,
    /// A genuine split.
    TwoPhase,
}

impl ReleasePhase {
    /// The spec's value for this phase.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::AllVapour => "all_vapour",
            Self::AllLiquid => "all_liquid",
            Self::TwoPhase => "two_phase",
        }
    }

    /// This port's word for one of the flash's phases, or `None` for a trivial solution.
    #[must_use]
    pub fn from_phase(phase: Phase) -> Option<Self> {
        match phase {
            Phase::AllVapour => Some(Self::AllVapour),
            Phase::AllLiquid => Some(Self::AllLiquid),
            Phase::TwoPhase => Some(Self::TwoPhase),
            Phase::Trivial => None,
        }
    }
}

/// One row of `stations`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HemReleaseStation {
    /// Which station this row is.
    pub station: ReleaseStation,
    /// Absolute pressure.
    #[serde(serialize_with = "scalar")]
    pub pressure: Pressure,
    /// Temperature.
    #[serde(serialize_with = "scalar")]
    pub temperature: ThermodynamicTemperature,
    /// Mass over EOS total volume.
    #[serde(serialize_with = "scalar")]
    pub density: MassDensity,
    /// Molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub molar_enthalpy: MolarEnergy,
    /// Molar entropy.
    #[serde(serialize_with = "scalar")]
    pub molar_entropy: MolarHeatCapacity,
    /// The station's axial velocity.
    #[serde(serialize_with = "scalar")]
    pub velocity: Velocity,
    /// The station's mass flux, before the discharge coefficient.
    #[serde(serialize_with = "scalar")]
    pub mass_flux: MassFlux,
    /// What the flash found.
    pub phase: ReleasePhase,
    /// The gas phase's share of the mass.
    pub gas_mass_fraction: f64,
}

/// One row of `diagnostics`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HemReleaseDiagnostic {
    /// The stable machine code.
    pub code: String,
    /// The engineering detail.
    pub message: String,
}

/// Result of `process.hem_release`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HemReleaseResult {
    /// `valid`, `valid_with_warnings`, `invalid` or `unsupported`.
    pub status: String,
    /// The release rate; absent when the result is not usable.
    #[serde(serialize_with = "optional_scalar")]
    pub mass_flow_rate: Option<MassRate>,
    /// Whether the accepted flow is choked.
    pub choked: bool,
    /// The throat's equilibrium sound speed, when the flow is choked and it converged.
    #[serde(serialize_with = "optional_scalar")]
    pub throat_sound_speed: Option<Velocity>,
    /// The upstream overall mole fractions.
    pub component_mole_fractions: Vec<f64>,
    /// The upstream overall mass fractions.
    pub component_mass_fractions: Vec<f64>,
    /// The resolved stations, upstream first.
    pub stations: Vec<HemReleaseStation>,
    /// The model's assumption and warning records.
    pub diagnostics: Vec<HemReleaseDiagnostic>,
    /// The diagnostics' machine codes alone.
    pub diagnostic_codes: Vec<String>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// One sample of the isentropic expansion.
#[derive(Debug, Clone)]
struct SampleState {
    state: EntropyState,
    velocity: f64,
    flux: f64,
}

/// A failure the class reports inside its result rather than as an exception.
#[derive(Debug, Clone)]
enum ReleaseFailure {
    /// The requested regime is outside the model.
    Unsupported(String),
    /// A required numerical or closure check failed.
    Invalid(String),
}

impl From<entropy_state::EntropyStateError> for ReleaseFailure {
    fn from(error: entropy_state::EntropyStateError) -> Self {
        Self::Invalid(error.0)
    }
}

/// The state one release sample reads.
struct Context<'a> {
    mixture: &'a Mixture,
    ideal_gas: &'a IdealGasModel,
    z: &'a [f64],
    names: &'a [String],
    molar_mass: f64,
    target: f64,
    upstream_h: f64,
    upstream_pressure: f64,
    start_temperature: f64,
}

impl Context<'_> {
    /// One isentropic sample at a pressure.
    fn sample(
        &self,
        pressure: f64,
        evaluations: &mut u32,
    ) -> std::result::Result<SampleState, ReleaseFailure> {
        let state = entropy_state::isentropic_state(
            self.mixture,
            self.ideal_gas,
            self.z,
            self.molar_mass,
            self.target,
            pressure,
            self.start_temperature,
            evaluations,
            RootPolicy::HEM,
        )?;
        check_closures(&state)?;
        check_support(&state, self.names, self.z)?;
        let drop = (self.upstream_h - state.molar_enthalpy) / self.molar_mass;
        if !drop.is_finite() || drop < ENTHALPY_FLOOR {
            return Err(ReleaseFailure::Invalid(
                "Expansion increased enthalpy or returned a nonfinite value".to_string(),
            ));
        }
        let velocity = (2.0 * drop.max(0.0)).sqrt();
        let flux = state.density * velocity;
        if !flux.is_finite() || flux < 0.0 {
            return Err(ReleaseFailure::Invalid(format!(
                "the mass flux at {pressure} Pa is {flux} kg/(m2 s)"
            )));
        }
        Ok(SampleState {
            state,
            velocity,
            flux,
        })
    }

    /// The golden section over a sampled local maximum, in linear pressure.
    fn refine(
        &self,
        mut low: f64,
        mut high: f64,
        evaluations: &mut u32,
    ) -> std::result::Result<SampleState, ReleaseFailure> {
        let ratio = (5.0_f64.sqrt() - 1.0) / 2.0;
        let mut left = self.sample(high - ratio * (high - low), evaluations)?;
        let mut right = self.sample(low + ratio * (high - low), evaluations)?;
        for _ in 0..REFINE_ITERATIONS {
            if high - low <= PRESSURE_TOLERANCE * self.upstream_pressure {
                return Ok(if left.flux > right.flux { left } else { right });
            }
            if left.flux < right.flux {
                low = left.state.pressure;
                left = right;
                right = self.sample(low + ratio * (high - low), evaluations)?;
            } else {
                high = right.state.pressure;
                right = left;
                left = self.sample(high - ratio * (high - low), evaluations)?;
            }
        }
        Err(ReleaseFailure::Invalid(
            "Critical pressure search did not converge".to_string(),
        ))
    }
}

/// The class's per-sample closure checks.
fn check_closures(state: &EntropyState) -> std::result::Result<(), ReleaseFailure> {
    if state.entropy_residual.abs() > ENTROPY_CLOSURE {
        return Err(ReleaseFailure::Invalid(format!(
            "Specific entropy closure failed: {} J/(kg K)",
            state.entropy_residual
        )));
    }
    if state.component_residual > BALANCE_CLOSURE {
        return Err(ReleaseFailure::Invalid(
            "Component inventory closure failed".to_string(),
        ));
    }
    if state.phase_composition_residual > BALANCE_CLOSURE {
        return Err(ReleaseFailure::Invalid(
            "Phase composition does not close".to_string(),
        ));
    }
    if state.fugacity_residual > FUGACITY_CLOSURE {
        return Err(ReleaseFailure::Invalid(
            "Interphase fugacity closure failed".to_string(),
        ));
    }
    if state.density <= 0.0
        || !state.density.is_finite()
        || state.temperature <= 0.0
        || !state.temperature.is_finite()
        || state.pressure <= 0.0
        || !state.pressure.is_finite()
    {
        return Err(ReleaseFailure::Invalid(format!(
            "the sample at {:.6} Pa is not a finite positive state",
            state.pressure
        )));
    }
    Ok(())
}

/// The class's upstream unsupported-regime check, before the first flash.
fn check_input_support(
    temperature: f64,
    names: &[String],
    z: &[f64],
) -> std::result::Result<(), ReleaseFailure> {
    for (name, fraction) in names.iter().zip(z) {
        if name.eq_ignore_ascii_case("CO2") && *fraction > 0.0 && temperature < CO2_SOLID_K {
            return Err(ReleaseFailure::Unsupported(
                "CO2 below 216.592 K requires a separately assessed solid-capable model"
                    .to_string(),
            ));
        }
    }
    Ok(())
}

/// The class's unsupported-regime checks this port can reach.
fn check_support(
    state: &EntropyState,
    names: &[String],
    z: &[f64],
) -> std::result::Result<(), ReleaseFailure> {
    if matches!(state.phase, Phase::Trivial) {
        return Err(ReleaseFailure::Invalid(format!(
            "the flash found a trivial solution at {:.6} Pa; the phase is not established",
            state.pressure
        )));
    }
    for (name, fraction) in names.iter().zip(z) {
        if name.eq_ignore_ascii_case("CO2") && *fraction > 0.0 && state.temperature < CO2_SOLID_K {
            return Err(ReleaseFailure::Unsupported(
                "CO2 below 216.592 K requires a separately assessed solid-capable model"
                    .to_string(),
            ));
        }
    }
    Ok(())
}

/// The release rate of a homogeneous-equilibrium short-opening orifice.
///
/// # Errors
/// `InvalidInput` for a state or opening the class refuses at construction, a composition
/// that is not one entry per component and sums to one, or whatever the databank refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seven
pub fn hem_release(
    components: &[String],
    upstream_z: &[f64],
    upstream_t: ThermodynamicTemperature,
    upstream_p: Pressure,
    diameter: Length,
    discharge_coefficient: f64,
    back_pressure: Pressure,
) -> Result<HemReleaseResult> {
    let spec = &model_gen::HEM_RELEASE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "upstream_t" => Some(upstream_t.value),
            "upstream_p" => Some(upstream_p.value),
            "diameter" => Some(diameter.value),
            "discharge_coefficient" => Some(discharge_coefficient),
            "back_pressure" => Some(back_pressure.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, ideal_gas) = databank::mixture_of(&names, Cubic::Pr, None)?;
    if upstream_z.len() != mixture.len() {
        return Err(AzothError::invalid_input(
            "upstream_z",
            format!(
                "a mixture of {} components needs that many mole fractions, and {upstream_z:?} \
                 has {}",
                mixture.len(),
                upstream_z.len()
            ),
        ));
    }
    check_composition(upstream_z, mixture.len(), "upstream_z")?;
    let molar_mass = entropy_state::molar_mass_of(&mixture, upstream_z)?;
    let mass_fractions = component_mass_fractions(&mixture, upstream_z)?;
    let mixture_names: Vec<String> = mixture
        .names()
        .map_or_else(|| components.to_vec(), <[String]>::to_vec);
    if let Err(failure) = check_input_support(upstream_t.value, &mixture_names, upstream_z) {
        return Ok(failure_result(failure, warnings));
    }

    let mut evaluations = 0_u32;
    let mut diagnostics = vec![
        HemReleaseDiagnostic {
            code: "MODEL_ASSUMPTIONS".to_string(),
            message: "Equilibrium, zero slip, isentropic acceleration; no pipe friction, heat \
                      transfer or solids"
                .to_string(),
        },
        HemReleaseDiagnostic {
            code: "EXIT_ALIASES_THROAT".to_string(),
            message: "Zero-length opening: exit equals accepted throat".to_string(),
        },
        HemReleaseDiagnostic {
            code: "ENTROPY_SOLVER".to_string(),
            message: if mixture.len() == 1 {
                "Pure-fluid PS flash retains saturation quality"
            } else {
                "Cold TP flashes in bracketed mixture entropy root"
            }
            .to_string(),
        },
    ];

    let upstream = match entropy_state::flashed_state(
        &mixture,
        &ideal_gas,
        upstream_z,
        molar_mass,
        upstream_p.value,
        upstream_t.value,
        &mut evaluations,
    ) {
        Ok(state) => state,
        Err(error) => {
            return Ok(failure("invalid", "CALCULATION_FAILED", error.0, warnings));
        }
    };
    if let Err(failure) = check_closures(&upstream)
        .and_then(|()| check_support(&upstream, &mixture_names, upstream_z))
    {
        return Ok(failure_result(failure, warnings));
    }
    let target = upstream.molar_entropy / molar_mass;
    let context = Context {
        mixture: &mixture,
        ideal_gas: &ideal_gas,
        z: upstream_z,
        names: &mixture_names,
        molar_mass,
        target,
        upstream_h: upstream.molar_enthalpy,
        upstream_pressure: upstream_p.value,
        start_temperature: upstream_t.value,
    };
    let upstream_sample = SampleState {
        state: upstream,
        velocity: 0.0,
        flux: 0.0,
    };
    let lower = back_pressure.value;
    let upper = upstream_p.value;

    if upper <= lower {
        diagnostics.push(HemReleaseDiagnostic {
            code: "NO_FORWARD_FLOW".to_string(),
            message: "Upstream pressure is at or below receiving pressure".to_string(),
        });
        let stations = vec![
            station_row(&upstream_sample, ReleaseStation::UpstreamStagnation),
            station_row(&upstream_sample, ReleaseStation::ThroatCritical),
            station_row(&upstream_sample, ReleaseStation::OrificeExit),
        ];
        return Ok(result(
            "valid",
            Some(0.0),
            false,
            None,
            upstream_z,
            mass_fractions,
            stations,
            diagnostics,
            warnings,
        ));
    }

    let mut samples = Vec::with_capacity(SAMPLES + 1);
    for index in 0..=SAMPLES {
        if index == SAMPLES {
            samples.push(upstream_sample.clone());
            continue;
        }
        let pressure = lower * ((upper / lower).ln() * index as f64 / SAMPLES as f64).exp();
        match context.sample(pressure, &mut evaluations) {
            Ok(sample) => samples.push(sample),
            Err(failure) => return Ok(failure_result(failure, warnings)),
        }
    }

    let mut best = samples[0].clone();
    for index in 1..SAMPLES {
        if samples[index].flux >= samples[index - 1].flux
            && samples[index].flux >= samples[index + 1].flux
        {
            match context.refine(
                samples[index - 1].state.pressure,
                samples[index + 1].state.pressure,
                &mut evaluations,
            ) {
                Ok(candidate) => {
                    if candidate.flux > best.flux {
                        best = candidate;
                    }
                }
                Err(failure) => return Ok(failure_result(failure, warnings)),
            }
        }
    }
    for sample in &samples {
        if sample.flux > best.flux {
            best = sample.clone();
        }
    }

    let choked = best.state.pressure > lower * (1.0 + CHOKED_MARGIN);
    let mut sound_speed = None;
    let mut warning = false;
    if choked {
        match equilibrium_speed_of(
            &mixture,
            &ideal_gas,
            kelvins(best.state.temperature),
            pascals(best.state.pressure),
            upstream_z,
            SOUND_SPEED_RELATIVE_STEP,
            Some(joules_per_kilogram_kelvin(target)),
        ) {
            Ok(acoustic) => {
                if acoustic.status == EquilibriumStatus::Converged {
                    if let Some(speed) = acoustic.sound_speed {
                        sound_speed = Some(speed);
                        let mach = best.velocity / speed;
                        if (mach - 1.0).abs() > MACH_BAND {
                            warning = true;
                            diagnostics.push(HemReleaseDiagnostic {
                                code: "THROAT_MACH_MISMATCH".to_string(),
                                message: format!(
                                    "Equilibrium Mach={mach}; phase-boundary maximum may be \
                                     nonsmooth"
                                ),
                            });
                        }
                    }
                }
                if sound_speed.is_none() {
                    warning = true;
                    diagnostics.push(HemReleaseDiagnostic {
                        code: "ACOUSTIC_UNAVAILABLE".to_string(),
                        message: format!("{}: {}", acoustic.status.as_str(), acoustic.message),
                    });
                }
            }
            Err(error) => {
                return Ok(failure(
                    "invalid",
                    "CALCULATION_FAILED",
                    format!("{error}"),
                    warnings,
                ));
            }
        }
    }

    let stations = vec![
        station_row(&upstream_sample, ReleaseStation::UpstreamStagnation),
        station_row(&best, ReleaseStation::ThroatCritical),
        station_row(&best, ReleaseStation::OrificeExit),
        station_row(&samples[0], ReleaseStation::AmbientExpanded),
    ];
    let rate =
        discharge_coefficient * std::f64::consts::PI * diameter.value.powi(2) / 4.0 * best.flux;
    let status = if warning {
        "valid_with_warnings"
    } else {
        "valid"
    };
    Ok(result(
        status,
        Some(rate),
        choked,
        sound_speed,
        upstream_z,
        mass_fractions,
        stations,
        diagnostics,
        warnings,
    ))
}

/// One station row from a sample.
fn station_row(sample: &SampleState, station: ReleaseStation) -> HemReleaseStation {
    HemReleaseStation {
        station,
        pressure: pascals(sample.state.pressure),
        temperature: kelvins(sample.state.temperature),
        density: kilograms_per_cubic_meter(sample.state.density),
        molar_enthalpy: joules_per_mole(sample.state.molar_enthalpy),
        molar_entropy: joules_per_mole_kelvin(sample.state.molar_entropy),
        velocity: meters_per_second(sample.velocity),
        mass_flux: kilograms_per_square_meter_second(sample.state.density * sample.velocity),
        phase: ReleasePhase::from_phase(sample.state.phase)
            .expect("the closures refused a trivial flash"),
        gas_mass_fraction: sample.state.gas_mass_fraction,
    }
}

/// The result of a successful calculation.
#[allow(clippy::too_many_arguments)] // one argument per result field
fn result(
    status: &str,
    flow: Option<f64>,
    choked: bool,
    sound_speed: Option<f64>,
    z: &[f64],
    mass_fractions: Vec<f64>,
    stations: Vec<HemReleaseStation>,
    diagnostics: Vec<HemReleaseDiagnostic>,
    warnings: Vec<Warning>,
) -> HemReleaseResult {
    HemReleaseResult {
        status: status.to_string(),
        mass_flow_rate: flow.map(kilograms_per_second),
        choked,
        throat_sound_speed: sound_speed.map(meters_per_second),
        component_mole_fractions: z.to_vec(),
        component_mass_fractions: mass_fractions,
        stations,
        diagnostic_codes: diagnostics.iter().map(|row| row.code.clone()).collect(),
        diagnostics,
        warnings,
    }
}

/// The class's failure result: no rate, no stations, one diagnostic.
fn failure(status: &str, code: &str, message: String, warnings: Vec<Warning>) -> HemReleaseResult {
    HemReleaseResult {
        status: status.to_string(),
        mass_flow_rate: None,
        choked: false,
        throat_sound_speed: None,
        component_mole_fractions: Vec::new(),
        component_mass_fractions: Vec::new(),
        stations: Vec::new(),
        diagnostics: vec![HemReleaseDiagnostic {
            code: code.to_string(),
            message,
        }],
        diagnostic_codes: vec![code.to_string()],
        warnings,
    }
}

/// The class's status for one failure.
fn failure_result(refusal: ReleaseFailure, warnings: Vec<Warning>) -> HemReleaseResult {
    match refusal {
        ReleaseFailure::Unsupported(message) => {
            failure("unsupported", "UNSUPPORTED_PHYSICS", message, warnings)
        }
        ReleaseFailure::Invalid(message) => {
            failure("invalid", "CALCULATION_FAILED", message, warnings)
        }
    }
}

/// The upstream overall mass fractions.
fn component_mass_fractions(mixture: &Mixture, z: &[f64]) -> Result<Vec<f64>> {
    let mut masses = Vec::with_capacity(z.len());
    let mut total = 0.0;
    for (index, (component, fraction)) in mixture.components().iter().zip(z).enumerate() {
        let mass = component.molar_mass.ok_or_else(|| {
            AzothError::invalid_input(
                "molar_mass",
                format!("component {index} carries no molar mass"),
            )
        })?;
        let mass = fraction * mass;
        total += mass;
        masses.push(mass);
    }
    if total <= 0.0 || !total.is_finite() {
        return Err(AzothError::invalid_input(
            "upstream_z",
            format!("the mixture's total mass is {total} kg/mol"),
        ));
    }
    Ok(masses.into_iter().map(|mass| mass / total).collect())
}
