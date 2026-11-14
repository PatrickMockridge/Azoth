//! `process.consequence_analysis` - the QRA composer as a registered id.
//!
//! Spec: `specs/models/process/consequence_analysis.toml`. The arithmetic is
//! `process/safety/qra/ConsequenceAnalysisEngine.java` at the pin over the ported
//! consequence kernels: the four dose formulas, the zero-dose answer, the per-outcome
//! frequency and the summed individual risk.

use azoth_core::units::{
    Length, Mass, MassFlux, MassRate, MolarMass, Pressure, ThermodynamicTemperature, Time,
    Velocity, meters,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::warnings as wire_warnings;
use crate::model_gen;
use crate::models::{gaussian_plume, jet_fire, pool_fire, probit, vce};

/// The gas constant the class's toxic conversion uses, J/(mol K).
const TOXIC_GAS_CONSTANT: f64 = 8.314;

/// One outcome branch, as the caller states it.
#[derive(Debug, Clone, PartialEq)]
pub struct ConsequenceOutcomeSpec {
    /// The branch type, which selects the fields the row must carry.
    pub outcome_type: String,
    /// The branch's probability conditional on the release.
    pub conditional_probability: f64,
    /// The probit intercept.
    pub probit_a: f64,
    /// The probit slope.
    pub probit_b: f64,
    /// The dose exponent.
    pub probit_n: f64,
    /// The receiver's exposure time.
    pub exposure_seconds: Option<Time>,
    /// The jet's fuel mass flow.
    pub mass_flow: Option<MassRate>,
    /// The lower heating value in J/kg.
    pub heat_of_combustion: Option<f64>,
    /// The jet's radiated fraction.
    pub radiative_fraction: Option<f64>,
    /// The jet's ambient humidity.
    pub relative_humidity: Option<f64>,
    /// The pool's equivalent circular diameter.
    pub pool_diameter: Option<Length>,
    /// The pool's mass-burning rate per unit area.
    pub burn_rate: Option<MassFlux>,
    /// The pool's ambient air density in kg/m3.
    pub air_density: Option<f64>,
    /// The VCE's flammable fuel mass.
    pub cloud_fuel_mass: Option<Mass>,
    /// The TNO strength class.
    pub strength_class: Option<f64>,
    /// The VCE's ambient pressure.
    pub ambient_pressure: Option<Pressure>,
    /// The toxic plume's source strength.
    pub emission_rate: Option<MassRate>,
    /// The toxic plume's effective release height.
    pub effective_stack_height: Option<Length>,
    /// The toxic plume's wind speed.
    pub wind_speed: Option<Velocity>,
    /// The toxic plume's Pasquill class.
    pub stability: Option<String>,
    /// The toxic plume's terrain category.
    pub terrain: Option<String>,
    /// The toxic component's molar mass.
    pub molar_mass: Option<MolarMass>,
    /// The toxic plume's temperature.
    pub temperature: Option<ThermodynamicTemperature>,
    /// The toxic plume's pressure.
    pub pressure: Option<Pressure>,
}

/// One outcome row of the result.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ConsequenceOutcomeResult {
    /// The class's fixed label.
    pub outcome: String,
    /// The branch type.
    pub outcome_type: String,
    /// The frequency, in 1/year.
    pub frequency_per_year: f64,
    /// The probit of the dose, or zero when the dose is not positive.
    pub fatality_probability: f64,
    /// The frequency times the fatality probability, in 1/year.
    pub fatality_frequency_per_year: f64,
}

/// Result of `process.consequence_analysis`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ConsequenceAnalysisResult {
    /// One row per input row, in caller order.
    pub outcomes: Vec<ConsequenceOutcomeResult>,
    /// The sum of the rows' fatality frequencies, in 1/year.
    pub individual_fatality_risk_per_year: f64,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// The class's fixed label for a branch type.
fn outcome_label(kind: &str) -> &'static str {
    match kind {
        "jet_fire" => "Jet fire",
        "pool_fire" => "Pool fire",
        "vce" => "VCE",
        "toxic_dispersion" => "Toxic dispersion",
        _ => "",
    }
}

/// The per-type fields a branch of this kind may carry.
fn allowed_fields(kind: &str) -> &'static [&'static str] {
    match kind {
        "jet_fire" => &[
            "mass_flow",
            "heat_of_combustion",
            "radiative_fraction",
            "relative_humidity",
            "exposure_seconds",
        ],
        "pool_fire" => &[
            "pool_diameter",
            "burn_rate",
            "air_density",
            "exposure_seconds",
        ],
        "vce" => &[
            "cloud_fuel_mass",
            "heat_of_combustion",
            "strength_class",
            "ambient_pressure",
        ],
        "toxic_dispersion" => &[
            "emission_rate",
            "effective_stack_height",
            "wind_speed",
            "stability",
            "terrain",
            "molar_mass",
            "temperature",
            "pressure",
            "exposure_seconds",
        ],
        _ => &[],
    }
}

/// A field the row's type requires.
fn require<T>(index: usize, kind: &str, field: &'static str, value: Option<T>) -> Result<T> {
    value.ok_or_else(|| {
        AzothError::invalid_input(
            format!("outcomes[{index}].{field}"),
            format!("a {kind} branch requires `{field}`"),
        )
    })
}

/// Refuse a non-finite scalar.
fn require_finite(name: String, value: f64) -> Result<()> {
    if value.is_finite() {
        Ok(())
    } else {
        Err(AzothError::invalid_input(
            name,
            "a non-finite value is not a number this calculation publishes",
        ))
    }
}

/// The fatality probability of a dose; zero where the class's probit answers `-Infinity`.
fn fatality(a: f64, b: f64, n: f64, dose: f64) -> Result<f64> {
    if dose <= 0.0 {
        return Ok(0.0);
    }
    Ok(probit::probit(a, b, n, Some(dose), None)?.probability)
}

/// The class's four consequence branches over one receiver, and the summed risk.
///
/// # Errors
/// `InvalidInput` for a probability outside `[0, 1]`, a missing or non-positive field its
/// type requires, a field belonging to another type, a non-positive receiver distance, a
/// negative release frequency, or a non-finite input.
pub fn consequence_analysis(
    release_frequency_per_year: f64,
    receiver_distance: Length,
    outcomes: &[ConsequenceOutcomeSpec],
) -> Result<ConsequenceAnalysisResult> {
    let spec = &model_gen::CONSEQUENCE_ANALYSIS_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "release_frequency_per_year" => Some(release_frequency_per_year),
            "receiver_distance" => Some(receiver_distance.value),
            _ => None,
        },
        &mut warnings,
    )?;
    require_finite(
        "release_frequency_per_year".to_string(),
        release_frequency_per_year,
    )?;
    require_finite("receiver_distance".to_string(), receiver_distance.value)?;

    let distance = receiver_distance.value;
    let mut rows = Vec::with_capacity(outcomes.len());
    for (index, row) in outcomes.iter().enumerate() {
        let kind = row.outcome_type.as_str();
        if !matches!(kind, "jet_fire" | "pool_fire" | "vce" | "toxic_dispersion") {
            return Err(AzothError::invalid_input(
                format!("outcomes[{index}].outcome_type"),
                format!("an outcome type is one of the class's four branches, and {kind:?} is not"),
            ));
        }
        require_finite(
            format!("outcomes[{index}].conditional_probability"),
            row.conditional_probability,
        )?;
        if !(0.0..=1.0).contains(&row.conditional_probability) {
            return Err(AzothError::invalid_input(
                format!("outcomes[{index}].conditional_probability"),
                format!(
                    "a branch probability is in [0, 1], and {} is not",
                    row.conditional_probability
                ),
            ));
        }
        for (field, value) in [
            ("probit_a", row.probit_a),
            ("probit_b", row.probit_b),
            ("probit_n", row.probit_n),
        ] {
            require_finite(format!("outcomes[{index}].{field}"), value)?;
        }

        let allowed = allowed_fields(kind);
        let present: [(&'static str, bool); 20] = [
            ("mass_flow", row.mass_flow.is_some()),
            ("heat_of_combustion", row.heat_of_combustion.is_some()),
            ("radiative_fraction", row.radiative_fraction.is_some()),
            ("relative_humidity", row.relative_humidity.is_some()),
            ("pool_diameter", row.pool_diameter.is_some()),
            ("burn_rate", row.burn_rate.is_some()),
            ("air_density", row.air_density.is_some()),
            ("cloud_fuel_mass", row.cloud_fuel_mass.is_some()),
            ("strength_class", row.strength_class.is_some()),
            ("ambient_pressure", row.ambient_pressure.is_some()),
            ("emission_rate", row.emission_rate.is_some()),
            (
                "effective_stack_height",
                row.effective_stack_height.is_some(),
            ),
            ("wind_speed", row.wind_speed.is_some()),
            ("stability", row.stability.is_some()),
            ("terrain", row.terrain.is_some()),
            ("molar_mass", row.molar_mass.is_some()),
            ("temperature", row.temperature.is_some()),
            ("pressure", row.pressure.is_some()),
            ("exposure_seconds", row.exposure_seconds.is_some()),
            ("outcome_type", false),
        ];
        for (field, is_present) in present {
            if is_present && !allowed.contains(&field) {
                return Err(AzothError::invalid_input(
                    format!("outcomes[{index}].{field}"),
                    format!("a {kind} branch does not carry `{field}`"),
                ));
            }
        }

        let dose = match kind {
            "jet_fire" => {
                let mass = require(index, kind, "mass_flow", row.mass_flow)?;
                let heat = require(index, kind, "heat_of_combustion", row.heat_of_combustion)?;
                let fraction = require(index, kind, "radiative_fraction", row.radiative_fraction)?;
                let exposure = require(index, kind, "exposure_seconds", row.exposure_seconds)?;
                let result = jet_fire::jet_fire(
                    mass,
                    heat,
                    fraction,
                    row.relative_humidity,
                    Some(meters(distance)),
                    None,
                )?;
                let flux = result.incident_heat_flux.map_or(0.0, |value| value.value);
                exposure.value * flux.powf(4.0 / 3.0) / 1.0e4
            }
            "pool_fire" => {
                let diameter = require(index, kind, "pool_diameter", row.pool_diameter)?;
                let burning = require(index, kind, "burn_rate", row.burn_rate)?;
                let exposure = require(index, kind, "exposure_seconds", row.exposure_seconds)?;
                let result = pool_fire::pool_fire(
                    diameter,
                    burning,
                    row.air_density,
                    Some(meters(distance)),
                    None,
                )?;
                let flux = result.incident_heat_flux.map_or(0.0, |value| value.value);
                exposure.value * flux.powf(4.0 / 3.0) / 1.0e4
            }
            "vce" => {
                let mass = require(index, kind, "cloud_fuel_mass", row.cloud_fuel_mass)?;
                let heat = require(index, kind, "heat_of_combustion", row.heat_of_combustion)?;
                let strength = require(index, kind, "strength_class", row.strength_class)?;
                let result = vce::vce(
                    mass,
                    heat,
                    strength,
                    row.ambient_pressure,
                    Some(meters(distance)),
                    None,
                )?;
                result.overpressure.map_or(0.0, |value| value.value)
            }
            "toxic_dispersion" => {
                let emission = require(index, kind, "emission_rate", row.emission_rate)?;
                let stack = require(
                    index,
                    kind,
                    "effective_stack_height",
                    row.effective_stack_height,
                )?;
                let wind = require(index, kind, "wind_speed", row.wind_speed)?;
                let stability = require(index, kind, "stability", row.stability.as_deref())?;
                let terrain = require(index, kind, "terrain", row.terrain.as_deref())?;
                let molar = require(index, kind, "molar_mass", row.molar_mass)?;
                let temperature = require(index, kind, "temperature", row.temperature)?;
                let pressure = require(index, kind, "pressure", row.pressure)?;
                let exposure = require(index, kind, "exposure_seconds", row.exposure_seconds)?;
                require_finite(format!("outcomes[{index}].emission_rate"), emission.value)?;
                require_finite(format!("outcomes[{index}].wind_speed"), wind.value)?;
                require_finite(format!("outcomes[{index}].molar_mass"), molar.value)?;
                require_finite(format!("outcomes[{index}].temperature"), temperature.value)?;
                require_finite(format!("outcomes[{index}].pressure"), pressure.value)?;
                if molar.value <= 0.0 {
                    return Err(AzothError::invalid_input(
                        format!("outcomes[{index}].molar_mass"),
                        "the ppm conversion divides by the molar mass",
                    ));
                }
                if temperature.value <= 0.0 || pressure.value <= 0.0 {
                    return Err(AzothError::invalid_input(
                        format!("outcomes[{index}].pressure"),
                        "the ppm conversion divides by the pressure and needs a positive temperature",
                    ));
                }
                let concentration = gaussian_plume::concentration(
                    emission.value,
                    stack.value,
                    wind.value,
                    gaussian_plume::Stability::parse(stability)?,
                    gaussian_plume::Terrain::parse(terrain)?,
                    distance,
                    0.0,
                    0.0,
                );
                let ppm = (concentration / molar.value) * TOXIC_GAS_CONSTANT * temperature.value
                    / pressure.value
                    * 1.0e6;
                ppm * exposure.value / 60.0
            }
            _ => {
                return Err(AzothError::invalid_input(
                    format!("outcomes[{index}].outcome_type"),
                    format!(
                        "an outcome type is one of the class's four branches, and {kind:?} is not"
                    ),
                ));
            }
        };

        let frequency = release_frequency_per_year * row.conditional_probability;
        let probability = fatality(row.probit_a, row.probit_b, row.probit_n, dose)?;
        rows.push(ConsequenceOutcomeResult {
            outcome: outcome_label(kind).to_string(),
            outcome_type: kind.to_string(),
            frequency_per_year: frequency,
            fatality_probability: probability,
            fatality_frequency_per_year: frequency * probability,
        });
    }

    let risk = rows.iter().map(|row| row.fatality_frequency_per_year).sum();
    Ok(ConsequenceAnalysisResult {
        outcomes: rows,
        individual_fatality_risk_per_year: risk,
        warnings,
    })
}
