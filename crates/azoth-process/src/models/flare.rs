//! `process.flare` - the flare's kernel as a registered id.
//!
//! Spec: `specs/models/process/flare.toml`. The arithmetic is [`crate::kernels::flare`]; what
//! is here is the boundary a case and a cross-impl test address.
//!
//! **The record passes through and the numbers are added.** `run` clones the inlet into the
//! outlet, so the five fields of the inlet port are the five fields of the outlet port - and
//! the class's own, the heat duty and the CO2 emission, are declared as outputs because they
//! are what the machine reports. **The radiation model, the dispersion surrogate and the
//! capacity check are the same kind of statement**: the setters are declared as optional inputs
//! and the numbers they feed as optional outputs, so a state without a distance or a design
//! limit reports no flux and no utilisation rather than a zero.

use azoth_core::units::{
    HeatFluxDensity, Length, MassRate, MolarEnergy, Power, Pressure, ThermodynamicTemperature,
    Velocity, VolumeRate, cubic_meters_per_second, joules_per_mole, meters, meters_per_second,
    pascals, watts_per_square_meter,
};
use azoth_core::{Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::flare as kernel;
use crate::kernels::flare::{
    DEFAULT_FLAME_HEIGHT, DEFAULT_RADIANT_FRACTION, DEFAULT_TIP_DIAMETER, FlareSettings,
    MIN_TIP_DIAMETER,
};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.flare`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FlareResult {
    /// Product molar flow, mol/s.
    pub product_n: f64,
    /// Product composition.
    pub product_z: Vec<f64>,
    /// Product pressure.
    #[serde(serialize_with = "scalar")]
    pub product_p: Pressure,
    /// Product temperature.
    #[serde(serialize_with = "scalar")]
    pub product_t: ThermodynamicTemperature,
    /// Product molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub product_h: MolarEnergy,
    /// The heat the flare releases, W.
    #[serde(serialize_with = "scalar")]
    pub heat_duty: Power,
    /// The carbon dioxide the combustion forms, kg/s.
    #[serde(serialize_with = "scalar")]
    pub co2_emission: MassRate,
    /// `estimateRadiationHeatFlux` at the stated distance; absent without one.
    #[serde(serialize_with = "optional_scalar")]
    pub radiation_flux: Option<HeatFluxDensity>,
    /// `radiationDistanceForFlux`; absent without a threshold.
    #[serde(serialize_with = "optional_scalar")]
    pub radiation_distance: Option<Length>,
    /// The dispersion surrogate's exit velocity.
    #[serde(serialize_with = "scalar")]
    pub exit_velocity: Velocity,
    /// The dispersion surrogate's momentum flux, `rho v^2`.
    #[serde(serialize_with = "scalar")]
    pub momentum_flux: Pressure,
    /// The surrogate's standard volume rate at the operating point.
    #[serde(serialize_with = "scalar")]
    pub standard_volume_rate: VolumeRate,
    /// The duty over the design heat release; absent without a positive design.
    pub heat_utilization: Option<f64>,
    /// The operating mass flow over the design mass flow; absent as above.
    pub mass_utilization: Option<f64>,
    /// The operating molar flow over the design molar flow; absent as above.
    pub molar_utilization: Option<f64>,
    /// `isOverloaded`; absent until a utilisation exists.
    pub overloaded: Option<bool>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl FlareResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's**: a flowsheet's are the *checker's* rather than
    /// [`apply_checks`]'.
    #[must_use]
    pub fn of(product: &Stream, numbers: &kernel::FlareNumbers, warnings: Vec<Warning>) -> Self {
        Self {
            product_n: product.n,
            product_z: product.z.clone(),
            product_p: product.p,
            product_t: product.t,
            product_h: joules_per_mole(product.h.value),
            heat_duty: numbers.heat_duty,
            co2_emission: numbers.co2_emission,
            radiation_flux: numbers.radiation_flux.map(watts_per_square_meter),
            radiation_distance: numbers.radiation_distance.map(meters),
            exit_velocity: meters_per_second(numbers.exit_velocity),
            momentum_flux: pascals(numbers.momentum_flux),
            standard_volume_rate: cubic_meters_per_second(numbers.standard_volumetric_flow),
            heat_utilization: numbers.heat_utilization,
            mass_utilization: numbers.mass_utilization,
            molar_utilization: numbers.molar_utilization,
            overloaded: numbers.overloaded,
            warnings,
        }
    }
}

/// A flare's steady state: the record through, and the numbers beside it.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if a component has no row in ISO 6976's table or
/// none in the element table, a stated distance or threshold is out of range, and whatever the
/// flashes refuse.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn flare(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    flame_height: Option<Length>,
    radiant_fraction: Option<f64>,
    tip_diameter: Option<Length>,
    ground_distance: Option<Length>,
    radiation_flux_threshold: Option<HeatFluxDensity>,
    design_heat_duty_capacity: Option<Power>,
    design_mass_flow_capacity: Option<MassRate>,
    design_molar_flow_capacity: Option<f64>,
) -> Result<FlareResult> {
    let spec = &model_gen::FLARE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_n" => Some(inlet_n),
            "ground_distance" => ground_distance.map(|value| value.value),
            "radiation_flux_threshold" => radiation_flux_threshold.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let inlet = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;
    // **The setters clamp rather than refuse**, and so does the boundary: a negative height is
    // zero, a fraction outside `[0, 1]` is the near end (a NaN keeps the default), and a tip
    // below `1e-4` m is `1e-4`.
    //
    // The class's two setters are `Math.max(0, value)` and `Math.min(Math.max(value, 0), 1)`,
    // and the port reproduces them rather than tightening into a refusal the class does not
    // make.
    // numerics-ok: `setRadiantFraction`'s own clamp, on the statement's first line.
    let settings = FlareSettings {
        flame_height: flame_height.map_or(DEFAULT_FLAME_HEIGHT, |value| value.value.max(0.0)),
        radiant_fraction: radiant_fraction.map_or(DEFAULT_RADIANT_FRACTION, |value| {
            if value.is_nan() {
                DEFAULT_RADIANT_FRACTION
            } else {
                // numerics-ok: `setRadiantFraction` is `Math.min(Math.max(value, 0), 1)`.
                value.clamp(0.0, 1.0)
            }
        }),
        tip_diameter: tip_diameter.map_or(DEFAULT_TIP_DIAMETER, |value| {
            value.value.max(MIN_TIP_DIAMETER)
        }),
        ground_distance: ground_distance.map(|value| value.value),
        radiation_flux_threshold: radiation_flux_threshold.map(|value| value.value),
        design_heat_duty_capacity: design_heat_duty_capacity.map(|value| value.value),
        design_mass_flow_capacity: design_mass_flow_capacity.map(|value| value.value),
        design_molar_flow_capacity,
    };
    let (product, numbers) = kernel::flare_with_settings(&inlet, &settings)?;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "heat_duty" => Some(numbers.heat_duty.value),
            "co2_emission" => Some(numbers.co2_emission.value),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(FlareResult::of(&product, &numbers, warnings))
}
