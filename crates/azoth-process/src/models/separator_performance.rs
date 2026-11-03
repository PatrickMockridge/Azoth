//! `process.separator_performance` - the droplet-based calculator as a registered id.
//!
//! Spec: `specs/models/process/separator_performance.toml`. The arithmetic is
//! [`crate::kernels::separator_performance`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **The class takes phase properties directly, and so does this boundary**: it does not flash.
//! The six fractions are the separator ids' own entrainment inputs, so a caller composes the
//! calculator with a vessel rather than enabling a flag on it.

use azoth_core::units::{
    DynamicViscosity, Length, MassDensity, Time, Velocity, meters, meters_per_second, seconds,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::separator_performance::{
    CurveSettings, DistributionKind, DropletDistribution, Internals, LookupPoint, Orientation,
    PerformanceSettings, separator_performance as kernel,
};
use crate::model_gen;
use crate::separator_design::inlet_device::InletDevice;

/// One row of a lookup table: a droplet diameter and its capture efficiency.
#[derive(Debug, Clone, PartialEq)]
pub struct EfficiencyPoint {
    /// The droplet diameter.
    pub diameter: Length,
    /// The capture efficiency, which the class clamps to `[0, 1]`.
    pub efficiency: f64,
}

/// Result of `process.separator_performance`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SeparatorPerformanceResult {
    /// The oil fraction carried out with the gas.
    pub oil_in_gas: f64,
    /// The water fraction carried out with the gas.
    pub water_in_gas: f64,
    /// The gas fraction carried under into the oil.
    pub gas_in_oil: f64,
    /// The gas fraction carried under into the water.
    pub gas_in_water: f64,
    /// The oil left in the water.
    pub oil_in_water: f64,
    /// The water left in the oil.
    pub water_in_oil: f64,
    /// The gravity stage's integrated efficiency.
    pub gravity_section_efficiency: f64,
    /// The mist eliminator's integrated efficiency.
    pub mist_eliminator_efficiency: f64,
    /// The two stages in series.
    pub overall_gas_liquid_efficiency: f64,
    /// The oil-in-water gravity efficiency.
    pub liquid_liquid_gravity_efficiency: f64,
    /// The gas section's critical diameter, m.
    #[serde(serialize_with = "scalar")]
    pub gravity_cut_diameter: Length,
    /// The predicted inlet regime's enum spelling; absent on the standard path.
    pub inlet_flow_regime: Option<String>,
    /// The inlet device's bulk efficiency.
    pub inlet_device_bulk_efficiency: Option<f64>,
    /// The inlet device's pressure drop, Pa.
    pub inlet_device_pressure_drop: Option<f64>,
    /// The inlet nozzle's velocity.
    #[serde(serialize_with = "optional_scalar")]
    pub inlet_nozzle_velocity: Option<Velocity>,
    /// The inlet nozzle's momentum flux, Pa.
    pub inlet_momentum_flux: Option<f64>,
    /// The actual Souders-Brown K.
    #[serde(serialize_with = "optional_scalar")]
    pub k_factor: Option<Velocity>,
    /// The K over its design value.
    pub k_factor_utilization: Option<f64>,
    /// Whether the mist eliminator is flooded.
    pub mist_eliminator_flooded: Option<bool>,
    /// The gas cross-section, m2.
    pub gas_area: Option<f64>,
    /// The liquid cross-section, m2.
    pub liquid_area: Option<f64>,
    /// The gas section's volume, m3.
    pub gas_volume: Option<f64>,
    /// The liquid section's volume, m3.
    pub liquid_volume: Option<f64>,
    /// The vessel's internal volume, m3.
    pub total_volume: Option<f64>,
    /// The gas section's settling height.
    #[serde(serialize_with = "optional_scalar")]
    pub effective_gas_settling_height: Option<Length>,
    /// The liquid section's settling height.
    #[serde(serialize_with = "optional_scalar")]
    pub effective_liquid_settling_height: Option<Length>,
    /// The gas's residence time.
    #[serde(serialize_with = "optional_scalar")]
    pub gas_residence_time: Option<Time>,
    /// The liquid's residence time.
    #[serde(serialize_with = "optional_scalar")]
    pub liquid_residence_time: Option<Time>,
    /// The gravity section's length.
    #[serde(serialize_with = "optional_scalar")]
    pub gravity_section_length: Option<Length>,
    /// The three-phase oil pad.
    #[serde(serialize_with = "optional_scalar")]
    pub oil_pad_thickness: Option<Length>,
    /// The three-phase water layer.
    #[serde(serialize_with = "optional_scalar")]
    pub water_layer_height: Option<Length>,
    /// API 12J's gas-section verdict.
    pub api_gas_liquid_compliant: Option<bool>,
    /// API 12J's liquid-section verdict.
    pub api_liquid_compliant: Option<bool>,
    /// Both verdicts.
    pub api_fully_compliant: Option<bool>,
    /// The class's gas-section sentence.
    pub api_gas_liquid_comment: Option<String>,
    /// The class's liquid-section sentence.
    pub api_liquid_comment: Option<String>,
    /// The cut diameter the API check read.
    #[serde(serialize_with = "optional_scalar")]
    pub api_gravity_cut_diameter: Option<Length>,
    /// The K over the API limit.
    pub api_k_factor_utilization: Option<f64>,
    /// The API limit on K.
    #[serde(serialize_with = "optional_scalar")]
    pub api_max_k_factor: Option<Velocity>,
    /// The API minimum liquid residence time.
    #[serde(serialize_with = "optional_scalar")]
    pub api_min_liquid_residence_time: Option<Time>,
    /// The downstream distribution's characteristic diameter.
    #[serde(serialize_with = "optional_scalar")]
    pub post_inlet_device_characteristic_diameter: Option<Length>,
    /// The downstream distribution's spread.
    pub post_inlet_device_spread: Option<f64>,
    /// The downstream distribution's volume-median diameter.
    #[serde(serialize_with = "optional_scalar")]
    pub post_inlet_d50: Option<Length>,
    /// The downstream distribution's Sauter mean diameter.
    #[serde(serialize_with = "optional_scalar")]
    pub post_inlet_d32: Option<Length>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// The class's standard separation calculation.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a value outside its declared range, an undeclared
/// orientation or distribution name, a distribution whose three fields are not stated
/// together, a curve missing what its type requires, or a lookup table on a type that does not
/// read one.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn separator_performance(
    gas_density: MassDensity,
    oil_density: MassDensity,
    water_density: MassDensity,
    gas_viscosity: DynamicViscosity,
    oil_viscosity: DynamicViscosity,
    water_viscosity: DynamicViscosity,
    gas_velocity: Velocity,
    vessel_diameter: Length,
    vessel_length: Length,
    orientation: &str,
    liquid_level_fraction: f64,
    mist_position: Option<f64>,
    weir_position: Option<f64>,
    gas_liquid_distribution: Option<&str>,
    gas_liquid_characteristic_diameter: Option<Length>,
    gas_liquid_spread: Option<f64>,
    gas_bubble_distribution: Option<&str>,
    gas_bubble_characteristic_diameter: Option<Length>,
    gas_bubble_spread: Option<f64>,
    oil_in_water_distribution: Option<&str>,
    oil_in_water_characteristic_diameter: Option<Length>,
    oil_in_water_spread: Option<f64>,
    water_in_oil_distribution: Option<&str>,
    water_in_oil_characteristic_diameter: Option<Length>,
    water_in_oil_spread: Option<f64>,
    mist_eliminator: Option<&str>,
    mist_eliminator_cut_diameter: Option<Length>,
    mist_eliminator_sharpness: Option<f64>,
    mist_eliminator_max_efficiency: Option<f64>,
    mist_eliminator_points: Option<&[EfficiencyPoint]>,
    oil_water_coalescer: Option<&str>,
    oil_water_coalescer_cut_diameter: Option<Length>,
    oil_water_coalescer_max_efficiency: Option<f64>,
    oil_water_coalescer_points: Option<&[EfficiencyPoint]>,
    include_gravity_section: Option<bool>,
    liquid_in_gas_calibration_factor: Option<f64>,
    gas_carry_under_calibration_factor: Option<f64>,
    liquid_liquid_calibration_factor: Option<f64>,
    enhanced: Option<bool>,
    inlet_pipe_diameter: Option<Length>,
    surface_tension: Option<f64>,
    oil_water_interfacial_tension: Option<f64>,
    oil_volume_fraction: Option<f64>,
    apply_turbulence_correction: Option<bool>,
    inlet_device: Option<&str>,
    override_bulk_efficiency: Option<f64>,
    override_dsd_multiplier: Option<f64>,
    design_k_factor: Option<Velocity>,
    liquid_liquid_residence_time: Option<Time>,
) -> Result<SeparatorPerformanceResult> {
    let spec = &model_gen::SEPARATOR_PERFORMANCE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "gas_density" => Some(gas_density.value),
            "oil_density" => Some(oil_density.value),
            "water_density" => Some(water_density.value),
            "gas_viscosity" => Some(gas_viscosity.value),
            "oil_viscosity" => Some(oil_viscosity.value),
            "water_viscosity" => Some(water_viscosity.value),
            "gas_velocity" => Some(gas_velocity.value),
            "vessel_diameter" => Some(vessel_diameter.value),
            "vessel_length" => Some(vessel_length.value),
            "liquid_level_fraction" => Some(liquid_level_fraction),
            "gas_liquid_characteristic_diameter" => {
                gas_liquid_characteristic_diameter.map(|value| value.value)
            }
            "gas_liquid_spread" => gas_liquid_spread,
            "gas_bubble_characteristic_diameter" => {
                gas_bubble_characteristic_diameter.map(|value| value.value)
            }
            "gas_bubble_spread" => gas_bubble_spread,
            "oil_in_water_characteristic_diameter" => {
                oil_in_water_characteristic_diameter.map(|value| value.value)
            }
            "oil_in_water_spread" => oil_in_water_spread,
            "water_in_oil_characteristic_diameter" => {
                water_in_oil_characteristic_diameter.map(|value| value.value)
            }
            "water_in_oil_spread" => water_in_oil_spread,
            "mist_eliminator_cut_diameter" => mist_eliminator_cut_diameter.map(|value| value.value),
            "mist_eliminator_sharpness" => mist_eliminator_sharpness,
            "mist_eliminator_max_efficiency" => mist_eliminator_max_efficiency,
            "oil_water_coalescer_cut_diameter" => {
                oil_water_coalescer_cut_diameter.map(|value| value.value)
            }
            "oil_water_coalescer_max_efficiency" => oil_water_coalescer_max_efficiency,
            "inlet_pipe_diameter" => inlet_pipe_diameter.map(|value| value.value),
            "surface_tension" => surface_tension,
            "oil_volume_fraction" => oil_volume_fraction,
            "design_k_factor" => design_k_factor.map(|value| value.value),
            "liquid_liquid_residence_time" => liquid_liquid_residence_time.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let settings = PerformanceSettings {
        gas_density: gas_density.value,
        oil_density: oil_density.value,
        water_density: water_density.value,
        gas_viscosity: gas_viscosity.value,
        oil_viscosity: oil_viscosity.value,
        water_viscosity: water_viscosity.value,
        gas_velocity: gas_velocity.value,
        vessel_diameter: vessel_diameter.value,
        vessel_length: vessel_length.value,
        orientation: orientation_of(orientation)?,
        liquid_level_fraction,
        // The class's `SeparatorGeometryCalculator` defaults.
        mist_position: mist_position.unwrap_or(0.90),
        weir_position: weir_position.unwrap_or(0.60),
        gas_liquid: optional_distribution(
            "gas_liquid",
            gas_liquid_distribution,
            gas_liquid_characteristic_diameter,
            gas_liquid_spread,
        )?,
        gas_bubble: optional_distribution(
            "gas_bubble",
            gas_bubble_distribution,
            gas_bubble_characteristic_diameter,
            gas_bubble_spread,
        )?,
        oil_in_water: optional_distribution(
            "oil_in_water",
            oil_in_water_distribution,
            oil_in_water_characteristic_diameter,
            oil_in_water_spread,
        )?,
        water_in_oil: optional_distribution(
            "water_in_oil",
            water_in_oil_distribution,
            water_in_oil_characteristic_diameter,
            water_in_oil_spread,
        )?,
        mist_eliminator: mist_curve(
            mist_eliminator,
            mist_eliminator_cut_diameter,
            mist_eliminator_sharpness,
            mist_eliminator_max_efficiency,
            mist_eliminator_points,
        )?,
        coalescer: coalescer_curve(
            oil_water_coalescer,
            oil_water_coalescer_cut_diameter,
            oil_water_coalescer_max_efficiency,
            oil_water_coalescer_points,
        )?,
        include_gravity_section: include_gravity_section.unwrap_or(true),
        // The class's setters are `Math.max(0.0, factor)`.
        liquid_in_gas_factor: liquid_in_gas_calibration_factor.map_or(1.0, |value| value.max(0.0)),
        gas_carry_under_factor: gas_carry_under_calibration_factor
            .map_or(1.0, |value| value.max(0.0)),
        liquid_liquid_factor: liquid_liquid_calibration_factor.map_or(1.0, |value| value.max(0.0)),
        enhanced: enhanced.unwrap_or(false),
        inlet_pipe_diameter: inlet_pipe_diameter.map_or(0.2, |value| value.value),
        surface_tension: surface_tension.unwrap_or(0.025),
        oil_water_interfacial_tension: oil_water_interfacial_tension.unwrap_or(0.030),
        oil_volume_fraction: oil_volume_fraction.unwrap_or(0.5),
        apply_turbulence_correction: apply_turbulence_correction.unwrap_or(true),
        inlet_device: inlet_device_of(inlet_device)?,
        // The class's setters treat zero as "no override".
        override_bulk_efficiency: override_bulk_efficiency.unwrap_or(0.0),
        override_dsd_multiplier: override_dsd_multiplier.unwrap_or(0.0),
        design_k_factor: design_k_factor.map(|value| value.value),
        liquid_liquid_residence_time: liquid_liquid_residence_time.map(|value| value.value),
    };
    let outcome = kernel(&settings)?;
    let geometry = outcome.geometry;
    let active = geometry.is_some();
    Ok(SeparatorPerformanceResult {
        oil_in_gas: outcome.oil_in_gas,
        water_in_gas: outcome.water_in_gas,
        gas_in_oil: outcome.gas_in_oil,
        gas_in_water: outcome.gas_in_water,
        oil_in_water: outcome.oil_in_water,
        water_in_oil: outcome.water_in_oil,
        gravity_section_efficiency: outcome.gravity_section_efficiency,
        mist_eliminator_efficiency: outcome.mist_eliminator_efficiency,
        overall_gas_liquid_efficiency: outcome.overall_gas_liquid_efficiency,
        liquid_liquid_gravity_efficiency: outcome.liquid_liquid_gravity_efficiency,
        gravity_cut_diameter: meters(outcome.gravity_cut_diameter),
        inlet_flow_regime: active.then(|| outcome.inlet_flow_regime.to_string()),
        inlet_device_bulk_efficiency: active.then_some(outcome.inlet_device_bulk_efficiency),
        inlet_device_pressure_drop: active.then_some(outcome.inlet_device_pressure_drop),
        inlet_nozzle_velocity: active.then(|| meters_per_second(outcome.inlet_nozzle_velocity)),
        inlet_momentum_flux: active.then_some(outcome.inlet_momentum_flux),
        k_factor: active.then(|| meters_per_second(outcome.k_factor)),
        k_factor_utilization: active.then_some(outcome.k_factor_utilization),
        mist_eliminator_flooded: active.then_some(outcome.mist_eliminator_flooded),
        gas_area: geometry.map(|value| value.gas_area),
        liquid_area: geometry.map(|value| value.liquid_area),
        gas_volume: geometry.map(|value| value.gas_volume),
        liquid_volume: geometry.map(|value| value.liquid_volume),
        total_volume: geometry.map(|value| value.total_volume),
        effective_gas_settling_height: geometry
            .map(|value| meters(value.effective_gas_settling_height)),
        effective_liquid_settling_height: geometry
            .map(|value| meters(value.effective_liquid_settling_height)),
        gas_residence_time: geometry.map(|value| seconds(value.gas_residence_time)),
        liquid_residence_time: geometry.map(|value| seconds(value.liquid_residence_time)),
        gravity_section_length: geometry.map(|value| meters(value.gravity_section_length)),
        oil_pad_thickness: geometry.map(|value| meters(value.oil_pad_thickness)),
        water_layer_height: geometry.map(|value| meters(value.water_layer_height)),
        api_gas_liquid_compliant: active.then_some(outcome.api_gas_liquid_compliant),
        api_liquid_compliant: active.then_some(outcome.api_liquid_compliant),
        api_fully_compliant: active.then_some(outcome.api_fully_compliant),
        api_gas_liquid_comment: active.then(|| outcome.api_gas_liquid_comment.clone()),
        api_liquid_comment: active.then(|| outcome.api_liquid_comment.clone()),
        api_gravity_cut_diameter: active
            .then(|| meters(outcome.api_gravity_cut_diameter_um / 1.0e6)),
        api_k_factor_utilization: active.then_some(outcome.api_k_factor_utilization),
        api_max_k_factor: active.then(|| meters_per_second(outcome.api_max_k_factor)),
        api_min_liquid_residence_time: active
            .then(|| seconds(outcome.api_min_liquid_residence_time)),
        post_inlet_device_characteristic_diameter: active
            .then(|| meters(outcome.post_inlet_device_characteristic_diameter)),
        post_inlet_device_spread: active.then_some(outcome.post_inlet_device_spread),
        post_inlet_d50: active.then(|| meters(outcome.post_inlet_d50)),
        post_inlet_d32: active.then(|| meters(outcome.post_inlet_d32)),
        warnings,
    })
}

/// `orientation`, refused by name rather than silently horizontal.
fn orientation_of(value: &str) -> Result<Orientation> {
    match value {
        "horizontal" => Ok(Orientation::Horizontal),
        "vertical" => Ok(Orientation::Vertical),
        other => Err(AzothError::invalid_input(
            "orientation",
            format!("'{other}' is not one of the two declared orientations"),
        )),
    }
}

/// The inlet device, by its spec spelling; absent is the class's `HALF_PIPE`.
fn inlet_device_of(value: Option<&str>) -> Result<InletDevice> {
    match value.unwrap_or("half_pipe") {
        "none" => Ok(InletDevice::None),
        "deflector_plate" => Ok(InletDevice::DeflectorPlate),
        "half_pipe" => Ok(InletDevice::HalfPipe),
        "inlet_vane" => Ok(InletDevice::InletVane),
        "inlet_cyclone" => Ok(InletDevice::InletCyclone),
        "schoepentoeter" => Ok(InletDevice::Schoepentoeter),
        "impingement_plate" => Ok(InletDevice::ImpingementPlate),
        other => Err(AzothError::invalid_input(
            "inlet_device",
            format!("'{other}' is not one of the seven declared devices"),
        )),
    }
}

/// `DistributionType`, refused by name.
fn distribution_kind(field: &str, value: &str) -> Result<DistributionKind> {
    match value {
        "rosin_rammler" => Ok(DistributionKind::RosinRammler),
        "log_normal" => Ok(DistributionKind::LogNormal),
        other => Err(AzothError::invalid_input(
            field,
            format!("'{other}' is not one of the two declared distributions"),
        )),
    }
}

/// One optional distribution, which is stated whole or not at all.
fn optional_distribution(
    role: &str,
    distribution: Option<&str>,
    characteristic_diameter: Option<Length>,
    spread: Option<f64>,
) -> Result<Option<DropletDistribution>> {
    match (distribution, characteristic_diameter, spread) {
        (None, None, None) => Ok(None),
        (Some(distribution), Some(characteristic_diameter), Some(spread)) => {
            Ok(Some(DropletDistribution {
                kind: distribution_kind(&format!("{role}_distribution"), distribution)?,
                characteristic_diameter: characteristic_diameter.value,
                spread,
            }))
        }
        _ => Err(AzothError::invalid_input(
            format!("{role}_distribution"),
            "a distribution needs its type, characteristic diameter and spread together",
        )),
    }
}

/// The mist eliminator's curve, with the class's per-type defaults resolved in the kernel.
fn mist_curve(
    internals: Option<&str>,
    cut_diameter: Option<Length>,
    sharpness: Option<f64>,
    max_efficiency: Option<f64>,
    points: Option<&[EfficiencyPoint]>,
) -> Result<Option<CurveSettings>> {
    let Some(internals) = internals else {
        return absent_curve(
            "mist_eliminator",
            cut_diameter,
            sharpness,
            max_efficiency,
            points,
        );
    };
    if internals == "none" {
        return absent_curve(
            "mist_eliminator",
            cut_diameter,
            sharpness,
            max_efficiency,
            points,
        );
    }
    let internals = match internals {
        "wire_mesh" => Internals::WireMesh,
        "vane_pack" => Internals::VanePack,
        "axial_cyclone" => Internals::AxialCyclone,
        "plate_pack" => Internals::PlatePack,
        "custom" => Internals::Custom,
        other => {
            return Err(AzothError::invalid_input(
                "mist_eliminator",
                format!("'{other}' is not one of the five declared internals"),
            ));
        }
    };
    Ok(Some(CurveSettings {
        internals,
        cut_diameter: cut_diameter.map(|value| value.value),
        sharpness,
        max_efficiency,
        points: lookup_points(points),
    }))
}

/// The liquid-liquid coalescer's curve; the class's `platePack` or a custom table.
fn coalescer_curve(
    internals: Option<&str>,
    cut_diameter: Option<Length>,
    max_efficiency: Option<f64>,
    points: Option<&[EfficiencyPoint]>,
) -> Result<Option<CurveSettings>> {
    let Some(internals) = internals else {
        return absent_curve(
            "oil_water_coalescer",
            cut_diameter,
            None,
            max_efficiency,
            points,
        );
    };
    if internals == "none" {
        return absent_curve(
            "oil_water_coalescer",
            cut_diameter,
            None,
            max_efficiency,
            points,
        );
    }
    let internals = match internals {
        "plate_pack" => Internals::PlatePack,
        "custom" => Internals::Custom,
        other => {
            return Err(AzothError::invalid_input(
                "oil_water_coalescer",
                format!("'{other}' is not one of the two declared coalescers"),
            ));
        }
    };
    Ok(Some(CurveSettings {
        internals,
        cut_diameter: cut_diameter.map(|value| value.value),
        sharpness: None,
        max_efficiency,
        points: lookup_points(points),
    }))
}

/// A curve that is not configured: its parameters must be absent too.
fn absent_curve(
    role: &str,
    cut_diameter: Option<Length>,
    sharpness: Option<f64>,
    max_efficiency: Option<f64>,
    points: Option<&[EfficiencyPoint]>,
) -> Result<Option<CurveSettings>> {
    if cut_diameter.is_some() || sharpness.is_some() || max_efficiency.is_some() || points.is_some()
    {
        return Err(AzothError::invalid_input(
            role,
            "curve parameters are stated without an internals type",
        ));
    }
    Ok(None)
}

/// The row list as the kernel's SI points.
fn lookup_points(points: Option<&[EfficiencyPoint]>) -> Option<Vec<LookupPoint>> {
    points.map(|rows| {
        rows.iter()
            .map(|row| LookupPoint {
                diameter: row.diameter.value,
                efficiency: row.efficiency,
            })
            .collect()
    })
}
