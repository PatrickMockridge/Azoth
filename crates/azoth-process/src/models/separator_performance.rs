//! `process.separator_performance` - the droplet-based calculator as a registered id.
//!
//! Spec: `specs/models/process/separator_performance.toml`. The arithmetic is
//! [`crate::kernels::separator_performance`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **The class takes phase properties directly, and so does this boundary**: it does not flash.
//! The six fractions are the separator ids' own entrainment inputs, so a caller composes the
//! calculator with a vessel rather than enabling a flag on it.

use azoth_core::units::{DynamicViscosity, Length, MassDensity, Velocity, meters};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::separator_performance::{
    CurveSettings, DistributionKind, DropletDistribution, Internals, LookupPoint, Orientation,
    PerformanceSettings, separator_performance as kernel,
};
use crate::model_gen;

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
    };
    let outcome = kernel(&settings)?;
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
