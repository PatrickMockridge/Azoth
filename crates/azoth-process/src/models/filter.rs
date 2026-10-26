//! `process.filter` - the filter's kernel as a registered id.
//!
//! Spec: `specs/models/process/filter.toml`. The arithmetic is
//! [`crate::kernels::filter`]; what is here is the boundary a case and a cross-impl test
//! address.

use azoth_core::units::{
    Area, DynamicViscosity, Length, MassRate, MolarEnergy, Pressure, ThermodynamicTemperature,
    VolumeRate, cubic_meter_per_hour, joules_per_mole,
};
use azoth_core::warning::{Warning, WarningCode};
use azoth_core::{AzothError, Result, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::filter::{
    self as kernel, BetaPoint, DropPoint, FilterSettings, Media, PressureDropModel,
};
use crate::model_gen;
use crate::stream::Stream;

/// One point of `pressure_drop_curve`.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterDropPoint {
    /// The actual volumetric flow.
    pub flow: VolumeRate,
    /// The clean drop there.
    pub drop: Pressure,
}

/// One point of `performance_curve`: a particle size and its beta ratio.
#[derive(Debug, Clone, PartialEq)]
pub struct FilterBetaPoint {
    /// The particle size.
    pub particle_size: Length,
    /// The beta ratio at that size.
    pub beta_ratio: f64,
}

/// Result of `process.filter`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct FilterResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure: `inlet_p - applied_drop`.
    #[serde(serialize_with = "scalar")]
    pub outlet_p: Pressure,
    /// Outlet temperature, which is the inlet's.
    #[serde(serialize_with = "scalar")]
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy, which is not the inlet's - the drop is isothermal.
    #[serde(serialize_with = "scalar")]
    pub outlet_h: MolarEnergy,
    /// The clean drop the selected model computed, before the clamp.
    #[serde(serialize_with = "scalar")]
    pub clean_drop: Pressure,
    /// The drop actually applied, which is the clean one unless it was clamped.
    #[serde(serialize_with = "scalar")]
    pub applied_drop: Pressure,
    /// The nominal removal efficiency, a fraction.
    pub nominal_removal_efficiency: f64,
    /// The outlet particle concentration as a mass fraction, zero when the model is off.
    pub outlet_particle_concentration: f64,
    /// The captured mass rate.
    #[serde(serialize_with = "scalar")]
    pub captured_rate: MassRate,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl FilterResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's**: a flowsheet's are the *checker's* rather than
    /// [`apply_checks`]'.
    #[must_use]
    pub fn of(outcome: &kernel::FilterOutcome, warnings: Vec<Warning>) -> Self {
        Self {
            outlet_n: outcome.outlet.n,
            outlet_z: outcome.outlet.z.clone(),
            outlet_p: outcome.outlet.p,
            outlet_t: outcome.outlet.t,
            outlet_h: joules_per_mole(outcome.outlet.h.value),
            clean_drop: outcome.clean_drop,
            applied_drop: outcome.applied_drop,
            nominal_removal_efficiency: outcome.nominal_removal_efficiency,
            outlet_particle_concentration: outcome.outlet_particle_concentration,
            captured_rate: outcome.captured_rate,
            warnings,
        }
    }
}

/// Drop a stream's pressure by the selected clean-drop model at a constant temperature.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, the selected model's own
/// settings are not a model, or the drop is not a number; and whatever the databank or the
/// flash refuses. A drop past the inlet pressure is **not** an error: the class clamps it and
/// so does this, saying so in a warning.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn filter(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    pressure_drop: Pressure,
    pressure_drop_model: Option<&str>,
    reference_flow: Option<VolumeRate>,
    flow_exponent: Option<f64>,
    pressure_drop_curve: Option<&[FilterDropPoint]>,
    media_area: Option<Area>,
    media_bed_depth: Option<Length>,
    media_particle_diameter: Option<Length>,
    media_void_fraction: Option<f64>,
    viscosity: Option<DynamicViscosity>,
    performance_curve: Option<&[FilterBetaPoint]>,
    particle_size: Option<Length>,
    nominal_removal_efficiency: Option<f64>,
    inlet_particle_concentration: Option<f64>,
) -> Result<FilterResult> {
    let spec = &model_gen::FILTER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "pressure_drop" => Some(pressure_drop.value),
            "flow_exponent" => flow_exponent,
            "media_void_fraction" => media_void_fraction,
            "nominal_removal_efficiency" => nominal_removal_efficiency,
            _ => None,
        },
        &mut warnings,
    )?;

    let model = match pressure_drop_model.unwrap_or("fixed") {
        "fixed" => PressureDropModel::Fixed,
        "flow_scaled" => PressureDropModel::FlowScaled,
        "tabulated" => PressureDropModel::Tabulated,
        "ergun" => PressureDropModel::Ergun,
        other => {
            return Err(AzothError::invalid_input(
                "pressure_drop_model",
                format!("{other} is not one of fixed, flow_scaled, tabulated or ergun"),
            ));
        }
    };
    let exponent = flow_exponent.unwrap_or(1.0);
    if exponent.is_nan() || exponent <= 0.0 {
        return Err(AzothError::invalid_input(
            "flow_exponent",
            format!("the scaling exponent is {exponent}, and a non-positive one is not a curve"),
        ));
    }
    if let Some(reference) = reference_flow {
        if reference.value.is_nan() || reference.value <= 0.0 {
            return Err(AzothError::invalid_input(
                "reference_flow",
                format!(
                    "the reference flow is {} m3/s, and the scaling divides by it",
                    reference.value
                ),
            ));
        }
    }
    let drop_curve: Vec<DropPoint> = pressure_drop_curve
        .unwrap_or(&[])
        .iter()
        .map(|point| DropPoint {
            flow_m3h: point.flow.get::<cubic_meter_per_hour>(),
            drop_pa: point.drop.value,
        })
        .collect();
    validate_drop_curve(&drop_curve)?;

    let media = Media {
        area_m2: media_area.map_or(Media::default().area_m2, |value| value.value),
        bed_depth_m: media_bed_depth.map_or(Media::default().bed_depth_m, |value| value.value),
        particle_diameter_m: media_particle_diameter
            .map_or(Media::default().particle_diameter_m, |value| value.value),
        void_fraction: media_void_fraction.unwrap_or(Media::default().void_fraction),
    };
    if model == PressureDropModel::Ergun {
        if viscosity.is_none() {
            return Err(AzothError::invalid_input(
                "viscosity",
                "the Ergun model reads the fluid's viscosity, and this filter was not given one",
            ));
        }
        if media.area_m2.is_nan()
            || media.area_m2 <= 0.0
            || media.bed_depth_m.is_nan()
            || media.bed_depth_m <= 0.0
            || media.particle_diameter_m.is_nan()
            || media.particle_diameter_m <= 0.0
        {
            return Err(AzothError::invalid_input(
                "media_area",
                "an Ergun bed needs a positive area, depth and particle diameter",
            ));
        }
        if media.void_fraction.is_nan() || media.void_fraction <= 0.0 || media.void_fraction >= 1.0
        {
            return Err(AzothError::invalid_input(
                "media_void_fraction",
                format!(
                    "the void fraction is {}, and an Ergun bed's is in (0, 1)",
                    media.void_fraction
                ),
            ));
        }
    }

    let beta_curve: Vec<BetaPoint> = performance_curve
        .unwrap_or(&[])
        .iter()
        .map(|point| BetaPoint {
            particle_size_m: point.particle_size.value,
            beta_ratio: point.beta_ratio,
        })
        .collect();
    validate_beta_curve(&beta_curve)?;
    let efficiency = nominal_removal_efficiency.unwrap_or(0.0);
    if !(0.0..=1.0).contains(&efficiency) {
        return Err(AzothError::invalid_input(
            "nominal_removal_efficiency",
            format!("the removal efficiency is {efficiency}, and a fraction is in [0, 1]"),
        ));
    }
    if let Some(concentration) = inlet_particle_concentration {
        if !concentration.is_finite() || concentration < 0.0 {
            return Err(AzothError::invalid_input(
                "inlet_particle_concentration",
                format!("the inlet concentration is {concentration}, not a non-negative fraction"),
            ));
        }
    }

    let settings = FilterSettings {
        model,
        reference_flow_m3h: reference_flow.map(|value| value.get::<cubic_meter_per_hour>()),
        flow_exponent: exponent,
        drop_curve: &drop_curve,
        media,
        viscosity_pa_s: viscosity.map(|value| value.value),
        beta_curve: &beta_curve,
        particle_size_m: particle_size.map_or(1.0e-5, |value| value.value),
        nominal_efficiency: efficiency,
    };

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;
    let outcome = kernel::filter_with(
        &feed,
        pressure_drop,
        &settings,
        inlet_particle_concentration,
    )?;

    // **The clamp is reported rather than hidden.** `Filter.run` logs a warning when the drop
    // it applies is below the one it was given, and the whole reason `applied_drop` is a
    // declared output is that the two can differ - so a result that carried the clamped value
    // and said nothing would be the silent-improvement failure this library is written
    // against. A negative `pressure_drop` is clamped by the same path and warned about by the
    // range check above.
    if (outcome.applied_drop.value - outcome.clean_drop.value).abs() > 0.0 {
        warnings.push(Warning::for_field(
            WarningCode::OutOfValidRange,
            "pressure_drop",
            format!(
                "a clean drop of {} Pa was computed and {} Pa was applied: the outlet is held a \
                 millionth of a bar above vacuum, which is the clamp `Filter.run` applies",
                outcome.clean_drop.value, outcome.applied_drop.value
            ),
        ));
    }

    Ok(FilterResult::of(&outcome, warnings))
}

fn validate_drop_curve(curve: &[DropPoint]) -> Result<()> {
    for (index, point) in curve.iter().enumerate() {
        if !point.flow_m3h.is_finite() || point.flow_m3h <= 0.0 {
            return Err(AzothError::invalid_input(
                "pressure_drop_curve",
                format!(
                    "point {index}'s flow is {} m3/hr, and the curve is over positive flows",
                    point.flow_m3h
                ),
            ));
        }
        if index > 0 && point.flow_m3h <= curve[index - 1].flow_m3h {
            return Err(AzothError::invalid_input(
                "pressure_drop_curve",
                "the curve's flows are not strictly increasing",
            ));
        }
        if !point.drop_pa.is_finite() || point.drop_pa < 0.0 {
            return Err(AzothError::invalid_input(
                "pressure_drop_curve",
                format!(
                    "point {index}'s drop is {} Pa, and a drop is non-negative",
                    point.drop_pa
                ),
            ));
        }
        if index > 0 && point.drop_pa < curve[index - 1].drop_pa {
            return Err(AzothError::invalid_input(
                "pressure_drop_curve",
                "the curve's drops are not non-decreasing",
            ));
        }
    }
    Ok(())
}

fn validate_beta_curve(curve: &[BetaPoint]) -> Result<()> {
    for (index, point) in curve.iter().enumerate() {
        if !point.particle_size_m.is_finite() || point.particle_size_m <= 0.0 {
            return Err(AzothError::invalid_input(
                "performance_curve",
                format!(
                    "point {index}'s particle size is {} m, and a size is positive",
                    point.particle_size_m
                ),
            ));
        }
        if index > 0 && point.particle_size_m <= curve[index - 1].particle_size_m {
            return Err(AzothError::invalid_input(
                "performance_curve",
                "the curve's particle sizes are not strictly increasing",
            ));
        }
        if !point.beta_ratio.is_finite() || point.beta_ratio < 1.0 {
            return Err(AzothError::invalid_input(
                "performance_curve",
                format!(
                    "point {index}'s beta ratio is {}, and a ratio is at least one",
                    point.beta_ratio
                ),
            ));
        }
    }
    Ok(())
}
