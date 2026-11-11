//! `process.filter` - the filter's kernel as a registered id.
//!
//! Spec: `specs/models/process/filter.toml`. The arithmetic is
//! [`crate::kernels::filter`]; what is here is the boundary a case and a cross-impl test
//! address.

use azoth_core::units::{
    Area, DynamicViscosity, Length, Mass, MassRate, MolarEnergy, Pressure, TemperatureInterval,
    ThermodynamicTemperature, Time, Velocity, Volume, VolumeRate, cubic_meter_per_hour,
    joules_per_mole, kelvins, kilogram_per_hour, kilograms, meters, meters_per_second, pascals,
    seconds, square_meters,
};
use azoth_core::warning::{Warning, WarningCode};
use azoth_core::{AzothError, Result, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::filter::{
    self as kernel, BetaPoint, DropPoint, FilterSettings, Media, PressureDropModel,
};
use crate::kernels::filter_design::{self as design_kernel, DesignRequest, FilterType};
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
    /// The clean drop plus the loading increase, before bypass and the clamp.
    #[serde(serialize_with = "scalar")]
    pub unrestricted_drop: Pressure,
    /// The bypassed flow fraction.
    pub bypass_fraction: f64,
    /// `solids_loading / loading_capacity`, or zero without a finite capacity.
    pub loading_fraction: f64,
    /// The fraction past the breakthrough start, capped at one.
    pub breakthrough_fraction: f64,
    /// The nominal efficiency retained by breakthrough and bypass.
    pub current_removal_efficiency: f64,
    /// The carried solids loading, kg.
    #[serde(serialize_with = "scalar")]
    pub solids_loading: Mass,
    /// `unrestricted_drop / terminal_pressure_drop`.
    pub differential_pressure_utilization: f64,
    /// Either the terminal drop or the loading capacity has been reached.
    pub replacement_required: bool,
    /// The holdup over the inlet flow, absent where the class returns positive infinity.
    #[serde(serialize_with = "optional_scalar")]
    pub holdup_residence_time: Option<Time>,
    /// The design pressure, absent unless `mechanical_design`.
    #[serde(serialize_with = "optional_scalar")]
    pub design_pressure: Option<Pressure>,
    /// The design temperature.
    #[serde(serialize_with = "optional_scalar")]
    pub design_temperature: Option<ThermodynamicTemperature>,
    /// The required element count.
    pub required_elements: Option<f64>,
    /// The vessel's inside diameter.
    #[serde(serialize_with = "optional_scalar")]
    pub vessel_inner_diameter: Option<Length>,
    /// The vessel's outside diameter.
    #[serde(serialize_with = "optional_scalar")]
    pub vessel_outer_diameter: Option<Length>,
    /// The tangent-to-tangent length.
    #[serde(serialize_with = "optional_scalar")]
    pub vessel_length: Option<Length>,
    /// The shell thickness.
    #[serde(serialize_with = "optional_scalar")]
    pub shell_thickness: Option<Length>,
    /// The head thickness.
    #[serde(serialize_with = "optional_scalar")]
    pub head_thickness: Option<Length>,
    /// The effective element area.
    #[serde(serialize_with = "optional_scalar")]
    pub effective_element_area: Option<Area>,
    /// The calculated face velocity.
    #[serde(serialize_with = "optional_scalar")]
    pub calculated_face_velocity: Option<Velocity>,
    /// The required nozzle diameter.
    #[serde(serialize_with = "optional_scalar")]
    pub required_nozzle_diameter: Option<Length>,
    /// The selected nominal nozzle diameter.
    #[serde(serialize_with = "optional_scalar")]
    pub selected_nozzle_diameter: Option<Length>,
    /// The velocity through the selected nozzle.
    #[serde(serialize_with = "optional_scalar")]
    pub calculated_nozzle_velocity: Option<Velocity>,
    /// The empty vessel weight.
    #[serde(serialize_with = "optional_scalar")]
    pub empty_vessel_weight: Option<Mass>,
    /// The total equipped weight.
    #[serde(serialize_with = "optional_scalar")]
    pub total_equipped_weight: Option<Mass>,
    /// The collapse-rating utilization.
    pub collapse_pressure_utilization: Option<f64>,
    /// Whether both differential-pressure checks pass.
    pub differential_pressure_design_acceptable: Option<bool>,
    /// The vessel cross-section a granular design writes to the filter.
    #[serde(serialize_with = "optional_scalar")]
    pub design_media_area: Option<Area>,
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
    pub fn of(
        outcome: &kernel::FilterOutcome,
        design: Option<&design_kernel::FilterDesign>,
        warnings: Vec<Warning>,
    ) -> Self {
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
            unrestricted_drop: outcome.unrestricted_drop,
            bypass_fraction: outcome.bypass_fraction,
            loading_fraction: outcome.loading_fraction,
            breakthrough_fraction: outcome.breakthrough_fraction,
            current_removal_efficiency: outcome.current_removal_efficiency,
            solids_loading: kilograms(outcome.solids_loading_kg),
            differential_pressure_utilization: outcome.differential_pressure_utilization,
            replacement_required: outcome.replacement_required,
            holdup_residence_time: outcome.holdup_residence_s.map(seconds),
            design_pressure: design.map(|value| pascals(value.design_pressure_pa)),
            design_temperature: design.map(|value| kelvins(value.design_temperature_k)),
            required_elements: design.map(|value| value.required_elements),
            vessel_inner_diameter: design.map(|value| meters(value.inner_diameter_m)),
            vessel_outer_diameter: design.map(|value| meters(value.outer_diameter_m)),
            vessel_length: design.map(|value| meters(value.vessel_length_m)),
            shell_thickness: design.map(|value| meters(value.shell_thickness_m)),
            head_thickness: design.map(|value| meters(value.head_thickness_m)),
            effective_element_area: design.map(|value| square_meters(value.element_area_m2)),
            calculated_face_velocity: design
                .map(|value| meters_per_second(value.calculated_face_velocity_mps)),
            required_nozzle_diameter: design.map(|value| meters(value.required_nozzle_diameter_m)),
            selected_nozzle_diameter: design.map(|value| meters(value.selected_nozzle_diameter_m)),
            calculated_nozzle_velocity: design
                .map(|value| meters_per_second(value.calculated_nozzle_velocity_mps)),
            empty_vessel_weight: design.map(|value| kilograms(value.empty_vessel_weight_kg)),
            total_equipped_weight: design.map(|value| kilograms(value.total_equipped_weight_kg)),
            collapse_pressure_utilization: design.map(|value| value.collapse_pressure_utilization),
            differential_pressure_design_acceptable: design
                .map(|value| value.differential_pressure_design_acceptable),
            design_media_area: design
                .and_then(|value| value.design_media_area_m2)
                .map(square_meters),
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
    solids_loading: Option<Mass>,
    loading_capacity: Option<Mass>,
    pressure_drop_increase_at_capacity: Option<Pressure>,
    breakthrough_start_fraction: Option<f64>,
    bypass_cracking_pressure: Option<Pressure>,
    terminal_pressure_drop: Option<Pressure>,
    time_step: Option<Time>,
    holdup_volume: Option<Volume>,
    solids_loading_rate: Option<MassRate>,
    backwash_removal_rate: Option<MassRate>,
    regeneration_removal_rate: Option<MassRate>,
    backwash_active: Option<bool>,
    regeneration_active: Option<bool>,
    mechanical_design: Option<bool>,
    filter_type: Option<&str>,
    number_of_elements: Option<f64>,
    max_face_velocity: Option<Velocity>,
    element_area: Option<Area>,
    max_nozzle_velocity: Option<Velocity>,
    min_ld_ratio: Option<f64>,
    max_ld_ratio: Option<f64>,
    design_pressure_margin: Option<f64>,
    minimum_design_pressure_margin: Option<Pressure>,
    design_temperature_margin: Option<TemperatureInterval>,
    allowable_stress: Option<Pressure>,
    joint_efficiency: Option<f64>,
    corrosion_allowance: Option<Length>,
    element_collapse_pressure: Option<Pressure>,
    element_integrity_verified: Option<bool>,
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
            "solids_loading" => solids_loading.map(|value| value.value),
            "pressure_drop_increase_at_capacity" => {
                pressure_drop_increase_at_capacity.map(|value| value.value)
            }
            "breakthrough_start_fraction" => breakthrough_start_fraction,
            "bypass_cracking_pressure" => bypass_cracking_pressure.map(|value| value.value),
            "terminal_pressure_drop" => terminal_pressure_drop.map(|value| value.value),
            "time_step" => time_step.map(|value| value.value),
            "holdup_volume" => holdup_volume.map(|value| value.value),
            "solids_loading_rate" => solids_loading_rate.map(|value| value.value),
            "backwash_removal_rate" => backwash_removal_rate.map(|value| value.value),
            "regeneration_removal_rate" => regeneration_removal_rate.map(|value| value.value),
            "number_of_elements" => number_of_elements,
            "max_face_velocity" => max_face_velocity.map(|value| value.value),
            "element_area" => element_area.map(|value| value.value),
            "max_nozzle_velocity" => max_nozzle_velocity.map(|value| value.value),
            "min_ld_ratio" => min_ld_ratio,
            "max_ld_ratio" => max_ld_ratio,
            "design_pressure_margin" => design_pressure_margin,
            "minimum_design_pressure_margin" => {
                minimum_design_pressure_margin.map(|value| value.value)
            }
            "corrosion_allowance" => corrosion_allowance.map(|value| value.value),
            "element_collapse_pressure" => element_collapse_pressure.map(|value| value.value),
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
    for (name, value) in [
        ("solids_loading", solids_loading.map(|value| value.value)),
        (
            "loading_capacity",
            loading_capacity.map(|value| value.value),
        ),
        (
            "pressure_drop_increase_at_capacity",
            pressure_drop_increase_at_capacity.map(|value| value.value),
        ),
        ("breakthrough_start_fraction", breakthrough_start_fraction),
        ("time_step", time_step.map(|value| value.value)),
        ("holdup_volume", holdup_volume.map(|value| value.value)),
        (
            "solids_loading_rate",
            solids_loading_rate.map(|value| value.get::<kilogram_per_hour>()),
        ),
        (
            "backwash_removal_rate",
            backwash_removal_rate.map(|value| value.get::<kilogram_per_hour>()),
        ),
        (
            "regeneration_removal_rate",
            regeneration_removal_rate.map(|value| value.get::<kilogram_per_hour>()),
        ),
    ] {
        if let Some(value) = value {
            if !value.is_finite() {
                return Err(AzothError::invalid_input(
                    name,
                    format!("{value} is not a finite number"),
                ));
            }
        }
    }
    if solids_loading_rate.is_some() && inlet_particle_concentration.is_some() {
        return Err(AzothError::invalid_input(
            "solids_loading_rate",
            "the generic captured rate and the concentration model each disable the other in \
             `Filter`, so the last setter decides; the port refuses both rather than choosing",
        ));
    }

    let design_requested = mechanical_design.unwrap_or(false);
    for (name, stated) in [
        ("filter_type", filter_type.is_some()),
        ("number_of_elements", number_of_elements.is_some()),
        ("max_face_velocity", max_face_velocity.is_some()),
        ("element_area", element_area.is_some()),
        ("max_nozzle_velocity", max_nozzle_velocity.is_some()),
        ("min_ld_ratio", min_ld_ratio.is_some()),
        ("max_ld_ratio", max_ld_ratio.is_some()),
        ("design_pressure_margin", design_pressure_margin.is_some()),
        (
            "minimum_design_pressure_margin",
            minimum_design_pressure_margin.is_some(),
        ),
        (
            "design_temperature_margin",
            design_temperature_margin.is_some(),
        ),
        ("allowable_stress", allowable_stress.is_some()),
        ("joint_efficiency", joint_efficiency.is_some()),
        ("corrosion_allowance", corrosion_allowance.is_some()),
        (
            "element_collapse_pressure",
            element_collapse_pressure.is_some(),
        ),
        (
            "element_integrity_verified",
            element_integrity_verified.is_some(),
        ),
    ] {
        if stated && !design_requested {
            return Err(AzothError::invalid_input(
                name,
                "a design input is stated without `mechanical_design`",
            ));
        }
    }
    for (name, value) in [
        ("number_of_elements", number_of_elements),
        (
            "max_face_velocity",
            max_face_velocity.map(|value| value.value),
        ),
        ("element_area", element_area.map(|value| value.value)),
        (
            "max_nozzle_velocity",
            max_nozzle_velocity.map(|value| value.value),
        ),
        ("min_ld_ratio", min_ld_ratio),
        ("max_ld_ratio", max_ld_ratio),
        ("design_pressure_margin", design_pressure_margin),
        (
            "minimum_design_pressure_margin",
            minimum_design_pressure_margin.map(|value| value.value),
        ),
        (
            "design_temperature_margin",
            design_temperature_margin.map(|value| value.value),
        ),
        (
            "allowable_stress",
            allowable_stress.map(|value| value.value),
        ),
        ("joint_efficiency", joint_efficiency),
        (
            "corrosion_allowance",
            corrosion_allowance.map(|value| value.value),
        ),
        (
            "element_collapse_pressure",
            element_collapse_pressure.map(|value| value.value),
        ),
    ] {
        if let Some(value) = value {
            if !value.is_finite() {
                return Err(AzothError::invalid_input(
                    name,
                    format!("{value} is not a finite number"),
                ));
            }
        }
    }
    let filter_type = match filter_type {
        Some(spelling) => FilterType::parse(spelling).ok_or_else(|| {
            AzothError::invalid_input(
                "filter_type",
                format!("{spelling} is not one of the eight filter constructions"),
            )
        })?,
        None => FilterType::Cartridge,
    };
    let number_of_elements = number_of_elements.unwrap_or(1.0);
    if number_of_elements < 1.0 || number_of_elements.fract() != 0.0 {
        return Err(AzothError::invalid_input(
            "number_of_elements",
            format!(
                "the element count is {number_of_elements}, and `setNumberOfElements` takes an \
                 integer of at least one"
            ),
        ));
    }
    let min_ld_ratio = min_ld_ratio.unwrap_or(2.0);
    let max_ld_ratio = max_ld_ratio.unwrap_or(5.0);
    if min_ld_ratio <= 0.0 || max_ld_ratio <= 0.0 || min_ld_ratio > max_ld_ratio {
        return Err(AzothError::invalid_input(
            "min_ld_ratio",
            format!(
                "the L/D window is {min_ld_ratio} to {max_ld_ratio}, and it must be positive with \
                 the minimum no greater than the maximum"
            ),
        ));
    }
    if let Some(cracking) = bypass_cracking_pressure {
        if !cracking.value.is_finite() || cracking.value <= 0.0 {
            return Err(AzothError::invalid_input(
                "bypass_cracking_pressure",
                format!(
                    "the cracking pressure is {} Pa, and `setBypassCrackingDeltaP` needs a \
                     finite positive one",
                    cracking.value
                ),
            ));
        }
    }
    if let Some(terminal) = terminal_pressure_drop {
        if !terminal.value.is_finite() || terminal.value <= 0.0 {
            return Err(AzothError::invalid_input(
                "terminal_pressure_drop",
                format!(
                    "the terminal drop is {} Pa, and `setTerminalDeltaP` needs a finite positive \
                     one",
                    terminal.value
                ),
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
    let loading = kernel::LoadingRequest {
        solids_loading_kg: solids_loading.map(|value| value.value),
        loading_capacity_kg: loading_capacity.map(|value| value.value),
        pressure_drop_increase_at_capacity_pa: pressure_drop_increase_at_capacity
            .map(|value| value.value),
        breakthrough_start_fraction,
        bypass_cracking_pressure_pa: bypass_cracking_pressure.map(|value| value.value),
        terminal_pressure_drop_pa: terminal_pressure_drop.map(|value| value.value),
        time_step_s: time_step.map(|value| value.value),
        holdup_volume_m3: holdup_volume.map(|value| value.value),
        solids_loading_rate_kg_h: solids_loading_rate.map(|value| value.get::<kilogram_per_hour>()),
        backwash_removal_rate_kg_h: backwash_removal_rate
            .map(|value| value.get::<kilogram_per_hour>()),
        regeneration_removal_rate_kg_h: regeneration_removal_rate
            .map(|value| value.get::<kilogram_per_hour>()),
        backwash_active: backwash_active.unwrap_or(false),
        regeneration_active: regeneration_active.unwrap_or(false),
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
        &loading,
    )?;

    // **The clamp is reported rather than hidden.** `Filter.run` logs a warning when the drop
    // it applies is below the one it was given, and the whole reason `applied_drop` is a
    // declared output is that the two can differ - so a result that carried the clamped value
    // and said nothing would be the silent-improvement failure this library is written
    // against. A negative `pressure_drop` is clamped by the same path and warned about by the
    // range check above.
    if (outcome.applied_drop.value - outcome.restricted_drop.value).abs() > 0.0 {
        warnings.push(Warning::for_field(
            WarningCode::OutOfValidRange,
            "pressure_drop",
            format!(
                "a restricted drop of {} Pa was computed and {} Pa was applied: the outlet is \
                 held a millionth of a bar above vacuum, which is the clamp `Filter.run` applies",
                outcome.restricted_drop.value, outcome.applied_drop.value
            ),
        ));
    }

    let design = if design_requested {
        let terminal = loading
            .terminal_pressure_drop_pa
            .unwrap_or(kernel::DEFAULT_TERMINAL_DROP_PA);
        let request = DesignRequest {
            filter_type,
            number_of_elements: number_of_elements.round() as u32,
            max_face_velocity_mps: max_face_velocity.map(|value| value.value),
            element_area_m2: element_area.map(|value| value.value),
            max_nozzle_velocity_mps: max_nozzle_velocity.map(|value| value.value),
            min_ld_ratio,
            max_ld_ratio,
            design_pressure_margin: design_pressure_margin.unwrap_or(1.10),
            minimum_design_pressure_margin_pa: minimum_design_pressure_margin
                .map_or(350_000.0, |value| value.value),
            design_temperature_margin_k: design_temperature_margin
                .map_or(25.0, |value| value.value),
            allowable_stress_pa: allowable_stress.map_or(138.0e6, |value| value.value),
            joint_efficiency: joint_efficiency.unwrap_or(0.85),
            corrosion_allowance_m: corrosion_allowance.map_or(0.003, |value| value.value),
            element_collapse_pressure_pa: element_collapse_pressure
                .map_or(1.0e6, |value| value.value),
            element_integrity_verified: element_integrity_verified.unwrap_or(false),
            performance_data: !beta_curve.is_empty() || efficiency > 0.0,
            media_bed_depth_m: media.bed_depth_m,
            terminal_pressure_drop_pa: terminal,
            unrestricted_drop_pa: outcome.unrestricted_drop.value,
            terminal_pressure_utilization: outcome.differential_pressure_utilization,
        };
        let design = design_kernel::design(&feed, &request)?;
        for warning in &design.warnings {
            warnings.push(Warning::for_field(
                WarningCode::OutOfValidRange,
                "mechanical_design",
                warning.message(),
            ));
        }
        Some(design)
    } else {
        None
    };

    Ok(FilterResult::of(&outcome, design.as_ref(), warnings))
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
