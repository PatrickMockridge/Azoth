//! `process.expander` - the Expander's kernel as a registered id.
//!
//! Spec: `specs/models/process/expander.toml`. The arithmetic is
//! [`crate::kernels::expander`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **The compressor's route with the efficiency multiplying rather than dividing**, which `Expander.run` is and which the physical rule requires: an expansion's isentropic difference is negative. **A `use_polytropic_calc` flag selects the stepped branch**, five increments by default.

use azoth_core::units::{
    Length, Mass, MolarEnergy, Power, Pressure, ThermodynamicTemperature, Velocity,
    joules_per_mole, watts,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::expander as kernel;
use crate::kernels::expander::expander_polytropic;
use crate::model_gen;
use crate::rotating_design::expander::{ExpanderDesign, design as design_expander};
use crate::stream::Stream;

/// Result of `process.expander`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ExpanderResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, which is the parameter the machine takes the stream to.
    #[serde(serialize_with = "scalar")]
    pub outlet_p: Pressure,
    /// Outlet temperature, solved from the shifted enthalpy.
    #[serde(serialize_with = "scalar")]
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy: the inlet's plus the step.
    #[serde(serialize_with = "scalar")]
    pub outlet_h: MolarEnergy,
    /// The recovered shaft power, W - negative, as the class's `getPower()` is.
    #[serde(serialize_with = "scalar")]
    pub power: Power,
    /// The `recoveredPower` constraint's utilization, when a positive rating is set.
    pub recovered_power_utilization: Option<f64>,
    /// `getMaxUtilization`, 0.0 without an enabled constraint.
    pub max_utilization: f64,
    /// `isCapacityExceeded`.
    pub capacity_exceeded: bool,
    /// `isHardLimitExceeded`.
    pub hard_limit_exceeded: bool,
    /// `getBottleneckConstraint`'s class-side name.
    pub bottleneck_constraint: Option<String>,
    /// The efficiency a stated power or an outlet-temperature solve produced.
    pub reported_isentropic_efficiency: Option<f64>,
    /// `ExpanderType`: `radial_inflow`, `axial` or `mixed_flow`.
    pub design_expander_type: Option<String>,
    /// `LoadType`: `generator`, `compressor` or `brake`.
    pub design_load_type: Option<String>,
    /// The stage count `selectExpanderType` assigns.
    pub design_number_of_stages: Option<f64>,
    /// The wheel's outer diameter, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_wheel_diameter: Option<Length>,
    /// The wheel's tip speed, m/s.
    #[serde(serialize_with = "optional_scalar")]
    pub design_tip_speed: Option<Velocity>,
    /// The shaft diameter, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_shaft_diameter: Option<Length>,
    /// The rated speed, rpm.
    pub design_rated_speed: Option<f64>,
    /// The maximum continuous speed, rpm.
    pub design_max_continuous_speed: Option<f64>,
    /// The trip speed, rpm.
    pub design_trip_speed: Option<f64>,
    /// The first critical speed, rpm.
    pub design_first_critical_speed: Option<f64>,
    /// The casing design pressure, Pa.
    #[serde(serialize_with = "optional_scalar")]
    pub design_casing_pressure: Option<Pressure>,
    /// The casing design temperature, K.
    #[serde(serialize_with = "optional_scalar")]
    pub design_casing_temperature: Option<ThermodynamicTemperature>,
    /// The casing wall thickness, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_casing_wall_thickness: Option<Length>,
    /// The bearing span, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_bearing_span: Option<Length>,
    /// The wheel weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub design_wheel_weight: Option<Mass>,
    /// The shaft weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub design_shaft_weight: Option<Mass>,
    /// The generator weight, kg; zero for a brake load.
    #[serde(serialize_with = "optional_scalar")]
    pub design_generator_weight: Option<Mass>,
    /// The seal type.
    pub design_seal_type: Option<String>,
    /// The bearing type.
    pub design_bearing_type: Option<String>,
    /// The casing inner diameter, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_inner_diameter: Option<Length>,
    /// The casing outer diameter, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_outer_diameter: Option<Length>,
    /// The tangent-to-tangent length, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_tantan_length: Option<Length>,
    /// The base class's vessel weight, kg: the casing.
    #[serde(serialize_with = "optional_scalar")]
    pub design_weight_vessel: Option<Mass>,
    /// The base class's internals weight, kg: the rotor bundle.
    #[serde(serialize_with = "optional_scalar")]
    pub design_weight_internals: Option<Mass>,
    /// The piping, lube and seal-system weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub design_weight_piping: Option<Mass>,
    /// The electrical and instrument weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub design_weight_electro_instrument: Option<Mass>,
    /// The structural steel and baseplate weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub design_weight_structural_steel: Option<Mass>,
    /// The total weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub design_weight_total: Option<Mass>,
    /// The module length, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_module_length: Option<Length>,
    /// The module width, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_module_width: Option<Length>,
    /// The module height, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_module_height: Option<Length>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl ExpanderResult {
    /// The record from the outlet and the recovered-power map.
    pub(crate) fn of(
        out: &Stream,
        power: Power,
        rated_recovered_power: Option<Power>,
        reported_isentropic_efficiency: Option<f64>,
        design: Option<ExpanderDesign>,
        warnings: Vec<Warning>,
    ) -> Self {
        let capacity = kernel::recovered_power_capacity(power, rated_recovered_power);
        let d = design.as_ref();
        Self {
            outlet_n: out.n,
            outlet_z: out.z.clone(),
            outlet_p: out.p,
            outlet_t: out.t,
            outlet_h: joules_per_mole(out.h.value),
            power,
            recovered_power_utilization: capacity.recovered_power_utilization,
            max_utilization: capacity.max_utilization,
            capacity_exceeded: capacity.capacity_exceeded,
            hard_limit_exceeded: capacity.hard_limit_exceeded,
            bottleneck_constraint: capacity.bottleneck_constraint.map(str::to_owned),
            reported_isentropic_efficiency,
            design_expander_type: d.map(|value| value.expander_type.to_owned()),
            design_load_type: d.map(|value| value.load_type.to_owned()),
            design_number_of_stages: d.map(|value| f64::from(value.number_of_stages)),
            design_wheel_diameter: d.map(|value| value.wheel_diameter),
            design_tip_speed: d.map(|value| value.tip_speed),
            design_shaft_diameter: d.map(|value| value.shaft_diameter),
            design_rated_speed: d.map(|value| value.rated_speed),
            design_max_continuous_speed: d.map(|value| value.max_continuous_speed),
            design_trip_speed: d.map(|value| value.trip_speed),
            design_first_critical_speed: d.map(|value| value.first_critical_speed),
            design_casing_pressure: d.map(|value| value.casing_pressure),
            design_casing_temperature: d.map(|value| value.casing_temperature),
            design_casing_wall_thickness: d.map(|value| value.casing_wall_thickness),
            design_bearing_span: d.map(|value| value.bearing_span),
            design_wheel_weight: d.map(|value| value.wheel_weight),
            design_shaft_weight: d.map(|value| value.shaft_weight),
            design_generator_weight: d.map(|value| value.generator_weight),
            design_seal_type: d.map(|value| value.seal_type.to_owned()),
            design_bearing_type: d.map(|value| value.bearing_type.to_owned()),
            design_inner_diameter: d.map(|value| value.inner_diameter),
            design_outer_diameter: d.map(|value| value.outer_diameter),
            design_tantan_length: d.map(|value| value.tantan_length),
            design_weight_vessel: d.map(|value| value.weight_vessel),
            design_weight_internals: d.map(|value| value.weight_internals),
            design_weight_piping: d.map(|value| value.weight_piping),
            design_weight_electro_instrument: d.map(|value| value.weight_electro_instrument),
            design_weight_structural_steel: d.map(|value| value.weight_structural_steel),
            design_weight_total: d.map(|value| value.weight_total),
            design_module_length: d.map(|value| value.module_length),
            design_module_width: d.map(|value| value.module_width),
            design_module_height: d.map(|value| value.module_height),
            warnings,
        }
    }
}

/// Drop a stream's pressure along an isentrope or in polytropic steps.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree or an efficiency is
/// outside `(0, 1]`, and whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are thirteen
pub fn expander(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_pressure: Pressure,
    isentropic_efficiency: f64,
    use_polytropic_calc: Option<bool>,
    polytropic_efficiency: Option<f64>,
    number_of_steps: Option<usize>,
    rated_recovered_power: Option<Power>,
    shaft_power: Option<Power>,
    outlet_temperature: Option<ThermodynamicTemperature>,
    mechanical_design: Option<bool>,
) -> Result<ExpanderResult> {
    let spec = &model_gen::EXPANDER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "isentropic_efficiency" => Some(isentropic_efficiency),
            "polytropic_efficiency" => polytropic_efficiency,
            _ => None,
        },
        &mut warnings,
    )?;

    if shaft_power.is_some() && outlet_temperature.is_some() {
        return Err(AzothError::invalid_input(
            "shaft_power",
            "a stated power and an outlet temperature are two modes; the class's `setPower` \
             and `solveEfficiency` are one at a time",
        ));
    }
    if (shaft_power.is_some() || outlet_temperature.is_some())
        && use_polytropic_calc.unwrap_or(false)
    {
        return Err(AzothError::invalid_input(
            "shaft_power",
            "the stated-power and outlet-temperature modes are the isentropic branch's",
        ));
    }
    for (name, value) in [
        ("shaft_power", shaft_power.map(|value| value.value)),
        (
            "outlet_temperature",
            outlet_temperature.map(|value| value.value),
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

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;
    let (out, reported) = if let Some(shaft_power) = shaft_power {
        let (out, efficiency) = kernel::expander_at_power(&feed, outlet_pressure, shaft_power)?;
        (out, efficiency)
    } else if let Some(target) = outlet_temperature {
        let (out, efficiency) =
            kernel::solve_efficiency(&feed, outlet_pressure, target, isentropic_efficiency)?;
        (out, Some(efficiency))
    } else if use_polytropic_calc.unwrap_or(false) {
        (
            expander_polytropic(
                &feed,
                outlet_pressure,
                polytropic_efficiency.unwrap_or(1.0),
                number_of_steps.unwrap_or(5),
            )?,
            None,
        )
    } else {
        (kernel(&feed, outlet_pressure, isentropic_efficiency)?, None)
    };
    let power = watts(feed.n * (out.h.value - feed.h.value));
    let design = if mechanical_design.unwrap_or(false) {
        Some(design_expander(&feed, out.t, power)?)
    } else {
        None
    };

    Ok(ExpanderResult::of(
        &out,
        power,
        rated_recovered_power,
        reported,
        design,
        warnings,
    ))
}
