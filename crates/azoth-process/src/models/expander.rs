//! `process.expander` - the Expander's kernel as a registered id.
//!
//! Spec: `specs/models/process/expander.toml`. The arithmetic is
//! [`crate::kernels::expander`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **The compressor's route with the efficiency multiplying rather than dividing**, which `Expander.run` is and which the physical rule requires: an expansion's isentropic difference is negative. **A `use_polytropic_calc` flag selects the stepped branch**, five increments by default.

use azoth_core::units::{
    MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole, watts,
};
use azoth_core::{Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::expander as kernel;
use crate::kernels::expander::expander_polytropic;
use crate::model_gen;
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
        warnings: Vec<Warning>,
    ) -> Self {
        let capacity = kernel::recovered_power_capacity(power, rated_recovered_power);
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
            warnings,
        }
    }
}

/// Drop a stream's pressure along an isentrope or in polytropic steps.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree or an efficiency is
/// outside `(0, 1]`, and whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eleven
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

    let feed = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;
    let out = if use_polytropic_calc.unwrap_or(false) {
        expander_polytropic(
            &feed,
            outlet_pressure,
            polytropic_efficiency.unwrap_or(1.0),
            number_of_steps.unwrap_or(5),
        )?
    } else {
        kernel(&feed, outlet_pressure, isentropic_efficiency)?
    };
    let power = watts(feed.n * (out.h.value - feed.h.value));

    Ok(ExpanderResult::of(
        &out,
        power,
        rated_recovered_power,
        warnings,
    ))
}
