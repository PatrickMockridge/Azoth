//! `process.heater` - the heater's kernel as a registered id.
//!
//! Spec: `specs/models/process/heater.toml`. The arithmetic is
//! [`crate::kernels::heater`]; what is here is the boundary a case and a cross-impl test
//! address.

use azoth_core::units::{
    MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole, watts,
};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::heater as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.heater`.
#[derive(Debug, Clone, PartialEq)]
pub struct HeaterResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure: `inlet_p - pressure_drop`.
    pub outlet_p: Pressure,
    /// Outlet temperature, which is the stated one or the one the shifted enthalpy reaches.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy.
    pub outlet_h: MolarEnergy,
    /// The duty moved, W: `inlet_n * (outlet_h - inlet_h)`, in every branch.
    pub outlet_duty: Power,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for HeaterResult {
    const CALC_ID: &'static str = "process.heater";
    const FIELDS: &'static [&'static str] = &[
        "outlet_n",
        "outlet_z",
        "outlet_p",
        "outlet_t",
        "outlet_h",
        "outlet_duty",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Heat or cool a stream to a stated temperature, or by a stated duty.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, if the pressure drop
/// leaves a non-positive pressure, or if a duty is asked of a stream carrying no flow, and
/// whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn heater(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
) -> Result<HeaterResult> {
    let spec = &model_gen::HEATER_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. The three optional inputs report `None` where they are omitted, which is a
    // skipped check and not a pass.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "outlet_temperature" => outlet_temperature.map(|t| t.value),
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
    let outcome = kernel(&feed, outlet_temperature, duty, pressure_drop)?;

    Ok(HeaterResult {
        outlet_n: outcome.outlet.n,
        outlet_z: outcome.outlet.z,
        outlet_p: outcome.outlet.p,
        outlet_t: outcome.outlet.t,
        outlet_h: joules_per_mole(outcome.outlet.h.value),
        outlet_duty: watts(outcome.duty.value),
        warnings,
    })
}
