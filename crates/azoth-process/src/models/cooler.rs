//! `process.cooler` - `Heater.run` reached through `Cooler`, as a registered id.
//!
//! Spec: `specs/models/process/cooler.toml`. The arithmetic is
//! [`crate::kernels::cooler`], which is [`crate::kernels::heater`]: `Cooler` overrides
//! `runTransient` and some getters and not `run`, and the probe's two captures are
//! byte-identical.

use azoth_core::units::{
    MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole, watts,
};
use azoth_core::{CalcResult, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::cooler as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.cooler`.
///
/// **The same fields as [`crate::models::HeaterResult`], and its own type.** A result class
/// is the registry's handle on an id - `result_fields` answers by id - so one class shared
/// by two ids would be one handle for two things, which is the shape every other pair in
/// this registry refuses too.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CoolerResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure: `inlet_p - pressure_drop`.
    #[serde(serialize_with = "scalar")]
    pub outlet_p: Pressure,
    /// Outlet temperature, which is the stated one or the one the shifted enthalpy reaches.
    #[serde(serialize_with = "scalar")]
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub outlet_h: MolarEnergy,
    /// The duty moved, W - negative when heat was removed, which is this entry's usual case.
    #[serde(serialize_with = "scalar")]
    pub outlet_duty: Power,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl CoolerResult {
    /// The result of one kernel call, for the reason [`crate::models::HeaterResult::of`] gives
    /// about the warnings.
    #[must_use]
    pub fn of(outcome: &crate::kernels::heater::HeaterOutcome, warnings: Vec<Warning>) -> Self {
        Self {
            outlet_n: outcome.outlet.n,
            outlet_z: outcome.outlet.z.clone(),
            outlet_p: outcome.outlet.p,
            outlet_t: outcome.outlet.t,
            outlet_h: joules_per_mole(outcome.outlet.h.value),
            outlet_duty: watts(outcome.duty.value),
            warnings,
        }
    }
}

impl CalcResult for CoolerResult {
    const CALC_ID: &'static str = "process.cooler";
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

/// Cool a stream to a stated temperature, or by a stated duty.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, if a temperature and a
/// duty are both given, if the pressure drop leaves a non-positive pressure, or if a duty is
/// asked of a stream carrying no flow, and whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn cooler(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
) -> Result<CoolerResult> {
    let spec = &model_gen::COOLER_SPEC;
    let mut warnings = Vec::new();
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

    Ok(CoolerResult::of(&outcome, warnings))
}
