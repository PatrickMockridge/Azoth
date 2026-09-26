//! `process.filter` - the filter's kernel as a registered id.
//!
//! Spec: `specs/models/process/filter.toml`. The arithmetic is
//! [`crate::kernels::filter`]; what is here is the boundary a case and a cross-impl test
//! address.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::warning::{Warning, WarningCode};
use azoth_core::{CalcResult, Result, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::filter as kernel;
use crate::model_gen;
use crate::stream::Stream;

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
    /// The drop actually applied, which is the requested one unless it was clamped.
    #[serde(serialize_with = "scalar")]
    pub applied_drop: Pressure,
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
            applied_drop: outcome.applied_drop,
            warnings,
        }
    }
}

impl CalcResult for FilterResult {
    const CALC_ID: &'static str = "process.filter";
    const FIELDS: &'static [&'static str] = &[
        "outlet_n",
        "outlet_z",
        "outlet_p",
        "outlet_t",
        "outlet_h",
        "applied_drop",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Drop a stream's pressure by a fixed amount at a constant temperature.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree, and whatever the databank
/// or the flash refuses. A drop past the inlet pressure is **not** an error: the class clamps
/// it and so does this, saying so in a warning.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are six
pub fn filter(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    pressure_drop: Pressure,
) -> Result<FilterResult> {
    let spec = &model_gen::FILTER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "pressure_drop" => Some(pressure_drop.value),
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
    let outcome = kernel(&feed, pressure_drop)?;

    // **The clamp is reported rather than hidden.** `Filter.run` logs a warning when the drop
    // it applies is below the one it was given, and the whole reason `applied_drop` is a
    // declared output is that the two can differ - so a result that carried the clamped value
    // and said nothing would be the silent-improvement failure this library is written
    // against. A negative `pressure_drop` is clamped by the same path and warned about by the
    // range check above.
    if (outcome.applied_drop.value - pressure_drop.value).abs() > 0.0 {
        warnings.push(Warning::for_field(
            WarningCode::OutOfValidRange,
            "pressure_drop",
            format!(
                "a drop of {} Pa was asked for and {} Pa was applied: the outlet is held a \
                 millionth of a bar above vacuum, which is the clamp `Filter.run` applies",
                pressure_drop.value, outcome.applied_drop.value
            ),
        ));
    }

    Ok(FilterResult::of(&outcome, warnings))
}
