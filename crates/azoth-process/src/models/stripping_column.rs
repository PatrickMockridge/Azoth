//! `process.stripping_column` - the tray stripper as a registered id.
//!
//! Spec: `specs/models/process/stripping_column.toml`. **The arithmetic is
//! [`crate::models::absorption_column`]'s**, because the class is: `StrippingColumn extends
//! AbsorptionColumn` and adds no equations - only the names of its two inlets and of its two
//! products. What is here is the boundary with those names, and the refusals the base carries.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning};
use serde::Serialize;

use crate::executor::json::{scalar, scalars, warnings as wire_warnings};
use crate::kernels::absorption_column::AbsorberOutcome;
use crate::models::absorption_column::absorption_column as absorber;

/// Result of `process.stripping_column`.
///
/// **The same record as the absorber's under the class's own getters**:
/// `getOverheadGasStream` is the stripped gas leaving the top tray and `getLeanLiquidStream` the
/// stripped liquid leaving the bottom, and there is no condenser and no reboiler, so no duties.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct StrippingColumnResult {
    /// Each tray's temperature, K, from the gas end at stage 0 up.
    #[serde(serialize_with = "scalars")]
    pub tray_temperature: Vec<ThermodynamicTemperature>,
    /// Each tray's pressure, Pa.
    #[serde(serialize_with = "scalars")]
    pub tray_pressure: Vec<Pressure>,
    /// Each tray's vapour traffic, mol/s.
    pub tray_gas_n: Vec<f64>,
    /// Each tray's liquid traffic, mol/s.
    pub tray_liquid_n: Vec<f64>,
    /// The stripped gas's molar flow, mol/s.
    pub overhead_gas_n: f64,
    /// The stripped gas's composition.
    pub overhead_gas_z: Vec<f64>,
    /// The stripped gas's pressure.
    #[serde(serialize_with = "scalar")]
    pub overhead_gas_p: Pressure,
    /// The stripped gas's temperature.
    #[serde(serialize_with = "scalar")]
    pub overhead_gas_t: ThermodynamicTemperature,
    /// The stripped gas's molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub overhead_gas_h: MolarEnergy,
    /// The stripped liquid's molar flow, mol/s.
    pub lean_liquid_n: f64,
    /// The stripped liquid's composition.
    pub lean_liquid_z: Vec<f64>,
    /// The stripped liquid's pressure.
    #[serde(serialize_with = "scalar")]
    pub lean_liquid_p: Pressure,
    /// The stripped liquid's temperature.
    #[serde(serialize_with = "scalar")]
    pub lean_liquid_t: ThermodynamicTemperature,
    /// The stripped liquid's molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub lean_liquid_h: MolarEnergy,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against both feeds, relative.
    pub mass_residual: f64,
    /// The enthalpy closure.
    pub energy_residual: f64,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl StrippingColumnResult {
    /// The result of one kernel call, which is the shared absorber's outcome under this id's names.
    ///
    /// **The warnings are the caller's**: a case's are `apply_checks`'s and a flowsheet's are the
    /// checker's, which report through the envelope rather than through a result.
    #[must_use]
    pub fn of(outcome: &AbsorberOutcome, warnings: Vec<Warning>) -> Self {
        Self {
            tray_temperature: outcome.trays.iter().map(|tray| tray.temperature).collect(),
            tray_pressure: outcome.trays.iter().map(|tray| tray.pressure).collect(),
            tray_gas_n: outcome.trays.iter().map(|tray| tray.gas_n).collect(),
            tray_liquid_n: outcome.trays.iter().map(|tray| tray.liquid_n).collect(),
            overhead_gas_n: outcome.gas_out.n,
            overhead_gas_z: outcome.gas_out.z.clone(),
            overhead_gas_p: outcome.gas_out.p,
            overhead_gas_t: outcome.gas_out.t,
            overhead_gas_h: joules_per_mole(outcome.gas_out.h.value),
            lean_liquid_n: outcome.liquid_out.n,
            lean_liquid_z: outcome.liquid_out.z.clone(),
            lean_liquid_p: outcome.liquid_out.p,
            lean_liquid_t: outcome.liquid_out.t,
            lean_liquid_h: joules_per_mole(outcome.liquid_out.h.value),
            iterations: outcome.iterations,
            temperature_residual: outcome.temperature_residual,
            mass_residual: outcome.mass_residual,
            energy_residual: outcome.energy_residual,
            warnings,
        }
    }
}

impl CalcResult for StrippingColumnResult {
    const CALC_ID: &'static str = "process.stripping_column";
    const FIELDS: &'static [&'static str] = &[
        "tray_temperature",
        "tray_pressure",
        "tray_gas_n",
        "tray_liquid_n",
        "overhead_gas_n",
        "overhead_gas_z",
        "overhead_gas_p",
        "overhead_gas_t",
        "overhead_gas_h",
        "lean_liquid_n",
        "lean_liquid_z",
        "lean_liquid_p",
        "lean_liquid_t",
        "lean_liquid_h",
        "iterations",
        "temperature_residual",
        "mass_residual",
        "energy_residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Solve a tray stripper.
///
/// **The base's own parameters under this class's names**, so the two ids read as the two
/// machines do: the stripping gas enters stage 0 and the rich liquid the top stage.
///
/// # Errors
/// Whatever the base refuses - the Murphree efficiencies of either kind, the strategies this
/// port does not carry, and its closure gate.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn stripping_column(
    stripping_gas_components: &[String],
    rich_liquid_components: &[String],
    stripping_gas_n: f64,
    stripping_gas_z: &[f64],
    stripping_gas_p: Pressure,
    stripping_gas_t: ThermodynamicTemperature,
    rich_liquid_n: f64,
    rich_liquid_z: &[f64],
    rich_liquid_p: Pressure,
    rich_liquid_t: ThermodynamicTemperature,
    number_of_stages: usize,
    top_pressure: Pressure,
    bottom_pressure: Pressure,
    temperature_tolerance: f64,
    max_iterations: usize,
    tray_temperatures: Option<&[f64]>,
    murphree_efficiency: Option<f64>,
    component_murphree_efficiency: Option<&[f64]>,
    max_allowable_gas_load_factor: Option<f64>,
    solver_type: Option<&str>,
) -> Result<StrippingColumnResult> {
    let out = absorber(
        stripping_gas_components,
        stripping_gas_n,
        stripping_gas_z,
        stripping_gas_p,
        stripping_gas_t,
        rich_liquid_components,
        rich_liquid_n,
        rich_liquid_z,
        rich_liquid_p,
        rich_liquid_t,
        number_of_stages,
        top_pressure,
        bottom_pressure,
        tray_temperatures,
        temperature_tolerance,
        max_iterations,
        murphree_efficiency,
        component_murphree_efficiency,
        max_allowable_gas_load_factor,
        solver_type,
    )?;

    Ok(StrippingColumnResult {
        tray_temperature: out.tray_temperature,
        tray_pressure: out.tray_pressure,
        tray_gas_n: out.tray_gas_n,
        tray_liquid_n: out.tray_liquid_n,
        overhead_gas_n: out.gas_out_n,
        overhead_gas_z: out.gas_out_z,
        overhead_gas_p: out.gas_out_p,
        overhead_gas_t: out.gas_out_t,
        overhead_gas_h: joules_per_mole(out.gas_out_h.value),
        lean_liquid_n: out.liquid_out_n,
        lean_liquid_z: out.liquid_out_z,
        lean_liquid_p: out.liquid_out_p,
        lean_liquid_t: out.liquid_out_t,
        lean_liquid_h: joules_per_mole(out.liquid_out_h.value),
        iterations: out.iterations,
        temperature_residual: out.temperature_residual,
        mass_residual: out.mass_residual,
        energy_residual: out.energy_residual,
        warnings: out.warnings,
    })
}
