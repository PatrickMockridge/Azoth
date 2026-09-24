//! `process.flare` - the flare's kernel as a registered id.
//!
//! Spec: `specs/models/process/flare.toml`. The arithmetic is [`crate::kernels::flare`]; what
//! is here is the boundary a case and a cross-impl test address.
//!
//! **The record passes through and two numbers are added.** `run` clones the inlet into the
//! outlet, so the five fields of the inlet port are the five fields of the outlet port - and
//! the class's own two, the heat duty and the CO2 emission, are declared as outputs because
//! they are what the machine reports.

use azoth_core::units::{
    MassRate, MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole,
};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::flare as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.flare`.
#[derive(Debug, Clone, PartialEq)]
pub struct FlareResult {
    /// Product molar flow, mol/s.
    pub product_n: f64,
    /// Product composition.
    pub product_z: Vec<f64>,
    /// Product pressure.
    pub product_p: Pressure,
    /// Product temperature.
    pub product_t: ThermodynamicTemperature,
    /// Product molar enthalpy.
    pub product_h: MolarEnergy,
    /// The heat the flare releases, W.
    pub heat_duty: Power,
    /// The carbon dioxide the combustion forms, kg/s.
    pub co2_emission: MassRate,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for FlareResult {
    const CALC_ID: &'static str = "process.flare";
    const FIELDS: &'static [&'static str] = &[
        "product_n",
        "product_z",
        "product_p",
        "product_t",
        "product_h",
        "heat_duty",
        "co2_emission",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// A flare's steady state: the record through, and the two numbers beside it.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if a component has no row in ISO 6976's table or
/// none in the element table, and whatever the flashes refuse.
pub fn flare(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
) -> Result<FlareResult> {
    let spec = &model_gen::FLARE_SPEC;
    let mut warnings = Vec::new();
    // The entry declares no parameters, so the only bound this can run is the one the family
    // puts on the flow - and `heat_duty` and `co2_emission` are derived checks rather than
    // input ones.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_n" => Some(inlet_n),
            _ => None,
        },
        &mut warnings,
    )?;

    let inlet = Stream::from_pt(
        components.to_vec(),
        inlet_z.to_vec(),
        inlet_n,
        inlet_p,
        inlet_t,
    )?;
    let (product, numbers) = kernel(&inlet)?;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "heat_duty" => Some(numbers.heat_duty.value),
            "co2_emission" => Some(numbers.co2_emission.value),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(FlareResult {
        product_n: product.n,
        product_z: product.z,
        product_p: product.p,
        product_t: product.t,
        product_h: joules_per_mole(product.h.value),
        heat_duty: numbers.heat_duty,
        co2_emission: numbers.co2_emission,
        warnings,
    })
}
