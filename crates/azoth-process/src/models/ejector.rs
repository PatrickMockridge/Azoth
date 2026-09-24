//! `process.ejector` - the ejector's kernel as a registered id.
//!
//! Spec: `specs/models/process/ejector.toml`. The arithmetic is [`crate::kernels::ejector`];
//! what is here is the boundary a case and a cross-impl test address.
//!
//! **Two inlets carrying two fluids, so ten fields of inlet and five of outlet.** A
//! two-inlet machine's inlets are not one port: the motive and the suction may be different
//! mixtures, and `_dispatch.fluid_inputs` derives each `components` prefix from the input's
//! own name - the arrangement `process.heat_exchanger` introduced.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::ejector::{EjectorSetup, ejector as kernel};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.ejector`.
#[derive(Debug, Clone, PartialEq)]
pub struct EjectorResult {
    /// Outlet molar flow, mol/s.
    pub outlet_n: f64,
    /// Outlet composition.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure.
    pub outlet_p: Pressure,
    /// Outlet temperature.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy.
    pub outlet_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for EjectorResult {
    const CALC_ID: &'static str = "process.ejector";
    const FIELDS: &'static [&'static str] = &[
        "outlet_n", "outlet_z", "outlet_p", "outlet_t", "outlet_h", "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Expand a motive stream, entrain a suction stream with it, and diffuse the mixture.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if an efficiency is outside `(0, 1]`, the
/// discharge pressure is not positive, or a stream's shapes disagree, and whatever the
/// flashes refuse.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are fourteen
pub fn ejector(
    motive_components: &[String],
    motive_n: f64,
    motive_z: &[f64],
    motive_p: Pressure,
    motive_t: ThermodynamicTemperature,
    suction_components: &[String],
    suction_n: f64,
    suction_z: &[f64],
    suction_p: Pressure,
    suction_t: ThermodynamicTemperature,
    discharge_pressure: Pressure,
    motive_nozzle_efficiency: f64,
    suction_nozzle_efficiency: f64,
    mixing_efficiency: f64,
    diffuser_efficiency: f64,
) -> Result<EjectorResult> {
    let spec = &model_gen::EJECTOR_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives: a bound on a pressure or an efficiency is resolvable here, and a bound this
    // leaves unrun would read as validation that did not happen.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "discharge_pressure" => Some(discharge_pressure.value),
            "motive_nozzle_efficiency" => Some(motive_nozzle_efficiency),
            "suction_nozzle_efficiency" => Some(suction_nozzle_efficiency),
            "mixing_efficiency" => Some(mixing_efficiency),
            "diffuser_efficiency" => Some(diffuser_efficiency),
            _ => None,
        },
        &mut warnings,
    )?;

    let motive = Stream::from_pt(
        motive_components.to_vec(),
        motive_z.to_vec(),
        motive_n,
        motive_p,
        motive_t,
    )?;
    let suction = Stream::from_pt(
        suction_components.to_vec(),
        suction_z.to_vec(),
        suction_n,
        suction_p,
        suction_t,
    )?;
    let outlet = kernel(
        &motive,
        &suction,
        EjectorSetup {
            discharge_pressure,
            motive_nozzle_efficiency,
            suction_nozzle_efficiency,
            mixing_efficiency,
            diffuser_efficiency,
        },
    )?;

    Ok(EjectorResult {
        outlet_n: outlet.n,
        outlet_z: outlet.z,
        outlet_p: outlet.p,
        outlet_t: outlet.t,
        outlet_h: joules_per_mole(outlet.h.value),
        warnings,
    })
}
