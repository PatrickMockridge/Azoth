//! `process.expander` - the Expander's kernel as a registered id.
//!
//! Spec: `specs/models/process/expander.toml`. The arithmetic is
//! [`crate::kernels::expander`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **The compressor's route with the efficiency multiplying rather than dividing**, which `Expander.run` is and which the physical rule requires: an expansion's isentropic difference is negative.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::expander as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.expander`.
#[derive(Debug, Clone, PartialEq)]
pub struct ExpanderResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, which is the parameter the machine takes the stream to.
    pub outlet_p: Pressure,
    /// Outlet temperature, solved from the shifted enthalpy.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy: the inlet's plus the isentropic step, over or times the
    /// efficiency.
    pub outlet_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ExpanderResult {
    const CALC_ID: &'static str = "process.expander";
    const FIELDS: &'static [&'static str] = &[
        "outlet_n", "outlet_z", "outlet_p", "outlet_t", "outlet_h", "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Expander a stream's pressure along an isentrope.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree or the efficiency is
/// outside `(0, 1]`, and whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seven
pub fn expander(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_pressure: Pressure,
    isentropic_efficiency: f64,
) -> Result<ExpanderResult> {
    let spec = &model_gen::EXPANDER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "isentropic_efficiency" => Some(isentropic_efficiency),
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
    let out = kernel(&feed, outlet_pressure, isentropic_efficiency)?;

    Ok(ExpanderResult {
        outlet_n: out.n,
        outlet_z: out.z,
        outlet_p: out.p,
        outlet_t: out.t,
        outlet_h: joules_per_mole(out.h.value),
        warnings,
    })
}
