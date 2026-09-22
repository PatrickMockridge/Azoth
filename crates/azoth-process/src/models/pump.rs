//! `process.pump` - the pump's kernel as a registered id.
//!
//! Spec: `specs/models/process/pump.toml`. The arithmetic is
//! [`crate::kernels::pump`]; what is here is the boundary a case and a cross-impl test
//! address, which is the same split every other namespace makes.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::pump as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.pump`.
#[derive(Debug, Clone, PartialEq)]
pub struct PumpResult {
    /// Molar flow out, mol/s. A pump moves the stream it is given, so this is the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, one entry per component.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, which is the parameter the pump raises the stream to.
    pub outlet_p: Pressure,
    /// Outlet temperature, solved from the shifted enthalpy.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy: the inlet's plus the shaft work per mole.
    pub outlet_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PumpResult {
    const CALC_ID: &'static str = "process.pump";
    const FIELDS: &'static [&'static str] = &[
        "outlet_n", "outlet_z", "outlet_p", "outlet_t", "outlet_h", "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Raise a stream's pressure, adding the pump's work as enthalpy.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the shapes disagree or the efficiency is
/// outside `(0, 1]`, and whatever the databank or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seven
pub fn pump(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_pressure: Pressure,
    isentropic_efficiency: f64,
) -> Result<PumpResult> {
    let spec = &model_gen::PUMP_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives: a closure that filters by name emits a warning for every check it does not
    // feed, which reads as a defect in the caller's state.
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

    Ok(PumpResult {
        outlet_n: out.n,
        outlet_z: out.z,
        outlet_p: out.p,
        outlet_t: out.t,
        outlet_h: joules_per_mole(out.h.value),
        warnings,
    })
}
