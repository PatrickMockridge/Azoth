//! `process.throttling_valve` - the valve's kernel as a registered id.
//!
//! Spec: `specs/models/process/throttling_valve.toml`. The arithmetic is
//! [`crate::kernels::throttling_valve`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **The kernel needed no change to become a port.** NeqSim's `ThrottlingValve.run`
//! flashes `PHflash` at the inlet's enthalpy, which is what `Stream::from_ph` does, and
//! `acceptNegativeDP` defaults to `true` - so a stated outlet pressure above the inlet is
//! taken at face value, exactly as this does. The entry was wrong about its *parameters*,
//! not about its arithmetic.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::throttling_valve as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.throttling_valve`.
#[derive(Debug, Clone, PartialEq)]
pub struct ThrottlingValveResult {
    /// Outlet molar flow, mol/s: the inlet's, since a valve adds no moles.
    pub outlet_n: f64,
    /// Outlet composition, which is the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, which is the parameter the valve drops the stream to.
    pub outlet_p: Pressure,
    /// Outlet temperature, solved from the unchanged enthalpy at the outlet pressure.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy, which is the inlet's: the drop is isenthalpic.
    pub outlet_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ThrottlingValveResult {
    const CALC_ID: &'static str = "process.throttling_valve";
    const FIELDS: &'static [&'static str] = &[
        "outlet_n", "outlet_z", "outlet_p", "outlet_t", "outlet_h", "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Drop a stream to a lower pressure without heat or work.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the inlet pressure is not positive and
/// whatever the databank refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are six
pub fn throttling_valve(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    outlet_pressure: Pressure,
) -> Result<ThrottlingValveResult> {
    let spec = &model_gen::THROTTLING_VALVE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "outlet_pressure" => Some(outlet_pressure.value),
            "inlet_t" => Some(inlet_t.value),
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
    let out = kernel(&feed, outlet_pressure)?;

    Ok(ThrottlingValveResult {
        outlet_n: out.n,
        outlet_z: out.z,
        outlet_p: out.p,
        outlet_t: out.t,
        outlet_h: joules_per_mole(out.h.value),
        warnings,
    })
}
