//! `process.pipe` - the line's kernel as a registered id.
//!
//! Spec: `specs/models/process/pipe.toml`. The arithmetic is [`crate::kernels::pipe`]; what
//! is here is the boundary a case and a cross-impl test address.

use azoth_core::units::{
    Length, MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole, pascals,
};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::pipe as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.pipe`.
#[derive(Debug, Clone, PartialEq)]
pub struct PipeResult {
    /// Molar flow out, mol/s: the inlet's.
    pub outlet_n: f64,
    /// Outlet composition, the inlet's.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, which the line *solves* rather than takes.
    pub outlet_p: Pressure,
    /// Outlet temperature, which is the inlet's.
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy, which is not the inlet's - the drop is isothermal.
    pub outlet_h: MolarEnergy,
    /// `inlet_p - outlet_p`, which is `AdiabaticPipe.getPressureDrop()`.
    pub pressure_drop: Pressure,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PipeResult {
    const CALC_ID: &'static str = "process.pipe";
    const FIELDS: &'static [&'static str] = &[
        "outlet_n",
        "outlet_z",
        "outlet_p",
        "outlet_t",
        "outlet_h",
        "pressure_drop",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Drop a stream's pressure along a line.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the geometry is not positive, and whatever the
/// flash, the density, the viscosity or the phase label refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn pipe(
    components: &[String],
    inlet_n: f64,
    inlet_z: &[f64],
    inlet_p: Pressure,
    inlet_t: ThermodynamicTemperature,
    length: Length,
    diameter: Length,
    roughness: Length,
) -> Result<PipeResult> {
    let spec = &model_gen::PIPE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "inlet_t" => Some(inlet_t.value),
            "length" => Some(length.value),
            "diameter" => Some(diameter.value),
            "roughness" => Some(roughness.value),
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
    let outcome = kernel(&feed, length, diameter, roughness)?;

    Ok(PipeResult {
        outlet_n: outcome.outlet.n,
        outlet_z: outcome.outlet.z,
        outlet_p: outcome.outlet.p,
        outlet_t: outcome.outlet.t,
        outlet_h: joules_per_mole(outcome.outlet.h.value),
        pressure_drop: pascals(inlet_p.value - outcome.outlet.p.value),
        warnings,
    })
}
