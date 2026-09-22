//! `process.heat_exchanger` - the exchanger's kernel as a registered id.
//!
//! Spec: `specs/models/process/heat_exchanger.toml`. The arithmetic is
//! [`crate::kernels::heat_exchanger`]; what is here is the boundary a case and a
//! cross-impl test address.
//!
//! **The first id whose two ports carry different fluids.** Every other unit operation
//! declares one `components`; this declares `hot_components` and `cold_components`, and
//! the boundary derives the two `Stream`s from them separately.

use azoth_core::units::{
    MolarEnergy, Pressure, ThermalConductance, ThermodynamicTemperature, joules_per_mole,
};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::heat_exchanger as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.heat_exchanger`.
///
/// Two single-multiplicity ports, so the record's five fields per port - the same shape
/// `process.separator` has, under this unit operation's own port names.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatExchangerResult {
    /// Hot outlet molar flow, mol/s.
    pub hot_out_n: f64,
    /// Hot outlet composition.
    pub hot_out_z: Vec<f64>,
    /// Hot outlet pressure.
    pub hot_out_p: Pressure,
    /// Hot outlet temperature.
    pub hot_out_t: ThermodynamicTemperature,
    /// Hot outlet molar enthalpy.
    pub hot_out_h: MolarEnergy,
    /// Cold outlet molar flow, mol/s.
    pub cold_out_n: f64,
    /// Cold outlet composition.
    pub cold_out_z: Vec<f64>,
    /// Cold outlet pressure.
    pub cold_out_p: Pressure,
    /// Cold outlet temperature.
    pub cold_out_t: ThermodynamicTemperature,
    /// Cold outlet molar enthalpy.
    pub cold_out_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for HeatExchangerResult {
    const CALC_ID: &'static str = "process.heat_exchanger";
    const FIELDS: &'static [&'static str] = &[
        "hot_out_n",
        "hot_out_z",
        "hot_out_p",
        "hot_out_t",
        "hot_out_h",
        "cold_out_n",
        "cold_out_z",
        "cold_out_p",
        "cold_out_t",
        "cold_out_h",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Exchange heat between two streams.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for the mode and arrangement refusals
/// [`crate::kernels::heat_exchanger`] lists, and whatever the databank refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are fourteen
pub fn heat_exchanger(
    hot_components: &[String],
    hot_in_n: f64,
    hot_in_z: &[f64],
    hot_in_p: Pressure,
    hot_in_t: ThermodynamicTemperature,
    cold_components: &[String],
    cold_in_n: f64,
    cold_in_z: &[f64],
    cold_in_p: Pressure,
    cold_in_t: ThermodynamicTemperature,
    ua: Option<ThermalConductance>,
    flow_arrangement: &str,
    hot_outlet_temperature: Option<ThermodynamicTemperature>,
    cold_outlet_temperature: Option<ThermodynamicTemperature>,
) -> Result<HeatExchangerResult> {
    let spec = &model_gen::HEAT_EXCHANGER_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. `ua` is optional, so an absent one is a *skipped check* rather than a pass -
    // which is the whole reason the bound is declared on it and not on the temperatures.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "ua" => ua.map(|value| value.value),
            "hot_in_t" => Some(hot_in_t.value),
            "cold_in_t" => Some(cold_in_t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let hot = Stream::from_pt(
        hot_components.to_vec(),
        hot_in_z.to_vec(),
        hot_in_n,
        hot_in_p,
        hot_in_t,
    )?;
    let cold = Stream::from_pt(
        cold_components.to_vec(),
        cold_in_z.to_vec(),
        cold_in_n,
        cold_in_p,
        cold_in_t,
    )?;
    let (hot_out, cold_out) = kernel(
        &hot,
        &cold,
        ua,
        flow_arrangement,
        hot_outlet_temperature,
        cold_outlet_temperature,
    )?;

    Ok(HeatExchangerResult {
        hot_out_n: hot_out.n,
        hot_out_z: hot_out.z,
        hot_out_p: hot_out.p,
        hot_out_t: hot_out.t,
        hot_out_h: joules_per_mole(hot_out.h.value),
        cold_out_n: cold_out.n,
        cold_out_z: cold_out.z,
        cold_out_p: cold_out.p,
        cold_out_t: cold_out.t,
        cold_out_h: joules_per_mole(cold_out.h.value),
        warnings,
    })
}
