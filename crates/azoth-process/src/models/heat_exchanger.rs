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
    MolarEnergy, Power, Pressure, ThermalConductance, ThermodynamicTemperature, joules_per_mole,
};
use azoth_core::{CalcResult, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};
use crate::kernels::heat_exchanger as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.heat_exchanger`.
///
/// Two single-multiplicity ports, so the record's five fields per port - the same shape
/// `process.separator` has, under this unit operation's own port names - and the duty and the
/// rating's own numbers beside them.
///
/// **Serialisable, so that a flowsheet's run publishes it.** The duty is on no outlet stream, and
/// neither are the four numbers a rating sized the exchanger by; the wire writes each through
/// [`crate::executor::json::scalar`] and writes `null` for the rating's where there was no rating
/// to size one, which is the same absence a flash's single phase is.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct HeatExchangerResult {
    /// Hot outlet molar flow, mol/s.
    pub hot_out_n: f64,
    /// Hot outlet composition.
    pub hot_out_z: Vec<f64>,
    #[serde(serialize_with = "scalar")]
    /// Hot outlet pressure.
    pub hot_out_p: Pressure,
    #[serde(serialize_with = "scalar")]
    /// Hot outlet temperature.
    pub hot_out_t: ThermodynamicTemperature,
    #[serde(serialize_with = "scalar")]
    /// Hot outlet molar enthalpy.
    pub hot_out_h: MolarEnergy,
    /// Cold outlet molar flow, mol/s.
    pub cold_out_n: f64,
    /// Cold outlet composition.
    pub cold_out_z: Vec<f64>,
    #[serde(serialize_with = "scalar")]
    /// Cold outlet pressure.
    pub cold_out_p: Pressure,
    #[serde(serialize_with = "scalar")]
    /// Cold outlet temperature.
    pub cold_out_t: ThermodynamicTemperature,
    #[serde(serialize_with = "scalar")]
    /// Cold outlet molar enthalpy.
    pub cold_out_h: MolarEnergy,
    /// The heat the hot side released, W: positive when the hot side cools, which is what the
    /// cold side's gain balances. The one figure that means the same thing in both modes.
    #[serde(serialize_with = "scalar")]
    pub duty: Power,
    /// `UA / C_min`, dimensionless, or `null` where one outlet was pinned instead of rated.
    pub ntu: Option<f64>,
    /// The swing the relation scaled, or `null` where one outlet was pinned.
    pub effectiveness: Option<f64>,
    /// The smaller of the two estimated capacities, W/K, or `null` where one outlet was pinned.
    #[serde(serialize_with = "optional_scalar")]
    pub c_min: Option<ThermalConductance>,
    /// The larger of them, W/K, or `null` where one outlet was pinned.
    #[serde(serialize_with = "optional_scalar")]
    pub c_max: Option<ThermalConductance>,
    /// `C_min / C_max`, dimensionless, or `null` where one outlet was pinned.
    pub capacity_ratio: Option<f64>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
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
        "duty",
        "ntu",
        "effectiveness",
        "c_min",
        "c_max",
        "capacity_ratio",
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
    let outcome = kernel(
        &hot,
        &cold,
        ua,
        flow_arrangement,
        hot_outlet_temperature,
        cold_outlet_temperature,
    )?;
    Ok(HeatExchangerResult::of(&outcome, warnings))
}

impl HeatExchangerResult {
    /// The result of one kernel call, as a flowsheet's dispatcher reads it.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case runs
    /// the spec's input checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    /// So an empty list here means "this kernel raised none" and not "nothing checked this call" -
    /// and this kernel raises none of its own, which is why the dispatcher hands it `Vec::new()`.
    #[must_use]
    pub fn of(outcome: &kernel::HeatExchangerOutcome, warnings: Vec<Warning>) -> Self {
        Self {
            hot_out_n: outcome.hot_out.n,
            hot_out_z: outcome.hot_out.z.clone(),
            hot_out_p: outcome.hot_out.p,
            hot_out_t: outcome.hot_out.t,
            hot_out_h: joules_per_mole(outcome.hot_out.h.value),
            cold_out_n: outcome.cold_out.n,
            cold_out_z: outcome.cold_out.z.clone(),
            cold_out_p: outcome.cold_out.p,
            cold_out_t: outcome.cold_out.t,
            cold_out_h: joules_per_mole(outcome.cold_out.h.value),
            duty: outcome.duty,
            ntu: outcome.rating.as_ref().map(|rating| rating.ntu),
            effectiveness: outcome.rating.as_ref().map(|rating| rating.effectiveness),
            c_min: outcome.rating.as_ref().map(|rating| rating.c_min),
            c_max: outcome.rating.as_ref().map(|rating| rating.c_max),
            capacity_ratio: outcome.rating.as_ref().map(|rating| rating.capacity_ratio),
            warnings,
        }
    }
}
