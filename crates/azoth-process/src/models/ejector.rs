//! `process.ejector` - the ejector's kernel as a registered id.
//!
//! Spec: `specs/models/process/ejector.toml`. The arithmetic is [`crate::kernels::ejector`];
//! what is here is the boundary a case and a cross-impl test address.
//!
//! **Two inlets carrying two fluids, so ten fields of inlet and five of outlet.** A
//! two-inlet machine's inlets are not one port: the motive and the suction may be different
//! mixtures, and `_dispatch.fluid_inputs` derives each `components` prefix from the input's
//! own name - the arrangement `process.heat_exchanger` introduced.

use azoth_core::units::{
    Area, Length, MolarEnergy, Pressure, ThermodynamicTemperature, Velocity, Volume, cubic_meters,
    joules_per_mole, meters, square_meters,
};
use azoth_core::{Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalar, warnings as wire_warnings};

use crate::kernels::ejector::{EjectorOutcome, EjectorSetup, ejector as kernel};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.ejector`.
///
/// **Serialisable, so that a flowsheet's run publishes it.** The mixing pressure and the four
/// velocities are on no outlet stream and are not inputs either - the class estimates the first
/// and derives the rest from the efficiencies - so without them a reader sees a discharge state
/// with no way to tell an ejector that barely drew from one at its limit.
///
/// The five are `Option`s for one degenerate case and not for a mode: a machine handed no flow
/// returns its motive stream unchanged, and a run on no fluid reached no pressure and no velocity.
///
/// The mechanical design that follows them is `updateDesign`'s container, and it is **plain
/// numbers**: the class calls `resetDesign` on that same branch, and its getters then answer
/// zero rather than nothing.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct EjectorResult {
    /// Outlet molar flow, mol/s.
    pub outlet_n: f64,
    /// Outlet composition.
    pub outlet_z: Vec<f64>,
    /// Outlet pressure.
    #[serde(serialize_with = "scalar")]
    pub outlet_p: Pressure,
    /// Outlet temperature.
    #[serde(serialize_with = "scalar")]
    pub outlet_t: ThermodynamicTemperature,
    /// Outlet molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub outlet_h: MolarEnergy,
    /// The pressure the two streams met at, Pa: the class's own estimate, not an input.
    #[serde(serialize_with = "optional_scalar")]
    pub mixing_pressure: Option<Pressure>,
    /// The motive nozzle's exit velocity, m/s.
    #[serde(serialize_with = "optional_scalar")]
    pub motive_nozzle_velocity: Option<Velocity>,
    /// The suction nozzle's, m/s.
    #[serde(serialize_with = "optional_scalar")]
    pub suction_nozzle_velocity: Option<Velocity>,
    /// The mixed stream's, m/s.
    #[serde(serialize_with = "optional_scalar")]
    pub mixing_velocity: Option<Velocity>,
    /// The diffuser's design velocity, m/s.
    #[serde(serialize_with = "optional_scalar")]
    pub diffuser_velocity: Option<Velocity>,
    /// The motive nozzle's throat area, m2.
    #[serde(serialize_with = "scalar")]
    pub motive_nozzle_throat_area: Area,
    /// Its throat diameter, m.
    #[serde(serialize_with = "scalar")]
    pub motive_nozzle_diameter: Length,
    /// The suction inlet's area, m2.
    #[serde(serialize_with = "scalar")]
    pub suction_inlet_area: Area,
    /// Its diameter, m.
    #[serde(serialize_with = "scalar")]
    pub suction_inlet_diameter: Length,
    /// The mixing chamber's area, m2.
    #[serde(serialize_with = "scalar")]
    pub mixing_chamber_area: Area,
    /// Its diameter, m.
    #[serde(serialize_with = "scalar")]
    pub mixing_chamber_diameter: Length,
    /// The diffuser outlet's area, m2.
    #[serde(serialize_with = "scalar")]
    pub diffuser_outlet_area: Area,
    /// Its diameter, m.
    #[serde(serialize_with = "scalar")]
    pub diffuser_outlet_diameter: Length,
    /// The motive nozzle's effective length, m.
    #[serde(serialize_with = "scalar")]
    pub motive_nozzle_effective_length: Length,
    /// The suction inlet's length, m.
    #[serde(serialize_with = "scalar")]
    pub suction_inlet_length: Length,
    /// The mixing chamber's length, m.
    #[serde(serialize_with = "scalar")]
    pub mixing_chamber_length: Length,
    /// The diffuser outlet's length, m.
    #[serde(serialize_with = "scalar")]
    pub diffuser_outlet_length: Length,
    /// The four stages' cylinder volumes, m3.
    #[serde(serialize_with = "scalar")]
    pub body_volume: Volume,
    /// The two connection pipes' cylinder volumes, m3.
    #[serde(serialize_with = "scalar")]
    pub connected_piping_volume: Volume,
    /// The two volumes' sum, m3.
    #[serde(serialize_with = "scalar")]
    pub total_volume: Volume,
    /// The suction connection's straight length, m.
    #[serde(serialize_with = "scalar")]
    pub suction_connection_length: Length,
    /// The discharge connection's straight length, m.
    #[serde(serialize_with = "scalar")]
    pub discharge_connection_length: Length,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl EjectorResult {
    /// The result of one kernel call, as a flowsheet's dispatcher reads it.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case runs
    /// the spec's input checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    #[must_use]
    pub fn of(outcome: &EjectorOutcome, warnings: Vec<Warning>) -> Self {
        let numbers = outcome.numbers.as_ref();
        Self {
            outlet_n: outcome.outlet.n,
            outlet_z: outcome.outlet.z.clone(),
            outlet_p: outcome.outlet.p,
            outlet_t: outcome.outlet.t,
            outlet_h: joules_per_mole(outcome.outlet.h.value),
            mixing_pressure: numbers.map(|n| n.mixing_pressure),
            motive_nozzle_velocity: numbers.map(|n| n.motive_nozzle_velocity),
            suction_nozzle_velocity: numbers.map(|n| n.suction_nozzle_velocity),
            mixing_velocity: numbers.map(|n| n.mixing_velocity),
            diffuser_velocity: numbers.map(|n| n.diffuser_velocity),
            motive_nozzle_throat_area: square_meters(outcome.design.motive_nozzle_throat_area),
            motive_nozzle_diameter: meters(outcome.design.motive_nozzle_diameter),
            suction_inlet_area: square_meters(outcome.design.suction_inlet_area),
            suction_inlet_diameter: meters(outcome.design.suction_inlet_diameter),
            mixing_chamber_area: square_meters(outcome.design.mixing_chamber_area),
            mixing_chamber_diameter: meters(outcome.design.mixing_chamber_diameter),
            diffuser_outlet_area: square_meters(outcome.design.diffuser_outlet_area),
            diffuser_outlet_diameter: meters(outcome.design.diffuser_outlet_diameter),
            motive_nozzle_effective_length: meters(outcome.design.motive_nozzle_effective_length),
            suction_inlet_length: meters(outcome.design.suction_inlet_length),
            mixing_chamber_length: meters(outcome.design.mixing_chamber_length),
            diffuser_outlet_length: meters(outcome.design.diffuser_outlet_length),
            body_volume: cubic_meters(outcome.design.body_volume),
            connected_piping_volume: cubic_meters(outcome.design.connected_piping_volume),
            total_volume: cubic_meters(outcome.design.total_volume),
            suction_connection_length: meters(outcome.design.suction_connection_length),
            discharge_connection_length: meters(outcome.design.discharge_connection_length),
            warnings,
        }
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
    let outcome = kernel(
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

    Ok(EjectorResult::of(&outcome, warnings))
}
