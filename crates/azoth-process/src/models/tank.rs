//! `process.tank` - the tank's kernel as a registered id.
//!
//! Spec: `specs/models/process/tank.toml`. The arithmetic is [`crate::kernels::tank`]; what
//! is here is the boundary a case and a cross-impl test address.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::tank as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.tank`.
#[derive(Debug, Clone, PartialEq)]
pub struct TankResult {
    /// Gas outlet molar flow, mol/s.
    pub gas_n: f64,
    /// Gas outlet composition.
    pub gas_z: Vec<f64>,
    /// Gas outlet pressure.
    pub gas_p: Pressure,
    /// Gas outlet temperature.
    pub gas_t: ThermodynamicTemperature,
    /// Gas outlet molar enthalpy.
    pub gas_h: MolarEnergy,
    /// Liquid outlet molar flow, mol/s.
    pub liquid_n: f64,
    /// Liquid outlet composition.
    pub liquid_z: Vec<f64>,
    /// Liquid outlet pressure.
    pub liquid_p: Pressure,
    /// Liquid outlet temperature.
    pub liquid_t: ThermodynamicTemperature,
    /// Liquid outlet molar enthalpy.
    pub liquid_h: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TankResult {
    const CALC_ID: &'static str = "process.tank";
    const FIELDS: &'static [&'static str] = &[
        "gas_n", "gas_z", "gas_p", "gas_t", "gas_h", "liquid_n", "liquid_z", "liquid_p",
        "liquid_t", "liquid_h", "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Join a tank's inlets and split the result into a gas and a liquid outlet.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the feeds' shapes disagree or carry no flow
/// in total, and whatever the databank or the flash refuses.
pub fn tank(
    components: &[String],
    feed_n: &[f64],
    feed_z: &[Vec<f64>],
    feed_p: &[Pressure],
    feed_t: &[ThermodynamicTemperature],
) -> Result<TankResult> {
    // **No parameters, so the one check it has is on a feed.** The entry declares nothing
    // between the feeds and the outlets - the design volume is the one thing a caller might
    // reach for and it is not a steady-state input - so the bound the family puts on a
    // temperature is the whole of it, and the feeds are one port, so it resolves the first.
    // The mixer's zero-total refusal and the flash's own are the rest.
    let mut warnings = Vec::new();
    apply_checks(
        model_gen::TANK_SPEC.input_checks(),
        |quantity| match quantity {
            "feed_t" => feed_t.first().map(|t| t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let feeds: Vec<Stream> = feed_n
        .iter()
        .zip(feed_z)
        .zip(feed_p)
        .zip(feed_t)
        .map(|(((n, z), p), t)| Stream::from_pt(components.to_vec(), z.clone(), *n, *p, *t))
        .collect::<Result<_>>()?;

    let (gas, liquid) = kernel(&feeds)?;

    Ok(TankResult {
        gas_n: gas.n,
        gas_z: gas.z,
        gas_p: gas.p,
        gas_t: gas.t,
        gas_h: joules_per_mole(gas.h.value),
        liquid_n: liquid.n,
        liquid_z: liquid.z,
        liquid_p: liquid.p,
        liquid_t: liquid.t,
        liquid_h: joules_per_mole(liquid.h.value),
        warnings,
    })
}
