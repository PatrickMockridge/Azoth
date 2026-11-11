//! `process.stirred_tank_reactor` - the reactor's kernel as a registered id.
//!
//! Spec: `specs/models/process/stirred_tank_reactor.toml`. The arithmetic is
//! [`crate::kernels::stirred_tank_reactor`]; what is here is the boundary a case and a
//! cross-impl test address.

use azoth_core::units::{
    MolarEnergy, Power, Pressure, ThermodynamicTemperature, Time, Volume, VolumetricPowerDensity,
    cubic_meters, hours, joules_per_mole, watts_per_cubic_meter,
};
use azoth_core::{Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::stirred_tank_reactor::{
    DEFAULT_AGITATOR_POWER_PER_VOLUME_W_PER_M3, DEFAULT_RESIDENCE_TIME_H, DEFAULT_VESSEL_VOLUME_M3,
    ReactorSetup, agitator_power, stirred_tank_reactor as kernel,
};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.stirred_tank_reactor`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct StirredTankReactorResult {
    /// Product molar flow, mol/s.
    pub product_n: f64,
    /// Product composition.
    pub product_z: Vec<f64>,
    /// Product pressure.
    #[serde(serialize_with = "scalar")]
    pub product_p: Pressure,
    /// Product temperature.
    #[serde(serialize_with = "scalar")]
    pub product_t: ThermodynamicTemperature,
    /// Product molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub product_h: MolarEnergy,
    /// The heat the vessel had to supply, W - nonzero only when it is isothermal.
    #[serde(serialize_with = "scalar")]
    pub heat_duty: Power,
    /// The agitator's power, kW: the stored specific power times the stored vessel volume.
    #[serde(serialize_with = "scalar")]
    pub agitator_power: Power,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl StirredTankReactorResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case
    /// runs the spec's checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    #[must_use]
    pub fn of(
        product: &Stream,
        duty: Power,
        agitator_power: Power,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            product_n: product.n,
            product_z: product.z.clone(),
            product_p: product.p,
            product_t: product.t,
            product_h: joules_per_mole(product.h.value),
            heat_duty: duty,
            agitator_power,
            warnings,
        }
    }
}

/// React a feed stoichiometrically and flash the product.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a reaction or a component the data does not
/// carry, a conversion outside `[0, 1]`, and whatever the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are fourteen
pub fn stirred_tank_reactor(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    reaction: &str,
    limiting_reactant: &str,
    conversion: f64,
    isothermal: bool,
    reactor_temperature: Option<ThermodynamicTemperature>,
    reactor_pressure: Option<Pressure>,
    pressure_drop: Pressure,
    vessel_volume: Option<Volume>,
    residence_time: Option<Time>,
    agitator_power_per_volume: Option<VolumetricPowerDensity>,
) -> Result<StirredTankReactorResult> {
    let spec = &model_gen::STIRRED_TANK_REACTOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "conversion" => Some(conversion),
            "pressure_drop" => Some(pressure_drop.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let (product, heat_duty) = kernel(
        &feed,
        &ReactorSetup {
            reaction: reaction.to_string(),
            limiting_reactant: limiting_reactant.to_string(),
            conversion,
            isothermal,
            reactor_temperature,
            reactor_pressure,
            pressure_drop,
        },
    )?;

    // The class's three stored fields, at its own initialisers when absent. **`run` reads
    // none of them**: the volume and the specific power are `toMap`'s echo and the agitator's
    // product, and the residence time is carried for a transient this port does not have.
    let vessel_volume = vessel_volume.unwrap_or_else(|| cubic_meters(DEFAULT_VESSEL_VOLUME_M3));
    let _residence_time = residence_time.unwrap_or_else(|| hours(DEFAULT_RESIDENCE_TIME_H));
    let per_volume = agitator_power_per_volume
        .unwrap_or_else(|| watts_per_cubic_meter(DEFAULT_AGITATOR_POWER_PER_VOLUME_W_PER_M3));

    Ok(StirredTankReactorResult::of(
        &product,
        heat_duty,
        agitator_power(per_volume, vessel_volume),
        warnings,
    ))
}
