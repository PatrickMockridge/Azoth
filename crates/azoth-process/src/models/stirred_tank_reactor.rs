//! `process.stirred_tank_reactor` - the reactor's kernel as a registered id.
//!
//! Spec: `specs/models/process/stirred_tank_reactor.toml`. The arithmetic is
//! [`crate::kernels::stirred_tank_reactor`]; what is here is the boundary a case and a
//! cross-impl test address.

use azoth_core::units::{MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{CalcResult, Result, Warning, apply_checks};

use crate::kernels::stirred_tank_reactor::{ReactorSetup, stirred_tank_reactor as kernel};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.stirred_tank_reactor`.
#[derive(Debug, Clone, PartialEq)]
pub struct StirredTankReactorResult {
    /// Product molar flow, mol/s.
    pub product_n: f64,
    /// Product composition.
    pub product_z: Vec<f64>,
    /// Product pressure.
    pub product_p: Pressure,
    /// Product temperature.
    pub product_t: ThermodynamicTemperature,
    /// Product molar enthalpy.
    pub product_h: MolarEnergy,
    /// The heat the vessel had to supply, W - nonzero only when it is isothermal.
    pub heat_duty: Power,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for StirredTankReactorResult {
    const CALC_ID: &'static str = "process.stirred_tank_reactor";
    const FIELDS: &'static [&'static str] = &[
        "product_n",
        "product_z",
        "product_p",
        "product_t",
        "product_h",
        "heat_duty",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// React a feed stoichiometrically and flash the product.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a reaction or a component the data does not
/// carry, a conversion outside `[0, 1]`, and whatever the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eleven
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

    Ok(StirredTankReactorResult {
        product_n: product.n,
        product_z: product.z,
        product_p: product.p,
        product_t: product.t,
        product_h: joules_per_mole(product.h.value),
        heat_duty,
        warnings,
    })
}
