//! `unit_ops.stirred_tank_reactor` - a stoichiometric conversion, then a flash.

use azoth_core::units::{
    Power, Pressure, ThermodynamicTemperature, joules_per_mole, pascals, watts,
};
use azoth_core::{AzothError, Result};
use azoth_eos::{ph_flash, pt_flash};
use azoth_reactions::databank::stoichiometry;

use crate::stream::Stream;

/// What a stirred-tank reactor is told: its reaction, its conversion, and how it is held.
///
/// **The palette entry declared none of this.** `unit_ops.stirred_tank_reactor` carries a feed
/// port and a product port and nothing between them, so the reaction the class applies, the
/// conversion it is applied at and the temperature the vessel is held at had no name in the
/// declaration - and a kernel with no parameters could only ever be a pass-through.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorSetup {
    /// The reaction's id in the reaction databank, which is what `addReaction` names.
    pub reaction: String,
    /// The limiting reactant: the component whose moles times `conversion` set the extent.
    pub limiting_reactant: String,
    /// The fractional conversion of the limiting reactant, in `[0, 1]`.
    pub conversion: f64,
    /// Whether the vessel holds its temperature, which is the class's `isothermal`.
    pub isothermal: bool,
    /// The temperature the vessel holds, which `run` reads only when isothermal.
    pub reactor_temperature: Option<ThermodynamicTemperature>,
    /// The pressure the vessel holds; absent means the feed's less `pressure_drop`.
    pub reactor_pressure: Option<Pressure>,
    /// The pressure drop the outlet takes where no reactor pressure is stated.
    pub pressure_drop: Pressure,
}

/// React a feed stoichiometrically, then flash it.
///
/// `run` clones the inlet, applies every reaction, sets the outlet pressure, holds the
/// temperature when it is isothermal, and flashes: a `TPflash` at the held temperature, or a
/// `PHflash` at the **inlet's own enthalpy** when it is adiabatic. There is no energy balance
/// over the reaction itself - the heat it releases is the difference between the two states,
/// which is what the flash answers.
///
/// **The reaction is a movement of moles.** `react` takes `moles(limiting) * conversion` and
/// adjusts every component by `that * coeff / |coeff(limiting)|`, so the extent is the
/// *limiting reactant's* and not the stoichiometry's. The per-second basis is the feed's: the
/// class's `getNumberOfmoles` is a rate here as everywhere in this library.
///
/// A stoichiometry row naming a substance the feed does not carry is skipped, which is the
/// class's own behaviour - `addComponent` throws for a component the fluid has no slot for and
/// the catch drops it - so a reaction whose products are not in the feed moves less material
/// than its stoichiometry says, and the port does not repair that.
///
/// # Errors
/// [`AzothError::InvalidInput`] if the reaction is not in the databank, if the limiting
/// reactant is not one of its rows or not a feed component, if the conversion is outside
/// `[0, 1]`, or if this conversion would take a component below zero. The flash's own refusals
/// pass through.
pub fn stirred_tank_reactor(feed: &Stream, setup: &ReactorSetup) -> Result<(Stream, Power)> {
    if !(0.0..=1.0).contains(&setup.conversion) {
        return Err(AzothError::invalid_input(
            "conversion",
            format!(
                "a conversion is a fraction of the limiting reactant, and {} is not in [0, 1]",
                setup.conversion
            ),
        ));
    }
    let rows = stoichiometry(&setup.reaction)?;
    if rows.is_empty() {
        return Err(AzothError::invalid_input(
            "reaction",
            format!(
                "`{}` is not a reaction the data carries, so there is no stoichiometry to apply",
                setup.reaction
            ),
        ));
    }
    let limiting = rows
        .iter()
        .find(|(component, _)| component.eq_ignore_ascii_case(&setup.limiting_reactant))
        .ok_or_else(|| {
            AzothError::invalid_input(
                "limiting_reactant",
                format!(
                    "`{}` is not a substance of reaction `{}`",
                    setup.limiting_reactant, setup.reaction
                ),
            )
        })?;
    let limiting_index = feed
        .components
        .iter()
        .position(|name| name.eq_ignore_ascii_case(&limiting.0))
        .ok_or_else(|| {
            AzothError::invalid_input(
                "components",
                format!(
                    "the limiting reactant `{}` is not one of the feed's substances",
                    limiting.0
                ),
            )
        })?;

    // `moles(limiting) * conversion`, then every coefficient scaled by the limiting one's
    // magnitude - `StoichiometricReaction.react`'s two lines.
    let limiting_moles = feed.n * feed.z[limiting_index];
    let reacted = limiting_moles * setup.conversion;
    let scale = reacted / limiting.1.abs();

    let mut amounts: Vec<f64> = feed.z.iter().map(|fraction| fraction * feed.n).collect();
    for (component, coefficient) in &rows {
        if let Some(index) = feed
            .components
            .iter()
            .position(|name| name.eq_ignore_ascii_case(component))
        {
            amounts[index] += scale * coefficient;
        }
    }
    if let Some(bad) = amounts.iter().position(|amount| *amount < 0.0) {
        return Err(AzothError::invalid_input(
            "conversion",
            format!(
                "this conversion takes {} below zero, so the reaction would consume more of it \
                 than the feed carries - the extent is the limiting reactant's, and the \
                 stoichiometry is what overdraws it",
                feed.components[bad]
            ),
        ));
    }
    let total: f64 = amounts.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "conversion",
            "the reaction leaves no moles in the reactor, so there is no state to flash",
        ));
    }
    let composition: Vec<f64> = amounts.iter().map(|amount| amount / total).collect();

    let p_out = match setup.reactor_pressure {
        Some(pressure) => pressure,
        None => pascals(feed.p.value - setup.pressure_drop.value),
    };
    if p_out.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "pressure_drop",
            format!(
                "a pressure drop of {} Pa takes a {} Pa feed to {} Pa, which is not a pressure",
                setup.pressure_drop.value, feed.p.value, p_out.value
            ),
        ));
    }

    let (mixture, ideal_gas) = feed.mixture()?;
    let temperature = if setup.isothermal {
        let held = setup.reactor_temperature.ok_or_else(|| {
            AzothError::invalid_input(
                "reactor_temperature",
                "an isothermal vessel holds a temperature, and none was stated",
            )
        })?;
        // The flash still runs, so the outlet's phases and its own enthalpy come from it -
        // which is what the class's `TPflash` does once it has set the temperature.
        let _ = pt_flash(&mixture, held, p_out, &composition)?;
        held.value
    } else {
        ph_flash::ph_flash(
            &mixture,
            &ideal_gas,
            p_out,
            joules_per_mole(feed.h.value),
            &composition,
        )?
        .temperature
        .value
    };

    let product = Stream::from_pt(
        feed.components.clone(),
        composition,
        total,
        p_out,
        azoth_core::units::kelvins(temperature),
    )?;

    // `heatDuty = outletEnthalpy - inletEnthalpy` when isothermal and not specified, and zero
    // otherwise: the class reports what it had to supply, and an adiabatic vessel supplies
    // nothing. Both are *totals* over one second of flow, which is a watt.
    let duty = if setup.isothermal {
        product.h.value * product.n - feed.h.value * feed.n
    } else {
        0.0
    };

    Ok((product, watts(duty)))
}
