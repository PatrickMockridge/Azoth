//! `unit_ops.separator` - a flash, split into vapour and liquid outlets.

use azoth_core::units::{Length, Power, Pressure, Velocity, joules_per_mole, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::hydrate_inhibitor_wt::PhaseLabel;
use azoth_eos::{Phase, RootSide, ph_flash, pt_flash};

use crate::segment::phase::{labelled_phases_at, mass_density, molar_mass_of};
use crate::stream::Stream;

/// `Separator.DEFAULT_LIQUID_DENSITY_FOR_SIZING`, which the sizing velocity substitutes where
/// the system carries no liquid.
const DEFAULT_LIQUID_DENSITY_FOR_SIZING: f64 = 1000.0;
/// `Separator.DEFAULT_LIQUID_DENSITY`, which the gas-load factor substitutes for a dry system
/// and for a second phase no denser than the gas.
const DEFAULT_LIQUID_DENSITY: f64 = 1000.0;
/// `Separator.MIN_LIQUID_GAS_DENSITY_DIFFERENCE`. **The load factor's rule and not the sizing
/// family's**: the two are deliberately separate.
const MIN_LIQUID_GAS_DENSITY_DIFFERENCE: f64 = 10.0;

/// `Separator`'s five steady capacity getters, on the system `separator` leaves behind.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Capacity {
    /// `getCapacityUtilization`.
    pub utilization: f64,
    /// `getMaxAllowableGasVelocity`; absent where the system carries no gas.
    pub max_velocity: Option<f64>,
    /// `getGasSuperficialVelocity`; absent as above.
    pub superficial_velocity: Option<f64>,
    /// `getGasLoadFactor()`; absent as above.
    pub load_factor: Option<f64>,
    /// `isOverloaded`: `utilization > 1.0`, exactly.
    pub overloaded: bool,
}

/// Flash a feed into vapour and liquid outlets.
///
/// The flash is at `P_in - pressure_drop` and at the **feed's own temperature**, which is
/// the vessel's: `Separator.run` reduces the pressure and runs a `TPflash`, so a pressure
/// drop is not a throttling and the outlets do not carry the inlet's enthalpy. A
/// `heat_input`, when given, moves the flash to the enthalpy it implies instead and the
/// outlet temperature is solved for.
///
/// `gas_in_liquid` carries that fraction of the vapour's moles into the liquid outlet,
/// which is NeqSim's `gasInLiquid` entrainment on its `feed`/`mole` basis. It moves
/// material between the phases and not energy: both outlets are then flashed at the same
/// temperature, so the balance across the transfer is not closed - which is what the class
/// does and is stated in the spec rather than corrected here. It is also the one thing
/// that re-runs the liquid outlet as a stream rather than leaving it as a phase, so on
/// that branch the liquid's enthalpy is the state its carried composition settles on.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the pressure drop would take the outlet to
/// a non-positive pressure or the entrainment fraction is outside `[0, 1]`, and whatever
/// the flash refuses.
pub fn separator(
    feed: &Stream,
    pressure_drop: Pressure,
    gas_in_liquid: f64,
    heat_input: Option<Power>,
) -> Result<(Stream, Stream)> {
    if !(0.0..=1.0).contains(&gas_in_liquid) {
        return Err(AzothError::invalid_input(
            "gas_in_liquid",
            format!("an entrainment fraction is a fraction, and {gas_in_liquid} is not in [0, 1]"),
        ));
    }
    let p_out = pascals(feed.p.value - pressure_drop.value);
    if p_out.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "pressure_drop",
            format!(
                "a pressure drop of {} Pa takes a {} Pa inlet to {} Pa, which is not a pressure",
                pressure_drop.value, feed.p.value, p_out.value
            ),
        ));
    }

    let (mixture, ideal_gas) = feed.mixture()?;

    let (temperature, x, y, beta) = match heat_input.filter(|power| power.value != 0.0) {
        Some(power) => {
            // One second of flow, for the reason `Heater.run` gives: the duty is W and the
            // enthalpy it moves is per mole, so it divides by the same flow a case states
            // in mol/s.
            let target = feed.h.value + power.value / feed.n;
            let flash = ph_flash::ph_flash(
                &mixture,
                &ideal_gas,
                p_out,
                joules_per_mole(target),
                &feed.z,
            )?;
            (
                flash.temperature,
                flash.x,
                flash.y,
                vapour_fraction(flash.phase, flash.vapour_fraction),
            )
        }
        None => {
            let flash = pt_flash(&mixture, feed.t, p_out, &feed.z)?;
            (
                feed.t,
                flash.x,
                flash.y,
                vapour_fraction(flash.phase, flash.vapour_fraction),
            )
        }
    };

    let n_vapour = feed.n * beta;
    let n_liquid = feed.n * (1.0 - beta);
    // NeqSim moves material between the phases only where both exist: `addPhaseFractionToPhase`
    // returns unchanged unless the from-phase and the to-phase are both present, so a
    // single-phase feed is untouched whatever the fraction says.
    let moved = if n_vapour > 0.0 && n_liquid > 0.0 {
        n_vapour * gas_in_liquid
    } else {
        0.0
    };

    let liquid_n = n_liquid + moved;
    let liquid_z: Vec<f64> = if moved > 0.0 {
        (0..x.len())
            .map(|i| (n_liquid * x[i] + moved * y[i]) / liquid_n)
            .collect()
    } else {
        x
    };

    // **The outlet is a phase, except where the class re-runs it as a stream, and the
    // difference between the two is not small.** `Separator.run` gives each outlet the
    // flash's own phase through `setThermoSystemFromPhase`, which keeps the phase's root -
    // and then calls `liquidOutStream.run(id)` under `gasInLiquid != 0.0`, where a stream
    // `run` is a `TPflash` of that composition at the outlet's own temperature. So the
    // vapour is always a phase root, and the liquid is one only while nothing was carried
    // into it: the capture's entrainment row reports `-13999.03` J/mol against the root's
    // `-14916.40`. The two fields that re-run the *vapour*, `oilInGas` and `waterInGas`,
    // are the three-phase entry's and are not carried here.
    let vapour = Stream::from_side(
        feed.components.clone(),
        y,
        n_vapour - moved,
        p_out,
        temperature,
        RootSide::Vapour,
    )?;
    let liquid = if gas_in_liquid != 0.0 {
        Stream::from_pt(
            feed.components.clone(),
            liquid_z,
            liquid_n,
            p_out,
            temperature,
        )?
    } else {
        Stream::from_side(
            feed.components.clone(),
            liquid_z,
            liquid_n,
            p_out,
            temperature,
            RootSide::Liquid,
        )?
    };
    Ok((vapour, liquid))
}

/// `Separator`'s steady capacity getters, on the system [`separator`] leaves behind.
///
/// **Two density rules, and they are not one.** The sizing family - the maximum velocity and
/// the utilisation over it - reads the liquid's own density and substitutes
/// [`DEFAULT_LIQUID_DENSITY_FOR_SIZING`] only where there is no liquid, with **no floor**.
/// The gas-load factor applies [`MIN_LIQUID_GAS_DENSITY_DIFFERENCE`] to the same density
/// first. A near-dry state makes the two disagree by construction, and neither is corrected
/// to agree.
///
/// The metric reads the system `run` left: `gas_outlet`'s temperature and pressure, the
/// vapour that is *left* after entrainment, and the liquid that absorbed it - whose density
/// `initPhysicalProperties` recomputes, which is why `liquid_outlet`'s composition is
/// re-flashed rather than the feed's `x` read.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] where the gas density is not positive, the liquid
/// is no denser than the gas (before or after the near-dry floor), or the gas area admits no
/// flow. NeqSim answers `NaN` from its getters in those states; this refuses the state.
pub fn capacity_states(
    feed: &Stream,
    gas_outlet: &Stream,
    liquid_outlet: &Stream,
    internal_diameter: Length,
    design_load_factor: Velocity,
    horizontal: bool,
    level_fraction: f64,
) -> Result<Capacity> {
    let (mixture, _ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, gas_outlet.t, gas_outlet.p, &feed.z)?;
    let labelled = labelled_phases_at(&mixture, gas_outlet.t, gas_outlet.p, &feed.z, &flash)?;

    let Some((_label, gas_z, gas_root)) = labelled
        .iter()
        .find(|(kind, _, _)| *kind == PhaseLabel::Gas)
    else {
        // The class's own first answer: a stream with no gas phase uses none of its capacity.
        return Ok(Capacity {
            utilization: 0.0,
            max_velocity: None,
            superficial_velocity: None,
            load_factor: None,
            overloaded: false,
        });
    };
    let gas_molar_mass = molar_mass_of(&mixture, gas_z)?;
    let gas_density = mass_density(
        &mixture,
        gas_outlet.t,
        gas_outlet.p,
        gas_z,
        *gas_root,
        gas_molar_mass,
    )?;
    if !gas_density.is_finite() || gas_density <= 0.0 {
        return Err(AzothError::invalid_input(
            "capacity_utilization",
            format!(
                "the gas density is {gas_density} kg/m3, and the Souders-Brown velocity divides \
                 by its square root"
            ),
        ));
    }

    let liquid_density = if labelled.iter().any(|(kind, _, _)| *kind != PhaseLabel::Gas) {
        // The root is the *entrained* composition's own flash: `initPhysicalProperties`
        // recomputes the density of the liquid that absorbed the carried moles.
        let entrained = pt_flash(&mixture, liquid_outlet.t, liquid_outlet.p, &liquid_outlet.z)?;
        mass_density(
            &mixture,
            liquid_outlet.t,
            liquid_outlet.p,
            &liquid_outlet.z,
            entrained.z_liquid,
            molar_mass_of(&mixture, &liquid_outlet.z)?,
        )?
    } else {
        DEFAULT_LIQUID_DENSITY_FOR_SIZING
    };

    let diameter = internal_diameter.value;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    let gas_area = if horizontal {
        area * (1.0 - level_fraction)
    } else {
        area
    };
    let volumetric_flow = gas_outlet.n * gas_molar_mass / gas_density;

    // The sizing rule's own positivity: a liquid denser than the gas is all it needs, and the
    // class's `getCapacityUtilization` answers `NaN` where this refuses.
    if !liquid_density.is_finite() || liquid_density <= gas_density {
        return Err(AzothError::invalid_input(
            "capacity_utilization",
            format!(
                "the liquid phase is {liquid_density} kg/m3 against the gas's {gas_density}, so \
                 `rho_l - rho_g` is not positive and the Souders-Brown velocity is the square \
                 root of a negative number. NeqSim answers NaN there; this refuses the state"
            ),
        ));
    }

    // The two refusals above are what make the radicand positive, and the diameter's own bound
    // what keeps the area off zero. guarded: separator_capacity.radicand_positive
    let max_velocity =
        design_load_factor.value * ((liquid_density - gas_density) / gas_density).sqrt();
    let max_flow = max_velocity * gas_area;
    if !max_flow.is_finite() || max_flow <= 0.0 {
        return Err(AzothError::invalid_input(
            "capacity_utilization",
            format!(
                "the maximum allowable gas flow is {max_flow} m3/s, and utilisation divides by it"
            ),
        ));
    }
    let utilization = volumetric_flow / max_flow;
    let superficial_velocity = volumetric_flow / gas_area;

    // **The load factor's own rule, applied after the sizing one.** The floor substitutes
    // `1000.0` for a second phase no denser than the gas, and the positivity is re-checked
    // against the substituted value.
    let floored_liquid_density = load_factor_liquid_density(liquid_density, gas_density);
    if !floored_liquid_density.is_finite() || floored_liquid_density <= gas_density {
        return Err(AzothError::invalid_input(
            "gas_load_factor",
            format!(
                "the liquid phase is {floored_liquid_density} kg/m3 after the near-dry floor \
                 against the gas's {gas_density}, and the load factor divides by their difference"
            ),
        ));
    }
    let term = (floored_liquid_density - gas_density) / gas_density;
    let load_factor = superficial_velocity * (1.0 / term).sqrt();
    if !load_factor.is_finite() {
        return Err(AzothError::invalid_input(
            "gas_load_factor",
            format!("the gas load factor is {load_factor} m/s"),
        ));
    }

    Ok(Capacity {
        utilization,
        max_velocity: Some(max_velocity),
        superficial_velocity: Some(superficial_velocity),
        load_factor: Some(load_factor),
        overloaded: utilization > 1.0,
    })
}

/// The gas-load factor's near-dry floor, kept apart from the sizing rule.
///
/// `getGasLoadFactor()` substitutes `1000.0` whenever `rho_l - rho_g` is under
/// [`MIN_LIQUID_GAS_DENSITY_DIFFERENCE`]; the sizing family does **not**. A state two
/// kilograms apart therefore has two liquid densities by design, and the pair is the reason
/// this is a function rather than an inline expression.
fn load_factor_liquid_density(liquid_density: f64, gas_density: f64) -> f64 {
    if liquid_density - gas_density < MIN_LIQUID_GAS_DENSITY_DIFFERENCE {
        DEFAULT_LIQUID_DENSITY
    } else {
        liquid_density
    }
}

/// The vapour fraction a flash's phase and `beta` imply.
///
/// `beta` is absent for a single-phase answer - the flash reports `None` rather than a
/// number outside `[0, 1]`, because there is no split to report - so the phase is what
/// says which side of the interval the outlet is on.
fn vapour_fraction(phase: Phase, beta: Option<f64>) -> f64 {
    match phase {
        Phase::TwoPhase => beta.unwrap_or(0.0),
        Phase::AllVapour => 1.0,
        Phase::AllLiquid | Phase::Trivial => 0.0,
    }
}

#[cfg(test)]
mod tests {
    use super::{DEFAULT_LIQUID_DENSITY, load_factor_liquid_density};

    /// **No captured state reaches the floor** - the pin's near-dry scrubber test is
    /// single-phase in this jar - so the branch is pinned on constructed densities.
    #[test]
    fn the_near_dry_floor_substitutes_and_exactly_ten_does_not() {
        assert_eq!(
            load_factor_liquid_density(1_009.0, 1_000.0),
            DEFAULT_LIQUID_DENSITY
        );
        assert_eq!(load_factor_liquid_density(1_010.0, 1_000.0), 1_010.0);
        assert_eq!(load_factor_liquid_density(542.32, 20.16), 542.32);
    }
}
