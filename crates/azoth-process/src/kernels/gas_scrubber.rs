//! `unit_ops.gas_scrubber` - `unit_ops.separator`'s arithmetic, under the other entry.

use azoth_core::units::{Length, Power, Pressure, Velocity};
use azoth_core::{AzothError, Result};
use azoth_eos::hydrate_inhibitor_wt::PhaseLabel;
use azoth_eos::pt_flash;

use crate::kernels::separator::separator;
use crate::segment::phase::{labelled_phases_at, mass_density, molar_mass_of};
use crate::stream::Stream;

/// `Separator.getMaxAllowableGasVelocity`'s own substitution where the system carries no liquid:
/// a bare `1000.0`, beside the `50.0` it puts on an absent gas phase. A stream with no gas never
/// reaches it, because `getCapacityUtilization` answers zero first.
const DEFAULT_LIQUID_DENSITY: f64 = 1000.0;

/// Flash a feed into vapour and liquid outlets, as a scrubber does.
///
/// **This is a delegation and the delegation is the port.** `GasScrubber extends Separator`
/// and **does not override `run`**: its own 112 lines are constructors, a mechanical design
/// and the capacity metric, and the steady state is the separator's flash exactly. So the
/// arithmetic is [`separator`]'s, and writing it out again would be a second implementation
/// of one vessel.
///
/// # Errors
/// The separator's own: a pressure drop that leaves a non-positive pressure, an entrainment
/// fraction outside `[0, 1]`, and whatever the flash refuses.
pub fn gas_scrubber(
    feed: &Stream,
    pressure_drop: Pressure,
    gas_in_liquid: f64,
    heat_input: Option<Power>,
) -> Result<(Stream, Stream)> {
    separator(feed, pressure_drop, gas_in_liquid, heat_input)
}

/// `Separator.getCapacityUtilization`, on the system [`gas_scrubber`] left behind.
///
/// `Q / (K sqrt((rho_l - rho_g)/rho_g) A)` with `A = pi d^2 / 4`. **Three details are the
/// class's and not the columns'.** `Q` is the vapour that is *left*: `Separator.run` applies the
/// entrainment to the system and only then does the getter read it, so the dry `rho_l` on the
/// capture is `542.3229333258249` and the entrained one `485.8516787302`. The liquid density is
/// the `oil` phase's, `aqueous`'s where there is no oil, and [`DEFAULT_LIQUID_DENSITY`] where
/// there is neither - no `10.0` floor, which is the columns' rule. And a system with no gas
/// phase answers `0.0`.
///
/// **The liquid's *type* is the flash's and its *density* is the entrained composition's.**
/// `hasPhaseType` reads a label `addPhaseFractionToPhase` never revisits, while
/// `initPhysicalProperties` does recompute the density of the oil that absorbed the carried
/// moles. `feed` is therefore the flashed feed's composition and `liquid_outlet` the state the
/// transfer left, and the two agree on a run with nothing carried.
///
/// # Errors
/// [`AzothError::InvalidInput`] where the liquid is no denser than the gas. NeqSim's
/// `getCapacityUtilization` answers `NaN` there; this refuses the state rather than answering it.
pub fn capacity_utilization(
    feed: &Stream,
    gas_outlet: &Stream,
    liquid_outlet: &Stream,
    internal_diameter: Length,
    design_gas_load_factor: Velocity,
) -> Result<f64> {
    let (mixture, _ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, gas_outlet.t, gas_outlet.p, &feed.z)?;
    let labelled = labelled_phases_at(&mixture, gas_outlet.t, gas_outlet.p, &feed.z, &flash)?;

    let Some((_label, gas_z, gas_root)) = labelled
        .iter()
        .find(|(kind, _, _)| *kind == PhaseLabel::Gas)
    else {
        // The class's own first answer: a stream with no gas phase uses none of its capacity.
        return Ok(0.0);
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

    let liquid_density = if labelled.iter().any(|(kind, _, _)| *kind != PhaseLabel::Gas) {
        // The root is the *entrained* composition's own, from its own flash: the oil that has
        // absorbed the carried moles sits on a different one.
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
        DEFAULT_LIQUID_DENSITY
    };

    // A `NaN` density is the same state and is named rather than left to the comparison's
    // negation, which reads as a mistyped `<` to anyone who has not thought about `NaN`.
    if liquid_density.is_nan() || gas_density.is_nan() || liquid_density <= gas_density {
        return Err(AzothError::invalid_input(
            "capacity_utilization",
            format!(
                "the liquid phase is {liquid_density} kg/m3 against the gas's {gas_density}, so \
                 `rho_l - rho_g` is not positive and the Souders-Brown velocity is the square \
                 root of a negative number. NeqSim's `getCapacityUtilization` answers NaN \
                 there; this refuses the state rather than answering it"
            ),
        ));
    }

    let diameter = internal_diameter.value;
    let area = std::f64::consts::PI * diameter * diameter / 4.0;
    // The refusal above is what makes the radicand positive, and the diameter's own bound what
    // keeps the area off zero. guarded: gas_scrubber_capacity.radicand_positive
    let velocity =
        design_gas_load_factor.value * ((liquid_density - gas_density) / gas_density).sqrt();
    let volumetric_flow = gas_outlet.n * gas_molar_mass / gas_density;
    Ok(volumetric_flow / (velocity * area))
}
