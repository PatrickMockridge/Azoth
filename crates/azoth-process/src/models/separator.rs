//! `process.separator` - the separator's kernel as a registered id.
//!
//! Spec: `specs/models/process/separator.toml`. The arithmetic is
//! [`crate::kernels::separator`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **Two outlets, so ten fields.** A single-multiplicity port crosses as the record's five
//! fields per port, and a separator has two of them, so a result is the pair written out
//! under the ports' own names rather than a list.

use azoth_core::units::{
    Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, Velocity, joules_per_mole,
    meters_per_second,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::kernels::separator as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.separator`.
#[derive(Debug, Clone, PartialEq)]
pub struct SeparatorResult {
    /// Vapour outlet molar flow, mol/s.
    pub vapour_n: f64,
    /// Vapour outlet composition.
    pub vapour_z: Vec<f64>,
    /// Vapour outlet pressure.
    pub vapour_p: Pressure,
    /// Vapour outlet temperature.
    pub vapour_t: ThermodynamicTemperature,
    /// Vapour outlet molar enthalpy.
    pub vapour_h: MolarEnergy,
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
    /// `Separator.getCapacityUtilization`, present only when the vessel's size and its design
    /// gas load factor were both stated.
    pub capacity_utilization: Option<f64>,
    /// `Separator.getMaxAllowableGasVelocity`, absent with no gas phase.
    pub max_allowable_gas_velocity: Option<Velocity>,
    /// `Separator.getGasSuperficialVelocity`, absent as above.
    pub gas_superficial_velocity: Option<Velocity>,
    /// `Separator.getGasLoadFactor()`, absent as above. Not the design factor.
    pub gas_load_factor: Option<Velocity>,
    /// `Separator.isOverloaded`.
    pub is_overloaded: Option<bool>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Flash a stream into vapour and liquid outlets.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the pressure drop takes the outlet below
/// zero, the entrainment fraction is outside `[0, 1]`, or the shapes disagree, and
/// whatever the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twelve
pub fn separator(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    pressure_drop: Pressure,
    gas_in_liquid: f64,
    heat_input: Option<Power>,
    internal_diameter: Option<Length>,
    design_gas_load_factor: Option<Velocity>,
    orientation: Option<&str>,
    design_liquid_level_fraction: Option<f64>,
) -> Result<SeparatorResult> {
    let spec = &model_gen::SEPARATOR_SPEC;
    let mut warnings = Vec::new();
    // One call over every declared input, for the reason `reactions.chemical_equilibrium`
    // gives. `heat_input` has no bound - a duty may heat or cool - so the `_` arm is the
    // whole of it. An absent optional input is `None` here, which `apply_checks` reports as a
    // skipped check rather than passing it silently.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "pressure_drop" => Some(pressure_drop.value),
            "gas_in_liquid" => Some(gas_in_liquid),
            "internal_diameter" => internal_diameter.map(|value| value.value),
            "design_gas_load_factor" => design_gas_load_factor.map(|value| value.value),
            "design_liquid_level_fraction" => design_liquid_level_fraction,
            _ => None,
        },
        &mut warnings,
    )?;

    // **Both or neither.** The capacity is a function of the vessel, and a caller who stated
    // half of it has stated nothing.
    if internal_diameter.is_some() != design_gas_load_factor.is_some() {
        let (named, missing) = if internal_diameter.is_some() {
            ("internal_diameter", "design_gas_load_factor")
        } else {
            ("design_gas_load_factor", "internal_diameter")
        };
        return Err(AzothError::invalid_input(
            missing,
            format!(
                "`{named}` was stated and `{missing}` was not, and the Souders-Brown capacity \
                 metric needs both - a vessel with a size but no design basis, or the reverse, \
                 has nothing to check"
            ),
        ));
    }
    // The two that shape the gas area have no meaning without the vessel they shape.
    if internal_diameter.is_none() {
        if orientation.is_some() {
            return Err(AzothError::invalid_input(
                "orientation",
                "`orientation` was stated without `internal_diameter` and \
                 `design_gas_load_factor`, and the gas area it shapes is not computed until the \
                 vessel is",
            ));
        }
        if design_liquid_level_fraction.is_some() {
            return Err(AzothError::invalid_input(
                "design_liquid_level_fraction",
                "`design_liquid_level_fraction` was stated without `internal_diameter` and \
                 `design_gas_load_factor`, and the gas area it shapes is not computed until the \
                 vessel is",
            ));
        }
    }
    let orientation = orientation.unwrap_or("horizontal");
    if orientation != "horizontal" && orientation != "vertical" {
        return Err(AzothError::invalid_input(
            "orientation",
            format!("a vessel is `horizontal` or `vertical`, and {orientation:?} is neither"),
        ));
    }
    let level_fraction = design_liquid_level_fraction.unwrap_or(0.8);

    // **A zero flow is not a state**, and the heat-per-mole is a division by it: refused
    // here with the input's own name rather than surfacing as a non-finite enthalpy.
    if feed_n <= 0.0 {
        return Err(azoth_core::AzothError::out_of_range(
            "feed_n",
            feed_n,
            "a unit operation needs a positive molar flow",
        ));
    }

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let (vapour, liquid) = kernel(&feed, pressure_drop, gas_in_liquid, heat_input)?;

    let capacity = match (internal_diameter, design_gas_load_factor) {
        (Some(diameter), Some(load_factor)) => Some(kernel::capacity_states(
            &feed,
            &vapour,
            &liquid,
            diameter,
            load_factor,
            orientation == "horizontal",
            level_fraction,
        )?),
        _ => None,
    };

    Ok(SeparatorResult {
        vapour_n: vapour.n,
        vapour_z: vapour.z,
        vapour_p: vapour.p,
        vapour_t: vapour.t,
        vapour_h: joules_per_mole(vapour.h.value),
        liquid_n: liquid.n,
        liquid_z: liquid.z,
        liquid_p: liquid.p,
        liquid_t: liquid.t,
        liquid_h: joules_per_mole(liquid.h.value),
        capacity_utilization: capacity.map(|states| states.utilization),
        max_allowable_gas_velocity: capacity
            .and_then(|states| states.max_velocity)
            .map(meters_per_second),
        gas_superficial_velocity: capacity
            .and_then(|states| states.superficial_velocity)
            .map(meters_per_second),
        gas_load_factor: capacity
            .and_then(|states| states.load_factor)
            .map(meters_per_second),
        is_overloaded: capacity.map(|states| states.overloaded),
        warnings,
    })
}
