//! `unit_ops.separator` - a flash, split into vapour and liquid outlets.

use azoth_core::units::{Length, Power, Pressure, Velocity, VolumeRate, joules_per_mole, pascals};
use azoth_core::{AzothError, Result};
use azoth_eos::hydrate_inhibitor_wt::PhaseLabel;
use azoth_eos::{Phase, RootSide, ph_flash, pt_flash};

use crate::segment::phase::{labelled_phases_at_from_feed, mass_density, molar_mass_of};
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
    /// The gas phase's mass density, kg/m3; absent with no gas.
    pub gas_density: Option<f64>,
    /// The sizing family's liquid density, kg/m3; absent with no gas.
    pub liquid_density: Option<f64>,
    /// The design gas area, m2; absent with no gas.
    pub gas_area: Option<f64>,
    /// The gas phase's volumetric flow, m3/s; absent with no gas.
    pub volumetric_flow: Option<f64>,
    /// The Souders-Brown maximum gas flow, m3/s; absent with no gas.
    pub max_flow: Option<f64>,
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

    // **A zero flow is not a state**, and the heat-per-mole divides by it. The Python
    // arithmetic refuses `feed_n <= 0` the same way; a NaN is left to the flash.
    if feed.n <= 0.0 {
        return Err(AzothError::out_of_range(
            "feed_n",
            feed.n,
            "a unit operation needs a positive molar flow",
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
            let (x, y) = phase_compositions(flash.phase, &feed.z, flash.x, flash.y);
            (
                flash.temperature,
                x,
                y,
                vapour_fraction(flash.phase, flash.vapour_fraction),
            )
        }
        None => {
            let flash = pt_flash(&mixture, feed.t, p_out, &feed.z)?;
            let (x, y) = phase_compositions(flash.phase, &feed.z, flash.x, flash.y);
            (
                feed.t,
                x,
                y,
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

// `Separator`'s steady capacity getters, on the system `separator` leaves behind.
//
// **Two density rules, and they are not one.** The sizing family - the maximum velocity and
// the utilisation over it - reads the liquid's own density and substitutes
// `DEFAULT_LIQUID_DENSITY_FOR_SIZING` only where there is no liquid, with **no floor**. The
// gas-load factor applies `MIN_LIQUID_GAS_DENSITY_DIFFERENCE` to the same density first. A
// near-dry state makes the two disagree by construction, and neither is corrected to agree.
//
// The metric reads the system `run` left: `gas_outlet`'s temperature and pressure, the vapour
// that is *left* after entrainment, and the liquid that absorbed it - whose density
// `initPhysicalProperties` recomputes, which is why `liquid_outlet`'s composition is re-flashed
// rather than the feed's `x` read.

/// The sizing family's evidence, before the gas-load factor is formed.
///
/// **The assessment reads this and not the load factor.** `SeparatorCapacityAssessment` folds
/// the sizing velocity and the two flows; `getGasLoadFactor()`'s near-dry floor is a getter it
/// never reports, so a state the floor refuses is not a state the assessment cannot evidence.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CapacitySizing {
    /// The gas phase's mass density, kg/m3.
    pub gas_density: f64,
    /// The sizing family's liquid density, kg/m3.
    pub liquid_density: f64,
    /// The design gas area, m2.
    pub gas_area: f64,
    /// The gas phase's volumetric flow, m3/s.
    pub volumetric_flow: f64,
    /// `getMaxAllowableGasVelocity`, m/s.
    pub max_velocity: f64,
    /// `getGasSuperficialVelocity`, m/s.
    pub superficial_velocity: f64,
    /// `getMaxAllowableGasFlowRate`, m3/s.
    pub max_flow: f64,
    /// `getCapacityUtilization`.
    pub utilization: f64,
}

/// The sizing half of [`capacity_states`]: densities, area, flows, velocity and utilisation.
///
/// `Ok(None)` is the class's own first answer - a stream with no gas phase uses none of its
/// capacity - and is not an error.
///
/// # Errors
/// The gas density, the liquid-minus-gas difference and the admissible flow are the class's
/// own evidence checks; NeqSim answers `NaN` where this refuses.
pub fn sizing_states(
    feed: &Stream,
    gas_outlet: &Stream,
    liquid_outlet: &Stream,
    internal_diameter: Length,
    design_load_factor: Velocity,
    horizontal: bool,
    level_fraction: f64,
) -> Result<Option<CapacitySizing>> {
    let (mixture, _ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, gas_outlet.t, gas_outlet.p, &feed.z)?;
    let labelled =
        labelled_phases_at_from_feed(&mixture, gas_outlet.t, gas_outlet.p, &feed.z, &flash)?;

    let Some((_label, gas_z, gas_root)) = labelled
        .iter()
        .find(|(kind, _, _)| *kind == PhaseLabel::Gas)
    else {
        return Ok(None);
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
    Ok(Some(CapacitySizing {
        gas_density,
        liquid_density,
        gas_area,
        volumetric_flow,
        max_velocity,
        superficial_velocity,
        max_flow,
        utilization,
    }))
}

pub fn capacity_states(
    feed: &Stream,
    gas_outlet: &Stream,
    liquid_outlet: &Stream,
    internal_diameter: Length,
    design_load_factor: Velocity,
    horizontal: bool,
    level_fraction: f64,
) -> Result<Capacity> {
    let Some(sizing) = sizing_states(
        feed,
        gas_outlet,
        liquid_outlet,
        internal_diameter,
        design_load_factor,
        horizontal,
        level_fraction,
    )?
    else {
        // The class's own first answer: a stream with no gas phase uses none of its capacity.
        return Ok(Capacity {
            utilization: 0.0,
            max_velocity: None,
            superficial_velocity: None,
            load_factor: None,
            overloaded: false,
            gas_density: None,
            liquid_density: None,
            gas_area: None,
            volumetric_flow: None,
            max_flow: None,
        });
    };

    // **The load factor's own rule, applied after the sizing one.** The floor substitutes
    // `1000.0` for a second phase no denser than the gas, and the positivity is re-checked
    // against the substituted value.
    let floored_liquid_density =
        load_factor_liquid_density(sizing.liquid_density, sizing.gas_density);
    if !floored_liquid_density.is_finite() || floored_liquid_density <= sizing.gas_density {
        return Err(AzothError::invalid_input(
            "gas_load_factor",
            format!(
                "the liquid phase is {floored_liquid_density} kg/m3 after the near-dry floor \
                 against the gas's {}, and the load factor divides by their difference",
                sizing.gas_density
            ),
        ));
    }
    let term = (floored_liquid_density - sizing.gas_density) / sizing.gas_density;
    let load_factor = sizing.superficial_velocity * (1.0 / term).sqrt();
    if !load_factor.is_finite() {
        return Err(AzothError::invalid_input(
            "gas_load_factor",
            format!("the gas load factor is {load_factor} m/s"),
        ));
    }

    Ok(Capacity {
        utilization: sizing.utilization,
        max_velocity: Some(sizing.max_velocity),
        superficial_velocity: Some(sizing.superficial_velocity),
        load_factor: Some(load_factor),
        overloaded: sizing.utilization > 1.0,
        gas_density: Some(sizing.gas_density),
        liquid_density: Some(sizing.liquid_density),
        gas_area: Some(sizing.gas_area),
        volumetric_flow: Some(sizing.volumetric_flow),
        max_flow: Some(sizing.max_flow),
    })
}

/// `Separator.getGasLoadFactor(int)` and `Separator.sizeFromFlow(safety)`, on the state
/// [`separator`] leaves behind.
#[derive(Debug, Clone, PartialEq)]
pub struct SizingGetters {
    /// `1 - liquidVolume / separatorVolume` horizontal, `1.0` vertical.
    pub gas_area_fraction: f64,
    /// The operating liquid volume, m3; the class's circular segment times the length.
    pub liquid_volume: f64,
    /// The vessel's total volume, m3; the cross-section times the length.
    pub separator_volume: f64,
    /// `getGasLoadFactor(int)` for each liquid phase, in the class's phase order.
    pub rated_by_phase: Vec<f64>,
    /// `sizeFromFlow(safety)`'s diameter, m; absent when the system carries no gas.
    pub required_diameter: Option<f64>,
}

/// The phase-indexed load factors and the Souders-Brown sizing.
///
/// **Two gas area fractions, and they are different fields.** The superficial velocity - and
/// therefore every load factor - divides by the **design** area, `1 - design_level_fraction`
/// horizontal. The indexed overload additionally divides its radicand by the **operating**
/// area fraction `1 - liquidVolume / separatorVolume`, which is the class's own holdup. A
/// vertical vessel is `1.0` on the second and the two overloads then agree.
///
/// `Ok(None)` is a system with no gas phase: the class's sizing returns without a diameter,
/// and the capacity family publishes nothing.
///
/// # Errors
/// The same density and radicand refusals [`sizing_states`] makes, and the liquid level's
/// zero-area case, which the class answers with a non-finite answer.
#[allow(clippy::too_many_arguments)] // one per declared input the fold reads
pub fn sizing_getters(
    feed: &Stream,
    gas_outlet: &Stream,
    liquid_outlet: &Stream,
    internal_diameter: Length,
    design_load_factor: Velocity,
    horizontal: bool,
    design_level_fraction: f64,
    separator_length: Length,
    liquid_level_fraction: f64,
    safety_factor: Option<f64>,
) -> Result<Option<SizingGetters>> {
    let Some(sizing) = sizing_states(
        feed,
        gas_outlet,
        liquid_outlet,
        internal_diameter,
        design_load_factor,
        horizontal,
        design_level_fraction,
    )?
    else {
        return Ok(None);
    };

    let diameter = internal_diameter.value;
    let length = separator_length.value;
    let cross_area = std::f64::consts::PI * diameter * diameter / 4.0;
    let separator_volume = cross_area * length;
    let level = if horizontal {
        liquid_level_fraction * diameter
    } else {
        liquid_level_fraction * length
    };
    let liquid_volume = if horizontal {
        liquid_area(level, diameter, true) * length
    } else {
        cross_area * level
    };
    let gas_area_fraction = if horizontal {
        1.0 - liquid_volume / separator_volume
    } else {
        1.0
    };
    if !gas_area_fraction.is_finite() || gas_area_fraction <= 0.0 {
        return Err(AzothError::invalid_input(
            "liquid_level_fraction",
            format!(
                "the operating liquid level of {liquid_level_fraction} leaves a gas area \
                 fraction of {gas_area_fraction}, and the phase-indexed load factor divides by it"
            ),
        ));
    }

    // **The phase densities are the system's own, after the run.** The class's
    // `getGasLoadFactor(int)` calls `initPhysicalProperties` and reads each phase, which is
    // the flash of each phase's composition - so the labelled phases at the outlet state are
    // the same states, and the entrained liquid is the liquid outlet's re-flash (which
    // `sizing_states` already measured for the no-argument form).
    let (mixture, _ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, gas_outlet.t, gas_outlet.p, &feed.z)?;
    let labelled =
        labelled_phases_at_from_feed(&mixture, gas_outlet.t, gas_outlet.p, &feed.z, &flash)?;

    let mut rated_by_phase = Vec::new();
    let mut oil_density = None;
    let mut aqueous_density = None;
    for (kind, z, root) in &labelled {
        if *kind == PhaseLabel::Gas {
            continue;
        }
        let density = mass_density(
            &mixture,
            gas_outlet.t,
            gas_outlet.p,
            z,
            *root,
            molar_mass_of(&mixture, z)?,
        )?;
        match kind {
            PhaseLabel::Oil => oil_density = Some(density),
            PhaseLabel::Aqueous => aqueous_density = Some(density),
            PhaseLabel::Gas => {}
        }
        let floored = load_factor_liquid_density(density, sizing.gas_density);
        if !floored.is_finite() || floored <= sizing.gas_density {
            return Err(AzothError::invalid_input(
                "gas_load_factor",
                format!(
                    "the liquid phase is {floored} kg/m3 after the near-dry floor against the \
                     gas's {}, and the indexed load factor divides by their difference",
                    sizing.gas_density
                ),
            ));
        }
        let term = (1.0 / gas_area_fraction) * (floored - sizing.gas_density) / sizing.gas_density;
        if !term.is_finite() || term <= 0.0 {
            return Err(AzothError::invalid_input(
                "gas_load_factor",
                format!("the phase-indexed load factor's radicand is {term}"),
            ));
        }
        rated_by_phase.push(sizing.superficial_velocity * (1.0 / term).sqrt());
    }

    // **`sizeFromFlow` uses the sizing family's liquid density, not the load factor's floor.**
    // The class reads the oil phase, then the aqueous phase, then `1000.0`.
    let required_diameter = safety_factor
        .map(|safety| {
            let liquid_density = oil_density
                .or(aqueous_density)
                .unwrap_or(DEFAULT_LIQUID_DENSITY_FOR_SIZING);
            if !liquid_density.is_finite() || liquid_density <= sizing.gas_density {
                return Err(AzothError::invalid_input(
                    "safety_factor",
                    format!(
                        "the liquid phase is {liquid_density} kg/m3 against the gas's {}, so the \
                         Souders-Brown velocity `sizeFromFlow` divides by is not positive",
                        sizing.gas_density
                    ),
                ));
            }
            let max_velocity = design_load_factor.value
                * ((liquid_density - sizing.gas_density) / sizing.gas_density).sqrt();
            let gas_area_share = if horizontal {
                1.0 - design_level_fraction
            } else {
                1.0
            };
            let required_gas_area = sizing.volumetric_flow * safety / max_velocity;
            let total_area = required_gas_area / gas_area_share;
            if !total_area.is_finite() || total_area <= 0.0 {
                return Err(AzothError::invalid_input(
                    "safety_factor",
                    format!("`sizeFromFlow` answers a total area of {total_area} m2"),
                ));
            }
            Ok((4.0 * total_area / std::f64::consts::PI).sqrt())
        })
        .transpose()?;

    Ok(Some(SizingGetters {
        gas_area_fraction,
        liquid_volume,
        separator_volume,
        rated_by_phase,
        required_diameter,
    }))
}

/// `Separator.liquidArea`: the circular segment below the liquid level, on a horizontal
/// vessel; the full cross-section on a vertical one.
fn liquid_area(level: f64, diameter: f64, horizontal: bool) -> f64 {
    let cross_area = std::f64::consts::PI * diameter * diameter / 4.0;
    if level <= 0.0 {
        return 0.0;
    }
    if level >= diameter {
        return cross_area;
    }
    if !horizontal {
        return cross_area;
    }
    let radius = diameter / 2.0;
    if level < radius {
        let d = radius - level;
        let theta = (d / radius).acos();
        let a = radius * theta.sin();
        theta * radius * radius - a * d
    } else if level > radius {
        let d = level - radius;
        let theta = (d / radius).acos();
        let a = radius * theta.sin();
        (std::f64::consts::PI - theta) * radius * radius + a * d
    } else {
        0.5 * std::f64::consts::PI * radius * radius
    }
}

/// `SeparatorCapacityAssessment`'s fixed diagnostics, verbatim.
const PROVENANCE_REQUIRED: &str = "Design-basis provenance is required";
const BAD_ORIENTATION: &str = "Separator orientation must be horizontal or vertical";
const GEOMETRY_REQUIRED: &str =
    "Finite geometry, design K-factor, and a calculated gas phase are required";
const BAD_LEVEL: &str = "Design liquid level fraction must be in [0, 1)";
const BAD_EVIDENCE: &str = "Separator capacity calculation returned non-finite evidence";
const FEASIBLE: &str = "Gas capacity is within the evidenced design limit";
const INFEASIBLE: &str = "Gas capacity exceeds the evidenced design limit";

/// `SeparatorCapacityAssessment`'s snapshot, in the class's field order.
#[derive(Debug, Clone, PartialEq)]
pub struct CapacityAssessment {
    /// `feasible`, `infeasible` or `unassessed`.
    pub status: &'static str,
    /// The class's sentence for the status.
    pub diagnostic: &'static str,
    /// `equipment` or `mechanicalDesign`; empty on `unassessed`.
    pub design_basis_source: &'static str,
    /// `soudersBrown` or `mechanicalDesignGasFlow`; empty on `unassessed`.
    pub governing_capacity_basis: &'static str,
    /// The design gas area, m2; absent on `unassessed`.
    pub gas_area: Option<f64>,
    /// The gas phase's density, kg/m3; absent on `unassessed`.
    pub gas_density: Option<f64>,
    /// The sizing family's liquid density, kg/m3; absent on `unassessed`.
    pub liquid_density: Option<f64>,
    /// `getMaxAllowableGasVelocity`, m/s; absent on `unassessed`.
    pub maximum_gas_velocity: Option<f64>,
    /// The gas phase's volumetric flow, m3/s; absent on `unassessed`.
    pub current_gas_flow: Option<f64>,
    /// The stated mechanical limit, m3/s, or `0.0` with none; absent on `unassessed`.
    pub mechanical_design_maximum_gas_flow: Option<f64>,
    /// `getMaxAllowableGasFlowRate`, m3/s; absent on `unassessed`.
    pub souders_brown_maximum_gas_flow: Option<f64>,
    /// The governed limit, m3/s; absent on `unassessed`.
    pub maximum_gas_flow: Option<f64>,
    /// Current flow over the governed limit; absent on `unassessed`.
    pub utilization: Option<f64>,
}

/// `SeparatorCapacityAssessment.from`, on the system [`separator`] leaves behind.
///
/// **Fail closed rather than refuse.** Blank provenance, an orientation that is neither of the
/// two, non-positive geometry or design K, no gas phase, a level outside `[0, 1)` and
/// non-finite evidence answer `unassessed` with the class's diagnostic. The checks run in the
/// class's order, so a state failing two carries the earlier sentence.
///
/// The mechanical limit is the class's `maxDesignGassVolumeFlow / 3600`: a positive limit at or
/// below the Souders-Brown flow governs, and a non-positive one is no mechanical basis.
///
/// # Errors
/// The separator's own where the flash cannot run, and [`azoth_core::AzothError::InvalidInput`]
/// for a stated mechanical limit that is not finite - an infinity is a value no snapshot
/// publishes. The assessment's other failing states are statuses, not errors.
#[allow(clippy::too_many_arguments)] // one per declared input the fold reads
pub fn capacity_assessment(
    feed: &Stream,
    gas_outlet: &Stream,
    liquid_outlet: &Stream,
    internal_diameter: Length,
    design_load_factor: Velocity,
    orientation: &str,
    level_fraction: f64,
    provenance: &str,
    mechanical_design_max_gas_flow: Option<VolumeRate>,
) -> Result<CapacityAssessment> {
    if provenance.trim().is_empty() {
        return Ok(unassessed(PROVENANCE_REQUIRED));
    }
    if orientation != "horizontal" && orientation != "vertical" {
        return Ok(unassessed(BAD_ORIENTATION));
    }
    if !finite_positive(internal_diameter.value) || !finite_positive(design_load_factor.value) {
        return Ok(unassessed(GEOMETRY_REQUIRED));
    }
    // **The class checks the gas phase before the level.** A liquid-only feed with a level of
    // one answers the geometry sentence, not the level's.
    if !has_gas_phase(feed, gas_outlet)? {
        return Ok(unassessed(GEOMETRY_REQUIRED));
    }
    // `contains` is the class's interval and also refuses a NaN, which is not in any range.
    if !(0.0..1.0).contains(&level_fraction) {
        return Ok(unassessed(BAD_LEVEL));
    }

    // **The sizing half only.** The assessment does not report `getGasLoadFactor`, so the
    // load factor's near-dry refusal is not evidence `unassessed` may rest on.
    let Some(sizing) = (match sizing_states(
        feed,
        gas_outlet,
        liquid_outlet,
        internal_diameter,
        design_load_factor,
        orientation == "horizontal",
        level_fraction,
    ) {
        Ok(sizing) => sizing,
        Err(_) => return Ok(unassessed(BAD_EVIDENCE)),
    }) else {
        return Ok(unassessed(GEOMETRY_REQUIRED));
    };

    let mechanical_flow = match mechanical_design_max_gas_flow {
        Some(limit) if !limit.value.is_finite() => {
            return Err(AzothError::invalid_input(
                "mechanical_design_max_gas_flow",
                format!(
                    "a stated mechanical gas limit is a finite flow, and {} is not",
                    limit.value
                ),
            ));
        }
        Some(limit) => limit.value,
        None => 0.0,
    };
    let mechanical_positive = finite_positive(mechanical_flow);
    let design_basis_source = if mechanical_positive {
        "mechanicalDesign"
    } else {
        "equipment"
    };
    let mechanical_governs = mechanical_positive && mechanical_flow <= sizing.max_flow;
    let maximum_flow = if mechanical_governs {
        mechanical_flow
    } else {
        sizing.max_flow
    };
    let governing_capacity_basis = if mechanical_governs {
        "mechanicalDesignGasFlow"
    } else {
        "soudersBrown"
    };
    if !finite_positive(maximum_flow) {
        return Ok(unassessed(BAD_EVIDENCE));
    }
    let utilization = sizing.volumetric_flow / maximum_flow;
    if !utilization.is_finite() || utilization < 0.0 {
        return Ok(unassessed(BAD_EVIDENCE));
    }
    let feasible = utilization <= 1.0;

    Ok(CapacityAssessment {
        status: if feasible { "feasible" } else { "infeasible" },
        diagnostic: if feasible { FEASIBLE } else { INFEASIBLE },
        design_basis_source,
        governing_capacity_basis,
        gas_area: Some(sizing.gas_area),
        gas_density: Some(sizing.gas_density),
        liquid_density: Some(sizing.liquid_density),
        maximum_gas_velocity: Some(sizing.max_velocity),
        current_gas_flow: Some(sizing.volumetric_flow),
        mechanical_design_maximum_gas_flow: Some(mechanical_flow),
        souders_brown_maximum_gas_flow: Some(sizing.max_flow),
        maximum_gas_flow: Some(maximum_flow),
        utilization: Some(utilization),
    })
}

/// Whether the system `run` left carries a gas phase - the class's `hasPhaseType("gas")`.
fn has_gas_phase(feed: &Stream, gas_outlet: &Stream) -> Result<bool> {
    let (mixture, _ideal_gas) = feed.mixture()?;
    let flash = pt_flash(&mixture, gas_outlet.t, gas_outlet.p, &feed.z)?;
    let labelled =
        labelled_phases_at_from_feed(&mixture, gas_outlet.t, gas_outlet.p, &feed.z, &flash)?;
    Ok(labelled.iter().any(|(kind, _, _)| *kind == PhaseLabel::Gas))
}

fn unassessed(diagnostic: &'static str) -> CapacityAssessment {
    CapacityAssessment {
        status: "unassessed",
        diagnostic,
        design_basis_source: "",
        governing_capacity_basis: "",
        gas_area: None,
        gas_density: None,
        liquid_density: None,
        maximum_gas_velocity: None,
        current_gas_flow: None,
        mechanical_design_maximum_gas_flow: None,
        souders_brown_maximum_gas_flow: None,
        maximum_gas_flow: None,
        utilization: None,
    }
}

fn finite_positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
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

/// The two compositions a flash's phase implies.
///
/// **A single phase has the feed's composition, and the flash's own `x` or `y` is not it.**
/// The flash answers a trial there - and on its negative-flash branch that trial is the other
/// side of the cubic, so a 95/5 methane/ethane gas at 50 bar would leave the vessel with the
/// trial composition and a density 2.6 times the feed's. The two-phase branch keeps the flash's
/// split.
fn phase_compositions(phase: Phase, z: &[f64], x: Vec<f64>, y: Vec<f64>) -> (Vec<f64>, Vec<f64>) {
    match phase {
        Phase::TwoPhase => (x, y),
        Phase::AllVapour | Phase::AllLiquid | Phase::Trivial => (z.to_vec(), z.to_vec()),
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
