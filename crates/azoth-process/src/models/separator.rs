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
    Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, Time, Velocity,
    joules_per_mole, meters, meters_per_second, pascals, seconds,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::kernels::separator as kernel;
use crate::model_gen;
use crate::separator_design::mechanical::{
    DesignInputs, VesselKind, mechanical_design as design_vessel,
};
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
    /// `1 - liquidVolume / separatorVolume` horizontal, `1.0` vertical.
    pub gas_area_fraction: Option<f64>,
    /// `getGasLoadFactor(int)` for each liquid phase, in the class's phase order.
    pub rated_gas_load_factor_by_phase: Option<Vec<f64>>,
    /// `sizeFromFlow(safety)`'s diameter; absent when no safety factor was stated.
    pub required_diameter: Option<Length>,
    /// `SeparatorMechanicalDesign.calcDesign`'s sized diameter; absent unless the design ran.
    pub design_inner_diameter: Option<Length>,
    /// The design's tangent-to-tangent length.
    pub design_tantan_length: Option<Length>,
    /// The design's outside diameter.
    pub design_outer_diameter: Option<Length>,
    /// The design's wall thickness.
    pub design_wall_thickness: Option<Length>,
    /// `getMaxDesignGassVolumeFlow`, m3/hr.
    pub max_design_gas_volume_flow: Option<f64>,
    /// The empty vessel shell's weight, kg.
    pub vessel_shell_weight: Option<f64>,
    /// The piping allowance, kg.
    pub piping_weight: Option<f64>,
    /// The structural allowance, kg.
    pub structural_weight: Option<f64>,
    /// The electrical allowance, kg.
    pub electrical_weight: Option<f64>,
    /// The internals' weight, kg; zero until the internals sections are ported.
    pub internals_weight: Option<f64>,
    /// The total skid weight, kg.
    pub total_weight: Option<f64>,
    /// The module's width.
    pub module_width: Option<Length>,
    /// The module's length.
    pub module_length: Option<Length>,
    /// The module's height.
    pub module_height: Option<Length>,
    /// The vessel's internal volume, m3.
    pub vessel_volume: Option<f64>,
    /// HHLL, as a fraction of the design diameter.
    pub hhll: Option<Length>,
    /// HLL.
    pub hll: Option<Length>,
    /// NLL.
    pub nll: Option<Length>,
    /// LLL.
    pub lll: Option<Length>,
    /// LLLL.
    pub llll: Option<Length>,
    /// The weir height.
    pub weir_height: Option<Length>,
    /// The normal interface level.
    pub nil: Option<Length>,
    /// The effective liquid separation length.
    pub effective_length_liquid: Option<Length>,
    /// The effective gas separation length.
    pub effective_length_gas: Option<Length>,
    /// The inlet nozzle's inside diameter.
    pub inlet_nozzle_id: Option<Length>,
    /// The gas outlet nozzle's inside diameter.
    pub gas_outlet_nozzle_id: Option<Length>,
    /// The oil outlet nozzle's inside diameter.
    pub oil_outlet_nozzle_id: Option<Length>,
    /// The surge volume between NLL and HLL, m3.
    pub surge_volume: Option<f64>,
    /// The surge time, seconds.
    pub surge_time: Option<Time>,
    /// The liquid retention time at NLL, seconds.
    pub liquid_retention_time: Option<Time>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Flash a stream into vapour and liquid outlets.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if the pressure drop takes the outlet below
/// zero, the entrainment fraction is outside `[0, 1]`, or the shapes disagree, and
/// whatever the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-three
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
    separator_length: Option<Length>,
    liquid_level_fraction: Option<f64>,
    safety_factor: Option<f64>,
    mechanical_design: Option<bool>,
    auto_size: Option<bool>,
    fg: Option<f64>,
    volume_safety_factor: Option<f64>,
    retention_time: Option<Time>,
    max_operation_pressure: Option<Pressure>,
    allowable_stress: Option<f64>,
    joint_efficiency: Option<f64>,
    corrosion_allowance: Option<f64>,
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
            "separator_length" => separator_length.map(|value| value.value),
            "liquid_level_fraction" => liquid_level_fraction,
            "safety_factor" => safety_factor,
            "fg" => fg,
            "volume_safety_factor" => volume_safety_factor,
            "retention_time" => retention_time.map(|value| value.value),
            "max_operation_pressure" => max_operation_pressure.map(|value| value.value),
            "allowable_stress" => allowable_stress,
            "joint_efficiency" => joint_efficiency,
            "corrosion_allowance" => corrosion_allowance,
            _ => None,
        },
        &mut warnings,
    )?;

    let design = mechanical_design.unwrap_or(false);
    let auto = auto_size.unwrap_or(false);
    // **A stated diameter needs its design basis; a design basis does not need a diameter.**
    // The capacity metric is the pair, but the mechanical design sizes its own diameter from
    // the flow and the design `K`, so `design_gas_load_factor` alone is a complete statement
    // there.
    if internal_diameter.is_some() && design_gas_load_factor.is_none() {
        return Err(AzothError::invalid_input(
            "design_gas_load_factor",
            "`internal_diameter` was stated and `design_gas_load_factor` was not, and the \
             Souders-Brown capacity metric needs both - a vessel with a size but no design basis \
             has nothing to check",
        ));
    }
    if design && auto {
        return Err(AzothError::invalid_input(
            "auto_size",
            "`mechanical_design` and `auto_size` are two ways to design the same vessel; state one",
        ));
    }
    if auto {
        if internal_diameter.is_some() {
            return Err(AzothError::invalid_input(
                "internal_diameter",
                "`auto_size` sizes the vessel itself, so a stated diameter would be overwritten",
            ));
        }
        if volume_safety_factor.is_some() {
            return Err(AzothError::invalid_input(
                "volume_safety_factor",
                "`auto_size` sets the volume factor from `safety_factor`; state it there",
            ));
        }
    }
    if !design && !auto && design_gas_load_factor.is_some() && internal_diameter.is_none() {
        return Err(AzothError::invalid_input(
            "internal_diameter",
            "`design_gas_load_factor` was stated without `internal_diameter` and without \
             `mechanical_design`, and the design basis has no vessel to be used on",
        ));
    }
    if design && design_gas_load_factor.is_none() {
        return Err(AzothError::invalid_input(
            "design_gas_load_factor",
            "`mechanical_design` was stated without `design_gas_load_factor`, and the design's \
             Souders-Brown basis is that number",
        ));
    }
    // The two that shape the gas area have no meaning without the vessel they shape.
    if internal_diameter.is_none() {
        if orientation.is_some() && !design && !auto {
            return Err(AzothError::invalid_input(
                "orientation",
                "`orientation` was stated without `internal_diameter`, `design_gas_load_factor` \
                 and `mechanical_design`, and the gas area it shapes is not computed until the \
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
        for (name, stated) in [
            ("separator_length", separator_length.is_some()),
            ("liquid_level_fraction", liquid_level_fraction.is_some()),
            // `safety_factor` is `autoSize`'s volume factor when `auto_size` is set.
            ("safety_factor", safety_factor.is_some() && !auto),
        ] {
            if stated {
                return Err(AzothError::invalid_input(
                    name,
                    format!(
                        "`{name}` was stated without `internal_diameter` and \
                         `design_gas_load_factor`, and the vessel it sizes is not computed until \
                         both are stated"
                    ),
                ));
            }
        }
    }
    if !design && !auto {
        for (name, stated) in [
            ("fg", fg.is_some()),
            ("volume_safety_factor", volume_safety_factor.is_some()),
            ("retention_time", retention_time.is_some()),
            ("max_operation_pressure", max_operation_pressure.is_some()),
            ("allowable_stress", allowable_stress.is_some()),
            ("joint_efficiency", joint_efficiency.is_some()),
            ("corrosion_allowance", corrosion_allowance.is_some()),
        ] {
            if stated {
                return Err(AzothError::invalid_input(
                    name,
                    format!("`{name}` was stated without `mechanical_design`"),
                ));
            }
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
    let sizing = match (internal_diameter, design_gas_load_factor) {
        (Some(diameter), Some(load_factor)) => kernel::sizing_getters(
            &feed,
            &vapour,
            &liquid,
            diameter,
            load_factor,
            orientation == "horizontal",
            level_fraction,
            separator_length.unwrap_or(meters(5.0)),
            liquid_level_fraction.unwrap_or(0.5),
            safety_factor,
        )?,
        _ => None,
    };
    let vessel = if design || auto {
        // `autoSize` reads the design standards' K when the separator has no stated one.
        let load_factor = match design_gas_load_factor {
            Some(value) => value.value,
            // The shipped design database's `GasLoadFactor` midpoint, which `autoSize` reads.
            None if auto => 0.137,
            None => {
                return Err(AzothError::invalid_input(
                    "design_gas_load_factor",
                    "`mechanical_design` needs the design's Souders-Brown basis",
                ));
            }
        };
        let inputs = DesignInputs {
            gas_load_factor: load_factor,
            // The database's `Fg` midpoint for the same standard.
            fg: fg.unwrap_or(if auto { 0.5 } else { 0.8 }),
            // `autoSize(safetyFactor)` sets the volume factor; its no-argument form is 1.2.
            volume_safety_factor: if auto {
                safety_factor.unwrap_or(1.2)
            } else {
                volume_safety_factor.unwrap_or(1.0)
            },
            // The database's `LiquidRetentionTime` is 150 s; `autoSize` keeps 120 s until
            // `sizeForLiquidHoldup` floors it by the minimum.
            retention_time_s: retention_time
                .map_or(if auto { 120.0 } else { 150.0 }, |value| value.value),
            max_operation_pressure: max_operation_pressure.unwrap_or(pascals(100.0e5)),
            // The spec declares pascals; the ASME formula is written in MPa.
            allowable_stress: allowable_stress.unwrap_or(95.147650602e6) / 1.0e6,
            joint_efficiency: joint_efficiency.unwrap_or(1.0),
            // The spec declares millimetres; the kernel's field is millimetres too, so a
            // converted value arrives in metres and is scaled back.
            corrosion_allowance_mm: corrosion_allowance.unwrap_or(0.0) * 1000.0,
            max_gas_velocity: 3.0,
            foam_allowance_factor: 1.35,
            min_oil_retention_s: 120.0,
            min_water_retention_s: 180.0,
            aqueous: false,
        };
        Some(design_vessel(
            &vapour,
            &liquid,
            orientation == "horizontal",
            VesselKind::Separator,
            &inputs,
            auto,
        )?)
    } else {
        None
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
        gas_area_fraction: sizing.as_ref().map(|states| states.gas_area_fraction),
        rated_gas_load_factor_by_phase: sizing
            .as_ref()
            .map(|states| states.rated_by_phase.clone())
            .filter(|values| !values.is_empty()),
        required_diameter: sizing
            .as_ref()
            .and_then(|states| states.required_diameter)
            .map(meters),
        design_inner_diameter: vessel.as_ref().map(|design| meters(design.inner_diameter)),
        design_tantan_length: vessel.as_ref().map(|design| meters(design.tantan_length)),
        design_outer_diameter: vessel.as_ref().map(|design| meters(design.outer_diameter)),
        design_wall_thickness: vessel.as_ref().map(|design| meters(design.wall_thickness)),
        max_design_gas_volume_flow: vessel
            .as_ref()
            .map(|design| design.max_design_gas_volume_flow_m3_hr / 3600.0),
        vessel_shell_weight: vessel.as_ref().map(|design| design.vessel_shell_weight_kg),
        piping_weight: vessel.as_ref().map(|design| design.piping_weight_kg),
        structural_weight: vessel.as_ref().map(|design| design.structural_weight_kg),
        electrical_weight: vessel.as_ref().map(|design| design.electrical_weight_kg),
        // No internals field exists in the port, and the class reports zero without the
        // `sectiontype` sections that would fill it.
        internals_weight: vessel.as_ref().map(|_| 0.0),
        total_weight: vessel.as_ref().map(|design| design.total_weight_kg),
        module_width: vessel.as_ref().map(|design| meters(design.module_width)),
        module_length: vessel.as_ref().map(|design| meters(design.module_length)),
        module_height: vessel.as_ref().map(|design| meters(design.module_height)),
        vessel_volume: vessel.as_ref().map(|design| design.vessel_volume_m3),
        hhll: vessel.as_ref().map(|design| meters(design.hhll)),
        hll: vessel.as_ref().map(|design| meters(design.hll)),
        nll: vessel.as_ref().map(|design| meters(design.nll)),
        lll: vessel.as_ref().map(|design| meters(design.lll)),
        llll: vessel.as_ref().map(|design| meters(design.llll)),
        weir_height: vessel.as_ref().map(|design| meters(design.weir_height)),
        nil: vessel.as_ref().map(|design| meters(design.nil)),
        effective_length_liquid: vessel
            .as_ref()
            .map(|design| meters(design.effective_length_liquid)),
        effective_length_gas: vessel
            .as_ref()
            .map(|design| meters(design.effective_length_gas)),
        inlet_nozzle_id: vessel.as_ref().map(|design| meters(design.inlet_nozzle_id)),
        gas_outlet_nozzle_id: vessel
            .as_ref()
            .map(|design| meters(design.gas_outlet_nozzle_id)),
        oil_outlet_nozzle_id: vessel
            .as_ref()
            .map(|design| meters(design.oil_outlet_nozzle_id)),
        surge_volume: vessel.as_ref().map(|design| design.surge_volume_m3),
        surge_time: vessel.as_ref().map(|design| seconds(design.surge_time_s)),
        liquid_retention_time: vessel
            .as_ref()
            .map(|design| seconds(design.liquid_retention_time_s)),
        warnings,
    })
}
