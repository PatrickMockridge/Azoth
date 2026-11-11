//! `unit_ops.splitter` - one stream, split into several with the same state.

use azoth_core::units::{Length, Pressure, Velocity, meters, meters_per_second, pascals};
use azoth_core::{AzothError, Result};

use crate::segment::phase::{system_mass_density_from_feed, system_volumetric_flow_from_feed};
use crate::stream::Stream;

/// `CapacityConstraint.MAX_UTILIZATION`: the ceiling on any utilization.
pub const MAX_UTILIZATION: f64 = 9.99;
/// `Splitter.designPressureDrop`, the class's `0.05` bar, in pascals.
pub const DEFAULT_DESIGN_PRESSURE_DROP_PA: f64 = 5000.0;
/// `Splitter.maxDesignVelocity`, m/s.
pub const DEFAULT_MAX_DESIGN_VELOCITY_M_S: f64 = 30.0;

/// Split a feed into several outlets with the same composition, pressure,
/// temperature and molar enthalpy, scaled by the given fractions.
///
/// The fractions are normalised, so `[0.3, 0.7]` and `[3.0, 7.0]` mean the same
/// split. Nothing else changes - a splitter does not flash, mix or lose pressure.
pub fn splitter(feed: &Stream, fractions: &[f64]) -> Result<Vec<Stream>> {
    if fractions.is_empty() {
        return Err(AzothError::invalid_input(
            "fractions",
            "a splitter needs at least one outlet",
        ));
    }
    let total: f64 = fractions.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "fractions",
            "the split fractions must sum to a positive value",
        ));
    }

    Ok(fractions
        .iter()
        .map(|&f| Stream {
            components: feed.components.clone(),
            z: feed.z.clone(),
            n: feed.n * f / total,
            p: feed.p,
            t: feed.t,
            h: feed.h,
            // A split moves no composition, so it moves no phase: each product is on the side its
            // feed was on, and a product of a feed nothing flashed is one nothing knows about.
            vapour_fraction: feed.vapour_fraction,
        })
        .collect())
}

/// `Splitter`'s capacity map and its three summaries, in SI.
///
/// The two constraints are the class's own: a `pressureDrop` whose value is always zero and a
/// `velocity` read from the inlet's system volume flow over the header area. The six per-constraint
/// fields are `None` when that constraint is not in the map, which is the class's lazy
/// initialization, not an absent answer.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SplitterCapacity {
    /// `pressureDrop`'s current value, Pa: always zero.
    pub pressure_drop_value: Option<Pressure>,
    /// `pressureDrop`'s design value, Pa.
    pub pressure_drop_design: Option<Pressure>,
    /// `pressureDrop`'s utilization.
    pub pressure_drop_utilization: Option<f64>,
    /// `velocity`'s current value, m/s.
    pub velocity_value: Option<Velocity>,
    /// `velocity`'s design value, m/s.
    pub velocity_design: Option<Velocity>,
    /// `velocity`'s utilization.
    pub velocity_utilization: Option<f64>,
    /// `getMaxUtilization`, 0.0 on an empty map.
    pub max_utilization: f64,
    /// `isCapacityExceeded`.
    pub capacity_exceeded: bool,
    /// `getBottleneckConstraint`'s class-side name.
    pub bottleneck_constraint: Option<&'static str>,
}

/// The capacity map `Splitter.splitterCapacityConstraints` holds, evaluated on `feed`.
///
/// **The map is lazy.** The class's constructor creates the mechanical design but no constraints;
/// `setDesignPressureDrop` and `setMaxDesignVelocity` each initialize it, and
/// `setCapacityAnalysisEnabled(true)` initializes it only when empty. An input that is `Some` is a
/// setter that was called, so it initializes the map; a non-positive limit suppresses its own
/// constraint but does not suppress the initialization.
///
/// # Errors
/// Whatever the feed's system volume flow raises - the capacity read flashes the feed for the
/// density `getFlowRate("m3/hr")` uses. The flash is skipped when the velocity constraint is
/// absent or the diameter is not positive, exactly the class's supplier guard.
pub fn capacity(
    feed: &Stream,
    design_pressure_drop: Option<Pressure>,
    max_design_velocity: Option<Velocity>,
    header_diameter: Option<Length>,
    analysis_enabled: bool,
) -> Result<SplitterCapacity> {
    if design_pressure_drop.is_none() && max_design_velocity.is_none() && !analysis_enabled {
        return Ok(SplitterCapacity::default());
    }

    let design_pressure_drop =
        design_pressure_drop.unwrap_or_else(|| pascals(DEFAULT_DESIGN_PRESSURE_DROP_PA));
    let max_design_velocity =
        max_design_velocity.unwrap_or_else(|| meters_per_second(DEFAULT_MAX_DESIGN_VELOCITY_M_S));

    let pressure_drop = (design_pressure_drop.value > 0.0).then(|| {
        let value = pascals(0.0);
        (
            value,
            design_pressure_drop,
            utilization(value.value, design_pressure_drop.value),
        )
    });

    let velocity = if max_design_velocity.value > 0.0 {
        let diameter = header_diameter.unwrap_or_else(|| meters(0.0));
        let value = if diameter.value > 0.0 {
            let volume_flow = system_volumetric_flow_from_feed(feed)?;
            let area = std::f64::consts::PI * diameter.value * diameter.value / 4.0;
            meters_per_second(volume_flow / area)
        } else {
            meters_per_second(0.0)
        };
        Some((
            value,
            max_design_velocity,
            utilization(value.value, max_design_velocity.value),
        ))
    } else {
        None
    };

    // `getBottleneckConstraint` walks the insertion-ordered map - pressureDrop first, then
    // velocity - and keeps a constraint only when its utilization is strictly above the running
    // maximum, which starts at zero. The always-zero pressure drop therefore never wins, and a
    // zero-utilization map has no bottleneck at all.
    let mut max_utilization = 0.0;
    let mut bottleneck_constraint = None;
    if let Some((_, _, value)) = pressure_drop {
        if value > max_utilization {
            max_utilization = value;
            bottleneck_constraint = Some("pressureDrop");
        }
    }
    if let Some((_, _, value)) = velocity {
        if value > max_utilization {
            max_utilization = value;
            bottleneck_constraint = Some("velocity");
        }
    }

    Ok(SplitterCapacity {
        pressure_drop_value: pressure_drop.map(|(value, _, _)| value),
        pressure_drop_design: pressure_drop.map(|(_, design, _)| design),
        pressure_drop_utilization: pressure_drop.map(|(_, _, value)| value),
        velocity_value: velocity.map(|(value, _, _)| value),
        velocity_design: velocity.map(|(_, design, _)| design),
        velocity_utilization: velocity.map(|(_, _, value)| value),
        max_utilization,
        capacity_exceeded: pressure_drop.is_some_and(|(_, _, value)| value > 1.0)
            || velocity.is_some_and(|(_, _, value)| value > 1.0),
        bottleneck_constraint,
    })
}

/// `CapacityConstraint.getUtilization` for one of the splitter's maximum constraints.
///
/// A design value that is not positive answers zero, as the class's does; the cap is the class's
/// own `MAX_UTILIZATION`.
fn utilization(value: f64, design: f64) -> f64 {
    if design <= 0.0 {
        0.0
    } else {
        (value / design).min(MAX_UTILIZATION)
    }
}

/// `SplitterMechanicalDesign.designVelocity`, m/s: the header and branch sizing velocity.
pub const DEFAULT_DESIGN_VELOCITY_M_S: f64 = 15.0;
/// `MechanicalDesign.maxOperationPressure`, the class's `100` bara, in pascals.
pub const DEFAULT_MAX_OPERATION_PRESSURE_PA: f64 = 100.0e5;
/// `MechanicalDesign.tensileStrength`, the class's `483` MPa, in pascals.
pub const DEFAULT_TENSILE_STRENGTH_PA: f64 = 483.0e6;
/// `MechanicalDesign.jointEfficiency`: read by the wall formula and then discarded.
pub const DEFAULT_JOINT_EFFICIENCY: f64 = 1.0;
/// `SplitterMechanicalDesign.calculateWeights`' steel density, kg/m3.
pub const STEEL_DENSITY_KG_M3: f64 = 7850.0;
/// The class's wall-thickness floor, mm.
pub const MIN_WALL_THICKNESS_MM: f64 = 3.0;
/// The class's branch-spacing floor, m.
pub const MIN_BRANCH_SPACING_M: f64 = 0.5;
/// The class's assumed branch length, m.
pub const BRANCH_LENGTH_M: f64 = 0.5;
/// `SplitterMechanicalDesign.calcDesign`'s splitting loss coefficient.
pub const SPLITTING_LOSS_COEFFICIENT: f64 = 0.5;
/// The B31.3 temperature coefficient the class hard-sets.
pub const WALL_Y_COEFFICIENT: f64 = 0.4;

/// The class's standard pipe sizes, in metres, from 1/2" to 48".
pub const STANDARD_PIPE_SIZES_M: [f64; 21] = [
    0.0127, 0.0254, 0.0381, 0.0508, 0.0635, 0.0762, 0.1016, 0.1270, 0.1524, 0.2032, 0.2540, 0.3048,
    0.3556, 0.4064, 0.4572, 0.5080, 0.6096, 0.7620, 0.9144, 1.0668, 1.2192,
];

/// The mechanical design's inputs, with `None` meaning the class's default.
///
/// The `Option`s are kept rather than resolved by the caller because their presence is itself
/// semantic: a design input stated without the flag is the separator's refused combination, and
/// the model id and the flowsheet dispatcher must reach that refusal through one function.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SplitterDesignRequest {
    /// Run `SplitterMechanicalDesign.calcDesign` at all.
    pub mechanical_design: bool,
    /// `setDesignVelocity`; absent is 15 m/s.
    pub design_velocity: Option<Velocity>,
    /// `setMaxOperationPressure`; absent is 100 bara.
    pub max_operation_pressure: Option<Pressure>,
    /// `setTensileStrength`; absent is 483 MPa.
    pub tensile_strength: Option<Pressure>,
    /// `setJointEfficiency`; absent is 1.0, and the wall formula discards it.
    pub joint_efficiency: Option<f64>,
    /// `setCorrosionAllowance`; absent is 0 mm.
    pub corrosion_allowance: Option<Length>,
}

/// `SplitterMechanicalDesign.calcDesign`'s answer, in SI.
///
/// The wall thickness is the class's mm converted to metres at this boundary; every other field
/// is the class's own SI quantity.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct SplitterDesign {
    /// `getHeaderDiameter`, m.
    pub header_diameter: f64,
    /// `getBranchDiameter`, m.
    pub branch_diameter: f64,
    /// `getHeaderLength`, m.
    pub header_length: f64,
    /// `getHeaderWallThickness` in mm, crossed in metres.
    pub wall_thickness: f64,
    /// `getManifoldVolume`, m3.
    pub manifold_volume: f64,
    /// `getTotalPressureDrop`, Pa.
    pub total_pressure_drop: f64,
    /// `weigthVesselShell`: the header's steel, kg.
    pub header_weight: f64,
    /// `weightPiping`: the branches' steel, kg.
    pub branch_weight: f64,
    /// `getWeightTotal`, kg.
    pub total_weight: f64,
}

/// `SplitterMechanicalDesign.calcDesign`, with the class's lazy call and its guards.
///
/// **The class does not call this from `run` and neither does the split.** The flag is the
/// caller's explicit request; with it off, any stated design input is refused by name, exactly
/// as the separator refuses a design parameter without its `mechanical_design` flag.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a design input without the flag, a stated header diameter
/// with it, a non-positive design velocity (the class divides by zero there) or a zero
/// wall-thickness denominator; and whatever the feed's system volume flow or density raises.
pub fn design(
    feed: &Stream,
    factors: &[f64],
    header_diameter: Option<Length>,
    request: &SplitterDesignRequest,
) -> Result<Option<SplitterDesign>> {
    if !request.mechanical_design {
        for (name, stated) in [
            ("design_velocity", request.design_velocity.is_some()),
            (
                "max_operation_pressure",
                request.max_operation_pressure.is_some(),
            ),
            ("tensile_strength", request.tensile_strength.is_some()),
            ("joint_efficiency", request.joint_efficiency.is_some()),
            ("corrosion_allowance", request.corrosion_allowance.is_some()),
        ] {
            if stated {
                return Err(AzothError::invalid_input(
                    name,
                    format!("`{name}` was stated without `mechanical_design`"),
                ));
            }
        }
        return Ok(None);
    }
    if header_diameter.is_some() {
        return Err(AzothError::invalid_input(
            "header_diameter",
            "`mechanical_design` sizes the header; state one or the other, not both",
        ));
    }

    let total: f64 = factors.iter().sum();
    if factors.is_empty() || total <= 0.0 {
        return Err(AzothError::invalid_input(
            "split_factors",
            "a splitter needs at least one outlet whose factors sum above zero",
        ));
    }

    let design_velocity = request
        .design_velocity
        .unwrap_or_else(|| meters_per_second(DEFAULT_DESIGN_VELOCITY_M_S));
    let max_operation_pressure = request
        .max_operation_pressure
        .unwrap_or_else(|| pascals(DEFAULT_MAX_OPERATION_PRESSURE_PA));
    let tensile_strength_mpa = request
        .tensile_strength
        .map_or(DEFAULT_TENSILE_STRENGTH_PA, |p| p.value)
        / 1.0e6;
    let joint_efficiency = request.joint_efficiency.unwrap_or(DEFAULT_JOINT_EFFICIENCY);
    let corrosion_allowance_mm = request
        .corrosion_allowance
        .map_or(0.0, |length| length.value)
        * 1000.0;

    let volume_flow = system_volumetric_flow_from_feed(feed)?;
    if volume_flow <= 0.0 || !volume_flow.is_finite() {
        // `calcDesign` returns before it divides or sizes anything; the class's fields stay at
        // their initial zeros and `toJson` reports them. They are answers, not absent ones.
        return Ok(Some(SplitterDesign::default()));
    }

    if design_velocity.value <= 0.0 || !design_velocity.value.is_finite() {
        return Err(AzothError::invalid_input(
            "design_velocity",
            "the header area divides by the design velocity, so it must be positive",
        ));
    }

    // The class's formula, in the class's units: P stays in bar while S is MPa. It is not a
    // dimensionally consistent B31.3 form, and reproducing it is the point.
    let design_pressure_bar = max_operation_pressure.value / 1.0e5 * 1.1;
    let allowable_stress_mpa = tensile_strength_mpa * 0.4;
    let denominator = 2.0 * allowable_stress_mpa + 2.0 * WALL_Y_COEFFICIENT * design_pressure_bar;
    if denominator <= 0.0 || !denominator.is_finite() {
        return Err(AzothError::invalid_input(
            "tensile_strength",
            "the wall-thickness denominator is zero when the tensile strength and the design \
             pressure are both zero",
        ));
    }
    // Read and discarded, as the class does; the binding keeps the read visible.
    let _ = joint_efficiency;

    let branches = factors.len() as f64;
    let max_ratio = factors.iter().map(|f| f / total).fold(0.0_f64, f64::max);

    let header_area = volume_flow / design_velocity.value;
    let header_diameter =
        round_to_standard_pipe_size((4.0 * header_area / std::f64::consts::PI).sqrt());
    let branch_area = volume_flow * max_ratio / design_velocity.value;
    let branch_diameter =
        round_to_standard_pipe_size((4.0 * branch_area / std::f64::consts::PI).sqrt());

    let branch_spacing = (2.0 * branch_diameter).max(MIN_BRANCH_SPACING_M);
    let header_length = (branches + 1.0) * branch_spacing;

    let wall_thickness_mm = ((design_pressure_bar * header_diameter * 1000.0) / denominator
        + corrosion_allowance_mm)
        .max(MIN_WALL_THICKNESS_MM);

    let header_section = std::f64::consts::PI * (header_diameter / 2.0).powi(2);
    let manifold_volume = header_section * header_length;

    let velocity = volume_flow / header_section;
    let density = system_mass_density_from_feed(feed)?;
    let total_pressure_drop = SPLITTING_LOSS_COEFFICIENT * 0.5 * density * velocity * velocity;

    let header_outer = header_diameter + 2.0 * wall_thickness_mm / 1000.0;
    let header_steel_area =
        std::f64::consts::PI / 4.0 * (header_outer.powi(2) - header_diameter.powi(2));
    let header_weight = header_steel_area * header_length * STEEL_DENSITY_KG_M3;

    let branch_outer = branch_diameter + 2.0 * wall_thickness_mm / 1000.0;
    let branch_steel_area =
        std::f64::consts::PI / 4.0 * (branch_outer.powi(2) - branch_diameter.powi(2));
    let branch_weight = branch_steel_area * BRANCH_LENGTH_M * STEEL_DENSITY_KG_M3 * branches;

    Ok(Some(SplitterDesign {
        header_diameter,
        branch_diameter,
        header_length,
        wall_thickness: wall_thickness_mm / 1000.0,
        manifold_volume,
        total_pressure_drop,
        header_weight,
        branch_weight,
        total_weight: header_weight + branch_weight,
    }))
}

/// `SplitterMechanicalDesign.roundToStandardPipeSize`: the first standard size at or above the
/// asked diameter, and 48" when it is above the whole list.
fn round_to_standard_pipe_size(diameter: f64) -> f64 {
    STANDARD_PIPE_SIZES_M
        .iter()
        .copied()
        .find(|&size| size >= diameter)
        .unwrap_or(STANDARD_PIPE_SIZES_M[STANDARD_PIPE_SIZES_M.len() - 1])
}

#[cfg(test)]
mod design_tests {
    use super::*;
    use azoth_core::units::kelvins;

    fn feed(n: f64) -> Stream {
        Stream::from_pt(
            vec!["n-butane".into(), "n-pentane".into()],
            vec![0.6, 0.4],
            n,
            pascals(1.0e6),
            kelvins(300.0),
        )
        .expect("the butane / pentane feed resolves and flashes")
    }

    fn requested() -> SplitterDesignRequest {
        SplitterDesignRequest {
            mechanical_design: true,
            ..SplitterDesignRequest::default()
        }
    }

    fn sized(n: f64, request: &SplitterDesignRequest) -> SplitterDesign {
        design(&feed(n), &[0.3, 0.7], None, request)
            .expect("the design runs")
            .expect("the flag is on")
    }

    #[test]
    fn the_nps_list_rounds_up_and_stops_at_48_inches() {
        assert_eq!(round_to_standard_pipe_size(0.001), 0.0127);
        assert_eq!(round_to_standard_pipe_size(0.0127), 0.0127);
        assert_eq!(round_to_standard_pipe_size(0.0128), 0.0254);
        assert_eq!(round_to_standard_pipe_size(1.3), 1.2192);
        assert_eq!(round_to_standard_pipe_size(f64::INFINITY), 1.2192);
    }

    #[test]
    fn the_default_design_sizes_the_header_and_the_branch() {
        let design = sized(1.0, &requested());
        assert_eq!(design.header_diameter, 0.0127);
        assert_eq!(design.branch_diameter, 0.0127);
        assert_eq!(design.header_length, 1.5);
        assert_eq!(design.wall_thickness, 0.003);
        assert!((design.manifold_volume - 1.9001e-4).abs() < 1e-6);
        assert!(design.total_pressure_drop > 0.0);
        assert!(design.header_weight > 0.0);
        assert_eq!(
            design.total_weight,
            design.header_weight + design.branch_weight
        );
    }

    #[test]
    fn a_larger_flow_grows_the_header_and_the_wall_past_the_floor() {
        let design = sized(20.0, &requested());
        assert_eq!(design.header_diameter, 0.0254);
        assert_eq!(design.branch_diameter, 0.0127);
        assert_eq!(design.header_length, 1.5);
        assert!((design.wall_thickness - 5.8894e-3).abs() < 1e-6);
    }

    #[test]
    fn the_pressure_tensile_and_corrosion_inputs_move_the_wall() {
        let pressure = SplitterDesignRequest {
            max_operation_pressure: Some(pascals(300.0e5)),
            ..requested()
        };
        assert!((sized(20.0, &pressure).wall_thickness - 12.8875e-3).abs() < 1e-6);

        let tensile = SplitterDesignRequest {
            tensile_strength: Some(pascals(300.0e6)),
            ..requested()
        };
        assert!((sized(20.0, &tensile).wall_thickness - 8.5183e-3).abs() < 1e-6);

        let corrosion = SplitterDesignRequest {
            corrosion_allowance: Some(meters(0.0015)),
            ..requested()
        };
        assert!((sized(20.0, &corrosion).wall_thickness - 7.3894e-3).abs() < 1e-6);
    }

    #[test]
    fn the_joint_efficiency_is_read_and_discarded() {
        let default = sized(20.0, &requested());
        let loose = SplitterDesignRequest {
            joint_efficiency: Some(0.7),
            ..requested()
        };
        assert_eq!(sized(20.0, &loose), default);
    }

    #[test]
    fn a_zero_flow_returns_the_class_s_zeros() {
        let design = sized(0.0, &requested());
        assert_eq!(design, SplitterDesign::default());
    }

    #[test]
    fn a_design_input_without_the_flag_is_refused() {
        let error = design(
            &feed(1.0),
            &[0.3, 0.7],
            None,
            &SplitterDesignRequest {
                mechanical_design: false,
                design_velocity: Some(meters_per_second(10.0)),
                ..SplitterDesignRequest::default()
            },
        )
        .expect_err("a stated design input without the flag");
        assert!(error.to_string().contains("design_velocity"));
    }

    #[test]
    fn a_stated_header_is_refused_with_the_flag() {
        let error = design(&feed(1.0), &[0.3, 0.7], Some(meters(0.02)), &requested())
            .expect_err("the design sizes the header");
        assert!(error.to_string().contains("header_diameter"));
    }

    #[test]
    fn a_non_positive_design_velocity_is_refused() {
        let error = design(
            &feed(1.0),
            &[0.3, 0.7],
            None,
            &SplitterDesignRequest {
                design_velocity: Some(meters_per_second(0.0)),
                ..requested()
            },
        )
        .expect_err("the header area divides by it");
        assert!(error.to_string().contains("design_velocity"));
    }

    #[test]
    fn a_zero_wall_denominator_is_refused() {
        let error = design(
            &feed(1.0),
            &[0.3, 0.7],
            None,
            &SplitterDesignRequest {
                max_operation_pressure: Some(pascals(0.0)),
                tensile_strength: Some(pascals(0.0)),
                ..requested()
            },
        )
        .expect_err("the wall formula divides by it");
        assert!(error.to_string().contains("tensile_strength"));
    }
}
