//! `process.splitter` - the splitter's kernel as a registered id.
//!
//! Spec: `specs/models/process/splitter.toml`. The arithmetic is
//! [`crate::kernels::splitter`]; what is here is the boundary a case and a cross-impl
//! test address.
//!
//! **This is the first id with a `many` outlet port**, so its result is the shape the
//! module doc describes: one field per record field, each carrying one entry per outlet,
//! and `z` a matrix with one row per outlet.

use azoth_core::units::{
    Length, Mass, MolarEnergy, Pressure, ThermodynamicTemperature, Velocity, Volume, cubic_meters,
    joules_per_mole, kilograms, meters, pascals,
};
use azoth_core::{Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, scalars, warnings as wire_warnings};
use crate::kernels::splitter as kernel;
use crate::kernels::splitter::{SplitterCapacity, SplitterDesign, SplitterDesignRequest};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.splitter`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SplitterResult {
    /// Molar flow of each outlet, mol/s.
    pub products_n: Vec<f64>,
    /// Composition of each outlet, one row per outlet.
    pub products_z: Vec<Vec<f64>>,
    /// Pressure of each outlet, which a splitter does not change.
    #[serde(serialize_with = "scalars")]
    pub products_p: Vec<Pressure>,
    /// Temperature of each outlet, which a splitter does not change.
    #[serde(serialize_with = "scalars")]
    pub products_t: Vec<ThermodynamicTemperature>,
    /// Molar enthalpy of each outlet, which is the feed's.
    #[serde(serialize_with = "scalars")]
    pub products_h: Vec<MolarEnergy>,
    /// The `pressureDrop` constraint's current value, Pa - always zero.
    #[serde(serialize_with = "optional_scalar")]
    pub pressure_drop_value: Option<Pressure>,
    /// The `pressureDrop` constraint's design value, Pa.
    #[serde(serialize_with = "optional_scalar")]
    pub pressure_drop_design: Option<Pressure>,
    /// The `pressureDrop` constraint's utilization.
    pub pressure_drop_utilization: Option<f64>,
    /// The `velocity` constraint's current value, m/s.
    #[serde(serialize_with = "optional_scalar")]
    pub velocity_value: Option<Velocity>,
    /// The `velocity` constraint's design value, m/s.
    #[serde(serialize_with = "optional_scalar")]
    pub velocity_design: Option<Velocity>,
    /// The `velocity` constraint's utilization.
    pub velocity_utilization: Option<f64>,
    /// `getMaxUtilization`, 0.0 on an empty map.
    pub max_utilization: f64,
    /// `isCapacityExceeded`.
    pub capacity_exceeded: bool,
    /// `getBottleneckConstraint`'s class-side name.
    pub bottleneck_constraint: Option<String>,
    /// `calcDesign`'s NPS-rounded header diameter, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_header_diameter: Option<Length>,
    /// `calcDesign`'s NPS-rounded branch diameter, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_branch_diameter: Option<Length>,
    /// `calcDesign`'s header length, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_header_length: Option<Length>,
    /// `calcDesign`'s wall thickness, m.
    #[serde(serialize_with = "optional_scalar")]
    pub design_wall_thickness: Option<Length>,
    /// `calcDesign`'s manifold volume, m3.
    #[serde(serialize_with = "optional_scalar")]
    pub manifold_volume: Option<Volume>,
    /// `calcDesign`'s total pressure drop, Pa.
    #[serde(serialize_with = "optional_scalar")]
    pub design_total_pressure_drop: Option<Pressure>,
    /// `calcDesign`'s header steel weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub header_weight: Option<Mass>,
    /// `calcDesign`'s branch steel weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub branch_weight: Option<Mass>,
    /// `calcDesign`'s total steel weight, kg.
    #[serde(serialize_with = "optional_scalar")]
    pub total_weight: Option<Mass>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl SplitterResult {
    /// The record from the outlets, the capacity report and the mechanical design.
    pub(crate) fn of(
        out: &[Stream],
        capacity: SplitterCapacity,
        design: Option<SplitterDesign>,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            products_n: out.iter().map(|stream| stream.n).collect(),
            products_z: out.iter().map(|stream| stream.z.clone()).collect(),
            products_p: out.iter().map(|stream| stream.p).collect(),
            products_t: out.iter().map(|stream| stream.t).collect(),
            products_h: out
                .iter()
                .map(|stream| joules_per_mole(stream.h.value))
                .collect(),
            pressure_drop_value: capacity.pressure_drop_value,
            pressure_drop_design: capacity.pressure_drop_design,
            pressure_drop_utilization: capacity.pressure_drop_utilization,
            velocity_value: capacity.velocity_value,
            velocity_design: capacity.velocity_design,
            velocity_utilization: capacity.velocity_utilization,
            max_utilization: capacity.max_utilization,
            capacity_exceeded: capacity.capacity_exceeded,
            bottleneck_constraint: capacity.bottleneck_constraint.map(str::to_owned),
            design_header_diameter: design.map(|d| meters(d.header_diameter)),
            design_branch_diameter: design.map(|d| meters(d.branch_diameter)),
            design_header_length: design.map(|d| meters(d.header_length)),
            design_wall_thickness: design.map(|d| meters(d.wall_thickness)),
            manifold_volume: design.map(|d| cubic_meters(d.manifold_volume)),
            design_total_pressure_drop: design.map(|d| pascals(d.total_pressure_drop)),
            header_weight: design.map(|d| kilograms(d.header_weight)),
            branch_weight: design.map(|d| kilograms(d.branch_weight)),
            total_weight: design.map(|d| kilograms(d.total_weight)),
            warnings,
        }
    }
}

/// Split a stream into several with the same state, in proportion to `split_factors`.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] if there is no outlet or the factors do not
/// sum to a positive value, and whatever the databank or the capacity read refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are sixteen
pub fn splitter(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    split_factors: &[f64],
    design_pressure_drop: Option<Pressure>,
    max_design_velocity: Option<Velocity>,
    header_diameter: Option<Length>,
    capacity_analysis_enabled: Option<bool>,
    mechanical_design: Option<bool>,
    design_velocity: Option<Velocity>,
    max_operation_pressure: Option<Pressure>,
    tensile_strength: Option<Pressure>,
    joint_efficiency: Option<f64>,
    corrosion_allowance: Option<Length>,
) -> Result<SplitterResult> {
    let spec = &model_gen::SPLITTER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_t" => Some(feed_t.value),
            "design_velocity" => design_velocity.map(|velocity| velocity.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let out = kernel(&feed, split_factors)?;
    let request = SplitterDesignRequest {
        mechanical_design: mechanical_design.unwrap_or(false),
        design_velocity,
        max_operation_pressure,
        tensile_strength,
        joint_efficiency,
        corrosion_allowance,
    };
    let design = kernel::design(&feed, split_factors, header_diameter, &request)?;
    // With the design on, the capacity map reads the sized diameter live, as the class's
    // `getHeaderDiameter` supplier does after `calcDesign`; with it off, the stated input stands.
    let effective_header_diameter =
        design.map_or(header_diameter, |d| Some(meters(d.header_diameter)));
    let capacity = kernel::capacity(
        &feed,
        design_pressure_drop,
        max_design_velocity,
        effective_header_diameter,
        capacity_analysis_enabled.unwrap_or(false),
    )?;

    Ok(SplitterResult::of(&out, capacity, design, warnings))
}
