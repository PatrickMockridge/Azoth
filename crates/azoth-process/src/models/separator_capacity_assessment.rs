//! `process.separator_capacity_assessment` - the separator's detached capacity report.
//!
//! Spec: `specs/models/process/separator_capacity_assessment.toml`. The arithmetic is
//! [`crate::kernels::separator::capacity_assessment`], a fold over the same run and capacity
//! surface `process.separator` carries.
//!
//! **A non-palette direct model.** It has no ports and no flowsheet entry: it is reached by
//! its case and its Python wrapper, like the release models, and sizes nothing.

use azoth_core::units::{
    Area, Length, MassDensity, Power, Pressure, ThermodynamicTemperature, Velocity, VolumeRate,
    cubic_meters_per_second, kilograms_per_cubic_meter, meters_per_second, square_meters,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{optional_scalar, warnings as wire_warnings};
use crate::kernels::separator as kernel;
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.separator_capacity_assessment`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SeparatorCapacityAssessmentResult {
    /// `feasible`, `infeasible` or `unassessed`.
    pub status: String,
    /// The class's sentence for the status.
    pub diagnostic: String,
    /// `equipment` or `mechanicalDesign`; empty on `unassessed`.
    pub design_basis_source: String,
    /// `soudersBrown` or `mechanicalDesignGasFlow`; empty on `unassessed`.
    pub governing_capacity_basis: String,
    /// The design gas area.
    #[serde(serialize_with = "optional_scalar")]
    pub gas_area: Option<Area>,
    /// The gas phase's mass density.
    #[serde(serialize_with = "optional_scalar")]
    pub gas_density: Option<MassDensity>,
    /// The sizing family's liquid density.
    #[serde(serialize_with = "optional_scalar")]
    pub liquid_density: Option<MassDensity>,
    /// `getMaxAllowableGasVelocity`.
    #[serde(serialize_with = "optional_scalar")]
    pub maximum_gas_velocity: Option<Velocity>,
    /// The gas phase's volumetric flow.
    #[serde(serialize_with = "optional_scalar")]
    pub current_gas_flow: Option<VolumeRate>,
    /// The stated mechanical limit, or zero with none.
    #[serde(serialize_with = "optional_scalar")]
    pub mechanical_design_maximum_gas_flow: Option<VolumeRate>,
    /// `getMaxAllowableGasFlowRate`.
    #[serde(serialize_with = "optional_scalar")]
    pub souders_brown_maximum_gas_flow: Option<VolumeRate>,
    /// The governed limit.
    #[serde(serialize_with = "optional_scalar")]
    pub maximum_gas_flow: Option<VolumeRate>,
    /// Current flow over the governed limit.
    pub utilization: Option<f64>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// Assess a separator's gas capacity, failing closed on incomplete evidence.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] where the pressure drop leaves a non-positive
/// pressure, the entrainment fraction is outside `[0, 1]`, or the flash refuses. The
/// assessment's own failing states are `unassessed` results, not errors.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are fourteen
pub fn separator_capacity_assessment(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    pressure_drop: Pressure,
    gas_in_liquid: f64,
    internal_diameter: Length,
    design_gas_load_factor: Velocity,
    orientation: &str,
    design_liquid_level_fraction: f64,
    design_basis_provenance: &str,
    heat_input: Option<Power>,
    mechanical_design_max_gas_flow: Option<VolumeRate>,
) -> Result<SeparatorCapacityAssessmentResult> {
    let spec = &model_gen::SEPARATOR_CAPACITY_ASSESSMENT_SPEC;
    let mut warnings = Vec::new();
    // The spec declares no `valid_range` rows: the class's failing states are statuses, not
    // typed refusals. The resolver is written out so a bound added to the spec is applied
    // rather than answered with a skipped-check warning.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "internal_diameter" => Some(internal_diameter.value),
            "design_gas_load_factor" => Some(design_gas_load_factor.value),
            "design_liquid_level_fraction" => Some(design_liquid_level_fraction),
            "mechanical_design_max_gas_flow" => {
                mechanical_design_max_gas_flow.map(|value| value.value)
            }
            _ => None,
        },
        &mut warnings,
    )?;

    if feed_n <= 0.0 {
        return Err(AzothError::out_of_range(
            "feed_n",
            feed_n,
            "a unit operation needs a positive molar flow",
        ));
    }
    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let (vapour, liquid) = kernel::separator(&feed, pressure_drop, gas_in_liquid, heat_input)?;
    let assessment = kernel::capacity_assessment(
        &feed,
        &vapour,
        &liquid,
        internal_diameter,
        design_gas_load_factor,
        orientation,
        design_liquid_level_fraction,
        design_basis_provenance,
        mechanical_design_max_gas_flow,
    )?;

    Ok(SeparatorCapacityAssessmentResult {
        status: assessment.status.to_string(),
        diagnostic: assessment.diagnostic.to_string(),
        design_basis_source: assessment.design_basis_source.to_string(),
        governing_capacity_basis: assessment.governing_capacity_basis.to_string(),
        gas_area: assessment.gas_area.map(square_meters),
        gas_density: assessment.gas_density.map(kilograms_per_cubic_meter),
        liquid_density: assessment.liquid_density.map(kilograms_per_cubic_meter),
        maximum_gas_velocity: assessment.maximum_gas_velocity.map(meters_per_second),
        current_gas_flow: assessment.current_gas_flow.map(cubic_meters_per_second),
        mechanical_design_maximum_gas_flow: assessment
            .mechanical_design_maximum_gas_flow
            .map(cubic_meters_per_second),
        souders_brown_maximum_gas_flow: assessment
            .souders_brown_maximum_gas_flow
            .map(cubic_meters_per_second),
        maximum_gas_flow: assessment.maximum_gas_flow.map(cubic_meters_per_second),
        utilization: assessment.utilization,
        warnings,
    })
}
