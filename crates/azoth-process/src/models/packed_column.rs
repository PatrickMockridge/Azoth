//! `process.packed_column` - the packed column as a registered id.
//!
//! Spec: `specs/models/process/packed_column.toml`. **The arithmetic is
//! [`crate::models::distillation_column`]'s**, because the class's is: `PackedColumn extends
//! `DistillationColumn`, its `run` is `super.run(id)` and the packing is read by a hydraulics
//! report afterwards. What this module adds is the stage count the class's constructor derives
//! from a packed height, and the refusal of the one packing parameter the class itself refuses.

use azoth_core::units::{MolarEnergy, Power, Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, CalcResult, Result, Warning};

use crate::kernels::packed_column::stage_count;
use crate::models::distillation_column::{DistillationColumnResult, distillation_column};

/// Result of `process.packed_column`.
///
/// **The same record as [`DistillationColumnResult`], and that is the id's claim rather than a
/// convenience**: the packing parameters reach the separation only through the stage count, and
/// the three quantities a caller would look for here instead - HETP, the theoretical stages and
/// the percent flood - are `ColumnInternalsDesigner`'s report on the far side of the solve and
/// are not ported.
#[derive(Debug, Clone, PartialEq)]
pub struct PackedColumnResult {
    /// Each tray's temperature, K.
    pub tray_temperature: Vec<ThermodynamicTemperature>,
    /// Each tray's pressure, Pa.
    pub tray_pressure: Vec<Pressure>,
    /// Each tray's vapour traffic, mol/s.
    pub tray_gas_n: Vec<f64>,
    /// Each tray's liquid traffic, mol/s.
    pub tray_liquid_n: Vec<f64>,
    /// Distillate molar flow, mol/s.
    pub distillate_n: f64,
    /// Distillate composition.
    pub distillate_z: Vec<f64>,
    /// Distillate pressure.
    pub distillate_p: Pressure,
    /// Distillate temperature.
    pub distillate_t: ThermodynamicTemperature,
    /// Distillate molar enthalpy.
    pub distillate_h: MolarEnergy,
    /// Bottoms molar flow, mol/s.
    pub bottoms_n: f64,
    /// Bottoms composition.
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    pub bottoms_p: Pressure,
    /// Bottoms temperature.
    pub bottoms_t: ThermodynamicTemperature,
    /// Bottoms molar enthalpy.
    pub bottoms_h: MolarEnergy,
    /// **The vapour each tray withdrew**, which is the base column's and empty for this id: the
    /// packing does not change a draw, and its fractions are not declared here.
    pub gas_side_draw_n: Vec<f64>,
    /// The liquid each tray withdrew as a liquid side draw.
    pub liquid_side_draw_n: Vec<f64>,
    /// The liquid each tray withdrew as a pumparound.
    pub pumparound_n: Vec<f64>,
    /// The condenser's duty, W.
    pub condenser_duty: Power,
    /// The reboiler's duty, W.
    pub reboiler_duty: Power,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against the feed, relative.
    pub mass_residual: f64,
    /// The enthalpy closure.
    pub energy_residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl From<DistillationColumnResult> for PackedColumnResult {
    fn from(out: DistillationColumnResult) -> Self {
        Self {
            tray_temperature: out.tray_temperature,
            tray_pressure: out.tray_pressure,
            tray_gas_n: out.tray_gas_n,
            tray_liquid_n: out.tray_liquid_n,
            distillate_n: out.distillate_n,
            distillate_z: out.distillate_z,
            distillate_p: out.distillate_p,
            distillate_t: out.distillate_t,
            distillate_h: out.distillate_h,
            bottoms_n: out.bottoms_n,
            bottoms_z: out.bottoms_z,
            bottoms_p: out.bottoms_p,
            bottoms_t: out.bottoms_t,
            bottoms_h: out.bottoms_h,
            gas_side_draw_n: out.gas_side_draw_n,
            liquid_side_draw_n: out.liquid_side_draw_n,
            pumparound_n: out.pumparound_n,
            condenser_duty: out.condenser_duty,
            reboiler_duty: out.reboiler_duty,
            iterations: out.iterations,
            temperature_residual: out.temperature_residual,
            mass_residual: out.mass_residual,
            energy_residual: out.energy_residual,
            warnings: out.warnings,
        }
    }
}

impl CalcResult for PackedColumnResult {
    const CALC_ID: &'static str = "process.packed_column";
    const FIELDS: &'static [&'static str] = &[
        "tray_temperature",
        "tray_pressure",
        "tray_gas_n",
        "tray_liquid_n",
        "distillate_n",
        "distillate_z",
        "distillate_p",
        "distillate_t",
        "distillate_h",
        "bottoms_n",
        "bottoms_z",
        "bottoms_p",
        "bottoms_t",
        "bottoms_h",
        "gas_side_draw_n",
        "liquid_side_draw_n",
        "pumparound_n",
        "condenser_duty",
        "reboiler_duty",
        "iterations",
        "temperature_residual",
        "mass_residual",
        "energy_residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Solve a packed column.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a `packing_hydraulic_capacity_factor` that is not positive
/// and finite, which is `setPackingHydraulicCapacityFactor`'s own check and the one parameter of
/// the packing group the class refuses rather than reads; and every error
/// [`distillation_column`] raises, which this id inherits because the solve is the base's.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-six
pub fn packed_column(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    packed_height: f64,
    feed_stage: usize,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Pressure,
    bottom_pressure: Pressure,
    reboiler_temperature: Option<ThermodynamicTemperature>,
    condenser_temperature: Option<ThermodynamicTemperature>,
    temperature_tolerance: f64,
    max_iterations: usize,
    packing_type: Option<&str>,
    structured_packing: Option<bool>,
    design_flood_fraction: Option<f64>,
    packing_hydraulic_capacity_factor: Option<f64>,
    column_diameter: Option<f64>,
    murphree_efficiency: Option<f64>,
    solver_type: Option<&str>,
    top_specification_type: Option<&str>,
    top_specification_target: Option<f64>,
    top_specification_component: Option<&str>,
    bottom_specification_type: Option<&str>,
    bottom_specification_target: Option<f64>,
    bottom_specification_component: Option<&str>,
) -> Result<PackedColumnResult> {
    // **The class's own check, on the one packing parameter it does not simply read.**
    // `setPackingHydraulicCapacityFactor` throws for a value that is not positive and finite, so
    // the class refuses a factor its own setter would reject rather than carrying it inertly as
    // it carries the other four.
    if let Some(factor) = packing_hydraulic_capacity_factor {
        if !factor.is_finite() || factor <= 0.0 {
            return Err(AzothError::invalid_input(
                "packing_hydraulic_capacity_factor",
                format!(
                    "a packing hydraulic capacity factor of {factor} is not positive and finite, \
                     which `PackedColumn.setPackingHydraulicCapacityFactor` refuses with an \
                     `IllegalArgumentException`"
                ),
            ));
        }
    }
    // **The other four are declarations the solve is indifferent to.** The class accepts them
    // and every one of them is read by `ColumnInternalsDesigner` after the column has
    // converged, so the model carries them rather than reading them - the shape
    // `process.absorption_column` gives `max_allowable_gas_load_factor`.
    let _ = (
        packing_type,
        structured_packing,
        design_flood_fraction,
        column_diameter,
    );

    let stages = stage_count(packed_height);
    let out = distillation_column(
        components,
        feed_n,
        feed_z,
        feed_p,
        feed_t,
        stages,
        feed_stage,
        has_reboiler,
        has_condenser,
        top_pressure,
        bottom_pressure,
        reboiler_temperature,
        condenser_temperature,
        temperature_tolerance,
        max_iterations,
        murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        // **A packed column's stages are equilibrium stages**: `PackedColumn` inherits
        // `setReactive`, and this id does not declare the section yet - the packing is a report
        // on the far side of the solve, and the reactive section is the base column's own.
        None,
        None,
        None,
        // **A packed column's draws are the base column's, and this id does not declare them.**
        None,
        None,
        None,
    )?;

    Ok(out.into())
}
