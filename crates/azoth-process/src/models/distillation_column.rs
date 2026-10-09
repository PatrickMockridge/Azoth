//! `process.distillation_column` - the column solve as a registered id.
//!
//! Spec: `specs/models/process/distillation_column.toml`. The arithmetic is
//! [`crate::kernels::distillation_column`]; what is here is the boundary a case and a
//! cross-impl test address, and the refusal of every parameter the palette declares and this
//! tranche does not implement.

use azoth_core::units::{
    Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole, meters,
    millimeters, pascals, watts,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::column::capacity::{
    DEFAULT_INTERNAL_DIAMETER_M, DEFAULT_MAX_ALLOWABLE_FS_FACTOR, FsLimits, fs_limits,
};
use crate::column::designer::{
    DesignerGeometry, DesignerReport, UNSIZED_COLUMN_DIAMETER_M, designer_report,
};
use crate::column::mechanical::{
    DEFAULT_CONTACTOR_INTERNALS_TYPE, DEFAULT_MATERIAL_GRADE, DEFAULT_MAX_FLOODING_FACTOR,
    DEFAULT_MAX_OPERATION_PRESSURE_BARA, DEFAULT_TRAY_EFFICIENCY, DEFAULT_TRAY_TYPE,
    MechanicalGeometry, MechanicalReport, mechanical_design,
};
use crate::column::murphree::Murphree;
use crate::executor::json::{scalar, scalars, warnings as wire_warnings};
use crate::kernels::distillation_column as kernel;
use crate::kernels::distillation_column::{Specification, SpecificationKind};
use crate::stream::Stream;
use crate::unported;

/// **The internals geometry's own defaults, re-stated from the designer's module** so a caller
/// of this id has one import: they are `ColumnInternalsDesigner`'s initialisers, none is read by
/// the solve, and [`crate::column::designer`] is where they are defined.
use crate::column::coupling::HydraulicCoupling;
pub use crate::column::designer::{
    DEFAULT_DESIGNER_FLOOD_FRACTION, DEFAULT_DOWNCOMMER_AREA_FRACTION, DEFAULT_HOLE_AREA_FRACTION,
    DEFAULT_HOLE_DIAMETER_MM, DEFAULT_INTERNALS_TYPE, DEFAULT_TRAY_SPACING_M,
    DEFAULT_WEIR_HEIGHT_M,
};

/// Result of `process.distillation_column`.
///
/// **The profile is a set of vectors, one entry per tray**, which is what a column's answer is
/// rather than a scalar: four quantities over the trays, then the two products' records, the
/// two duties and the three residuals. The vector lengths are the tray count, which the inputs
/// decide - a shape `eos.pt_phase_envelope` already exercises with its trace points.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct DistillationColumnResult {
    /// Each tray's temperature, K.
    #[serde(serialize_with = "scalars")]
    pub tray_temperature: Vec<ThermodynamicTemperature>,
    /// Each tray's pressure, Pa.
    #[serde(serialize_with = "scalars")]
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
    #[serde(serialize_with = "scalar")]
    pub distillate_p: Pressure,
    /// Distillate temperature.
    #[serde(serialize_with = "scalar")]
    pub distillate_t: ThermodynamicTemperature,
    /// Distillate molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub distillate_h: MolarEnergy,
    /// Bottoms molar flow, mol/s.
    pub bottoms_n: f64,
    /// Bottoms composition.
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    #[serde(serialize_with = "scalar")]
    pub bottoms_p: Pressure,
    /// Bottoms temperature.
    #[serde(serialize_with = "scalar")]
    pub bottoms_t: ThermodynamicTemperature,
    /// Bottoms molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub bottoms_h: MolarEnergy,
    /// **The vapour each tray withdrew**, one entry per tray and zero where it drew none.
    pub gas_side_draw_n: Vec<f64>,
    /// The liquid each tray withdrew as a liquid side draw.
    pub liquid_side_draw_n: Vec<f64>,
    /// The liquid each tray withdrew as a pumparound.
    pub pumparound_n: Vec<f64>,
    /// The condenser's duty, W.
    #[serde(serialize_with = "scalar")]
    pub condenser_duty: Power,
    /// The reboiler's duty, W.
    #[serde(serialize_with = "scalar")]
    pub reboiler_duty: Power,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against the feed, relative.
    pub mass_residual: f64,
    /// The enthalpy closure.
    pub energy_residual: f64,
    /// `getFsFactor`: the gas load factor over the **total** cross-section, `u·sqrt(rho_g)`.
    pub fs_factor: f64,
    /// `getFsFactorUtilization`: the factor over `max_allowable_fs_factor`.
    pub fs_factor_utilization: f64,
    /// `isFsFactorWithinDesignLimit`.
    pub fs_factor_within_design_limit: bool,
    /// `getMinimumDiameterForFsLimit`: the diameter at which the factor would reach the limit.
    #[serde(serialize_with = "scalar")]
    pub minimum_diameter_for_fs_limit: Length,
    /// **`ColumnInternalsDesigner`'s own report**, which is the trayed internals tree: the
    /// diameter it resolved, the tray it sized from, and every tray's load and efficiency at that
    /// diameter. The first three are read *after* the solve, from the trays it left behind.
    #[serde(serialize_with = "scalar")]
    pub required_diameter: Length,
    /// The tray the diameter was sized from: the largest vapour mass flow, the first of equals.
    pub controlling_tray_index: usize,
    /// Every tray's own verdict, and-ed together.
    pub internals_design_ok: bool,
    /// The largest per-tray load, in per cent of flood.
    pub max_percent_flood: f64,
    /// The smallest load above zero, which is not the smallest: an exact zero is skipped.
    pub min_percent_flood: f64,
    /// The mean of the per-tray efficiencies, or the class's own `0.65` with no tray walked.
    pub average_tray_efficiency: f64,
    /// The trays' pressure drops, summed.
    #[serde(serialize_with = "scalar")]
    pub total_pressure_drop: Pressure,
    /// The same sum in millibars, which the class publishes beside it.
    pub total_pressure_drop_mbar: f64,
    /// Each tray's load, in per cent of flood, at the diameter the sizing resolved.
    pub tray_percent_flood: Vec<f64>,
    /// Each tray's own total pressure drop, Pa - the entries `total_pressure_drop` sums.
    pub tray_pressure_drop: Vec<f64>,
    /// Each tray's efficiency - the entries `average_tray_efficiency` means.
    pub tray_efficiency: Vec<f64>,
    /// **`DistillationColumnMechanicalDesign`'s own report**, the vessel around the internals
    /// tree. `getColumnDiameter()`, which is the *designer's* answer rather than the
    /// Souders-Brown one it was driven at.
    #[serde(serialize_with = "scalar")]
    pub vessel_diameter: Length,
    /// `getColumnHeight()`: `actual_trays * spacing + 1 + 2 + 2 * 0.5`.
    #[serde(serialize_with = "scalar")]
    pub vessel_height: Length,
    /// `getColumnWallThickness()` - recomputed at the final diameter. **The magnitude is SI and
    /// the declared unit is how it is shown**: a dimensioned scalar crosses as a base value the
    /// boundary presents in the spec's unit, which is what `from_si` is the inverse of.
    #[serde(serialize_with = "scalar")]
    pub vessel_wall_thickness: Length,
    /// `getActualTrays()`: `ceil(trays / tray_efficiency)`.
    pub actual_trays: usize,
    /// `getFloodingFactor()`, at the *first pass's* diameter.
    pub flooding_factor: f64,
    /// `getWeirLoading()`, m3/hr per metre of weir, at the same diameter.
    pub weir_loading: f64,
    /// `getTrayPressureDrop()`, mbar/tray - the estimate the class computes and discards.
    pub tray_pressure_drop_mbar: f64,
    /// `getTotalPressureDrop()` - the designer's sum, which overwrote that estimate. SI, shown
    /// in bar.
    #[serde(serialize_with = "scalar")]
    pub total_pressure_drop_bar: Pressure,
    /// `getReboilerDuty()`. SI, shown in kW.
    #[serde(serialize_with = "scalar")]
    pub reboiler_duty_kw: Power,
    /// `getCondenserDuty()`, **absolute**. SI, shown in kW.
    #[serde(serialize_with = "scalar")]
    pub condenser_duty_kw: Power,
    /// `getMaterialGrade()`, carried and reported.
    pub material_grade: String,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl DistillationColumnResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's**: a case's are the spec's `apply_checks` and a flowsheet's
    /// are the checker's, which report through the envelope rather than through a result.
    ///
    /// **The capacity family is the caller's too**, because it reads two inputs the solve does
    /// not: the diameter and the limit. [`fs_limits`] is what computes it, and both callers - this
    /// id's own function and the dispatcher - hand the same value in.
    #[must_use]
    pub fn of(
        outcome: &kernel::ColumnOutcome,
        limits: &FsLimits,
        internals: &DesignerReport,
        mechanical: &MechanicalReport,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            tray_temperature: outcome.trays.iter().map(|t| t.temperature).collect(),
            tray_pressure: outcome.trays.iter().map(|t| t.pressure).collect(),
            tray_gas_n: outcome.trays.iter().map(|t| t.gas_n).collect(),
            tray_liquid_n: outcome.trays.iter().map(|t| t.liquid_n).collect(),
            distillate_n: outcome.distillate.n,
            distillate_z: outcome.distillate.z.clone(),
            distillate_p: outcome.distillate.p,
            distillate_t: outcome.distillate.t,
            distillate_h: joules_per_mole(outcome.distillate.h.value),
            bottoms_n: outcome.bottoms.n,
            bottoms_z: outcome.bottoms.z.clone(),
            bottoms_p: outcome.bottoms.p,
            bottoms_t: outcome.bottoms.t,
            bottoms_h: joules_per_mole(outcome.bottoms.h.value),
            gas_side_draw_n: outcome
                .gas_side_draws
                .iter()
                .map(|draw| draw.as_ref().map_or(0.0, |draw| draw.n))
                .collect(),
            liquid_side_draw_n: outcome
                .liquid_side_draws
                .iter()
                .map(|draw| draw.as_ref().map_or(0.0, |draw| draw.n))
                .collect(),
            pumparound_n: outcome
                .pumparounds
                .iter()
                .map(|draw| draw.as_ref().map_or(0.0, |draw| draw.n))
                .collect(),
            condenser_duty: outcome.condenser_duty,
            reboiler_duty: outcome.reboiler_duty,
            iterations: outcome.iterations,
            temperature_residual: outcome.temperature_residual,
            mass_residual: outcome.mass_residual,
            energy_residual: outcome.energy_residual,
            fs_factor: limits.fs_factor,
            fs_factor_utilization: limits.fs_factor_utilization,
            fs_factor_within_design_limit: limits.fs_factor_within_design_limit,
            minimum_diameter_for_fs_limit: limits.minimum_diameter_for_fs_limit,
            required_diameter: internals.required_diameter,
            controlling_tray_index: internals.controlling_tray_index,
            internals_design_ok: internals.design_ok,
            max_percent_flood: internals.max_percent_flood,
            min_percent_flood: internals.min_percent_flood,
            average_tray_efficiency: internals.average_tray_efficiency,
            total_pressure_drop: internals.total_pressure_drop,
            total_pressure_drop_mbar: internals.total_pressure_drop_mbar,
            tray_percent_flood: internals.trays.iter().map(|t| t.percent_flood).collect(),
            tray_pressure_drop: internals
                .trays
                .iter()
                .map(|t| t.total_pressure_drop.value)
                .collect(),
            tray_efficiency: internals.trays.iter().map(|t| t.tray_efficiency).collect(),
            vessel_diameter: mechanical.vessel_diameter,
            vessel_height: mechanical.vessel_height,
            vessel_wall_thickness: meters(mechanical.vessel_wall_thickness_mm / 1000.0),
            actual_trays: mechanical.actual_trays,
            flooding_factor: mechanical.flooding_factor,
            weir_loading: mechanical.weir_loading,
            tray_pressure_drop_mbar: mechanical.tray_pressure_drop_mbar,
            total_pressure_drop_bar: pascals(mechanical.total_pressure_drop_bar * 1.0e5),
            reboiler_duty_kw: watts(mechanical.reboiler_duty_kw * 1000.0),
            condenser_duty_kw: watts(mechanical.condenser_duty_kw * 1000.0),
            material_grade: mechanical.material_grade.clone(),
            warnings,
        }
    }
}

/// **The two efficiency fields, resolved and clamped, or `None` where neither is stated.**
///
/// `DistillationColumn` clamps rather than refusing - `clampMurphreeEfficiency` is
/// `max(0, min(1, x))` on both setters, so a request outside `[0, 1]` arrives at the solve as its
/// nearest end. The per-stage vector's *length* is the one thing refused rather than clamped,
/// because a length is a statement about the column and not a value.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a per-stage vector whose length is not the **tray** count,
/// which `DistillationColumn.setMurphreeEfficiencies` refuses in its own words. A tray is a
/// stage *or* an end - `numberOfTrays` counts both - so an override on the reboiler or the
/// condenser is a legitimate entry and the vector is one longer per end.
pub fn build_murphree(
    column_wide: Option<f64>,
    per_stage: Option<&[f64]>,
    tray_count: usize,
) -> Result<Option<Murphree>> {
    if column_wide.is_none() && per_stage.is_none() {
        return Ok(None);
    }
    let efficiency = Murphree {
        column_wide: Murphree::clamp(column_wide.unwrap_or(1.0)),
        per_stage: per_stage.map(|values| values.iter().copied().map(Murphree::clamp).collect()),
    };
    Ok(Some(efficiency.checked(tray_count)?.clone()))
}

/// **Solve a distillation column and hand back the kernel's own outcome**, which is what a
/// subclass reads after `super.run(id)`.
///
/// `PackedColumn.run` is `super.run(id)` and then `calcPackingHydraulics()`, and the second step
/// reads the **trays** - the middle one's flows, densities, viscosity and surface tension. A
/// record carries each tray's temperature, pressure and two flows but not its compositions, so
/// the id that inherits this one needs the outcome rather than the record. [`distillation_column`]
/// is the record this returns to every other caller.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a stage count or a feed stage outside the
/// column, for a specification missing the component it constrains, **for every declared
/// parameter whose arithmetic this tranche does not port** - the Murphree efficiency and the
/// nine other solving strategies, each refused by naming the class that would close it - and
/// [`azoth_core::AzothError::SolverNotConverged`] when the solve or a specification misses
/// its gate.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-two
pub(crate) fn distillation_column_outcome(
    coupling: Option<&crate::column::coupling::HydraulicCoupling>,
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    number_of_stages: usize,
    feed_stage: usize,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Pressure,
    bottom_pressure: Pressure,
    reboiler_temperature: Option<ThermodynamicTemperature>,
    condenser_temperature: Option<ThermodynamicTemperature>,
    temperature_tolerance: f64,
    max_iterations: usize,
    murphree_efficiency: Option<f64>,
    tray_murphree_efficiency: Option<&[f64]>,
    solver_type: Option<&str>,
    top_specification_type: Option<&str>,
    top_specification_target: Option<f64>,
    top_specification_component: Option<&str>,
    bottom_specification_type: Option<&str>,
    bottom_specification_target: Option<f64>,
    bottom_specification_component: Option<&str>,
    reactive: Option<bool>,
    reactive_start_tray: Option<usize>,
    reactive_end_tray: Option<usize>,
    gas_side_draw_fractions: Option<&[f64]>,
    liquid_side_draw_fractions: Option<&[f64]>,
    pumparound_fractions: Option<&[f64]>,
    side_draw_flow_tray: Option<usize>,
    side_draw_flow_phase: Option<&str>,
    side_draw_flow_target: Option<f64>,
    side_draw_flow_tolerance: Option<f64>,
    side_draw_flow_max_iterations: Option<usize>,
    pumparound_return_tray: Option<usize>,
    pumparound_draw_tray: Option<usize>,
    pumparound_draw_fraction: Option<f64>,
    pumparound_temperature_drop: Option<f64>,
    pumparound_tolerance: Option<f64>,
    pumparound_max_iterations: Option<usize>,
) -> Result<(kernel::ColumnOutcome, Vec<Warning>)> {
    let pumparounds = build_pumparound_returns(
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        pumparound_temperature_drop,
        pumparound_tolerance,
        pumparound_max_iterations,
    )?;
    let side_draw_flows = build_side_draw_flow(
        side_draw_flow_tray,
        side_draw_flow_phase,
        side_draw_flow_target,
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
    )?;
    // **The length is the tray count, not the stage count.** `setMurphreeEfficiencies` is held
    // to `numberOfTrays`, which counts the ends too - so a six-tray column takes six overrides
    // whatever its middle stage count is.
    let murphree_efficiency = build_murphree(
        murphree_efficiency,
        tray_murphree_efficiency,
        number_of_stages + usize::from(has_reboiler) + usize::from(has_condenser),
    )?;
    // **The end is the parameter's name, not a value.** `ColumnSpecification` carries a
    // location - TOP or BOTTOM - and the class holds exactly two of them, so a declaration
    // that names its parameters `top_*` and `bottom_*` has already said which is which;
    // `setTopSpecification` refusing a BOTTOM location is the same statement from the other
    // side.
    let top_specification = build_specification(
        top_specification_type,
        top_specification_target,
        top_specification_component,
        "top",
    )?;
    let bottom_specification = build_specification(
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        "bottom",
    )?;

    let spec = &crate::model_gen::DISTILLATION_COLUMN_SPEC;
    let mut warnings = Vec::new();
    // **Only the bounds this solve resolves.** It is `process.packed_column`'s solve as well as
    // this id's - `PackedColumn.run` is `super.run(id)` - so a bound the *distillation* spec
    // declares and this solve does not take would otherwise arrive as a `RangeCheckSkipped` on
    // every packed-column case. `tray_efficiency` is exactly that one, and [`distillation_column`]
    // applies it where the value is known.
    const RESOLVED: &[&str] = &[
        "number_of_stages",
        "top_pressure",
        "bottom_pressure",
        "temperature_tolerance",
        "feed_t",
    ];
    apply_checks(
        spec.input_checks()
            .filter(|check| RESOLVED.contains(&check.quantity)),
        |quantity| match quantity {
            "number_of_stages" => Some(number_of_stages as f64),
            "top_pressure" => Some(top_pressure.value),
            "bottom_pressure" => Some(bottom_pressure.value),
            "temperature_tolerance" => Some(temperature_tolerance),
            "feed_t" => Some(feed_t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // **The eight this port does not carry are refused from the declaration, one arm each.**
    // The arm is what makes the pair countable: the key is a literal here and in the Python
    // half, so `tools/check_unported.py` reads two sets of the same shape rather than two
    // sentences. The class that would close each gap comes from the row, so the map that used
    // to be hand-written in both languages is gone with the sentences that needed it.
    let solver = match solver_type.unwrap_or("direct_substitution") {
        "direct_substitution" => kernel::SolverType::DirectSubstitution,
        "naphtali_sandholm" => kernel::SolverType::NaphtaliSandholm,
        "damped_substitution" => return Err(unported::refuse("solver_type=damped_substitution")),
        "inside_out" => return Err(unported::refuse("solver_type=inside_out")),
        "matrix_inside_out" => return Err(unported::refuse("solver_type=matrix_inside_out")),
        "wegstein" => return Err(unported::refuse("solver_type=wegstein")),
        "sum_rates" => return Err(unported::refuse("solver_type=sum_rates")),
        "newton" => return Err(unported::refuse("solver_type=newton")),
        "mesh_residual" => return Err(unported::refuse("solver_type=mesh_residual")),
        "auto" => return Err(unported::refuse("solver_type=auto")),
        other => {
            return Err(AzothError::invalid_input(
                "solver_type",
                format!(
                    "`{other}` is not one of `ColumnSolverFactory`'s ten strategies, which are \
                     `direct_substitution`, `damped_substitution`, `inside_out`, \
                     `matrix_inside_out`, `wegstein`, `sum_rates`, `newton`, \
                     `naphtali_sandholm`, `mesh_residual` and `auto`"
                ),
            ));
        }
    };

    let reactive = kernel::reactive_section(reactive, reactive_start_tray, reactive_end_tray)?;

    // **A reactive section under the simultaneous solve is refused rather than ignored.** This
    // port's mesh takes its fugacities from the mesh's own mixture, where NeqSim's trays take
    // theirs from the tray's own flash - so here the flag could only be silently non-reactive,
    // which is the failure mode this port refuses everywhere else.
    if reactive != kernel::ReactiveSection::None && solver == kernel::SolverType::NaphtaliSandholm {
        return Err(unported::refuse("reactive@solver_type=naphtali_sandholm"));
    }

    // **Side draws under the simultaneous solve are refused for the same reason the reactive
    // section is**: this port's mesh takes its fugacities from its own mixture rather than from a
    // tray's flash, so a draw has no tray outlet there to split.
    let draws_stated = gas_side_draw_fractions.is_some()
        || liquid_side_draw_fractions.is_some()
        || pumparound_fractions.is_some();
    if draws_stated && solver == kernel::SolverType::NaphtaliSandholm {
        return Err(unported::refuse(
            "gas_side_draw_fractions@solver_type=naphtali_sandholm",
        ));
    }

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let setup = kernel::ColumnSetup {
        feed,
        feed_stage,
        number_of_stages,
        has_reboiler,
        has_condenser,
        top_pressure,
        bottom_pressure,
        condenser_temperature,
        reboiler_temperature,
        temperature_tolerance,
        max_iterations,
        murphree_efficiency: murphree_efficiency.clone(),
        // This id is the base column; `AbsorptionColumn`'s override is the absorber's.
        absorber_murphree: None,
        // A case and a capture take the class's cold seed.
        initial_state: None,
        top_specification,
        bottom_specification,
        top_feed: None,
        tray_temperatures: None,
        solver_type: solver,
        gas_side_draw_fractions: gas_side_draw_fractions.map(|v| v.to_vec()),
        side_draw_flows,
        pumparound_returns: pumparounds.0,
        pumparound_inlets: Vec::new(),
        pumparound_tolerance: pumparounds.1,
        pumparound_max_iterations: pumparounds.2,
        liquid_side_draw_fractions: liquid_side_draw_fractions.map(|v| v.to_vec()),
        pumparound_fractions: pumparound_fractions.map(|v| v.to_vec()),
        reactive,
    };
    // **The coupling is the outer loop and the kernel call is its inner solve**, so the flag
    // decides which entry runs rather than adding a step inside either.
    let out = match coupling {
        Some(coupling) => {
            let (mut out, tear) =
                crate::column::coupling::distillation_column_coupled(&setup, coupling)?;
            out.tear = Some(kernel::TearDiagnostics {
                iterations: tear.iterations as usize,
                residual: tear.residual,
                converged: tear.converged,
                rejected_candidates: 0,
                rollbacks: 0,
                inner_iterations: tear.inner_iterations,
                // **The coupling's tear has no candidate search and no specified flow**, which
                // the capture's every coupled row states: zero rejections, zero rollbacks, an
                // empty candidate history and no draw fraction.
                fraction: 0.0,
                actual_flow: 0.0,
                history: String::new(),
            });
            out
        }
        None => kernel(&setup)?,
    };

    // **What the stages had to fall back on joins what the checks found.** A reactive tray's
    // fallback is the kernel's own caveat, and a caller that asked for reactive trays learns
    // from the result whether they ran reactively.
    warnings.extend(out.warnings.iter().cloned());

    Ok((out, warnings))
}

/// Solve a distillation column, as its record.
///
/// # Errors
/// Every error [`distillation_column_outcome`] raises.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-two
pub fn distillation_column(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    number_of_stages: usize,
    feed_stage: usize,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Pressure,
    bottom_pressure: Pressure,
    reboiler_temperature: Option<ThermodynamicTemperature>,
    condenser_temperature: Option<ThermodynamicTemperature>,
    temperature_tolerance: f64,
    max_iterations: usize,
    murphree_efficiency: Option<f64>,
    tray_murphree_efficiency: Option<&[f64]>,
    solver_type: Option<&str>,
    top_specification_type: Option<&str>,
    top_specification_target: Option<f64>,
    top_specification_component: Option<&str>,
    bottom_specification_type: Option<&str>,
    bottom_specification_target: Option<f64>,
    bottom_specification_component: Option<&str>,
    reactive: Option<bool>,
    reactive_start_tray: Option<usize>,
    reactive_end_tray: Option<usize>,
    gas_side_draw_fractions: Option<&[f64]>,
    liquid_side_draw_fractions: Option<&[f64]>,
    pumparound_fractions: Option<&[f64]>,
    side_draw_flow_tray: Option<usize>,
    side_draw_flow_phase: Option<&str>,
    side_draw_flow_target: Option<f64>,
    side_draw_flow_tolerance: Option<f64>,
    side_draw_flow_max_iterations: Option<usize>,
    pumparound_return_tray: Option<usize>,
    pumparound_draw_tray: Option<usize>,
    pumparound_draw_fraction: Option<f64>,
    pumparound_temperature_drop: Option<f64>,
    pumparound_tolerance: Option<f64>,
    pumparound_max_iterations: Option<usize>,
    column_diameter: Option<f64>,
    max_allowable_fs_factor: Option<f64>,
    internals_type: Option<&str>,
    tray_spacing: Option<f64>,
    weir_height: Option<f64>,
    hole_diameter: Option<f64>,
    hole_area_fraction: Option<f64>,
    downcommer_area_fraction: Option<f64>,
    design_flood_fraction: Option<f64>,
    column_diameter_override: Option<f64>,
    hydraulic_pressure_drop_coupling: Option<bool>,
    hydraulic_pressure_drop_internals_type: Option<&str>,
    tray_efficiency: Option<f64>,
    max_flooding_factor: Option<f64>,
    tray_type: Option<&str>,
    contactor_internals_type: Option<&str>,
    material_grade: Option<&str>,
    max_operation_pressure: Option<f64>,
) -> Result<DistillationColumnResult> {
    // **The coupling's own settings, resolved from the two declared inputs.** The flag decides
    // whether the outer loop runs at all; the type is `hydraulicPressureDropInternalsType`, whose
    // class default is `sieve` and which `calcColumnInternals` passes to the designer alone.
    let coupling = hydraulic_pressure_drop_coupling
        .unwrap_or(false)
        .then(|| HydraulicCoupling {
            internals_type: hydraulic_pressure_drop_internals_type
                .unwrap_or(DEFAULT_INTERNALS_TYPE)
                .to_string(),
            ..HydraulicCoupling::default()
        });
    let (out, mut warnings) = distillation_column_outcome(
        coupling.as_ref(),
        components,
        feed_n,
        feed_z,
        feed_p,
        feed_t,
        number_of_stages,
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
        tray_murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        reactive,
        reactive_start_tray,
        reactive_end_tray,
        gas_side_draw_fractions,
        liquid_side_draw_fractions,
        pumparound_fractions,
        side_draw_flow_tray,
        side_draw_flow_phase,
        side_draw_flow_target,
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        pumparound_temperature_drop,
        pumparound_tolerance,
        pumparound_max_iterations,
    )?;
    // **The declared value this id refuses**, from its own `[[unported]]` row: the class's packed
    // arm is a second machine - a preset, a bed height and a hydraulic capacity factor this entry
    // carries none of - so it is refused by name rather than half-carried.
    if contactor_internals_type.is_some_and(|kind| kind.eq_ignore_ascii_case("packed")) {
        return Err(unported::refuse("contactor_internals_type=packed"));
    }
    // **The vessel sizing's own bound, applied where the value is known.** The shared solve
    // resolves the five *it* takes; this one is the record's, and an omitted efficiency is the
    // default the sizing runs at rather than a bound nothing could check.
    apply_checks(
        crate::model_gen::DISTILLATION_COLUMN_SPEC
            .input_checks()
            .filter(|check| check.quantity == "tray_efficiency"),
        |_| Some(tray_efficiency.unwrap_or(DEFAULT_TRAY_EFFICIENCY)),
        &mut warnings,
    )?;
    // **The internal diameter the class carries, which its own constructor defaults to `1.0`.**
    // A stated one that is not positive is the class's own "not stated" value, so it is the
    // default that answers rather than a refusal.
    let diameter = meters(column_diameter.unwrap_or(DEFAULT_INTERNAL_DIAMETER_M));
    let limits = fs_limits(
        &out.distillate,
        diameter,
        max_allowable_fs_factor.unwrap_or(DEFAULT_MAX_ALLOWABLE_FS_FACTOR),
    )?;
    // **The internals tree, read after the solve and from the trays it left.** Its own geometry
    // inputs reach nothing on the run path, which is the class's own split: `calculateTrayed`
    // walks `getTrays()` and reads the profile the solve published.
    let internals = designer_report(
        &out.trays,
        components,
        &DesignerGeometry {
            internals_type: internals_type.unwrap_or(DEFAULT_INTERNALS_TYPE).to_string(),
            tray_spacing: meters(tray_spacing.unwrap_or(DEFAULT_TRAY_SPACING_M)),
            weir_height: meters(weir_height.unwrap_or(DEFAULT_WEIR_HEIGHT_M)),
            hole_diameter: millimeters(hole_diameter.unwrap_or(DEFAULT_HOLE_DIAMETER_MM)),
            hole_area_fraction: hole_area_fraction.unwrap_or(DEFAULT_HOLE_AREA_FRACTION),
            downcommer_area_fraction: downcommer_area_fraction
                .unwrap_or(DEFAULT_DOWNCOMMER_AREA_FRACTION),
            design_flood_fraction: design_flood_fraction.unwrap_or(DEFAULT_DESIGNER_FLOOD_FRACTION),
            column_diameter_override: meters(
                column_diameter_override.unwrap_or(UNSIZED_COLUMN_DIAMETER_M),
            ),
        },
    )?;
    // **The vessel around the internals, on the same trays and the same geometry.** The two
    // duties and the tray type are the only things here the designer above does not also take;
    // everything else is stated once and read by both.
    let mechanical = mechanical_design(
        &out.trays,
        components,
        &out.distillate,
        &out.bottoms,
        has_reboiler.then_some(out.reboiler_duty),
        has_condenser.then_some(out.condenser_duty),
        &MechanicalGeometry {
            tray_type: tray_type.unwrap_or(DEFAULT_TRAY_TYPE).to_string(),
            contactor_internals_type: contactor_internals_type
                .unwrap_or(DEFAULT_CONTACTOR_INTERNALS_TYPE)
                .to_string(),
            tray_efficiency: tray_efficiency.unwrap_or(DEFAULT_TRAY_EFFICIENCY),
            max_flooding_factor: max_flooding_factor.unwrap_or(DEFAULT_MAX_FLOODING_FACTOR),
            material_grade: material_grade.unwrap_or(DEFAULT_MATERIAL_GRADE).to_string(),
            max_operation_pressure_bara: max_operation_pressure
                .unwrap_or(DEFAULT_MAX_OPERATION_PRESSURE_BARA),
            tray_spacing: meters(tray_spacing.unwrap_or(DEFAULT_TRAY_SPACING_M)),
            weir_height: meters(weir_height.unwrap_or(DEFAULT_WEIR_HEIGHT_M)),
            hole_diameter: millimeters(hole_diameter.unwrap_or(DEFAULT_HOLE_DIAMETER_MM)),
            hole_area_fraction: hole_area_fraction.unwrap_or(DEFAULT_HOLE_AREA_FRACTION),
            downcommer_area_fraction: downcommer_area_fraction
                .unwrap_or(DEFAULT_DOWNCOMMER_AREA_FRACTION),
            column_diameter_override: meters(
                column_diameter_override.unwrap_or(UNSIZED_COLUMN_DIAMETER_M),
            ),
        },
    )?;
    Ok(DistillationColumnResult::of(
        &out,
        &limits,
        &internals,
        &mechanical,
        warnings,
    ))
}

/// **The one side-draw flow specification this model declares**, or an empty vector.
///
/// **One and not a list, because one is what the port implements.** `addSideDrawFlowSpecification`
/// appends to a list and the class solves a list of them as *coordinated* tear variables - a
/// different search, refused by name. So the declaration here is the single-specification case
/// the kernel carries, expressed as five scalars rather than as a record the spec schema has no
/// input type for.
///
/// **A target without a tray or a phase is refused**, the same judgement the two end
/// specifications make: a declaration that cannot be read is not a default worth guessing.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a partial specification, for a phase the enum does not
/// carry, and for a target that is not finite and non-negative - the class's own
/// `ColumnSideDrawSpecification` constructor checks.
pub fn build_side_draw_flow(
    tray: Option<usize>,
    phase: Option<&str>,
    target: Option<f64>,
    tolerance: Option<f64>,
    max_iterations: Option<usize>,
) -> Result<Vec<kernel::SideDrawFlow>> {
    let (Some(tray), Some(phase), Some(target)) = (tray, phase, target) else {
        if tray.is_some() || phase.is_some() || target.is_some() {
            return Err(AzothError::invalid_input(
                "side_draw_flow_tray",
                "a side-draw flow specification is stated by its tray, its phase and its target \
                 together: `addSideDrawFlowSpecification(tray, phase, flow, unit)` takes all \
                 three, and a declaration that states one or two of them says nothing",
            ));
        }
        return Ok(Vec::new());
    };
    if !target.is_finite() || target < 0.0 {
        return Err(AzothError::invalid_input(
            "side_draw_flow_target",
            format!(
                "a side-draw target flow of {target} is not finite and non-negative, which \
                 `ColumnSideDrawSpecification`'s own constructor refuses"
            ),
        ));
    }
    let phase = match phase {
        "gas" => kernel::SideDrawPhase::Gas,
        "liquid" => kernel::SideDrawPhase::Liquid,
        other => {
            return Err(AzothError::invalid_input(
                "side_draw_flow_phase",
                format!(
                    "`{other}` is not a side-draw phase: `SideDrawPhase` carries `gas` and \
                     `liquid`"
                ),
            ));
        }
    };
    Ok(vec![kernel::SideDrawFlow {
        tray,
        phase,
        target,
        tolerance: tolerance.unwrap_or(1.0e-4),
        max_iterations: max_iterations.unwrap_or(12),
    }])
}

/// **The one pumparound with a return this model declares**, or nothing.
///
/// **One and not a list**, for the reason the side-draw flow specification gives:
/// `addLiquidPumparound` owns one draw tray each and refuses a second pumparound on the same
/// tray, so the case this port carries is one pumparound stated as scalars. A partial statement
/// is refused, which is the judgement the end specifications make.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a partial statement; the kernel checks the trays, the
/// fraction and the drop again, by tray, because that is where `addLiquidPumparound` does.
pub fn build_pumparound_returns(
    return_tray: Option<usize>,
    draw_tray: Option<usize>,
    draw_fraction: Option<f64>,
    temperature_drop: Option<f64>,
    tolerance: Option<f64>,
    max_iterations: Option<usize>,
) -> Result<(
    Vec<crate::column::pumparound::PumparoundReturn>,
    Option<f64>,
    Option<usize>,
)> {
    let (Some(return_tray), Some(draw_tray), Some(draw_fraction), Some(temperature_drop)) =
        (return_tray, draw_tray, draw_fraction, temperature_drop)
    else {
        let none_stated = return_tray.is_none()
            && draw_tray.is_none()
            && draw_fraction.is_none()
            && temperature_drop.is_none();
        if none_stated {
            return Ok((Vec::new(), tolerance, max_iterations));
        }
        return Err(AzothError::invalid_input(
            "pumparound_return_tray",
            "a pumparound with a return is stated by its draw tray, its return tray, its \
             fraction and its temperature drop together: `addLiquidPumparound(name, drawTray, \
             returnTray, fraction, drop)` takes all four, and a declaration that states some of \
             them says nothing",
        ));
    };
    Ok((
        vec![crate::column::pumparound::PumparoundReturn {
            draw_tray,
            return_tray,
            fraction: draw_fraction,
            temperature_drop,
        }],
        tolerance,
        max_iterations,
    ))
}

/// One end's specification from its three declared parameters.
///
/// `None` where no type was stated, which is the temperature-pinned route the model takes by
/// default. A target without a type, or a type without a target, is a declaration that cannot
/// be read rather than a default worth guessing.
fn build_specification(
    kind: Option<&str>,
    target: Option<f64>,
    component: Option<&str>,
    which: &str,
) -> Result<Option<Specification>> {
    let Some(kind) = kind else {
        if target.is_some() || component.is_some() {
            return Err(AzothError::invalid_input(
                "specification",
                format!(
                    "the {which} has a target or a component and no type, so nothing says what \
                     it constrains"
                ),
            ));
        }
        return Ok(None);
    };
    let kind = match kind {
        "product_purity" => SpecificationKind::ProductPurity,
        "component_recovery" => SpecificationKind::ComponentRecovery,
        "product_flow_rate" => SpecificationKind::ProductFlowRate,
        "reflux_ratio" => SpecificationKind::RefluxRatio,
        "duty" => SpecificationKind::Duty,
        other => {
            return Err(AzothError::invalid_input(
                "specification",
                format!(
                    "{other} is not one of `ColumnSpecification`'s five types: product_purity, \
                     reflux_ratio, component_recovery, product_flow_rate or duty"
                ),
            ));
        }
    };
    let target = target.ok_or_else(|| {
        AzothError::invalid_input(
            "specification",
            format!("the {which} specification states a type and no target"),
        )
    })?;
    // A purity or a recovery constrains a named component, and the class's own `validate`
    // refuses one that does not name it; a flow rate, a ratio and a duty do not read the name.
    if matches!(
        kind,
        SpecificationKind::ProductPurity | SpecificationKind::ComponentRecovery
    ) && component.is_none()
    {
        return Err(AzothError::invalid_input(
            "specification",
            format!("a {which} purity or recovery constrains a component, and none was stated"),
        ));
    }
    Ok(Some(Specification {
        kind,
        target,
        component: component.map(String::from),
    }))
}
