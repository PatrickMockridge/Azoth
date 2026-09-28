//! `process.distillation_column` - the column solve as a registered id.
//!
//! Spec: `specs/models/process/distillation_column.toml`. The arithmetic is
//! [`crate::kernels::distillation_column`]; what is here is the boundary a case and a
//! cross-impl test address, and the refusal of every parameter the palette declares and this
//! tranche does not implement.

use azoth_core::units::{MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};
use serde::Serialize;

use crate::column::murphree::Murphree;
use crate::executor::json::{scalar, scalars, warnings as wire_warnings};
use crate::kernels::distillation_column as kernel;
use crate::kernels::distillation_column::{Specification, SpecificationKind};
use crate::stream::Stream;

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
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl DistillationColumnResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's**: a case's are the spec's `apply_checks` and a flowsheet's
    /// are the checker's, which report through the envelope rather than through a result.
    #[must_use]
    pub fn of(outcome: &kernel::ColumnOutcome, warnings: Vec<Warning>) -> Self {
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
            warnings,
        }
    }
}

impl CalcResult for DistillationColumnResult {
    const CALC_ID: &'static str = "process.distillation_column";
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

/// Solve a distillation column.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a stage count or a feed stage outside the
/// column, for a specification missing the component it constrains, **for every declared
/// parameter whose arithmetic this tranche does not port** - the Murphree efficiency and the
/// nine other solving strategies, each refused by naming the class that would close it - and
/// [`azoth_core::AzothError::SolverNotConverged`] when the solve or a specification misses
/// its gate.
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
) -> Result<DistillationColumnResult> {
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
    apply_checks(
        spec.input_checks(),
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

    let solver = match solver_type.unwrap_or("direct_substitution") {
        "direct_substitution" => kernel::SolverType::DirectSubstitution,
        "naphtali_sandholm" => kernel::SolverType::NaphtaliSandholm,
        other => {
            return Err(AzothError::invalid_input(
                "solver_type",
                format!(
                    "solver_type = {other} is not ported: `ColumnSolverFactory.{}` is the class \
                     that would close it. **The capture measures why it is refused rather than \
                     ported**: `validation/neqsim/captures/process_column_solvers.tsv` runs the \
                     binary column under all ten strategies and puts every one of them within \
                     `2.5e-6` K of every other on tray 1 and within `1.1e-7` relative on the \
                     distillate, so they are path variants rather than different physics - the \
                     ten land on three bit-identical states, and each state is where a solve \
                     stopped.",
                    unported_class(other)
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
        return Err(AzothError::invalid_input(
            "reactive",
            "a reactive section under `naphtali_sandholm` is not ported: `NaphtaliSandholmSolver` \
             reads its fugacities from the MESH equations, and this port's mesh does not route a \
             tray's flash at all, so the flag would be ignored rather than honoured",
        ));
    }

    // **Side draws under the simultaneous solve are refused for the same reason the reactive
    // section is**: this port's mesh takes its fugacities from its own mixture rather than from a
    // tray's flash, so a draw has no tray outlet there to split.
    let draws_stated = gas_side_draw_fractions.is_some()
        || liquid_side_draw_fractions.is_some()
        || pumparound_fractions.is_some();
    if draws_stated && solver == kernel::SolverType::NaphtaliSandholm {
        return Err(AzothError::invalid_input(
            "gas_side_draw_fractions",
            "side draws under `naphtali_sandholm` are not ported: this port's mesh solves the \
             MESH equations together and never forms a tray's own outlet streams, so there is \
             nothing for a fraction to split",
        ));
    }

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let out = kernel(&kernel::ColumnSetup {
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
    })?;

    // **What the stages had to fall back on joins what the checks found.** A reactive tray's
    // fallback is the kernel's own caveat, and a caller that asked for reactive trays learns
    // from the result whether they ran reactively.
    warnings.extend(out.warnings.iter().cloned());

    Ok(DistillationColumnResult::of(&out, warnings))
}

/// The `ColumnSolverFactory` class behind each strategy this port does not carry.
///
/// **Named per strategy, because that is what a refusal owes a caller.** `columnSolver` hands
/// back one of these for each `SolverType`, and `AutoSolver` is the ladder rather than a
/// method: `candidateSolvers` returns `NAPHTALI_SANDHOLM` first, then `MATRIX_INSIDE_OUT`,
/// `INSIDE_OUT` and `DAMPED_SUBSTITUTION`, and this port has the first and the last of those.
///
/// **The eight are refused as path variants rather than as physics**, which the capture
/// measures: on the binary column every one of the ten lands within `2.5e-6` K of every other
/// on tray 1 and within `1.1e-7` relative on the distillate, on three bit-identical states -
/// and `inside_out`, `matrix_inside_out`, `mesh_residual` and a fallen-back `wegstein` land on
/// the substitution core's own state exactly. That is a measured non-port, the shape `P10` used
/// for the adaptive-derivative refinement that converges on nothing.
pub(crate) fn unported_class(strategy: &str) -> &'static str {
    match strategy {
        "damped_substitution" => "DampedSubstitutionSolver",
        "inside_out" => "InsideOutSolver",
        "matrix_inside_out" => "MatrixInsideOutSolver",
        "wegstein" => "WegsteinSolver",
        "sum_rates" => "SumRatesSolver",
        "newton" => "TemperatureNewtonSolver",
        "mesh_residual" => "MeshResidualSolver",
        _ => "AutoSolver",
    }
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
