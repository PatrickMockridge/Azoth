//! `process.distillation_column` - the column solve as a registered id.
//!
//! Spec: `specs/models/process/distillation_column.toml`. The arithmetic is
//! [`crate::kernels::distillation_column`]; what is here is the boundary a case and a
//! cross-impl test address, and the refusal of every parameter the palette declares and this
//! tranche does not implement.

use azoth_core::units::{MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::kernels::distillation_column as kernel;
use crate::kernels::distillation_column::{Specification, SpecificationKind};
use crate::stream::Stream;

/// Result of `process.distillation_column`.
///
/// **The profile is a set of vectors, one entry per tray**, which is what a column's answer is
/// rather than a scalar: four quantities over the trays, then the two products' records, the
/// two duties and the three residuals. The vector lengths are the tray count, which the inputs
/// decide - a shape `eos.pt_phase_envelope` already exercises with its trace points.
#[derive(Debug, Clone, PartialEq)]
pub struct DistillationColumnResult {
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
    solver_type: Option<&str>,
    top_specification_type: Option<&str>,
    top_specification_target: Option<f64>,
    top_specification_component: Option<&str>,
    bottom_specification_type: Option<&str>,
    bottom_specification_target: Option<f64>,
    bottom_specification_component: Option<&str>,
) -> Result<DistillationColumnResult> {
    refuse_unported(murphree_efficiency)?;

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
                    "solver_type = {other} is not ported. `direct_substitution` is the \
                     class's own default and `naphtali_sandholm` is `NaphtaliSandholmSolver`; \
                     the rest are `ColumnSolverFactory`'s inside-out family and damping \
                     strategies, where `auto` is a ladder over them rather than one method"
                ),
            ));
        }
    };

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
        top_specification,
        bottom_specification,
        solver_type: solver,
    })?;

    Ok(DistillationColumnResult {
        tray_temperature: out.trays.iter().map(|t| t.temperature).collect(),
        tray_pressure: out.trays.iter().map(|t| t.pressure).collect(),
        tray_gas_n: out.trays.iter().map(|t| t.gas_n).collect(),
        tray_liquid_n: out.trays.iter().map(|t| t.liquid_n).collect(),
        distillate_n: out.distillate.n,
        distillate_z: out.distillate.z,
        distillate_p: out.distillate.p,
        distillate_t: out.distillate.t,
        distillate_h: joules_per_mole(out.distillate.h.value),
        bottoms_n: out.bottoms.n,
        bottoms_z: out.bottoms.z,
        bottoms_p: out.bottoms.p,
        bottoms_t: out.bottoms.t,
        bottoms_h: joules_per_mole(out.bottoms.h.value),
        condenser_duty: out.condenser_duty,
        reboiler_duty: out.reboiler_duty,
        iterations: out.iterations,
        temperature_residual: out.temperature_residual,
        mass_residual: out.mass_residual,
        energy_residual: out.energy_residual,
        warnings,
    })
}

/// Refuse every parameter the palette declares and this tranche does not implement.
///
/// **Declared and refused, rather than withdrawn.** `unit_ops.distillation_column`'s palette
/// entry declares these because they are the machine's own form fields - the same reason
/// `unit_ops.throttling_valve`'s `valve_opening` was *withdrawn*, and the opposite conclusion:
/// there nothing read it and nothing was coming, and here the class that reads each one is
/// named and its stage is known. A form field that errors with the reason beats a form field
/// that is silently absent, and beats one that answers with the ideal stage.
fn refuse_unported(murphree_efficiency: Option<f64>) -> Result<()> {
    if let Some(efficiency) = murphree_efficiency {
        return Err(AzothError::invalid_input(
            "murphree_efficiency",
            format!(
                "a Murphree efficiency of {efficiency} is not ported: \
                 `SimpleTray.setMurphreeEfficiency` and the per-tray correction the column \
                 solver applies after each run are the classes that would close it. Omitted \
                 means the ideal stage, which is the class's own default of one"
            ),
        ));
    }
    Ok(())
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
