//! `process.absorption_column` - the tray absorber as a registered id.
//!
//! Spec: `specs/models/process/absorption_column.toml`. The arithmetic is
//! [`crate::kernels::absorption_column`], which is the column's with two feeds and no ends;
//! what is here is the boundary a case and a cross-impl test address, and the refusal of every
//! parameter the palette declares and this stage does not implement.

use azoth_core::units::{
    Length, MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole, meters,
};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};
use serde::Serialize;

use crate::column::absorber_murphree::AbsorberMurphree;
use crate::column::capacity::{
    DEFAULT_INTERNAL_DIAMETER_M, DEFAULT_MAX_ALLOWABLE_FS_FACTOR_ABSORBER,
    DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR, FsLimits, GasLoadLimits, fs_limits, gas_load_limits,
};
use crate::column::murphree::Murphree;
use crate::executor::json::{scalar, scalars, warnings as wire_warnings};
use crate::kernels::absorption_column::AbsorberOutcome;
use crate::kernels::absorption_column::AbsorberSetup;
use crate::kernels::absorption_column::absorption_column as kernel;
use crate::kernels::distillation_column::{SolverType, tray_streams};
use crate::models::distillation_column::{build_pumparound_returns, build_side_draw_flow};
use crate::stream::Stream;
use crate::unported;

/// Result of `process.absorption_column`.
///
/// **The profile is the answer, as it is for the column, and the products are the class's own
/// getters**: `getGasOutStream` is the treated gas overhead and `getLiquidOutStream` the
/// loaded solvent. There is no condenser and no reboiler, so there are no duties.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct AbsorptionColumnResult {
    /// Each tray's temperature, K, from the gas end at stage 0 up.
    #[serde(serialize_with = "scalars")]
    pub tray_temperature: Vec<ThermodynamicTemperature>,
    /// Each tray's pressure, Pa.
    #[serde(serialize_with = "scalars")]
    pub tray_pressure: Vec<Pressure>,
    /// Each tray's vapour traffic, mol/s.
    pub tray_gas_n: Vec<f64>,
    /// Each tray's liquid traffic, mol/s.
    pub tray_liquid_n: Vec<f64>,
    /// The treated gas's molar flow, mol/s.
    pub gas_out_n: f64,
    /// The treated gas's composition.
    pub gas_out_z: Vec<f64>,
    /// The treated gas's pressure.
    #[serde(serialize_with = "scalar")]
    pub gas_out_p: Pressure,
    /// The treated gas's temperature.
    #[serde(serialize_with = "scalar")]
    pub gas_out_t: ThermodynamicTemperature,
    /// The treated gas's molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub gas_out_h: MolarEnergy,
    /// The loaded solvent's molar flow, mol/s.
    pub liquid_out_n: f64,
    /// The loaded solvent's composition.
    pub liquid_out_z: Vec<f64>,
    /// The loaded solvent's pressure.
    #[serde(serialize_with = "scalar")]
    pub liquid_out_p: Pressure,
    /// The loaded solvent's temperature.
    #[serde(serialize_with = "scalar")]
    pub liquid_out_t: ThermodynamicTemperature,
    /// The loaded solvent's molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub liquid_out_h: MolarEnergy,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against both feeds, relative.
    pub mass_residual: f64,
    /// The enthalpy closure.
    pub energy_residual: f64,
    /// `getFsFactor`: the Fs family is the base class's, so an absorber carries it too.
    pub fs_factor: f64,
    /// `getFsFactorUtilization`, against this class's own `3.0` limit.
    pub fs_factor_utilization: f64,
    /// `isFsFactorWithinDesignLimit`.
    pub fs_factor_within_design_limit: bool,
    /// `getMinimumDiameterForFsLimit`.
    #[serde(serialize_with = "scalar")]
    pub minimum_diameter_for_fs_limit: Length,
    /// `getGasLoadFactor`: the Souders-Brown `Ks`, which is `AbsorptionColumn`'s own method.
    pub gas_load_factor: f64,
    /// `getGasLoadFactorUtilization`.
    pub gas_load_factor_utilization: f64,
    /// `isGasLoadFactorWithinDesignLimit`.
    pub gas_load_factor_within_design_limit: bool,
    /// `getMinimumDiameterForGasLoadLimit`.
    #[serde(serialize_with = "scalar")]
    pub minimum_diameter_for_gas_load_limit: Length,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl AbsorptionColumnResult {
    /// The result of one kernel call.
    ///
    /// **The warnings are the caller's**: a case's are the spec's `apply_checks` and a flowsheet's
    /// are the checker's, which report through the envelope rather than through a result.
    ///
    /// **Both capacity families are the caller's too**, because they read two inputs the solve
    /// does not: the diameter and the two limits.
    #[must_use]
    pub fn of(
        outcome: &AbsorberOutcome,
        fs: &FsLimits,
        gas_load: &GasLoadLimits,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            tray_temperature: outcome.trays.iter().map(|tray| tray.temperature).collect(),
            tray_pressure: outcome.trays.iter().map(|tray| tray.pressure).collect(),
            tray_gas_n: outcome.trays.iter().map(|tray| tray.gas_n).collect(),
            tray_liquid_n: outcome.trays.iter().map(|tray| tray.liquid_n).collect(),
            gas_out_n: outcome.gas_out.n,
            gas_out_z: outcome.gas_out.z.clone(),
            gas_out_p: outcome.gas_out.p,
            gas_out_t: outcome.gas_out.t,
            gas_out_h: joules_per_mole(outcome.gas_out.h.value),
            liquid_out_n: outcome.liquid_out.n,
            liquid_out_z: outcome.liquid_out.z.clone(),
            liquid_out_p: outcome.liquid_out.p,
            liquid_out_t: outcome.liquid_out.t,
            liquid_out_h: joules_per_mole(outcome.liquid_out.h.value),
            iterations: outcome.iterations,
            temperature_residual: outcome.temperature_residual,
            mass_residual: outcome.mass_residual,
            energy_residual: outcome.energy_residual,
            fs_factor: fs.fs_factor,
            fs_factor_utilization: fs.fs_factor_utilization,
            fs_factor_within_design_limit: fs.fs_factor_within_design_limit,
            minimum_diameter_for_fs_limit: fs.minimum_diameter_for_fs_limit,
            gas_load_factor: gas_load.gas_load_factor,
            gas_load_factor_utilization: gas_load.gas_load_factor_utilization,
            gas_load_factor_within_design_limit: gas_load.gas_load_factor_within_design_limit,
            minimum_diameter_for_gas_load_limit: gas_load.minimum_diameter_for_gas_load_limit,
            warnings,
        }
    }
}

impl CalcResult for AbsorptionColumnResult {
    const CALC_ID: &'static str = "process.absorption_column";
    const FIELDS: &'static [&'static str] = &[
        "tray_temperature",
        "tray_pressure",
        "tray_gas_n",
        "tray_liquid_n",
        "gas_out_n",
        "gas_out_z",
        "gas_out_p",
        "gas_out_t",
        "gas_out_h",
        "liquid_out_n",
        "liquid_out_z",
        "liquid_out_p",
        "liquid_out_t",
        "liquid_out_h",
        "iterations",
        "temperature_residual",
        "mass_residual",
        "energy_residual",
        "fs_factor",
        "fs_factor_utilization",
        "fs_factor_within_design_limit",
        "minimum_diameter_for_fs_limit",
        "gas_load_factor",
        "gas_load_factor_utilization",
        "gas_load_factor_within_design_limit",
        "minimum_diameter_for_gas_load_limit",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Solve a tray absorber.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a declared parameter whose arithmetic this
/// stage does not port - each refusal names the class that would close it - and
/// [`azoth_core::AzothError::SolverNotConverged`] when the solve misses its gate.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn absorption_column(
    gas_components: &[String],
    gas_n: f64,
    gas_z: &[f64],
    gas_p: Pressure,
    gas_t: ThermodynamicTemperature,
    solvent_components: &[String],
    solvent_n: f64,
    solvent_z: &[f64],
    solvent_p: Pressure,
    solvent_t: ThermodynamicTemperature,
    number_of_stages: usize,
    top_pressure: Pressure,
    bottom_pressure: Pressure,
    tray_temperatures: Option<&[f64]>,
    temperature_tolerance: f64,
    max_iterations: usize,
    murphree_efficiency: Option<f64>,
    component_murphree_efficiency: Option<&[f64]>,
    max_allowable_gas_load_factor: Option<f64>,
    reactive: Option<bool>,
    reactive_start_tray: Option<usize>,
    reactive_end_tray: Option<usize>,
    solver_type: Option<&str>,
    tray_murphree_efficiency: Option<&[f64]>,
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
) -> Result<AbsorptionColumnResult> {
    // **The two efficiency fields, resolved and clamped.** `AbsorptionColumn` inherits the
    // base's two setters and adds `setComponentMurphreeEfficiency`'s map, so the correction
    // reads component → per-tray → column-wide in that order; the per-tray-per-component map
    // `setComponentMurphreeEfficiency(int, String, double)` writes has no palette spelling and
    // is not carried.
    let murphree = if murphree_efficiency.is_none() && component_murphree_efficiency.is_none() {
        None
    } else {
        Some(
            AbsorberMurphree {
                base: Murphree {
                    column_wide: Murphree::clamp(murphree_efficiency.unwrap_or(1.0)),
                    per_stage: tray_murphree_efficiency
                        .map(|values| values.iter().copied().map(Murphree::clamp).collect()),
                },
                per_component: component_murphree_efficiency
                    .map(|values| values.iter().copied().map(Murphree::clamp).collect()),
            }
            .checked(gas_components.len())?
            .clone(),
        )
    };
    // The same eight `process.distillation_column` refuses, from this id's own declaration.
    let solver = match solver_type.unwrap_or("direct_substitution") {
        "direct_substitution" => SolverType::DirectSubstitution,
        "naphtali_sandholm" => SolverType::NaphtaliSandholm,
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

    let spec = &crate::model_gen::ABSORPTION_COLUMN_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "number_of_stages" => Some(number_of_stages as f64),
            "top_pressure" => Some(top_pressure.value),
            "bottom_pressure" => Some(bottom_pressure.value),
            "temperature_tolerance" => Some(temperature_tolerance),
            "gas_t" => Some(gas_t.value),
            "solvent_t" => Some(solvent_t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let gas = Stream::from_pt(gas_components.to_vec(), gas_z.to_vec(), gas_n, gas_p, gas_t)?;
    let solvent = Stream::from_pt(
        solvent_components.to_vec(),
        solvent_z.to_vec(),
        solvent_n,
        solvent_p,
        solvent_t,
    )?;
    if gas.components != solvent.components {
        return Err(AzothError::invalid_input(
            "solvent_components",
            "the two inlets of one column carry the same substances in the same order: the \
             class's own feed handling indexes both by the same component list, and a solvent \
             that names different substances in a different order is a different fluid, not an \
             inlet",
        ));
    }

    let out = kernel(&AbsorberSetup {
        gas,
        solvent,
        number_of_stages,
        top_pressure,
        bottom_pressure,
        tray_temperatures: tray_temperatures.map(<[f64]>::to_vec),
        murphree,
        // **An absorber has no ends, so a draw on tray 0 or the top tray is legitimate** - the
        // base kernel's `refuse_end_draws` keys its refusal on the two end flags, both false here.
        gas_side_draw_fractions: gas_side_draw_fractions.map(<[f64]>::to_vec),
        liquid_side_draw_fractions: liquid_side_draw_fractions.map(<[f64]>::to_vec),
        pumparound_fractions: pumparound_fractions.map(<[f64]>::to_vec),
        side_draw_flows: build_side_draw_flow(
            side_draw_flow_tray,
            side_draw_flow_phase,
            side_draw_flow_target,
            side_draw_flow_tolerance,
            side_draw_flow_max_iterations,
        )?,
        pumparound_returns: build_pumparound_returns(
            pumparound_return_tray,
            pumparound_draw_tray,
            pumparound_draw_fraction,
            pumparound_temperature_drop,
            pumparound_tolerance,
            pumparound_max_iterations,
        )?
        .0,
        pumparound_tolerance,
        pumparound_max_iterations,
        temperature_tolerance,
        max_iterations,
        solver_type: solver,
        reactive: crate::kernels::distillation_column::reactive_section(
            reactive,
            reactive_start_tray,
            reactive_end_tray,
        )?,
    })?;

    warnings.extend(out.warnings.iter().cloned());

    // **The two limits are read after the solve, and one of them used to be discarded here.**
    // `max_allowable_gas_load_factor` was declared by the palette and accepted inertly, because
    // nothing on the run path reads it - which is true and is not a reason to answer nothing
    // with it: `getGasLoadFactor` is the quantity it is the limit *of*.
    let diameter = meters(column_diameter.unwrap_or(DEFAULT_INTERNAL_DIAMETER_M));
    let fs = fs_limits(
        &out.gas_out,
        diameter,
        max_allowable_fs_factor.unwrap_or(DEFAULT_MAX_ALLOWABLE_FS_FACTOR_ABSORBER),
    )?;
    // **The liquid outlet's own system, which the class reads phase 0 of.** `getLiquidOutStream`
    // answers a stream carrying the *tray's* two-phase system, so its phase 0 is the vapour the
    // bottom tray leaves with - not this port's liquid product, which is that tray's liquid phase
    // alone. `column/tray_streams` rebuilds the vapour from the profile, which is what makes the
    // reading available at all.
    let (bottom_vapour, _) = tray_streams(&out.trays[0], gas_components)?;
    let gas_load = gas_load_limits(
        &out.gas_out,
        &bottom_vapour,
        diameter,
        max_allowable_gas_load_factor.unwrap_or(DEFAULT_MAX_ALLOWABLE_GAS_LOAD_FACTOR),
    )?;

    Ok(AbsorptionColumnResult::of(&out, &fs, &gas_load, warnings))
}
