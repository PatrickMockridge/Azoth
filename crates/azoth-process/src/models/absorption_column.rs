//! `process.absorption_column` - the tray absorber as a registered id.
//!
//! Spec: `specs/models/process/absorption_column.toml`. The arithmetic is
//! [`crate::kernels::absorption_column`], which is the column's with two feeds and no ends;
//! what is here is the boundary a case and a cross-impl test address, and the refusal of every
//! parameter the palette declares and this stage does not implement.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::kernels::absorption_column::AbsorberSetup;
use crate::kernels::absorption_column::absorption_column as kernel;
use crate::kernels::distillation_column::SolverType;
use crate::models::distillation_column::unported_class;
use crate::stream::Stream;

/// Result of `process.absorption_column`.
///
/// **The profile is the answer, as it is for the column, and the products are the class's own
/// getters**: `getGasOutStream` is the treated gas overhead and `getLiquidOutStream` the
/// loaded solvent. There is no condenser and no reboiler, so there are no duties.
#[derive(Debug, Clone, PartialEq)]
pub struct AbsorptionColumnResult {
    /// Each tray's temperature, K, from the gas end at stage 0 up.
    pub tray_temperature: Vec<ThermodynamicTemperature>,
    /// Each tray's pressure, Pa.
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
    pub gas_out_p: Pressure,
    /// The treated gas's temperature.
    pub gas_out_t: ThermodynamicTemperature,
    /// The treated gas's molar enthalpy.
    pub gas_out_h: MolarEnergy,
    /// The loaded solvent's molar flow, mol/s.
    pub liquid_out_n: f64,
    /// The loaded solvent's composition.
    pub liquid_out_z: Vec<f64>,
    /// The loaded solvent's pressure.
    pub liquid_out_p: Pressure,
    /// The loaded solvent's temperature.
    pub liquid_out_t: ThermodynamicTemperature,
    /// The loaded solvent's molar enthalpy.
    pub liquid_out_h: MolarEnergy,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against both feeds, relative.
    pub mass_residual: f64,
    /// The enthalpy closure.
    pub energy_residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
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
    solver_type: Option<&str>,
) -> Result<AbsorptionColumnResult> {
    refuse_unported(
        murphree_efficiency,
        component_murphree_efficiency,
        max_allowable_gas_load_factor,
    )?;
    let solver = match solver_type.unwrap_or("direct_substitution") {
        "direct_substitution" => SolverType::DirectSubstitution,
        "naphtali_sandholm" => SolverType::NaphtaliSandholm,
        other => {
            return Err(AzothError::invalid_input(
                "solver_type",
                format!(
                    "solver_type = {other} is not ported: `ColumnSolverFactory.{}` is the class \
                     that would close it, and the strategies are path variants of one another \
                     rather than different physics - see `process.distillation_column`'s own \
                     `solver_type`, whose capture measures it.",
                    unported_class(other)
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
        temperature_tolerance,
        max_iterations,
        solver_type: solver,
    })?;

    Ok(AbsorptionColumnResult {
        tray_temperature: out.trays.iter().map(|tray| tray.temperature).collect(),
        tray_pressure: out.trays.iter().map(|tray| tray.pressure).collect(),
        tray_gas_n: out.trays.iter().map(|tray| tray.gas_n).collect(),
        tray_liquid_n: out.trays.iter().map(|tray| tray.liquid_n).collect(),
        gas_out_n: out.gas_out.n,
        gas_out_z: out.gas_out.z,
        gas_out_p: out.gas_out.p,
        gas_out_t: out.gas_out.t,
        gas_out_h: joules_per_mole(out.gas_out.h.value),
        liquid_out_n: out.liquid_out.n,
        liquid_out_z: out.liquid_out.z,
        liquid_out_p: out.liquid_out.p,
        liquid_out_t: out.liquid_out.t,
        liquid_out_h: joules_per_mole(out.liquid_out.h.value),
        iterations: out.iterations,
        temperature_residual: out.temperature_residual,
        mass_residual: out.mass_residual,
        energy_residual: out.energy_residual,
        warnings,
    })
}

/// Refuse every parameter the palette declares and this stage does not implement.
///
/// **Declared and refused, rather than withdrawn**, as `process.distillation_column`'s own
/// refusals are: a form field that errors with the class that would close it beats one that is
/// silently absent. The gas-load factor is the third shape: the class *accepts* it and no part
/// of `run` reads it, so the model carries it as a declaration the solve is indifferent to.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a Murphree efficiency of either kind.
fn refuse_unported(
    murphree_efficiency: Option<f64>,
    component_murphree_efficiency: Option<&[f64]>,
    max_allowable_gas_load_factor: Option<f64>,
) -> Result<()> {
    if let Some(efficiency) = murphree_efficiency {
        return Err(AzothError::invalid_input(
            "murphree_efficiency",
            format!(
                "a Murphree efficiency of {efficiency} is not ported: \
                 `SimpleTray.setMurphreeEfficiency` and the per-tray correction `applyMurphreeCorrection` \
                 applies are the classes that would close it. Omitted means the ideal stage, \
                 which is the class's own default of one"
            ),
        ));
    }
    if let Some(efficiencies) = component_murphree_efficiency {
        return Err(AzothError::invalid_input(
            "component_murphree_efficiency",
            format!(
                "{} component Murphree efficiencies are not ported: \
                 `AbsorptionColumn.setComponentMurphreeEfficiency(int, String, double)` and the \
                 `applyMurphreeCorrection` override it feeds are the classes that would close \
                 them",
                efficiencies.len()
            ),
        ));
    }
    // The design limit: read by `isGasLoadFactorWithinDesignLimit`, \
    // `getGasLoadFactorUtilization` and `getMinimumDiameterForGasLoadLimit`, and by nothing on
    // the run path - so it is accepted and the separation is indifferent to it.
    let _ = max_allowable_gas_load_factor;
    Ok(())
}
