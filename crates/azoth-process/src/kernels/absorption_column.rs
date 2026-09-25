//! `unit_ops.absorption_column` - the tray absorber, which is the column's kernel twice fed.
//!
//! `AbsorptionColumn extends DistillationColumn` and overrides **no `run`**: its constructor
//! passes `false, false`, so it is the base's counter-current stage equations with no condenser
//! and no reboiler, and its two inlets are fixed - the gas at stage 0 and the solvent at the top
//! stage. So this kernel is the column's with a shape, and nothing here is arithmetic the
//! column's kernel does not already carry.
//!
//! **The class's own tests solve an isothermal column, and that is a degenerate one.** They pin
//! every stage with `SimpleTray.setOutletTemperature`, which makes `solveSequential`'s gate - the
//! mean tray-temperature change - zero, so the solve stops after its first sweep with no liquid
//! traffic at all. `validation/neqsim/captures/process_absorber.tsv` holds those rows as
//! evidence and the unpinned columns as the oracles; `tray_temperatures` is the mechanism they
//! use, and a caller who states it gets the same one-iteration state the class does.

use azoth_core::units::Pressure;
use azoth_core::{AzothError, Result};

use super::distillation_column::{ColumnOutcome, ColumnSetup, SolverType, TrayProfile};
use crate::stream::Stream;

/// What an absorber's solve hands back: the profile and the two products.
#[derive(Debug, Clone)]
pub struct AbsorberOutcome {
    /// The trays, numbered 0 at the bottom where the gas enters.
    pub trays: Vec<TrayProfile>,
    /// The gas leaving the top - `getGasOutStream`, or the stripper's overhead gas.
    pub gas_out: Stream,
    /// The liquid leaving the bottom - `getLiquidOutStream`, or the stripper's lean liquid.
    pub liquid_out: Stream,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against both feeds, relative.
    pub mass_residual: f64,
    /// What the stages had to fall back on, one entry per kind - see `ColumnOutcome`'s own.
    pub warnings: Vec<azoth_core::warning::Warning>,
    /// The enthalpy closure, which is open here - see the spec's assumptions.
    pub energy_residual: f64,
}

/// Everything an absorber's solve takes.
#[derive(Debug, Clone)]
pub struct AbsorberSetup {
    /// The gas, which enters stage 0.
    pub gas: Stream,
    /// The solvent, which enters the top stage.
    pub solvent: Stream,
    /// The trays, which the constructor requires to be positive.
    pub number_of_stages: usize,
    /// The pressure at the top stage.
    pub top_pressure: Pressure,
    /// The pressure at the bottom stage.
    pub bottom_pressure: Pressure,
    /// One outlet-temperature pin per tray, `NaN` where a tray has none.
    pub tray_temperatures: Option<Vec<f64>>,
    /// The convergence tolerance on the mean tray-temperature change.
    pub temperature_tolerance: f64,
    /// The iteration cap.
    pub max_iterations: usize,
    /// Which of the base's strategies to run.
    pub solver_type: SolverType,
}

/// Solve a tray absorber, or a stripper - the same equations with the inlets named otherwise.
///
/// # Errors
/// Whatever the column's kernel refuses, and a stage count of zero.
pub fn absorption_column(setup: &AbsorberSetup) -> Result<AbsorberOutcome> {
    if setup.number_of_stages == 0 {
        return Err(AzothError::invalid_input(
            "number_of_stages",
            "an absorber with no trays is not an absorber",
        ));
    }
    // **The two inlets are the class's own positions**: the gas at stage 0 and the solvent at
    // `getNumberOfTrays() - 1`, which is what `addGasInStream` and `addSolventInStream` call.
    let out: ColumnOutcome = super::distillation_column::distillation_column(&ColumnSetup {
        feed: setup.gas.clone(),
        feed_stage: 0,
        number_of_stages: setup.number_of_stages,
        has_reboiler: false,
        has_condenser: false,
        top_pressure: setup.top_pressure,
        bottom_pressure: setup.bottom_pressure,
        condenser_temperature: None,
        reboiler_temperature: None,
        temperature_tolerance: setup.temperature_tolerance,
        max_iterations: setup.max_iterations,
        top_specification: None,
        bottom_specification: None,
        top_feed: Some(setup.solvent.clone()),
        tray_temperatures: setup.tray_temperatures.clone(),
        solver_type: setup.solver_type,
        // **The absorber's own section is the class's**: `AbsorptionColumn` inherits
        // `setReactive`, and its model does not declare the input yet, so it is `None` here.
        reactive: super::distillation_column::ReactiveSection::None,
        // The absorber's own fractions are the class's too, and its model does not declare them.
        gas_side_draw_fractions: None,
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
    })?;

    Ok(AbsorberOutcome {
        trays: out.trays,
        warnings: out.warnings,
        gas_out: out.distillate,
        liquid_out: out.bottoms,
        iterations: out.iterations,
        temperature_residual: out.temperature_residual,
        mass_residual: out.mass_residual,
        energy_residual: out.energy_residual,
    })
}
