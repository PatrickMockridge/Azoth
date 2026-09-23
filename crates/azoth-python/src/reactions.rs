//! The reactions namespace's extension surface.
//!
//! One function per registered `reactions.*` calc, thin: extract the SI magnitude,
//! call the kernel, wrap the result. The spec-driven behaviour lives in
//! `azoth-reactions`; what is here is the boundary.
//!
//! **The reaction's name crosses unresolved**, as a string. The Rust side reads the
//! table itself - which source, which row, which coefficients - so the Python caller
//! never handles a coefficient and the two languages cannot disagree about which row
//! answered.

use pyo3::prelude::*;

use azoth_reactions::chemical_equilibrium::ConcentrationBasis;
use azoth_reactions::databank::ReactionDataSource;
use azoth_reactions::reactive_phase_equilibrium::ReactionSeed;

use crate::errors::to_pyerr;
use crate::results::{
    PyEquilibriumConstantResult, PyReactivePhaseEquilibriumResult, PyReactiveTpFlashResult,
};

/// One reaction's equilibrium constant, its derivative and its heat of reaction.
///
/// `source` is parsed rather than matched here, so a source this library does not carry
/// is refused by the same code path the Rust tests exercise.
#[pyfunction]
#[pyo3(signature = (source, reaction, T))]
#[pyo3(text_signature = "(source, reaction, T)")]
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn equilibrium_constant(
    py: Python<'_>,
    source: &str,
    reaction: &str,
    T: f64,
) -> PyResult<PyEquilibriumConstantResult> {
    let parsed: ReactionDataSource = source.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::equilibrium_constant::equilibrium_constant(
        parsed,
        reaction,
        azoth_core::units::kelvins(T),
    )
    .map(|r| PyEquilibriumConstantResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The standard-state reference potentials of a fluid's reactive components.
///
/// `components` crosses as names and unresolved, for the same reason the reaction name
/// does: the Rust side reads the tables, chooses the basis and propagates, so no
/// stoichiometric coefficient and no rank decision reaches Python.
#[pyfunction]
#[pyo3(signature = (components, source, T))]
#[pyo3(text_signature = "(components, source, T)")]
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn reference_potentials(
    py: Python<'_>,
    components: Vec<String>,
    source: &str,
    T: f64,
) -> PyResult<crate::results::PyReferencePotentialsResult> {
    let parsed: ReactionDataSource = source.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::reference_potentials::reference_potentials(
        &components,
        parsed,
        azoth_core::units::kelvins(T),
    )
    .map(|r| crate::results::PyReferencePotentialsResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The reactive equilibrium composition of a phase.
///
/// **The matrix crosses nested**, so the parameter list is exactly the spec's declared
/// inputs - a flat vector plus a row count would be two parameters where the spec has one,
/// and the generated stub would disagree with this signature.
///
/// `whole_system` is NeqSim's `getNumberOfPhases() == 1`, which decides whether the solve
/// corrects its conservation coupling to the phase's own element amounts.
#[pyfunction]
#[pyo3(signature = (a_matrix, b, whole_system, moles, chem_ref, log_activity, T, max_iterations, tolerance, concentration_basis, solvent_weight, solvent_mask, phase_moles))]
#[pyo3(
    text_signature = "(a_matrix, b, whole_system, moles, chem_ref, log_activity, T, max_iterations, tolerance, concentration_basis, solvent_weight, solvent_mask, phase_moles)"
)]
#[allow(non_snake_case)] // `T` is the symbol in the chemistry
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twelve
pub fn chemical_equilibrium(
    py: Python<'_>,
    a_matrix: Vec<Vec<f64>>,
    b: Vec<f64>,
    whole_system: bool,
    moles: Vec<f64>,
    chem_ref: Vec<f64>,
    log_activity: Vec<f64>,
    T: f64,
    max_iterations: u32,
    tolerance: f64,
    concentration_basis: &str,
    solvent_weight: f64,
    solvent_mask: Vec<f64>,
    phase_moles: f64,
) -> PyResult<crate::results::PyChemicalEquilibriumResult> {
    let parsed_basis: ConcentrationBasis =
        concentration_basis.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::chemical_equilibrium::chemical_equilibrium(
        &a_matrix,
        &b,
        whole_system,
        &moles,
        &chem_ref,
        &log_activity,
        T,
        max_iterations,
        tolerance,
        parsed_basis,
        solvent_weight,
        &solvent_mask,
        phase_moles,
    )
    .map(|r| crate::results::PyChemicalEquilibriumResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The phase's reactive equilibrium, as the operation the facade is.
///
/// **`skipped` is part of the result and not an error.** A phase that is neither aqueous
/// nor liquid nor oil has no phase for a water-based equilibrium to be solved in, and the
/// composition comes back as it went in - so the flag is what separates that from a solve
/// that ran and did not converge.
#[pyfunction]
#[pyo3(
    signature = (components, source, phase, moles, phase_charge, phase_moles, whole_system, log_activity, T, max_iterations, tolerance, seed)
)]
#[pyo3(
    text_signature = "(components, source, phase, moles, phase_charge, phase_moles, whole_system, log_activity, T, max_iterations, tolerance, seed)"
)]
#[allow(non_snake_case)] // `T` is the symbol in the chemistry
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twelve
pub fn reactive_phase_equilibrium(
    py: Python<'_>,
    components: Vec<String>,
    source: &str,
    phase: &str,
    moles: Vec<f64>,
    phase_charge: f64,
    phase_moles: f64,
    whole_system: bool,
    log_activity: Vec<f64>,
    T: f64,
    max_iterations: u32,
    tolerance: f64,
    seed: &str,
) -> PyResult<PyReactivePhaseEquilibriumResult> {
    let parsed: ReactionDataSource = source.parse().map_err(|e| to_pyerr(py, e))?;
    let parsed_seed: ReactionSeed = seed.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::reactive_phase_equilibrium::reactive_phase_equilibrium(
        &components,
        parsed,
        phase,
        &moles,
        phase_charge,
        phase_moles,
        whole_system,
        &log_activity,
        T,
        max_iterations,
        tolerance,
        parsed_seed,
    )
    .map(|r| PyReactivePhaseEquilibriumResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The reactive flash: simultaneous chemical and phase equilibrium at fixed `T` and `P`.
///
/// `moles` is the **overall component amounts and not a composition**, which is what the
/// driver's `getOverallMoles` reads and what its frozen element inventory is built from.
/// A charged component is refused: the RAND solve's ionic branch is not ported.
#[pyfunction]
#[pyo3(signature = (components, T, P, moles, max_phases))]
#[pyo3(text_signature = "(components, T, P, moles, max_phases)")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the flash's own name
pub fn reactive_tp_flash(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    moles: Vec<f64>,
    max_phases: f64,
) -> PyResult<PyReactiveTpFlashResult> {
    azoth_reactions::reactive_tp_flash::reactive_tp_flash(
        &components,
        T,
        P,
        &moles,
        max_phases as usize,
    )
    .map(|r| PyReactiveTpFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
