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

use azoth_reactions::databank::ReactionDataSource;

use crate::errors::to_pyerr;
use crate::transport_gen::{
    PyEquilibriumConstantResult, PyReactivePhFlashResult, PyReactiveTpFlashResult,
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
) -> PyResult<crate::transport_gen::PyReferencePotentialsResult> {
    let parsed: ReactionDataSource = source.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::reference_potentials::reference_potentials(
        &components,
        parsed,
        azoth_core::units::kelvins(T),
    )
    .map(|r| crate::transport_gen::PyReferencePotentialsResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The reactive flash: simultaneous chemical and phase equilibrium at fixed `T` and `P`.
///
/// `moles` is the **overall component amounts and not a composition**, which is what the
/// driver's `getOverallMoles` reads and what its frozen element inventory is built from.
/// A charged component is refused: the RAND solve's ionic branch is not ported.
#[pyfunction]
#[pyo3(signature = (components, T, P, moles, max_phases, cubic = "srk"))]
#[pyo3(text_signature = "(components, T, P, moles, max_phases, cubic = \"srk\")")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the flash's own name
pub fn reactive_tp_flash(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    moles: Vec<f64>,
    max_phases: f64,
    cubic: &str,
) -> PyResult<PyReactiveTpFlashResult> {
    azoth_reactions::reactive_tp_flash::reactive_tp_flash(
        &components,
        cubic
            .parse()
            .map_err(pyo3::exceptions::PyValueError::new_err)?,
        T,
        P,
        &moles,
        max_phases as usize,
    )
    .map(|r| PyReactiveTpFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The reactive PH flash: the temperature a specified enthalpy asks for.
///
/// `enthalpy` is the **thermochemical** specification - the fluid's sensible enthalpy plus the
/// formation inventory - and `T` is where the secant search starts.
#[pyfunction]
#[pyo3(signature = (components, T, P, moles, enthalpy, max_phases, cubic = "srk"))]
#[pyo3(text_signature = "(components, T, P, moles, enthalpy, max_phases, cubic = \"srk\")")]
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seven
#[allow(non_snake_case)] // `T` and `P` are the symbols in the flash's own name
pub fn reactive_ph_flash(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    moles: Vec<f64>,
    enthalpy: f64,
    max_phases: f64,
    cubic: &str,
) -> PyResult<PyReactivePhFlashResult> {
    azoth_reactions::reactive_ph_flash::reactive_ph_flash(
        &components,
        cubic
            .parse()
            .map_err(pyo3::exceptions::PyValueError::new_err)?,
        T,
        P,
        &moles,
        enthalpy,
        max_phases as usize,
    )
    .map(|r| PyReactivePhFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
