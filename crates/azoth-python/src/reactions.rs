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
use crate::results::PyEquilibriumConstantResult;

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
