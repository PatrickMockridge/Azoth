//! The equations-of-state calculations, exposed to Python.
//!
//! Same rule as `hydraulics.rs` and `thermal.rs`: every parameter below is an SI
//! magnitude, not a `uom` quantity. Here that rule costs nothing, because the
//! quantities in this namespace are dimensionless by construction - the reduced
//! variables and constitutive coefficients an equation of state is written in -
//! so there is no conversion to do at the boundary and no unit string to keep in
//! step with the spec.

use azoth_eos as eos;
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::PyPrKappaResult;

/// The Peng-Robinson alpha-function coefficient.
///
/// `omega` is the acentric factor. A `kappa` below zero is returned carrying an
/// `OUT_OF_VALID_RANGE` warning rather than refused - see the calc's own
/// documentation for why that is the honest answer rather than a refusal.
#[pyfunction]
#[pyo3(signature = (omega))]
#[pyo3(text_signature = "(omega)")]
pub fn pr_kappa(py: Python<'_>, omega: f64) -> PyResult<PyPrKappaResult> {
    eos::pr_kappa(omega)
        .map(|r| PyPrKappaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}
