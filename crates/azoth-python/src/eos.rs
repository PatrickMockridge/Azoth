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
use crate::results::{PyPrAlphaAbResult, PyPrKappaResult, PyPrZFactorResult, PyPrsvKappaResult};

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

/// The Peng-Robinson alpha function and the reduced attraction parameters.
///
/// `kappa` comes from `pr_kappa`; `Tr` and `Pr` are reduced against the caller's
/// critical point. All six quantities are dimensionless, so nothing here touches
/// units in either direction.
#[pyfunction]
#[pyo3(signature = (kappa, Tr, Pr))]
#[pyo3(text_signature = "(kappa, Tr, Pr)")]
#[allow(non_snake_case)] // `Tr` and `Pr` are the symbols in the published equation
pub fn pr_alpha_ab(py: Python<'_>, kappa: f64, Tr: f64, Pr: f64) -> PyResult<PyPrAlphaAbResult> {
    eos::pr_alpha_ab(kappa, Tr, Pr)
        .map(|r| PyPrAlphaAbResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Peng-Robinson compressibility factor.
///
/// Both arguments are the cubic's dimensionless parameters, so nothing here
/// touches units. `root_structure` crosses as the spec's spelling and the bridge
/// rebuilds the enum - the same arrangement `reynolds_number`'s `regime` uses.
#[pyfunction]
#[pyo3(signature = (a_reduced, b_reduced))]
#[pyo3(text_signature = "(a_reduced, b_reduced)")]
pub fn pr_z_factor(py: Python<'_>, a_reduced: f64, b_reduced: f64) -> PyResult<PyPrZFactorResult> {
    eos::pr_z_factor(a_reduced, b_reduced)
        .map(|r| PyPrZFactorResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The PRSV alpha-function coefficient.
///
/// All three arguments are dimensionless, and `kappa1` is the caller's - this
/// library ships no fitted values for it, which is the one thing about PRSV a
/// caller is most likely to expect and not get.
#[pyfunction]
#[pyo3(signature = (omega, Tr, kappa1))]
#[pyo3(text_signature = "(omega, Tr, kappa1)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn prsv_kappa(py: Python<'_>, omega: f64, Tr: f64, kappa1: f64) -> PyResult<PyPrsvKappaResult> {
    eos::prsv_kappa(omega, Tr, kappa1)
        .map(|r| PyPrsvKappaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}
