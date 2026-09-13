//! The equations-of-state calculations, exposed to Python.
//!
//! Same rule as `hydraulics.rs` and `thermal.rs`: every parameter below is an SI
//! magnitude, not a `uom` quantity. Here that rule costs nothing, because the
//! quantities in this namespace are dimensionless by construction - the reduced
//! variables and constitutive coefficients an equation of state is written in -
//! so there is no conversion to do at the boundary and no unit string to keep in
//! step with the spec.

use azoth_core::units::{cubic_meters_per_mole, kelvins, kilograms_per_mole, pascals};
use azoth_eos as eos;
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::{
    PyPrAlphaAbResult, PyPrDepartureResult, PyPrKappaResult, PyPrMassDensityResult,
    PyPrMolarVolumeResult, PyPrZFactorResult, PyPrsvKappaResult, PyPureSaturationResult,
    PyRachfordRiceBinaryResult, PyVdw1fMixBinaryResult,
};

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

/// The Peng-Robinson fugacity coefficient and departure functions.
///
/// Five dimensionless arguments and three dimensionless outputs. `kappa` may come
/// from either `pr_kappa` or `prsv_kappa` - this calc uses it only through the
/// logarithmic derivative of the alpha function, which both correlations feed.
#[pyfunction]
#[pyo3(signature = (a_reduced, b_reduced, z, kappa, Tr))]
#[pyo3(text_signature = "(a_reduced, b_reduced, z, kappa, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_departure(
    py: Python<'_>,
    a_reduced: f64,
    b_reduced: f64,
    z: f64,
    kappa: f64,
    Tr: f64,
) -> PyResult<PyPrDepartureResult> {
    eos::pr_departure(a_reduced, b_reduced, z, kappa, Tr)
        .map(|r| PyPrDepartureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The van der Waals one-fluid mixture parameters for a binary.
///
/// Six dimensionless arguments and two dimensionless outputs. `k12` is the
/// caller's - this library ships no fitted binary parameters.
#[pyfunction]
#[pyo3(signature = (z1, a1, a2, b1, b2, k12))]
#[pyo3(text_signature = "(z1, a1, a2, b1, b2, k12)")]
pub fn vdw1f_mix_binary(
    py: Python<'_>,
    z1: f64,
    a1: f64,
    a2: f64,
    b1: f64,
    b2: f64,
    k12: f64,
) -> PyResult<PyVdw1fMixBinaryResult> {
    eos::vdw1f_mix_binary(z1, a1, a2, b1, b2, k12)
        .map(|r| PyVdw1fMixBinaryResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The binary Rachford-Rice vapour fraction.
///
/// A `beta` outside `[0, 1]` crosses as a value carrying an `OUT_OF_VALID_RANGE`
/// warning rather than as an error - the feed is single phase, which is a real
/// answer and not a failure.
#[pyfunction]
#[pyo3(signature = (z1, K1, K2))]
#[pyo3(text_signature = "(z1, K1, K2)")]
#[allow(non_snake_case)] // `K1` and `K2` are the symbols in the published equation
pub fn rachford_rice_binary(
    py: Python<'_>,
    z1: f64,
    K1: f64,
    K2: f64,
) -> PyResult<PyRachfordRiceBinaryResult> {
    eos::rachford_rice_binary(z1, K1, K2)
        .map(|r| PyRachfordRiceBinaryResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Molar volume from a compressibility factor.
///
/// The one dimensional calc in this namespace, so unlike its neighbours this takes
/// and returns SI magnitudes with units attached - `T` in kelvin and `P` in pascals
/// in, a molar volume out.
#[pyfunction]
#[pyo3(signature = (z, T, P))]
#[pyo3(text_signature = "(z, T, P)")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the published equation
pub fn pr_molar_volume(py: Python<'_>, z: f64, T: f64, P: f64) -> PyResult<PyPrMolarVolumeResult> {
    eos::pr_molar_volume(z, kelvins(T), pascals(P))
        .map(|r| PyPrMolarVolumeResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Mass density from a molar mass and a molar volume.
///
/// `M` is in kg/mol - not the g/mol a table would quote - and `v` in m**3/mol.
/// Both are SI magnitudes crossing this boundary, converted once in Python.
#[pyfunction]
#[pyo3(signature = (M, v))]
#[pyo3(text_signature = "(M, v)")]
#[allow(non_snake_case)] // `M` is the symbol in the equation
pub fn pr_mass_density(py: Python<'_>, M: f64, v: f64) -> PyResult<PyPrMassDensityResult> {
    eos::pr_mass_density(kilograms_per_mole(M), cubic_meters_per_mole(v))
        .map(|r| PyPrMassDensityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The saturation pressure of a pure component.
///
/// A *model* rather than a calculation - its spec fixes a procedure and it composes
/// the kernels above rather than adding arithmetic of its own. It crosses the
/// boundary like any other function: SI magnitudes in, a quantity out.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, T))]
#[pyo3(text_signature = "(Tc, Pc, omega, T)")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the chemistry
pub fn pure_saturation(
    py: Python<'_>,
    Tc: f64,
    Pc: f64,
    omega: f64,
    T: f64,
) -> PyResult<PyPureSaturationResult> {
    eos::pure_saturation(kelvins(Tc), pascals(Pc), omega, kelvins(T))
        .map(|r| PyPureSaturationResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Every model id `azoth-eos` implements.
///
/// Separate from `calc_ids()` on purpose: the calc registry's id list is asserted to
/// be *exactly* the specs under `specs/calcs/`, so a model appearing there would
/// break that contract rather than extend it. Models have their own list, generated
/// from their own tree.
#[pyfunction]
#[must_use]
pub fn model_ids() -> Vec<String> {
    eos::model_gen::models()
        .iter()
        .map(|m| m.id.to_string())
        .collect()
}

/// The algorithm scheme each model runs, by model id.
///
/// The third leg of the same contract `solver_kinds()` provides for solvers: the
/// schema's `algorithm.scheme`, the implementations' own names, and this list must
/// agree, and `test_model_contract.py` asserts it rather than trusting three
/// hand-edited lists to stay in step.
#[pyfunction]
#[must_use]
pub fn model_schemes(model_id: &str) -> Vec<String> {
    match eos::model_gen::model(model_id) {
        Some(spec) => vec![spec.algorithm.scheme.to_string()],
        None => Vec::new(),
    }
}
