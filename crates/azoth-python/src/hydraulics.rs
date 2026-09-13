//! The five calculations, exposed to Python.
//!
//! # Why these take plain floats
//!
//! Every parameter below is an SI magnitude, not a `uom` quantity. Unit handling
//! happens once, in Python, in `azoth.hydraulics`, before the call crosses this
//! boundary.
//!
//! That is deliberate rather than a shortcut. If the extension did its own unit
//! conversion, the two backends would each have a unit path, and the failure they
//! could then produce - agreeing on the numbers while disagreeing about what the
//! numbers are in - is precisely the kind of bug this project is built to
//! prevent. One conversion site, used by both backends, cannot disagree with
//! itself.
//!
//! The cost is that `azoth._core` is not units-safe on its own. It is a private
//! module: the public API is `azoth.hydraulics`, which is.

use azoth_core::units::{
    DynamicViscosity, cubic_meters_per_second, kilograms_per_cubic_meter, meters,
    meters_per_second, pascal_seconds,
};
use azoth_hydraulics as hyd;
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::{
    PyColebrookResult, PyDarcyWeisbachResult, PyHaalandResult, PyKFactorsResult, PyPumpPowerResult,
    PyReynoldsNumberResult, PySwameeJainResult,
};

/// Reynolds number for flow in a circular pipe.
///
/// All arguments are SI magnitudes. See the module documentation for why.
#[pyfunction]
#[pyo3(signature = (rho, v, D, mu))]
#[pyo3(text_signature = "(rho, v, D, mu)")]
#[allow(non_snake_case)] // symbols from the published equation
pub fn reynolds_number(
    py: Python<'_>,
    rho: f64,
    v: f64,
    D: f64,
    mu: f64,
) -> PyResult<PyReynoldsNumberResult> {
    hyd::reynolds_number(
        kilograms_per_cubic_meter(rho),
        meters_per_second(v),
        meters(D),
        pascal_seconds(mu),
    )
    .map(|r| PyReynoldsNumberResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// Solve the Colebrook-White equation for the Darcy friction factor.
#[pyfunction]
#[pyo3(signature = (re, relative_roughness))]
#[pyo3(text_signature = "(re, relative_roughness)")]
pub fn friction_factor_colebrook(
    py: Python<'_>,
    re: f64,
    relative_roughness: f64,
) -> PyResult<PyColebrookResult> {
    hyd::friction_factor_colebrook(re, relative_roughness)
        .map(|r| PyColebrookResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Explicit Swamee-Jain approximation to the Colebrook friction factor.
#[pyfunction]
#[pyo3(signature = (re, relative_roughness))]
#[pyo3(text_signature = "(re, relative_roughness)")]
pub fn friction_factor_swamee_jain(
    py: Python<'_>,
    re: f64,
    relative_roughness: f64,
) -> PyResult<PySwameeJainResult> {
    hyd::friction_factor_swamee_jain(re, relative_roughness)
        .map(|r| PySwameeJainResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Explicit Haaland approximation to the Colebrook friction factor.
#[pyfunction]
#[pyo3(signature = (re, relative_roughness))]
#[pyo3(text_signature = "(re, relative_roughness)")]
pub fn friction_factor_haaland(
    py: Python<'_>,
    re: f64,
    relative_roughness: f64,
) -> PyResult<PyHaalandResult> {
    hyd::friction_factor_haaland(re, relative_roughness)
        .map(|r| PyHaalandResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Shaft power a pump must be supplied with.
///
/// `eta` is dimensionless and arrives as a plain float. See the module
/// documentation for why every other argument is an SI magnitude.
#[pyfunction]
#[pyo3(signature = (rho, q, H, eta))]
#[pyo3(text_signature = "(rho, q, H, eta)")]
#[allow(non_snake_case)] // `H` is the symbol in the published equation
pub fn pump_power(
    py: Python<'_>,
    rho: f64,
    q: f64,
    H: f64,
    eta: f64,
) -> PyResult<PyPumpPowerResult> {
    hyd::pump_power(
        kilograms_per_cubic_meter(rho),
        cubic_meters_per_second(q),
        meters(H),
        eta,
    )
    .map(|r| PyPumpPowerResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// Total resistance coefficient for a list of fittings.
#[pyfunction]
#[pyo3(signature = (fittings, f_t))]
#[pyo3(text_signature = "(fittings, f_t)")]
pub fn crane_k_factors(
    py: Python<'_>,
    fittings: Vec<String>,
    f_t: f64,
) -> PyResult<PyKFactorsResult> {
    let ids: Vec<&str> = fittings.iter().map(String::as_str).collect();
    hyd::crane_k_factors(&ids, f_t)
        .map(|r| PyKFactorsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Pressure drop over a length of straight pipe.
///
/// `mu` is optional: it is needed only to check the flow regime. Omit it and the
/// result says the regime went unchecked rather than implying it passed.
#[pyfunction]
#[pyo3(signature = (f, L, D, rho, v, mu=None))]
#[pyo3(text_signature = "(f, L, D, rho, v, mu=None)")]
#[allow(non_snake_case)] // `L` and `D` are the symbols in the published equation
pub fn darcy_weisbach(
    py: Python<'_>,
    f: f64,
    L: f64,
    D: f64,
    rho: f64,
    v: f64,
    mu: Option<f64>,
) -> PyResult<PyDarcyWeisbachResult> {
    let viscosity: Option<DynamicViscosity> = mu.map(pascal_seconds);
    hyd::darcy_weisbach(
        f,
        meters(L),
        meters(D),
        kilograms_per_cubic_meter(rho),
        meters_per_second(v),
        viscosity,
    )
    .map(|r| PyDarcyWeisbachResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
