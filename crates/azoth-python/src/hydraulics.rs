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
    DynamicViscosity, cubic_meters_per_second, kilograms_per_cubic_meter, kilograms_per_second,
    meters, meters_per_second, newtons_per_meter, pascal_seconds, pascals,
};
use azoth_hydraulics as hyd;
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::{
    PyChokedFlowAreaResult, PyColebrookResult, PyControlValveCvResult, PyDarcyWeisbachResult,
    PyHaalandResult, PyKFactorsResult, PyOrificeFlowResult, PyPumpPowerResult,
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

/// Volumetric flow through an orifice.
///
/// `d` arrives as the SI **base** magnitude, so a 50 mm bore crosses as `0.05` - the
/// spec declares millimetres, but unit handling happens once in Python and what
/// crosses is metres. `Cd` is dimensionless and arrives as a plain float.
#[pyfunction]
#[pyo3(signature = (d, dP, rho, Cd))]
#[pyo3(text_signature = "(d, dP, rho, Cd)")]
#[allow(non_snake_case)] // `dP` and `Cd` are the symbols in the published equation
pub fn orifice_flow(
    py: Python<'_>,
    d: f64,
    dP: f64,
    rho: f64,
    Cd: f64,
) -> PyResult<PyOrificeFlowResult> {
    hyd::orifice_flow(meters(d), pascals(dP), kilograms_per_cubic_meter(rho), Cd)
        .map(|r| PyOrificeFlowResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Liquid flow through a control valve.
///
/// `Cv` and `SG` are dimensionless and arrive as plain floats. `dP` is the SI base
/// magnitude, so it crosses in pascals whatever unit the spec declares.
#[pyfunction]
#[pyo3(signature = (Cv, dP, SG))]
#[pyo3(text_signature = "(Cv, dP, SG)")]
#[allow(non_snake_case)] // `Cv` and `SG` are the symbols in the published relation
pub fn control_valve_cv(
    py: Python<'_>,
    Cv: f64,
    dP: f64,
    SG: f64,
) -> PyResult<PyControlValveCvResult> {
    hyd::control_valve_cv(Cv, pascals(dP), SG)
        .map(|r| PyControlValveCvResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Throat area required for a choked gas flow.
///
/// All arguments are SI magnitudes and `k` is dimensionless, so it arrives as a plain
/// float.
#[pyfunction]
#[pyo3(signature = (m_dot, P0, rho0, k))]
#[pyo3(text_signature = "(m_dot, P0, rho0, k)")]
#[allow(non_snake_case)] // `P0` and `rho0` are the symbols in the published relation
pub fn choked_flow_area(
    py: Python<'_>,
    m_dot: f64,
    P0: f64,
    rho0: f64,
    k: f64,
) -> PyResult<PyChokedFlowAreaResult> {
    hyd::choked_flow_area(
        kilograms_per_second(m_dot),
        pascals(P0),
        kilograms_per_cubic_meter(rho0),
        k,
    )
    .map(|r| PyChokedFlowAreaResult::from(&r))
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
/// A packed bed's hydraulics, from the packing's name and the state it runs at.
#[pyfunction]
#[pyo3(signature = (packing, column_diameter, packed_height, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, vapor_viscosity, liquid_viscosity, surface_tension, vapor_diffusivity, liquid_diffusivity, hydraulic_capacity_factor))]
#[pyo3(
    text_signature = "(packing, column_diameter, packed_height, vapor_mass_flow, liquid_mass_flow, vapor_density, liquid_density, vapor_viscosity, liquid_viscosity, surface_tension, vapor_diffusivity, liquid_diffusivity, hydraulic_capacity_factor)"
)]
#[allow(clippy::too_many_arguments)] // One per declared input, and a bed reads many.
pub fn packing_hydraulics(
    py: Python<'_>,
    packing: &str,
    column_diameter: f64,
    packed_height: f64,
    vapor_mass_flow: f64,
    liquid_mass_flow: f64,
    vapor_density: f64,
    liquid_density: f64,
    vapor_viscosity: f64,
    liquid_viscosity: f64,
    surface_tension: f64,
    vapor_diffusivity: f64,
    liquid_diffusivity: f64,
    hydraulic_capacity_factor: f64,
) -> PyResult<crate::results::PyPackingHydraulicsResult> {
    hyd::packing_hydraulics::packing_hydraulics(
        packing,
        hyd::packing_hydraulics::PackingState {
            column_diameter: meters(column_diameter),
            packed_height,
            vapor_mass_flow: kilograms_per_second(vapor_mass_flow),
            liquid_mass_flow: kilograms_per_second(liquid_mass_flow),
            vapor_density: kilograms_per_cubic_meter(vapor_density),
            liquid_density: kilograms_per_cubic_meter(liquid_density),
            vapor_viscosity: pascal_seconds(vapor_viscosity),
            liquid_viscosity: pascal_seconds(liquid_viscosity),
            surface_tension: newtons_per_meter(surface_tension),
            vapor_diffusivity,
            liquid_diffusivity,
            hydraulic_capacity_factor,
        },
    )
    .map(|r| crate::results::PyPackingHydraulicsResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

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
