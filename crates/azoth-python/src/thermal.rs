//! The thermal calculations, exposed to Python.
//!
//! Same rule as `hydraulics.rs`: every parameter below is an SI magnitude, not a
//! `uom` quantity. Unit handling happens once, in Python, before the call crosses
//! this boundary. One conversion site used by both backends cannot disagree with
//! itself about what a number is in, which is the failure this arrangement exists
//! to make impossible - and which the `mm` entry in the unit vocabulary did
//! produce, silently, for a unit nothing was using.

use azoth_core::units::{kelvin_intervals, meters, square_meters, watts_per_meter_kelvin};
use azoth_thermal as therm;
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::PyConductionPlaneWallResult;

/// Steady conduction through a plane wall.
///
/// All arguments are SI magnitudes. `dT` is a temperature *difference* rather than
/// an absolute temperature, which is why it is built as a `TemperatureInterval` -
/// conflating the two would add 273.15 without saying so.
#[pyfunction]
#[pyo3(signature = (k, A, dT, L))]
#[pyo3(text_signature = "(k, A, dT, L)")]
#[allow(non_snake_case)] // `A`, `dT` and `L` are the symbols in the published equation
pub fn conduction_plane_wall(
    py: Python<'_>,
    k: f64,
    A: f64,
    dT: f64,
    L: f64,
) -> PyResult<PyConductionPlaneWallResult> {
    therm::conduction_plane_wall(
        watts_per_meter_kelvin(k),
        square_meters(A),
        kelvin_intervals(dT),
        meters(L),
    )
    .map(|r| PyConductionPlaneWallResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
