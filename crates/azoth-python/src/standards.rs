//! The standards namespace, exposed to Python.
//!
//! Thin by design: this is a binding to `azoth-standards`, not a second implementation.

use azoth_core::units::kelvins;
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::PyIso6976Result;

/// `standards.iso6976` - the calorific values and density of a natural gas.
///
/// **The component names cross unresolved**, as the process models' do: the Rust side
/// resolves each against the standard's own table rather than the component databank, so the
/// two languages cannot disagree about which row answered.
#[pyfunction]
#[pyo3(signature = (components, z, volumetric_reference_temperature, energy_reference_temperature))]
#[pyo3(
    text_signature = "(components, z, volumetric_reference_temperature, energy_reference_temperature)"
)]
pub fn iso6976(
    py: Python<'_>,
    components: Vec<String>,
    z: Vec<f64>,
    volumetric_reference_temperature: f64,
    energy_reference_temperature: f64,
) -> PyResult<PyIso6976Result> {
    azoth_standards::iso6976(
        &components,
        &z,
        kelvins(volumetric_reference_temperature),
        kelvins(energy_reference_temperature),
    )
    .map(|r| PyIso6976Result::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
