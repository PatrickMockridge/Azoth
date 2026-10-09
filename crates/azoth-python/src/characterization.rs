//! The characterisation calculations, exposed to Python.
//!
//! Same rule as the other bindings: every parameter below is an SI magnitude, not a `uom`
//! quantity. Unit handling happens once, in Python, before the call crosses this boundary - one
//! conversion site used by both backends cannot disagree with itself about what a number is in.

use azoth_characterization::{TbpModel, tbp_cut_properties as kernel};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole, kelvins};
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::transport_gen::PyTbpCutPropertiesResult;

/// A TBP cut's critical properties, by any of NeqSim's ten models.
///
/// All arguments are SI magnitudes. `model` is one of the names the spec's enum declares;
/// `None` takes NeqSim's own default. An unrecognised name falls back to that default, which is
/// what the NeqSim class does - the spec's vocabulary check is what refuses.
#[pyfunction]
#[pyo3(signature = (model, molar_mass, density, boiling_point=None))]
#[pyo3(text_signature = "(model, molar_mass, density, boiling_point=None)")]
pub fn tbp_cut_properties(
    py: Python<'_>,
    model: Option<&str>,
    molar_mass: f64,
    density: f64,
    boiling_point: Option<f64>,
) -> PyResult<PyTbpCutPropertiesResult> {
    let parsed: Option<TbpModel> = model.map(|name| name.parse().unwrap_or_default());
    kernel(
        parsed,
        kilograms_per_mole(molar_mass),
        kilograms_per_cubic_meter(density),
        boiling_point.map(kelvins),
    )
    .map(|result| PyTbpCutPropertiesResult::from(&result))
    .map_err(|error| to_pyerr(py, error))
}
