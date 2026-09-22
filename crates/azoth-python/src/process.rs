//! The process layer, exposed to Python.
//!
//! Thin by design: this is a binding to `azoth-process`, not a second
//! implementation. The `Stream` value and the kernels run in Rust; Python supplies
//! SI magnitudes at the boundary and reads the streams back.

use std::path::Path;

use azoth_core::units::{joules_per_mole, kelvins, pascals, watts, watts_per_kelvin};
use azoth_process::{self, Stream};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

use crate::errors::to_pyerr;

/// A material stream, as the process layer carries it. All quantities are SI
/// magnitudes: `p` in pascals, `t` in kelvin, `h` in J/mol, `n` in mol/s.
#[pyclass(frozen, module = "azoth._core", name = "Stream")]
pub struct PyStream {
    #[pyo3(get)]
    pub components: Vec<String>,
    #[pyo3(get)]
    pub z: Vec<f64>,
    #[pyo3(get)]
    pub n: f64,
    #[pyo3(get)]
    pub p: f64,
    #[pyo3(get)]
    pub t: f64,
    #[pyo3(get)]
    pub h: f64,
}

impl PyStream {
    fn from_inner(inner: Stream) -> Self {
        PyStream {
            components: inner.components,
            z: inner.z,
            n: inner.n,
            p: inner.p.value,
            t: inner.t.value,
            h: inner.h.value,
        }
    }

    fn to_stream(&self) -> Stream {
        Stream {
            components: self.components.clone(),
            z: self.z.clone(),
            n: self.n,
            p: pascals(self.p),
            t: kelvins(self.t),
            h: joules_per_mole(self.h),
        }
    }
}

#[pymethods]
impl PyStream {
    /// A stream at a known pressure and temperature, its enthalpy computed.
    #[new]
    #[pyo3(signature = (components, z, n, p, t))]
    fn new(
        py: Python<'_>,
        components: Vec<String>,
        z: Vec<f64>,
        n: f64,
        p: f64,
        t: f64,
    ) -> PyResult<Self> {
        Stream::from_pt(components, z, n, pascals(p), kelvins(t))
            .map(PyStream::from_inner)
            .map_err(|e| to_pyerr(py, e))
    }

    /// A stream at a known pressure and molar enthalpy, its temperature solved.
    #[staticmethod]
    #[pyo3(signature = (components, z, n, p, h))]
    fn from_ph(
        py: Python<'_>,
        components: Vec<String>,
        z: Vec<f64>,
        n: f64,
        p: f64,
        h: f64,
    ) -> PyResult<Self> {
        Stream::from_ph(components, z, n, pascals(p), joules_per_mole(h))
            .map(PyStream::from_inner)
            .map_err(|e| to_pyerr(py, e))
    }
}

fn wrap_streams(streams: Vec<Stream>) -> Vec<PyStream> {
    streams.into_iter().map(PyStream::from_inner).collect()
}

/// Split a stream into several with the same state, scaled by `fractions`.
#[pyfunction]
pub fn splitter_stream(
    py: Python<'_>,
    feed: &PyStream,
    fractions: Vec<f64>,
) -> PyResult<Vec<PyStream>> {
    azoth_process::kernels::splitter(&feed.to_stream(), &fractions)
        .map(wrap_streams)
        .map_err(|e| to_pyerr(py, e))
}

/// Join several inlets into one, conserving molar flow and enthalpy.
#[pyfunction]
#[pyo3(signature = (inlets, outlet_pressure = None))]
pub fn mixer_stream(
    py: Python<'_>,
    inlets: Vec<Py<PyStream>>,
    outlet_pressure: Option<f64>,
) -> PyResult<PyStream> {
    let streams: Vec<Stream> = inlets.iter().map(|s| s.borrow(py).to_stream()).collect();
    azoth_process::kernels::mixer(&streams, outlet_pressure.map(pascals))
        .map(PyStream::from_inner)
        .map_err(|e| to_pyerr(py, e))
}

/// Flash a stream into vapour and liquid outlets (SI: Pa, W, dimensionless).
#[pyfunction]
#[pyo3(signature = (feed, pressure_drop, gas_in_liquid = 0.0, heat_input = None))]
pub fn separator_stream(
    py: Python<'_>,
    feed: &PyStream,
    pressure_drop: f64,
    gas_in_liquid: f64,
    heat_input: Option<f64>,
) -> PyResult<(PyStream, PyStream)> {
    azoth_process::kernels::separator(
        &feed.to_stream(),
        pascals(pressure_drop),
        gas_in_liquid,
        heat_input.map(watts),
    )
    .map(|(v, l)| (PyStream::from_inner(v), PyStream::from_inner(l)))
    .map_err(|e| to_pyerr(py, e))
}

/// Drop a stream to `outlet_pressure` without heat or work (SI, Pa).
#[pyfunction]
pub fn throttling_valve_stream(
    py: Python<'_>,
    feed: &PyStream,
    outlet_pressure: f64,
) -> PyResult<PyStream> {
    azoth_process::kernels::throttling_valve(&feed.to_stream(), pascals(outlet_pressure))
        .map(PyStream::from_inner)
        .map_err(|e| to_pyerr(py, e))
}

/// Exchange heat between two streams (SI: W/K, K).
#[pyfunction]
#[pyo3(signature = (hot, cold, ua = None, flow_arrangement = "counterflow", hot_outlet_temperature = None, cold_outlet_temperature = None))]
pub fn heat_exchanger_stream(
    py: Python<'_>,
    hot: &PyStream,
    cold: &PyStream,
    ua: Option<f64>,
    flow_arrangement: &str,
    hot_outlet_temperature: Option<f64>,
    cold_outlet_temperature: Option<f64>,
) -> PyResult<(PyStream, PyStream)> {
    azoth_process::kernels::heat_exchanger(
        &hot.to_stream(),
        &cold.to_stream(),
        ua.map(watts_per_kelvin),
        flow_arrangement,
        hot_outlet_temperature.map(kelvins),
        cold_outlet_temperature.map(kelvins),
    )
    .map(|(h, c)| (PyStream::from_inner(h), PyStream::from_inner(c)))
    .map_err(|e| to_pyerr(py, e))
}

/// Raise a liquid stream to `outlet_pressure` (SI, Pa) at `efficiency`.
///
/// **The Stream-level kernel, not the `process.pump` id.** The two are one arithmetic under
/// two call shapes: this takes a `Stream` and returns one, which is what a flowsheet's
/// connection carries, while the id takes the record field by field so that a case, a
/// cross-impl test and a NeqSim capture can address it. The suffix says which is which.
#[pyfunction]
pub fn pump_stream(
    py: Python<'_>,
    feed: &PyStream,
    outlet_pressure: f64,
    efficiency: f64,
) -> PyResult<PyStream> {
    azoth_process::kernels::pump(&feed.to_stream(), pascals(outlet_pressure), efficiency)
        .map(PyStream::from_inner)
        .map_err(|e| to_pyerr(py, e))
}

/// `process.pump` - the pump's kernel as a registered id.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure, isentropic_efficiency))]
#[pyo3(
    text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure, isentropic_efficiency)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seven
pub fn pump(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    outlet_pressure: f64,
    isentropic_efficiency: f64,
) -> PyResult<crate::results::PyPumpResult> {
    azoth_process::pump(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        pascals(outlet_pressure),
        isentropic_efficiency,
    )
    .map(|r| crate::results::PyPumpResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.throttling_valve` - the valve's kernel as a registered id.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure))]
#[pyo3(text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure)")]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are six
pub fn throttling_valve(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    outlet_pressure: f64,
) -> PyResult<crate::results::PyThrottlingValveResult> {
    azoth_process::throttling_valve(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        pascals(outlet_pressure),
    )
    .map(|r| crate::results::PyThrottlingValveResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.separator` - the separator's kernel as a registered id.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, pressure_drop, gas_in_liquid, heat_input = None))]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, pressure_drop, gas_in_liquid, heat_input=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn separator(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    pressure_drop: f64,
    gas_in_liquid: f64,
    heat_input: Option<f64>,
) -> PyResult<crate::results::PySeparatorResult> {
    azoth_process::separator(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        pascals(pressure_drop),
        gas_in_liquid,
        heat_input.map(watts),
    )
    .map(|r| crate::results::PySeparatorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.mixer` - the mixer's kernel as a registered id.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, outlet_pressure = None))]
#[pyo3(text_signature = "(components, feed_n, feed_z, feed_p, feed_t, outlet_pressure=None)")]
#[allow(non_snake_case)] // the record's own field names
pub fn mixer(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: Vec<f64>,
    feed_z: Vec<Vec<f64>>,
    feed_p: Vec<f64>,
    feed_t: Vec<f64>,
    outlet_pressure: Option<f64>,
) -> PyResult<crate::results::PyMixerResult> {
    let pressures: Vec<azoth_core::units::Pressure> = feed_p.into_iter().map(pascals).collect();
    let temperatures: Vec<azoth_core::units::ThermodynamicTemperature> =
        feed_t.into_iter().map(kelvins).collect();
    azoth_process::mixer(
        &components,
        &feed_n,
        &feed_z,
        &pressures,
        &temperatures,
        outlet_pressure.map(pascals),
    )
    .map(|r| crate::results::PyMixerResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.splitter` - the splitter's kernel as a registered id.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, split_factors))]
#[pyo3(text_signature = "(components, feed_n, feed_z, feed_p, feed_t, split_factors)")]
#[allow(non_snake_case)] // the record's own field names
pub fn splitter(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    split_factors: Vec<f64>,
) -> PyResult<crate::results::PySplitterResult> {
    azoth_process::splitter(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        &split_factors,
    )
    .map(|r| crate::results::PySplitterResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// Validate a flowsheet's TOML against a palette directory, returning the
/// diagnostic lines (empty when the flowsheet is clean).
#[pyfunction]
pub fn validate_flowsheet(flowsheet: &str, palette_dir: &str) -> PyResult<Vec<String>> {
    let palette =
        azoth_process::load_palette(Path::new(palette_dir)).map_err(PyValueError::new_err)?;
    let sheet = azoth_process::parse_flowsheet(flowsheet).map_err(PyValueError::new_err)?;

    let mut lines = Vec::new();
    for d in azoth_process::validate_palette(&palette) {
        lines.push(format!("palette: {d:?}"));
    }
    for d in azoth_process::validate(&sheet, &palette) {
        lines.push(format!("{d:?}"));
    }
    Ok(lines)
}
