//! The process layer, exposed to Python.
//!
//! Thin by design: this is a binding to `azoth-process`, not a second
//! implementation. The `Stream` value and the kernels run in Rust; Python supplies
//! SI magnitudes at the boundary and reads the streams back.

use std::path::Path;

use azoth_core::units::{joules_per_mole, kelvins, meters, pascals, watts, watts_per_kelvin};
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

/// `process.compressor` - the isentropic route as a registered id.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure, isentropic_efficiency))]
#[pyo3(
    text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure, isentropic_efficiency)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seven
pub fn compressor(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    outlet_pressure: f64,
    isentropic_efficiency: f64,
) -> PyResult<crate::results::PyCompressorResult> {
    azoth_process::compressor(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        pascals(outlet_pressure),
        isentropic_efficiency,
    )
    .map(|r| crate::results::PyCompressorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.expander` - the same route with the efficiency multiplying.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure, isentropic_efficiency))]
#[pyo3(
    text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_pressure, isentropic_efficiency)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are seven
pub fn expander(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    outlet_pressure: f64,
    isentropic_efficiency: f64,
) -> PyResult<crate::results::PyExpanderResult> {
    azoth_process::expander(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        pascals(outlet_pressure),
        isentropic_efficiency,
    )
    .map(|r| crate::results::PyExpanderResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.filter` - the filter's kernel as a registered id.
///
/// The drop is required: the class defaults `deltaP` to `0.01` bar, which the palette entry
/// does not declare, and inventing a default here would be answering a question the caller
/// did not ask.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, pressure_drop))]
#[pyo3(text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, pressure_drop)")]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are six
pub fn filter(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    pressure_drop: f64,
) -> PyResult<crate::results::PyFilterResult> {
    azoth_process::filter(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        pascals(pressure_drop),
    )
    .map(|r| crate::results::PyFilterResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.pipe` - the line's kernel as a registered id.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, length, diameter, roughness))]
#[pyo3(
    text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, length, diameter, roughness)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn pipe(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    length: f64,
    diameter: f64,
    roughness: f64,
) -> PyResult<crate::results::PyPipeResult> {
    azoth_process::pipe(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        meters(length),
        meters(diameter),
        meters(roughness),
    )
    .map(|r| crate::results::PyPipeResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.cooler` - `Heater.run` reached through `Cooler`, as a registered id.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_temperature = None, duty = None, pressure_drop = None))]
#[pyo3(
    text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_temperature=None, duty=None, pressure_drop=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn cooler(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    outlet_temperature: Option<f64>,
    duty: Option<f64>,
    pressure_drop: Option<f64>,
) -> PyResult<crate::results::PyCoolerResult> {
    azoth_process::cooler(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        outlet_temperature.map(kelvins),
        duty.map(watts),
        pressure_drop.map(pascals),
    )
    .map(|r| crate::results::PyCoolerResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.heater` - the heater's kernel as a registered id.
///
/// **The three optional parameters are the class's three setters, and two of them are
/// exclusive.** A stated outlet temperature and a stated duty clear each other's flags in
/// `Heater`, so the pair is refused here rather than resolved; `pressure_drop` is
/// independent of both.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_temperature = None, duty = None, pressure_drop = None))]
#[pyo3(
    text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t, outlet_temperature=None, duty=None, pressure_drop=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn heater(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
    outlet_temperature: Option<f64>,
    duty: Option<f64>,
    pressure_drop: Option<f64>,
) -> PyResult<crate::results::PyHeaterResult> {
    azoth_process::heater(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
        outlet_temperature.map(kelvins),
        duty.map(watts),
        pressure_drop.map(pascals),
    )
    .map(|r| crate::results::PyHeaterResult::from(&r))
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

/// `process.heat_exchanger` - the exchanger's kernel as a registered id.
///
/// **Two component lists, because this is the only unit operation whose ports do not
/// share one fluid.** The names cross unresolved, as every other process model's do, and
/// the Rust side resolves each through the databank `Stream::mixture()` uses.
#[pyfunction]
// **The two component lists lead, and the rest follow in the spec's order**, which is
// what `tools/gen_stub.py` emits: it writes every fluid the model names before the
// declared inputs, so a signature that interleaved them would be a stub that lies.
#[pyo3(signature = (hot_components, cold_components, hot_in_n, hot_in_z, hot_in_p, hot_in_t, cold_in_n, cold_in_z, cold_in_p, cold_in_t, flow_arrangement, ua = None, hot_outlet_temperature = None, cold_outlet_temperature = None))]
#[pyo3(
    text_signature = "(hot_components, cold_components, hot_in_n, hot_in_z, hot_in_p, hot_in_t, cold_in_n, cold_in_z, cold_in_p, cold_in_t, flow_arrangement, ua=None, hot_outlet_temperature=None, cold_outlet_temperature=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are fourteen
pub fn heat_exchanger(
    py: Python<'_>,
    hot_components: Vec<String>,
    cold_components: Vec<String>,
    hot_in_n: f64,
    hot_in_z: Vec<f64>,
    hot_in_p: f64,
    hot_in_t: f64,
    cold_in_n: f64,
    cold_in_z: Vec<f64>,
    cold_in_p: f64,
    cold_in_t: f64,
    flow_arrangement: &str,
    ua: Option<f64>,
    hot_outlet_temperature: Option<f64>,
    cold_outlet_temperature: Option<f64>,
) -> PyResult<crate::results::PyHeatExchangerResult> {
    azoth_process::heat_exchanger(
        &hot_components,
        hot_in_n,
        &hot_in_z,
        pascals(hot_in_p),
        kelvins(hot_in_t),
        &cold_components,
        cold_in_n,
        &cold_in_z,
        pascals(cold_in_p),
        kelvins(cold_in_t),
        ua.map(watts_per_kelvin),
        flow_arrangement,
        hot_outlet_temperature.map(kelvins),
        cold_outlet_temperature.map(kelvins),
    )
    .map(|r| crate::results::PyHeatExchangerResult::from(&r))
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

/// `process.gas_scrubber` - the separator's kernel under the other entry.
///
/// `GasScrubber` does not override `run`, so this is `process.separator`'s arithmetic and
/// the same arguments.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, pressure_drop, gas_in_liquid, heat_input = None))]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, pressure_drop, gas_in_liquid, heat_input=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eight
pub fn gas_scrubber(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    pressure_drop: f64,
    gas_in_liquid: f64,
    heat_input: Option<f64>,
) -> PyResult<crate::results::PyGasScrubberResult> {
    azoth_process::gas_scrubber(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        pascals(pressure_drop),
        gas_in_liquid,
        heat_input.map(watts),
    )
    .map(|r| crate::results::PyGasScrubberResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.shortcut_distillation_column` - the FUG column as a registered id.
///
/// **The first `procedure` in this namespace.** Its `[algorithm]` block is the Underwood
/// bisection the two implementations both run; nothing in Rust or Python reads it at
/// runtime, and the case is what holds both to it.
#[pyfunction]
#[pyo3(
    signature = (components, feed_n, feed_z, feed_p, feed_t, light_key, heavy_key, light_key_recovery_distillate, heavy_key_recovery_bottoms, reflux_ratio_multiplier, condenser_pressure = None, reboiler_pressure = None)
)]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, light_key, heavy_key, light_key_recovery_distillate, heavy_key_recovery_bottoms, reflux_ratio_multiplier, condenser_pressure=None, reboiler_pressure=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are eleven
pub fn shortcut_distillation_column(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    light_key: &str,
    heavy_key: &str,
    light_key_recovery_distillate: f64,
    heavy_key_recovery_bottoms: f64,
    reflux_ratio_multiplier: f64,
    condenser_pressure: Option<f64>,
    reboiler_pressure: Option<f64>,
) -> PyResult<crate::results::PyShortcutDistillationColumnResult> {
    azoth_process::shortcut_distillation_column(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        light_key,
        heavy_key,
        light_key_recovery_distillate,
        heavy_key_recovery_bottoms,
        reflux_ratio_multiplier,
        condenser_pressure.map(pascals),
        reboiler_pressure.map(pascals),
    )
    .map(|r| crate::results::PyShortcutDistillationColumnResult::from(&r))
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

/// `process.manifold` - the manifold's kernel as a registered id.
///
/// **`many` at both ends**: the feeds cross as `process.mixer`'s vectors do and the outlets
/// as `process.splitter`'s do, because the manifold is those two composed.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, split_factors))]
#[pyo3(text_signature = "(components, feed_n, feed_z, feed_p, feed_t, split_factors)")]
#[allow(non_snake_case)] // the record's own field names
pub fn manifold(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: Vec<f64>,
    feed_z: Vec<Vec<f64>>,
    feed_p: Vec<f64>,
    feed_t: Vec<f64>,
    split_factors: Vec<f64>,
) -> PyResult<crate::results::PyManifoldResult> {
    let pressures: Vec<azoth_core::units::Pressure> = feed_p.into_iter().map(pascals).collect();
    let temperatures: Vec<azoth_core::units::ThermodynamicTemperature> =
        feed_t.into_iter().map(kelvins).collect();
    azoth_process::manifold(
        &components,
        &feed_n,
        &feed_z,
        &pressures,
        &temperatures,
        &split_factors,
    )
    .map(|r| crate::results::PyManifoldResult::from(&r))
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
