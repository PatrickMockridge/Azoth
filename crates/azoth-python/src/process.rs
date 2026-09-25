//! The process layer, exposed to Python.
//!
//! Thin by design: this is a binding to `azoth-process`, not a second
//! implementation. The `Stream` value and the kernels run in Rust; Python supplies
//! SI magnitudes at the boundary and reads the streams back.

use std::path::Path;

use azoth_core::units::{
    joules_per_mole, kelvins, kilograms_per_cubic_meter, meters, pascals, square_meters_per_second,
    watts, watts_per_kelvin, watts_per_square_meter_kelvin,
};
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

/// `process.distillation_column` - the column solve as a registered id.
///
/// **The first process id whose outputs are vectors of a *computed* length** - one entry per
/// tray - and the first that refuses declared parameters by name. Both are the boundary's
/// business rather than the kernel's.
#[pyfunction]
#[pyo3(
    signature = (components, feed_n, feed_z, feed_p, feed_t, number_of_stages, feed_stage, has_reboiler, has_condenser, top_pressure, bottom_pressure, temperature_tolerance = 1e-6, max_iterations = 200, reboiler_temperature = None, condenser_temperature = None, murphree_efficiency = None, solver_type = None, top_specification_type = None, top_specification_target = None, top_specification_component = None, bottom_specification_type = None, bottom_specification_target = None, bottom_specification_component = None, reactive = None, reactive_start_tray = None, reactive_end_tray = None)
)]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, number_of_stages, feed_stage, has_reboiler, has_condenser, top_pressure, bottom_pressure, temperature_tolerance, max_iterations, reboiler_temperature=None, condenser_temperature=None, murphree_efficiency=None, solver_type=None, top_specification_type=None, top_specification_target=None, top_specification_component=None, bottom_specification_type=None, bottom_specification_target=None, bottom_specification_component=None, reactive=None, reactive_start_tray=None, reactive_end_tray=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-two
pub fn distillation_column(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    number_of_stages: usize,
    feed_stage: usize,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: f64,
    bottom_pressure: f64,
    temperature_tolerance: f64,
    max_iterations: usize,
    reboiler_temperature: Option<f64>,
    condenser_temperature: Option<f64>,
    murphree_efficiency: Option<f64>,
    solver_type: Option<&str>,
    top_specification_type: Option<&str>,
    top_specification_target: Option<f64>,
    top_specification_component: Option<&str>,
    bottom_specification_type: Option<&str>,
    bottom_specification_target: Option<f64>,
    bottom_specification_component: Option<&str>,
    reactive: Option<bool>,
    reactive_start_tray: Option<usize>,
    reactive_end_tray: Option<usize>,
) -> PyResult<crate::results::PyDistillationColumnResult> {
    azoth_process::distillation_column(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        number_of_stages,
        feed_stage,
        has_reboiler,
        has_condenser,
        pascals(top_pressure),
        pascals(bottom_pressure),
        reboiler_temperature.map(kelvins),
        condenser_temperature.map(kelvins),
        temperature_tolerance,
        max_iterations,
        murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        reactive,
        reactive_start_tray,
        reactive_end_tray,
    )
    .map(|r| crate::results::PyDistillationColumnResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.absorption_column` - the tray absorber as a registered id.
///
/// **The column's own profile shape with the class's own product names**, and two inlets where
/// the column has one: the gas at stage 0 through `feed` and the solvent at the top stage
/// through `top_feed`, which is `addGasInStream` and `addSolventInStream`.
#[pyfunction]
#[pyo3(
    signature = (gas_components, solvent_components, gas_n, gas_z, gas_p, gas_t, solvent_n, solvent_z, solvent_p, solvent_t, number_of_stages, top_pressure, bottom_pressure, temperature_tolerance, max_iterations, tray_temperatures = None, murphree_efficiency = None, component_murphree_efficiency = None, max_allowable_gas_load_factor = None, solver_type = None)
)]
#[pyo3(
    text_signature = "(gas_components, solvent_components, gas_n, gas_z, gas_p, gas_t, solvent_n, solvent_z, solvent_p, solvent_t, number_of_stages, top_pressure, bottom_pressure, temperature_tolerance, max_iterations, tray_temperatures=None, murphree_efficiency=None, component_murphree_efficiency=None, max_allowable_gas_load_factor=None, solver_type=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn absorption_column(
    py: Python<'_>,
    // **The two component lists come first**, which is the order the generated stub states:
    // `gen_stub` puts every `type = "components"` input ahead of the scalars, and a binding
    // whose positional order disagreed with its own stub would be a trap for a positional
    // caller - see `process.distillation_column`'s own note on the optional inputs.
    gas_components: Vec<String>,
    solvent_components: Vec<String>,
    gas_n: f64,
    gas_z: Vec<f64>,
    gas_p: f64,
    gas_t: f64,
    solvent_n: f64,
    solvent_z: Vec<f64>,
    solvent_p: f64,
    solvent_t: f64,
    number_of_stages: usize,
    top_pressure: f64,
    bottom_pressure: f64,
    temperature_tolerance: f64,
    max_iterations: usize,
    tray_temperatures: Option<Vec<f64>>,
    murphree_efficiency: Option<f64>,
    component_murphree_efficiency: Option<Vec<f64>>,
    max_allowable_gas_load_factor: Option<f64>,
    solver_type: Option<&str>,
) -> PyResult<crate::results::PyAbsorptionColumnResult> {
    azoth_process::absorption_column(
        &gas_components,
        gas_n,
        &gas_z,
        pascals(gas_p),
        kelvins(gas_t),
        &solvent_components,
        solvent_n,
        &solvent_z,
        pascals(solvent_p),
        kelvins(solvent_t),
        number_of_stages,
        pascals(top_pressure),
        pascals(bottom_pressure),
        tray_temperatures.as_deref(),
        temperature_tolerance,
        max_iterations,
        murphree_efficiency,
        component_murphree_efficiency.as_deref(),
        max_allowable_gas_load_factor,
        solver_type,
    )
    .map(|r| crate::results::PyAbsorptionColumnResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.packed_column` - the packed column as a registered id.
///
/// **The base column's own machine at a derived stage count**: `PackedColumn extends
/// DistillationColumn`, its `run` is `super.run(id)` and the packing is read by a hydraulics
/// report afterwards, so `packed_height` reaches the separation only through `estimateStages`.
#[pyfunction]
#[pyo3(
    signature = (components, feed_n, feed_z, feed_p, feed_t, packed_height, feed_stage, has_reboiler, has_condenser, top_pressure, bottom_pressure, temperature_tolerance, max_iterations, reboiler_temperature = None, condenser_temperature = None, packing_type = None, structured_packing = None, design_flood_fraction = None, packing_hydraulic_capacity_factor = None, column_diameter = None, murphree_efficiency = None, solver_type = None, top_specification_type = None, top_specification_target = None, top_specification_component = None, bottom_specification_type = None, bottom_specification_target = None, bottom_specification_component = None)
)]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, packed_height, feed_stage, has_reboiler, has_condenser, top_pressure, bottom_pressure, temperature_tolerance, max_iterations, reboiler_temperature=None, condenser_temperature=None, packing_type=None, structured_packing=None, design_flood_fraction=None, packing_hydraulic_capacity_factor=None, column_diameter=None, murphree_efficiency=None, solver_type=None, top_specification_type=None, top_specification_target=None, top_specification_component=None, bottom_specification_type=None, bottom_specification_target=None, bottom_specification_component=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn packed_column(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    packed_height: f64,
    feed_stage: usize,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: f64,
    bottom_pressure: f64,
    temperature_tolerance: f64,
    max_iterations: usize,
    reboiler_temperature: Option<f64>,
    condenser_temperature: Option<f64>,
    packing_type: Option<&str>,
    structured_packing: Option<bool>,
    design_flood_fraction: Option<f64>,
    packing_hydraulic_capacity_factor: Option<f64>,
    column_diameter: Option<f64>,
    murphree_efficiency: Option<f64>,
    solver_type: Option<&str>,
    top_specification_type: Option<&str>,
    top_specification_target: Option<f64>,
    top_specification_component: Option<&str>,
    bottom_specification_type: Option<&str>,
    bottom_specification_target: Option<f64>,
    bottom_specification_component: Option<&str>,
) -> PyResult<crate::results::PyPackedColumnResult> {
    azoth_process::packed_column(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        packed_height,
        feed_stage,
        has_reboiler,
        has_condenser,
        pascals(top_pressure),
        pascals(bottom_pressure),
        reboiler_temperature.map(kelvins),
        condenser_temperature.map(kelvins),
        temperature_tolerance,
        max_iterations,
        packing_type,
        structured_packing,
        design_flood_fraction,
        packing_hydraulic_capacity_factor,
        column_diameter,
        murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
    )
    .map(|r| crate::results::PyPackedColumnResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.stripping_column` - the tray stripper as a registered id.
///
/// **The absorber's own machine under the class's own names**: `StrippingColumn extends
/// AbsorptionColumn` and renames its two inlets and its two products, so the binding is the same
/// call with the other labels.
#[pyfunction]
#[pyo3(
    signature = (stripping_gas_components, rich_liquid_components, stripping_gas_n, stripping_gas_z, stripping_gas_p, stripping_gas_t, rich_liquid_n, rich_liquid_z, rich_liquid_p, rich_liquid_t, number_of_stages, top_pressure, bottom_pressure, temperature_tolerance, max_iterations, tray_temperatures = None, murphree_efficiency = None, component_murphree_efficiency = None, max_allowable_gas_load_factor = None, solver_type = None)
)]
#[pyo3(
    text_signature = "(stripping_gas_components, rich_liquid_components, stripping_gas_n, stripping_gas_z, stripping_gas_p, stripping_gas_t, rich_liquid_n, rich_liquid_z, rich_liquid_p, rich_liquid_t, number_of_stages, top_pressure, bottom_pressure, temperature_tolerance, max_iterations, tray_temperatures=None, murphree_efficiency=None, component_murphree_efficiency=None, max_allowable_gas_load_factor=None, solver_type=None)"
)]
#[allow(non_snake_case)] // the record's own field names
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn stripping_column(
    py: Python<'_>,
    stripping_gas_components: Vec<String>,
    rich_liquid_components: Vec<String>,
    stripping_gas_n: f64,
    stripping_gas_z: Vec<f64>,
    stripping_gas_p: f64,
    stripping_gas_t: f64,
    rich_liquid_n: f64,
    rich_liquid_z: Vec<f64>,
    rich_liquid_p: f64,
    rich_liquid_t: f64,
    number_of_stages: usize,
    top_pressure: f64,
    bottom_pressure: f64,
    temperature_tolerance: f64,
    max_iterations: usize,
    tray_temperatures: Option<Vec<f64>>,
    murphree_efficiency: Option<f64>,
    component_murphree_efficiency: Option<Vec<f64>>,
    max_allowable_gas_load_factor: Option<f64>,
    solver_type: Option<&str>,
) -> PyResult<crate::results::PyStrippingColumnResult> {
    azoth_process::stripping_column(
        &stripping_gas_components,
        &rich_liquid_components,
        stripping_gas_n,
        &stripping_gas_z,
        pascals(stripping_gas_p),
        kelvins(stripping_gas_t),
        rich_liquid_n,
        &rich_liquid_z,
        pascals(rich_liquid_p),
        kelvins(rich_liquid_t),
        number_of_stages,
        pascals(top_pressure),
        pascals(bottom_pressure),
        temperature_tolerance,
        max_iterations,
        tray_temperatures.as_deref(),
        murphree_efficiency,
        component_murphree_efficiency.as_deref(),
        max_allowable_gas_load_factor,
        solver_type,
    )
    .map(|r| crate::results::PyStrippingColumnResult::from(&r))
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

/// `process.tank` - the tank's kernel as a registered id.
///
/// **No parameters at all.** A tank's steady state is the flash it holds, so the only
/// inputs are the feeds: `many` at the inlet, as `process.mixer`'s are.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t))]
#[pyo3(text_signature = "(components, feed_n, feed_z, feed_p, feed_t)")]
#[allow(non_snake_case)] // the record's own field names
pub fn tank(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: Vec<f64>,
    feed_z: Vec<Vec<f64>>,
    feed_p: Vec<f64>,
    feed_t: Vec<f64>,
) -> PyResult<crate::results::PyTankResult> {
    let pressures: Vec<azoth_core::units::Pressure> = feed_p.into_iter().map(pascals).collect();
    let temperatures: Vec<azoth_core::units::ThermodynamicTemperature> =
        feed_t.into_iter().map(kelvins).collect();
    azoth_process::tank(&components, &feed_n, &feed_z, &pressures, &temperatures)
        .map(|r| crate::results::PyTankResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// `process.stirred_tank_reactor` - the reactor's kernel as a registered id.
///
/// **Eight parameters where the palette entry declared none.** The reaction, the limiting
/// reactant, the conversion, the isothermal flag, the held temperature and pressure and the
/// pressure drop are all `run`'s own, and a kernel that took none of them could only be a
/// pass-through.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, reaction, limiting_reactant, conversion, isothermal, reactor_temperature = None, reactor_pressure = None, pressure_drop = None))]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, reaction, limiting_reactant, conversion, isothermal, reactor_temperature=None, reactor_pressure=None, pressure_drop=None)"
)]
#[allow(clippy::too_many_arguments)] // one argument per declared input, and there are twelve
pub fn stirred_tank_reactor(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    reaction: &str,
    limiting_reactant: &str,
    conversion: f64,
    isothermal: bool,
    reactor_temperature: Option<f64>,
    reactor_pressure: Option<f64>,
    pressure_drop: Option<f64>,
) -> PyResult<crate::results::PyStirredTankReactorResult> {
    azoth_process::stirred_tank_reactor(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        reaction,
        limiting_reactant,
        conversion,
        isothermal,
        reactor_temperature.map(kelvins),
        reactor_pressure.map(pascals),
        pascals(pressure_drop.unwrap_or(0.0)),
    )
    .map(|r| crate::results::PyStirredTankReactorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.plug_flow_reactor` - the reactor's kernel as a registered id.
///
/// **Twenty-eight arguments where the palette entry declared two.** The entry declared `length`
/// and `diameter`, and `run` reads the geometry, the energy mode and its coolant, the march's
/// controls, the catalyst bed and the whole rate law; each is here.
///
/// The optional ones are last, which is `gen_stub`'s rule: a `.pyi` with a defaulted parameter
/// before a non-defaulted one does not parse.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, length, diameter, number_of_tubes, energy_mode, coolant_temperature, overall_heat_transfer_coefficient, number_of_steps, integration_method, property_update_frequency, thermodynamic_coupling, reaction, reaction_orders, rate_type, pre_exponential_factor, activation_energy, temperature_exponent, heat_of_reaction, catalyst_bulk_density = None, catalyst_activity_factor = None, catalyst_particle_diameter = None, catalyst_void_fraction = None, catalyst_molecular_diffusivity = None, catalyst_effectiveness_enabled = None, key_component = None))]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, length, diameter, number_of_tubes, energy_mode, coolant_temperature, overall_heat_transfer_coefficient, number_of_steps, integration_method, property_update_frequency, thermodynamic_coupling, reaction, reaction_orders, rate_type, pre_exponential_factor, activation_energy, temperature_exponent, heat_of_reaction, catalyst_bulk_density=None, catalyst_activity_factor=None, catalyst_particle_diameter=None, catalyst_void_fraction=None, catalyst_molecular_diffusivity=None, catalyst_effectiveness_enabled=None, key_component=None)"
)]
#[allow(clippy::too_many_arguments)] // one argument per declared input, and there are twenty-nine
pub fn plug_flow_reactor(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    length: f64,
    diameter: f64,
    number_of_tubes: f64,
    energy_mode: &str,
    coolant_temperature: f64,
    overall_heat_transfer_coefficient: f64,
    number_of_steps: f64,
    integration_method: &str,
    property_update_frequency: f64,
    thermodynamic_coupling: &str,
    reaction: &str,
    reaction_orders: Vec<f64>,
    rate_type: &str,
    pre_exponential_factor: f64,
    activation_energy: f64,
    temperature_exponent: f64,
    heat_of_reaction: f64,
    catalyst_bulk_density: Option<f64>,
    catalyst_activity_factor: Option<f64>,
    catalyst_particle_diameter: Option<f64>,
    catalyst_void_fraction: Option<f64>,
    catalyst_molecular_diffusivity: Option<f64>,
    catalyst_effectiveness_enabled: Option<bool>,
    key_component: Option<String>,
) -> PyResult<crate::results::PyPlugFlowReactorResult> {
    azoth_process::plug_flow_reactor(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        length,
        diameter,
        number_of_tubes,
        energy_mode,
        kelvins(coolant_temperature),
        watts_per_square_meter_kelvin(overall_heat_transfer_coefficient),
        number_of_steps,
        integration_method,
        property_update_frequency,
        thermodynamic_coupling,
        reaction,
        &reaction_orders,
        rate_type,
        pre_exponential_factor,
        activation_energy,
        temperature_exponent,
        heat_of_reaction,
        catalyst_bulk_density.map(kilograms_per_cubic_meter),
        catalyst_activity_factor,
        catalyst_particle_diameter,
        catalyst_void_fraction,
        catalyst_molecular_diffusivity.map(square_meters_per_second),
        catalyst_effectiveness_enabled,
        key_component,
    )
    .map(|r| crate::results::PyPlugFlowReactorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.gibbs_reactor` - the equilibrium composition as a registered id.
///
/// **The equilibrium temperature is the feed's.** `GibbsReactor` has no temperature setter: it
/// reads `system.getTemperature()` from the fluid it is handed, so `feed_t` is the only one, and
/// the palette entry that used to declare a separate `temperature` parameter was wrong.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, energy_mode, damping_composition, max_iterations, convergence_tolerance, min_iterations))]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, energy_mode, damping_composition, max_iterations, convergence_tolerance, min_iterations)"
)]
#[allow(clippy::too_many_arguments)] // one argument per declared input
pub fn gibbs_reactor(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    energy_mode: &str,
    damping_composition: f64,
    max_iterations: f64,
    convergence_tolerance: f64,
    min_iterations: f64,
) -> PyResult<crate::results::PyGibbsReactorResult> {
    azoth_process::gibbs_reactor(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        energy_mode,
        damping_composition,
        max_iterations,
        convergence_tolerance,
        min_iterations,
    )
    .map(|r| crate::results::PyGibbsReactorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.flare` - the flare's kernel as a registered id.
///
/// **The record through and two numbers beside it.** `Flare.run` clones the inlet into the
/// outlet, so the five fields cross unchanged; the duty and the CO2 emission are the class's
/// own, and they are outputs here because they are what the machine reports.
#[pyfunction]
#[pyo3(signature = (components, inlet_n, inlet_z, inlet_p, inlet_t))]
#[pyo3(text_signature = "(components, inlet_n, inlet_z, inlet_p, inlet_t)")]
pub fn flare(
    py: Python<'_>,
    components: Vec<String>,
    inlet_n: f64,
    inlet_z: Vec<f64>,
    inlet_p: f64,
    inlet_t: f64,
) -> PyResult<crate::results::PyFlareResult> {
    azoth_process::flare(
        &components,
        inlet_n,
        &inlet_z,
        pascals(inlet_p),
        kelvins(inlet_t),
    )
    .map(|r| crate::results::PyFlareResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.ejector` - the ejector's kernel as a registered id.
///
/// **Two inlets carrying two fluids and one outlet.** The four efficiencies are the class's
/// five parameters less the mixing pressure, which the palette does not declare and `run`
/// therefore estimates.
#[pyfunction]
#[pyo3(signature = (motive_components, suction_components, motive_n, motive_z, motive_p, motive_t, suction_n, suction_z, suction_p, suction_t, discharge_pressure, motive_nozzle_efficiency, suction_nozzle_efficiency, mixing_efficiency, diffuser_efficiency))]
#[pyo3(
    text_signature = "(motive_components, suction_components, motive_n, motive_z, motive_p, motive_t, suction_n, suction_z, suction_p, suction_t, discharge_pressure, motive_nozzle_efficiency, suction_nozzle_efficiency, mixing_efficiency, diffuser_efficiency)"
)]
#[allow(clippy::too_many_arguments)] // one argument per declared input, and there are fifteen
pub fn ejector(
    py: Python<'_>,
    motive_components: Vec<String>,
    suction_components: Vec<String>,
    motive_n: f64,
    motive_z: Vec<f64>,
    motive_p: f64,
    motive_t: f64,
    suction_n: f64,
    suction_z: Vec<f64>,
    suction_p: f64,
    suction_t: f64,
    discharge_pressure: f64,
    motive_nozzle_efficiency: f64,
    suction_nozzle_efficiency: f64,
    mixing_efficiency: f64,
    diffuser_efficiency: f64,
) -> PyResult<crate::results::PyEjectorResult> {
    azoth_process::ejector(
        &motive_components,
        motive_n,
        &motive_z,
        pascals(motive_p),
        kelvins(motive_t),
        &suction_components,
        suction_n,
        &suction_z,
        pascals(suction_p),
        kelvins(suction_t),
        pascals(discharge_pressure),
        motive_nozzle_efficiency,
        suction_nozzle_efficiency,
        mixing_efficiency,
        diffuser_efficiency,
    )
    .map(|r| crate::results::PyEjectorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// `process.three_phase_separator` - the separator's kernel as a registered id.
///
/// **Six entrainment fractions and three outlets.** The fractions cross in `run`'s own
/// order, which is the order they are applied in, and the outlets are the record's five
/// fields under `vapour`, `light_liquid` and `heavy_liquid`.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, pressure_drop, gas_in_aqueous, gas_in_oil, oil_in_aqueous, oil_in_gas, aqueous_in_gas, aqueous_in_oil, heat_input = None))]
#[pyo3(
    text_signature = "(components, feed_n, feed_z, feed_p, feed_t, pressure_drop, gas_in_aqueous, gas_in_oil, oil_in_aqueous, oil_in_gas, aqueous_in_gas, aqueous_in_oil, heat_input=None)"
)]
#[allow(clippy::too_many_arguments)] // one argument per declared input, and there are thirteen
pub fn three_phase_separator(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    pressure_drop: f64,
    gas_in_aqueous: f64,
    gas_in_oil: f64,
    oil_in_aqueous: f64,
    oil_in_gas: f64,
    aqueous_in_gas: f64,
    aqueous_in_oil: f64,
    heat_input: Option<f64>,
) -> PyResult<crate::results::PyThreePhaseSeparatorResult> {
    azoth_process::three_phase_separator(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        pascals(pressure_drop),
        gas_in_aqueous,
        gas_in_oil,
        oil_in_aqueous,
        oil_in_gas,
        aqueous_in_gas,
        aqueous_in_oil,
        heat_input.map(watts),
    )
    .map(|r| crate::results::PyThreePhaseSeparatorResult::from(&r))
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

/// `process.component_splitter` - the component splitter's kernel as a registered id.
///
/// The factors are one per **component**, not one per outlet.
#[pyfunction]
#[pyo3(signature = (components, feed_n, feed_z, feed_p, feed_t, split_factors))]
#[pyo3(text_signature = "(components, feed_n, feed_z, feed_p, feed_t, split_factors)")]
#[allow(non_snake_case)] // the record's own field names
pub fn component_splitter(
    py: Python<'_>,
    components: Vec<String>,
    feed_n: f64,
    feed_z: Vec<f64>,
    feed_p: f64,
    feed_t: f64,
    split_factors: Vec<f64>,
) -> PyResult<crate::results::PyComponentSplitterResult> {
    azoth_process::component_splitter(
        &components,
        feed_n,
        &feed_z,
        pascals(feed_p),
        kelvins(feed_t),
        &split_factors,
    )
    .map(|r| crate::results::PyComponentSplitterResult::from(&r))
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
