//! `azoth._core` - the compiled core, exposed to Python.
//!
//! # What this module is and is not
//!
//! It is a numerical core with a private name. The public Python API is
//! `azoth.hydraulics`, which decides whether to call this module or the pure
//! Python reference implementation and presents identical results either way.
//!
//! Two consequences worth stating, because they look like omissions:
//!
//! * **Arguments are SI magnitudes, not quantities.** Unit handling happens once,
//!   in Python, so the two backends cannot disagree about what a number is in.
//!   See `hydraulics.rs` and `thermal.rs`.
//! * **Exception classes are imported from `azoth.core.errors`, not defined
//!   here.** Both backends then raise the *same* class objects, so
//!   ``except OutOfRangeError`` works regardless of which one answered. See
//!   `errors.rs`.
//!
//! # Introspection
//!
//! `warning_codes`, `unit_names`, `solver_kinds`, `result_fields`, `calc_ids` and
//! `version` exist so the test suite can assert cross-language agreement without
//! parsing Rust source. They are the mechanism behind the claims that the two
//! implementations share a warning vocabulary, a unit vocabulary, a solver
//! vocabulary and a result shape.
//!
//! `data_files`, `fittings_rows` and `fluid_rows` do the same job for the data the
//! calcs are built from. See `data.rs` - the claim that both languages read the same
//! tables was made in four docstrings before anything checked it.

use pyo3::prelude::*;
use pyo3::types::PyModule;

mod batch;
mod data;
mod eos;
mod errors;
mod hydraulics;
mod results;
mod thermal;

use results::{
    PyChokedFlowAreaResult, PyColebrookResult, PyConductionPlaneWallResult, PyControlValveCvResult,
    PyDarcyWeisbachResult, PyHaalandResult, PyKComponent, PyKFactorsResult, PyOrificeFlowResult,
    PyPrAlphaAbResult, PyPrKappaResult, PyPumpPowerResult, PyQty, PyReynoldsNumberResult,
    PySwameeJainResult, PyWarning,
};

#[pymodule]
fn _core(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Compiled core for azoth. Use azoth.hydraulics instead.",
    )?;

    // Transport and result types.
    m.add_class::<PyQty>()?;
    m.add_class::<PyWarning>()?;
    m.add_class::<PyKComponent>()?;
    m.add_class::<PyReynoldsNumberResult>()?;
    m.add_class::<PyColebrookResult>()?;
    m.add_class::<PySwameeJainResult>()?;
    m.add_class::<PyHaalandResult>()?;
    m.add_class::<PyConductionPlaneWallResult>()?;
    m.add_class::<PyPrKappaResult>()?;
    m.add_class::<PyPrAlphaAbResult>()?;
    m.add_class::<PyPumpPowerResult>()?;
    m.add_class::<PyOrificeFlowResult>()?;
    m.add_class::<PyControlValveCvResult>()?;
    m.add_class::<PyChokedFlowAreaResult>()?;

    // Data transport, for the cross-language data comparison.
    m.add_class::<batch::PyBatchColumn>()?;
    m.add_class::<batch::PyBatchResult>()?;
    m.add_class::<data::PyDataFile>()?;
    m.add_class::<data::PyFittingRow>()?;
    m.add_class::<data::PyFluidRow>()?;
    m.add_class::<PyKFactorsResult>()?;
    m.add_class::<PyDarcyWeisbachResult>()?;

    // Calculations.
    m.add_function(wrap_pyfunction!(hydraulics::reynolds_number, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::friction_factor_colebrook, m)?)?;
    m.add_function(wrap_pyfunction!(
        hydraulics::friction_factor_swamee_jain,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(hydraulics::friction_factor_haaland, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::crane_k_factors, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::darcy_weisbach, m)?)?;

    m.add_function(wrap_pyfunction!(hydraulics::pump_power, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::orifice_flow, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::control_valve_cv, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::choked_flow_area, m)?)?;

    // Thermal calculations.
    m.add_function(wrap_pyfunction!(thermal::conduction_plane_wall, m)?)?;

    // Equations of state.
    m.add_function(wrap_pyfunction!(eos::pr_kappa, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_alpha_ab, m)?)?;

    // Introspection.
    m.add_function(wrap_pyfunction!(batch::batch_run, m)?)?;
    m.add_function(wrap_pyfunction!(data::data_files, m)?)?;
    m.add_function(wrap_pyfunction!(data::fittings_rows, m)?)?;
    m.add_function(wrap_pyfunction!(data::fluid_rows, m)?)?;
    m.add_function(wrap_pyfunction!(results::warning_codes, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_names, m)?)?;
    m.add_function(wrap_pyfunction!(results::solver_kinds, m)?)?;
    m.add_function(wrap_pyfunction!(results::result_fields, m)?)?;
    m.add_function(wrap_pyfunction!(results::calc_ids, m)?)?;
    m.add_function(wrap_pyfunction!(results::version, m)?)?;

    // Re-export the Python exception classes so both backends raise the same
    // objects rather than two lookalike hierarchies.
    errors::register(py, m)?;

    Ok(())
}
