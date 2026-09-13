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
//!   See `hydraulics.rs`.
//! * **Exception classes are imported from `azoth.core.errors`, not defined
//!   here.** Both backends then raise the *same* class objects, so
//!   ``except OutOfRangeError`` works regardless of which one answered. See
//!   `errors.rs`.
//!
//! # Introspection
//!
//! `warning_codes`, `result_fields`, `calc_ids` and `version` exist so the test
//! suite can assert cross-language agreement without parsing Rust source. They
//! are the mechanism behind the claims that the two implementations share a
//! warning vocabulary and a result shape.

use pyo3::prelude::*;
use pyo3::types::PyModule;

mod errors;
mod hydraulics;
mod results;

use results::{
    PyColebrookResult, PyDarcyWeisbachResult, PyKComponent, PyKFactorsResult, PyQty,
    PyReynoldsNumberResult, PySwameeJainResult, PyWarning,
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
    m.add_class::<PyKFactorsResult>()?;
    m.add_class::<PyDarcyWeisbachResult>()?;

    // Calculations.
    m.add_function(wrap_pyfunction!(hydraulics::reynolds_number, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::friction_factor_colebrook, m)?)?;
    m.add_function(wrap_pyfunction!(
        hydraulics::friction_factor_swamee_jain,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(hydraulics::crane_k_factors, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::darcy_weisbach, m)?)?;

    // Introspection.
    m.add_function(wrap_pyfunction!(results::warning_codes, m)?)?;
    m.add_function(wrap_pyfunction!(results::result_fields, m)?)?;
    m.add_function(wrap_pyfunction!(results::calc_ids, m)?)?;
    m.add_function(wrap_pyfunction!(results::version, m)?)?;

    // Re-export the Python exception classes so both backends raise the same
    // objects rather than two lookalike hierarchies.
    errors::register(py, m)?;

    Ok(())
}
