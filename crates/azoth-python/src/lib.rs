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
//! calcs are built from: the claim that both languages read the same tables is
//! asserted by comparing parsed rows and embedded bytes. See `data.rs`.

use pyo3::prelude::*;
use pyo3::types::PyModule;

mod batch;
mod data;
mod eos;
mod errors;
mod hydraulics;
mod overlay;
mod process;
mod provenance;
mod reactions;
mod register_gen;
mod registry_tables_gen;
mod results;
mod standards;
mod thermal;
mod transport_gen;
mod wrappers_gen;

use results::{PyKComponent, PyQty, PyWarning};

#[pymodule]
fn _core(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Compiled core for azoth. Use azoth.hydraulics instead.",
    )?;

    register_gen::register(m)?;

    // Transport and result types.
    m.add_class::<PyQty>()?;
    m.add_class::<PyWarning>()?;
    m.add_class::<PyKComponent>()?;
    // An associating mixture's parameters, which every model whose Python side takes a
    // `Mixture` requires. A record rather than a defaulted argument: a caller that forgets
    // must get a `TypeError`, not a plausible answer for a different fluid.
    m.add_class::<eos::PyAssociationSpec>()?;

    // Data transport, for the cross-language data comparison.
    m.add_class::<batch::PyBatchColumn>()?;
    m.add_class::<batch::PyBatchResult>()?;
    m.add_class::<data::PyDataFile>()?;
    m.add_class::<data::PyFittingRow>()?;
    m.add_class::<data::PyFluidRow>()?;
    m.add_class::<data::PyComponentRow>()?;
    m.add_class::<data::PyKijRow>()?;

    // Calculations.

    // Thermal calculations.

    // Equations of state.

    // Models: the same shape, a different spec tree and generator.

    // Introspection.
    m.add_function(wrap_pyfunction!(batch::batch_run, m)?)?;
    m.add_function(wrap_pyfunction!(data::data_files, m)?)?;
    m.add_function(wrap_pyfunction!(data::fittings_rows, m)?)?;
    m.add_function(wrap_pyfunction!(data::fluid_rows, m)?)?;
    m.add_function(wrap_pyfunction!(data::component_rows, m)?)?;
    m.add_function(wrap_pyfunction!(data::kij_rows, m)?)?;
    m.add_function(wrap_pyfunction!(results::warning_codes, m)?)?;
    m.add_class::<overlay::PyOverlay>()?;
    m.add_function(wrap_pyfunction!(overlay::overlay, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::card_overlay, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::card_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::card_model_components, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_entry_row, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_origins, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_component_rows, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_kij_rows, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_cpa_kij_rows, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_names, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_dimensions, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_si_factor, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_slots, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_sets, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_affine_si, m)?)?;
    m.add_function(wrap_pyfunction!(results::solver_kinds, m)?)?;
    m.add_function(wrap_pyfunction!(eos::model_ids, m)?)?;
    m.add_function(wrap_pyfunction!(eos::model_schemes, m)?)?;
    m.add_function(wrap_pyfunction!(eos::model_kind, m)?)?;
    m.add_function(wrap_pyfunction!(registry_tables_gen::result_fields, m)?)?;
    m.add_function(wrap_pyfunction!(registry_tables_gen::calc_ids, m)?)?;
    m.add_function(wrap_pyfunction!(provenance::provenance_block, m)?)?;
    m.add_function(wrap_pyfunction!(results::version, m)?)?;

    // The process layer: a stream value and the unit-operation kernels, plus the
    // flowsheet checker. A binding, not a second implementation.
    m.add_class::<process::PyStream>()?;
    m.add_function(wrap_pyfunction!(process::splitter_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::mixer_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::separator_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::throttling_valve_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::heat_exchanger_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::pump_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::validate_flowsheet, m)?)?;
    m.add_function(wrap_pyfunction!(process::run_flowsheet, m)?)?;
    m.add_function(wrap_pyfunction!(process::catalogue, m)?)?;
    m.add_class::<process::PySession>()?;

    // The standards namespace.

    // Re-export the Python exception classes so both backends raise the same
    // objects rather than two lookalike hierarchies.
    errors::register(py, m)?;

    Ok(())
}
