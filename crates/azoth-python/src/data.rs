//! The data tables, exposed to Python so the two sides can be compared.
//!
//! # Why this exists
//!
//! `crates/azoth-hydraulics/src/fluids.rs`, `fittings.rs` and
//! `python/src/azoth/_data.py` all used to claim a cross-language test that compared
//! the parsed data values. No such test existed, and it could not be written: the
//! extension exposed the calculations and not the tables they were built from, so a
//! Python test had no Rust-parsed value to compare against. Four docstrings asserting
//! a check that did not run is exactly the failure this project is organised against,
//! and it was found in its own source.
//!
//! So this module is the missing half. It exposes three things:
//!
//! * the parsed rows, so every field can be compared row by row;
//! * the **raw embedded text**, so the comparison can be made on bytes rather than on
//!   values. Parsed values agreeing is a weaker claim: two different files can parse
//!   to the same values, and a stale copy bundled into a wheel would look precisely
//!   like that - an agreement that proves nothing about which file was read;
//! * the repo-relative path of each file, so Python knows which file on disk is
//!   supposed to be the same one, rather than hardcoding the mapping in a test.
//!
//! # What this is not
//!
//! Not a data API. Nothing a calculation needs comes from here - the calcs read their
//! own tables inside the crate, and the Python reference reads the files. This is
//! introspection in the same sense `warning_codes` and `result_fields` are: it exists
//! so a test can assert a cross-language claim instead of asserting it in prose.

use azoth_hydraulics::{fittings, fluids};
use pyo3::prelude::*;

/// One file this build embeds, with its bytes.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "DataFile")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyDataFile {
    /// The name the file is addressed by, e.g. `water` or `fittings`.
    #[pyo3(get)]
    pub name: String,
    /// Repo-relative path, e.g. `data/fluids/water.csv`.
    #[pyo3(get)]
    pub path: String,
    /// The exact text this build embedded. UTF-8, because `include_str!` requires it.
    #[pyo3(get)]
    pub text: String,
}

#[pymethods]
impl PyDataFile {
    fn __repr__(&self) -> String {
        format!("DataFile({}, {} bytes)", self.path, self.text.len())
    }
}

/// One row of the fittings registry, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "FittingRow"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFittingRow {
    /// Stable identifier.
    #[pyo3(get)]
    pub id: String,
    /// Coarse category.
    #[pyo3(get)]
    pub family: String,
    /// Human-readable name.
    #[pyo3(get)]
    pub name: String,
    /// Equivalent length ratio, `L_eq / D`.
    #[pyo3(get)]
    pub n_ld: f64,
    /// What the coefficient is multiplied by, conventionally `f_t`.
    #[pyo3(get)]
    pub f_t_basis: String,
    /// Where the value came from, or a statement that it came from nowhere.
    #[pyo3(get)]
    pub citation: String,
    /// `verified`, `unverified` or `estimated_dummy`.
    #[pyo3(get)]
    pub verify_status: String,
    /// The document the value was read from, if any.
    #[pyo3(get)]
    pub source_ref: Option<String>,
    /// Where inside that document, if any.
    #[pyo3(get)]
    pub source_locator: Option<String>,
}

/// One row of a fluid property table, transported.
///
/// Carries the fluid's name on every row rather than once per table, so a comparison
/// in Python can flatten every row of every fluid into one list and still know which
/// fluid each came from. A table-per-call shape would need the test to iterate fluids
/// and tables separately, which is one more place to get the loops wrong.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "FluidRow")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFluidRow {
    /// The fluid this row belongs to.
    #[pyo3(get)]
    pub fluid: String,
    /// Temperature in degrees Celsius.
    #[pyo3(get)]
    pub temperature_c: f64,
    /// Density in kg/m^3.
    #[pyo3(get)]
    pub density_kg_m3: f64,
    /// Dynamic viscosity in Pa*s.
    #[pyo3(get)]
    pub dynamic_viscosity_pa_s: f64,
    /// Where the value came from.
    #[pyo3(get)]
    pub citation: String,
    /// `verified`, `unverified` or `estimated_dummy`.
    #[pyo3(get)]
    pub verify_status: String,
}

/// Every data file this build embeds, with its bytes.
///
/// # Errors
/// Returns an error if the embedded text is malformed, which is a build-time
/// invariant rather than a user condition.
#[pyfunction]
#[must_use]
pub fn data_files() -> Vec<PyDataFile> {
    let mut out = vec![PyDataFile {
        name: "fittings".to_string(),
        path: fittings::embedded_path().to_string(),
        text: fittings::embedded_csv().to_string(),
    }];
    for fluid in fluids::available_fluids() {
        if let (Some(path), Some(text)) =
            (fluids::embedded_path(fluid), fluids::embedded_csv(fluid))
        {
            out.push(PyDataFile {
                name: fluid.to_string(),
                path: path.to_string(),
                text: text.to_string(),
            });
        }
    }
    out
}

/// Every row of the fittings registry, as this crate parsed it.
///
/// # Errors
/// Returns an error if the embedded registry is malformed.
#[pyfunction]
pub fn fittings_rows(py: Python<'_>) -> PyResult<Vec<PyFittingRow>> {
    let rows = fittings::registry().map_err(|e| crate::errors::to_pyerr(py, e))?;
    Ok(rows
        .iter()
        .map(|row| PyFittingRow {
            id: row.id.clone(),
            family: row.family.clone(),
            name: row.name.clone(),
            n_ld: row.n_ld,
            f_t_basis: row.f_t_basis.clone(),
            citation: row.citation.clone(),
            verify_status: row.status.as_str().to_string(),
            source_ref: row.source_ref.clone(),
            source_locator: row.source_locator.clone(),
        })
        .collect())
}

/// Every row of a built-in fluid's table, as this crate parsed it.
///
/// # Errors
/// Returns an error for a fluid this build does not carry, rather than an empty list:
/// an unknown fluid and a fluid with no rows are different problems, and returning
/// `[]` for both would make a typo look like an empty table.
#[pyfunction]
pub fn fluid_rows(py: Python<'_>, name: &str) -> PyResult<Vec<PyFluidRow>> {
    let table = fluids::provider_for(name).map_err(|e| crate::errors::to_pyerr(py, e))?;
    Ok(table
        .points
        .iter()
        .map(|point| PyFluidRow {
            fluid: table.name.clone(),
            temperature_c: point.temperature_c,
            density_kg_m3: point.density_kg_m3,
            dynamic_viscosity_pa_s: point.dynamic_viscosity_pa_s,
            citation: point.citation.clone(),
            verify_status: point.status.as_str().to_string(),
        })
        .collect())
}
