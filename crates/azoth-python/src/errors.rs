//! Mapping Rust errors onto the Python exception hierarchy.
//!
//! The extension does **not** define its own exception classes. It imports the
//! ones in `azoth.core.errors` and raises those, so a caller writing
//! ``except OutOfRangeError`` catches errors from either backend. Two classes
//! with the same name and different identities would be a trap: the code would
//! look right, pass on whichever backend the author tested, and fail on the
//! other.
//!
//! The lookup happens at raise time rather than being cached at module init.
//! Errors are rare, a dict lookup on a already-imported module is nearly free,
//! and caching `Py<PyType>` values would need a `GILOnceCell` and a decision
//! about what to do when the import fails - which is a lot of machinery for a
//! path that runs when something has already gone wrong.

use azoth_core::AzothError;
use pyo3::prelude::*;
use pyo3::types::{PyModule, PyTuple};

/// The Python module holding the exception hierarchy.
const ERRORS_MODULE: &str = "azoth.core.errors";

/// Convert a [`AzothError`] into the matching Python exception.
///
/// Every variant maps to a class carrying the same information, so a caller can
/// inspect `err.field()` and get the offending input without parsing a message.
/// The `_` arm is deliberate: `AzothError` is `#[non_exhaustive]`, and a
/// variant added on the Rust side without a Python counterpart should surface as
/// the base class rather than failing to convert at all.
pub fn to_pyerr(py: Python<'_>, error: AzothError) -> PyErr {
    let class_name = match &error {
        AzothError::InvalidInput { .. } => "InvalidInputError",
        AzothError::OutOfRange { .. } => "OutOfRangeError",
        AzothError::SolverNotConverged { .. } => "SolverNotConvergedError",
        AzothError::UnknownFitting { .. } => "UnknownFittingError",
        AzothError::UnverifiedCalculation { .. } => "UnverifiedCalculationError",
        _ => "AzothError",
    };
    let message = error.to_string();

    let built: PyResult<PyErr> = (|| {
        let module = PyModule::import(py, ERRORS_MODULE)?;
        let class = module.getattr(class_name)?;

        // Construct through the class where it takes structured arguments, so
        // the Python exception carries the same fields the Rust error does and
        // `err.field()` works identically on both sides.
        let instance = match &error {
            AzothError::InvalidInput { field, reason } => {
                class.call1((field.as_str(), reason.as_str()))?
            }
            AzothError::OutOfRange {
                field,
                value,
                detail,
            } => class.call1((field.as_str(), *value, detail.as_str()))?,
            AzothError::SolverNotConverged {
                iterations,
                residual,
                tolerance,
            } => class.call1((*iterations, *residual, *tolerance))?,
            AzothError::UnknownFitting { id } => class.call1((id.as_str(),))?,
            AzothError::UnverifiedCalculation { id } => class.call1((id.as_str(),))?,
            _ => {
                let args = PyTuple::new(py, [message.as_str()])?;
                class.call1(args)?
            }
        };
        Ok(PyErr::from_value(instance))
    })();

    // If the import itself fails - a broken installation, or the extension
    // imported before the package - fall back to a plain RuntimeError carrying
    // the original message rather than masking the real error with an
    // ImportError.
    built.unwrap_or_else(|fallback| {
        PyErr::new::<pyo3::exceptions::PyRuntimeError, _>(format!(
            "azoth: {message} (could not import {ERRORS_MODULE} to build the \
             matching exception: {fallback})"
        ))
    })
}

/// Re-export the Python exception classes from the extension module, so
/// ``azoth._core.OutOfRangeError is azoth.core.errors.OutOfRangeError``.
///
/// The identity test in the suite is what makes "both backends raise the same
/// class" checkable rather than a claim. It asserts the two name sets are equal
/// rather than checking a few, because a hand-written list here is exactly the
/// thing that drifts: `KeycardError` was absent from this one while
/// `azoth.core.errors` declared it.
pub fn register(py: Python<'_>, module: &Bound<'_, PyModule>) -> PyResult<()> {
    let errors = PyModule::import(py, ERRORS_MODULE)?;
    for name in [
        "AzothError",
        "InvalidInputError",
        "KeycardError",
        "OutOfRangeError",
        "PropertyUnavailableError",
        "SolverNotConvergedError",
        "UnitMismatchError",
        "UnknownFittingError",
        "UnverifiedCalculationError",
    ] {
        module.add(name, errors.getattr(name)?)?;
    }
    Ok(())
}
