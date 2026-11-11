//! Run one registered calculation by id.
//!
//! **The generated half does the dispatch; this half owns the contract.** A quantity crosses as
//! an SI magnitude and an output quantity crosses back as `{magnitude_si, unit}`; a kernel's
//! `AzothError` is an outcome rather than an exception, because a refused calculation is a state
//! a sheet draws.

use serde_json::Value;

use crate::workspace::CalcOutcome;

/// What one calculation answered.
#[derive(Debug, Clone, PartialEq)]
pub enum Outcome {
    Ok(CalcOutcome),
    Error {
        message: String,
        field: Option<String>,
    },
}

/// Run `id` with the inputs a sheet holds.
///
/// An error for a missing input or a kernel's own refusal is `Outcome::Error`, not a panic: the
/// sheet shows the sentence beside the fields that caused it.
#[must_use]
pub fn run(id: &str, inputs: &Value) -> Outcome {
    crate::dispatch_gen::run(id, inputs)
}
