//! Result types for the equations-of-state calculations.
//!
//! The same contract as every other namespace's results: one struct per
//! calculation, field names identical to the Python result dataclass and listed in
//! [`CalcResult::FIELDS`], with a test asserting the three agree. See
//! `crates/azoth-core/src/result.rs` for why the duplication is deliberate.

use azoth_core::{CalcResult, Warning};

/// Result of `eos.pr_kappa`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrKappaResult {
    /// The Peng-Robinson alpha-function coefficient. Dimensionless, and a property
    /// of the substance alone - it carries no temperature or pressure dependence.
    pub kappa: f64,
    /// Caveats. A negative `kappa` is returned rather than refused, carrying
    /// `OutOfValidRange`: the arithmetic is well defined and inspecting the limit
    /// is a legitimate thing for a caller to be doing.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PrKappaResult {
    const CALC_ID: &'static str = "eos.pr_kappa";
    const FIELDS: &'static [&'static str] = &["kappa", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
