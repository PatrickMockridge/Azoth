//! Errors: conditions that stop a calculation.
//!
//! The split between an error and a warning is deliberate and load-bearing, and
//! this crate applies it consistently:
//!
//! * **Error** - the calculation is mathematically undefined, or the input is
//!   meaningless. `Re = 0` in the Colebrook equation divides by zero; a negative
//!   diameter is not a state anything can be computed from.
//! * **Warning** (`crate::warning`) - the calculation is well defined but lies
//!   outside the range where it has been validated. The number is real; its
//!   trustworthiness is the issue.
//!
//! Getting that boundary wrong in either direction is a real failure. Treating a
//! range violation as an error makes the library unusable for the exploratory
//! work engineers actually do. Treating a division by zero as a warning returns
//! `inf` dressed up as a result.
//!
//! Library code does not panic. The only `panic!` permitted is for genuinely
//! unreachable invariants, and there are none in this crate today.

use crate::warning::Warning;

/// Anything that can stop a calculation.
///
/// Every variant names the offending field, because "calculation failed" is not
/// an actionable message and an engineer chasing a bad number needs to know
/// which input to look at.
///
/// `PartialEq` but not `Eq`: several variants carry an `f64`, and `Eq` would be
/// a lie about NaN.
#[derive(Debug, Clone, PartialEq, thiserror::Error)]
#[non_exhaustive]
pub enum ChemEngError {
    /// An input is the wrong shape or value for any computation to proceed.
    #[error("invalid input `{field}`: {reason}")]
    InvalidInput {
        /// Name of the offending input.
        field: String,
        /// Why it cannot be used.
        reason: String,
    },

    /// An input violates a range check marked `severity: error` in the spec,
    /// meaning the calculation is undefined rather than merely out of range.
    #[error("`{field}` = {value} violates a hard limit: {detail}")]
    OutOfRange {
        /// Name of the offending input.
        field: String,
        /// The value supplied.
        value: f64,
        /// The limit it violated and why that limit is hard.
        detail: String,
    },

    /// An iterative solver hit its cap without meeting tolerance.
    #[error(
        "solver did not converge after {iterations} iterations (residual {residual:.3e}, \
         tolerance {tolerance:.3e})"
    )]
    SolverNotConverged {
        /// Iterations performed.
        iterations: u32,
        /// Final residual.
        residual: f64,
        /// Tolerance that was not met.
        tolerance: f64,
    },

    /// A fitting id could not be resolved against the fittings registry.
    ///
    /// An error rather than a skip: silently treating an unknown fitting as
    /// zero loss would under-report pressure drop, which is the dangerous
    /// direction to be wrong in.
    #[error("unknown fitting `{id}`; known ids are in the fittings registry")]
    UnknownFitting {
        /// The id that was not found.
        id: String,
    },

    /// The calculation has no verified source, so it must not be presented as a
    /// validated result.
    #[error("calculation `{id}` has no verified source")]
    UnverifiedCalculation {
        /// The calc id.
        id: String,
    },
}

impl ChemEngError {
    /// Construct an [`ChemEngError::InvalidInput`].
    pub fn invalid_input(field: impl Into<String>, reason: impl Into<String>) -> Self {
        Self::InvalidInput {
            field: field.into(),
            reason: reason.into(),
        }
    }

    /// Construct an [`ChemEngError::OutOfRange`].
    pub fn out_of_range(field: impl Into<String>, value: f64, detail: impl Into<String>) -> Self {
        Self::OutOfRange {
            field: field.into(),
            value,
            detail: detail.into(),
        }
    }

    /// The field this error is about, when it is about one.
    #[must_use]
    pub fn field(&self) -> Option<&str> {
        match self {
            Self::InvalidInput { field, .. } | Self::OutOfRange { field, .. } => Some(field),
            Self::UnknownFitting { id } | Self::UnverifiedCalculation { id } => Some(id),
            Self::SolverNotConverged { .. } => None,
        }
    }
}

/// Result alias used throughout the workspace.
pub type Result<T> = std::result::Result<T, ChemEngError>;

/// A value plus the warnings accumulated while producing it.
///
/// Used internally by the calc implementations to thread warnings through the
/// steps of a multi-stage calculation - `pipe` composes Reynolds number,
/// friction factor and fitting losses, and each stage can contribute warnings
/// that must survive to the final result rather than being dropped at the
/// boundary between stages.
#[derive(Debug, Clone, PartialEq)]
pub struct Warned<T> {
    /// The computed value.
    pub value: T,
    /// Warnings accumulated so far.
    pub warnings: Vec<Warning>,
}

impl<T> Warned<T> {
    /// Wrap a value with no warnings.
    pub fn clean(value: T) -> Self {
        Self {
            value,
            warnings: Vec::new(),
        }
    }

    /// Attach a warning.
    #[must_use]
    pub fn with(mut self, warning: Warning) -> Self {
        self.warnings.push(warning);
        self
    }

    /// Consume, returning the value and its warnings separately.
    #[must_use]
    pub fn split(self) -> (T, Vec<Warning>) {
        (self.value, self.warnings)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::warning::WarningCode;

    #[test]
    fn errors_name_their_field() {
        let e = ChemEngError::invalid_input("mu", "must be positive, got 0");
        assert_eq!(e.field(), Some("mu"));
        assert!(e.to_string().contains("mu"));
        assert!(e.to_string().contains("must be positive"));

        let e = ChemEngError::out_of_range("D", 0.0, "L/D is singular at D = 0");
        assert_eq!(e.field(), Some("D"));
        assert!(e.to_string().contains("D"));

        let e = ChemEngError::UnknownFitting {
            id: "no_such_fitting".into(),
        };
        assert_eq!(e.field(), Some("no_such_fitting"));
        assert!(e.to_string().contains("no_such_fitting"));
    }

    #[test]
    fn solver_error_reports_its_numbers() {
        let e = ChemEngError::SolverNotConverged {
            iterations: 100,
            residual: 1e-7,
            tolerance: 1e-12,
        };
        let msg = e.to_string();
        assert!(msg.contains("100"), "{msg}");
        assert!(msg.contains("1.000e-7"), "{msg}");
        assert!(msg.contains("1.000e-12"), "{msg}");
        assert_eq!(e.field(), None, "a solver failure is not about one input");
    }

    #[test]
    fn warned_threads_warnings() {
        let w = Warned::clean(1.0).with(Warning::new(WarningCode::TransitionalFlow, "band"));
        assert_eq!(w.value, 1.0);
        assert_eq!(w.warnings.len(), 1);
        let (value, warnings) = w.split();
        assert_eq!(value, 1.0);
        assert_eq!(warnings[0].code, WarningCode::TransitionalFlow);
    }
}
