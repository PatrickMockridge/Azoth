//! Warnings: conditions that do not stop a calculation but do change how much
//! its answer should be trusted.
//!
//! Warnings are deliberately *not* errors. A pressure drop computed at
//! transitional Reynolds number is real arithmetic applied to a friction factor
//! that is genuinely indeterminate; refusing to return it would be less useful
//! than returning it with a warning, because the caller often needs the number to
//! decide what to do next. The one thing the library must never do is return that
//! number silently - a result that cannot be told apart from a validated one is
//! the failure this type exists to prevent.

/// Machine-readable warning identifiers.
///
/// This list is part of the crate's public contract: the Python implementation
/// mirrors it exactly, and a test asserts the two sets are equal, so a warning
/// added on one side and not the other fails CI rather than quietly diverging.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum WarningCode {
    /// The value is outside the range in which the calculation is validated.
    OutOfValidRange,
    /// A range check could not be evaluated because an optional input was not
    /// supplied. Distinct from `OutOfValidRange`: the value was never checked,
    /// which is not the same as having passed.
    RangeCheckSkipped,
    /// Flow is in the transitional band, where the friction factor is
    /// indeterminate rather than merely uncertain.
    TransitionalFlow,
    /// An iterative solver hit its iteration cap without meeting tolerance. The
    /// returned value is the last iterate, not a converged result.
    SolverNotConverged,
    /// The calculation's source has not been verified against a primary
    /// reference, so the equation's provenance is unconfirmed.
    UnverifiedSource,
    /// The calculation depends on data that is a placeholder, not engineering
    /// data. Nothing computed from it should be used for design.
    EstimatedData,
}

impl WarningCode {
    /// Stable string form, used in generated docs, error messages and the
    /// cross-language parity test. Never change an existing value; add a new
    /// variant instead.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::OutOfValidRange => "OUT_OF_VALID_RANGE",
            Self::RangeCheckSkipped => "RANGE_CHECK_SKIPPED",
            Self::TransitionalFlow => "TRANSITIONAL_FLOW",
            Self::SolverNotConverged => "SOLVER_NOT_CONVERGED",
            Self::UnverifiedSource => "UNVERIFIED_SOURCE",
            Self::EstimatedData => "ESTIMATED_DATA",
        }
    }

    /// Every variant, for the parity test against the Python side.
    #[must_use]
    pub const fn all() -> &'static [Self] {
        &[
            Self::OutOfValidRange,
            Self::RangeCheckSkipped,
            Self::TransitionalFlow,
            Self::SolverNotConverged,
            Self::UnverifiedSource,
            Self::EstimatedData,
        ]
    }
}

impl std::fmt::Display for WarningCode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// A single warning attached to a result.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Warning {
    /// What kind of problem this is.
    pub code: WarningCode,
    /// Human-readable explanation. Should say what the consequence is, not just
    /// restate the code.
    pub message: String,
    /// The input or output the warning is about, when it is about one.
    pub field: Option<String>,
}

impl Warning {
    /// A warning not attached to a specific field.
    #[must_use]
    pub fn new(code: WarningCode, message: impl Into<String>) -> Self {
        Self {
            code,
            message: message.into(),
            field: None,
        }
    }

    /// A warning about a specific input or output.
    #[must_use]
    pub fn for_field(
        code: WarningCode,
        field: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            code,
            message: message.into(),
            field: Some(field.into()),
        }
    }
}

impl std::fmt::Display for Warning {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match &self.field {
            Some(field) => write!(f, "[{}] {}: {}", self.code, field, self.message),
            None => write!(f, "[{}] {}", self.code, self.message),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn codes_are_unique() {
        let mut seen = std::collections::HashSet::new();
        for code in WarningCode::all() {
            assert!(seen.insert(code.as_str()), "duplicate warning code {code}");
        }
    }

    #[test]
    fn as_str_is_screaming_snake_case() {
        // The Python side uses the same strings as a str-valued enum, so the
        // format is a cross-language contract, not a style preference.
        for code in WarningCode::all() {
            let s = code.as_str();
            assert!(
                s.chars().all(|c| c.is_ascii_uppercase() || c == '_'),
                "{s} is not SCREAMING_SNAKE_CASE"
            );
        }
    }

    #[test]
    fn warnings_are_comparable_for_the_parity_test() {
        let a = Warning::new(WarningCode::TransitionalFlow, "in the 2000-4000 band");
        let b = Warning::for_field(WarningCode::OutOfValidRange, "re", "below 4000");
        assert_ne!(a, b);
        assert_eq!(a, a.clone());
        assert_eq!(b.field.as_deref(), Some("re"));
        assert_eq!(a.field, None);
    }
}
