//! How much a data value can be trusted.
//!
//! Shared by the fittings registry and the fluid property tables, because the
//! distinction those tables need to make is the same one:
//!
//! * **Estimated dummy** - a placeholder invented so the software has something
//!   to run against. Not engineering data. The Crane coefficients are currently
//!   all in this state.
//! * **Unverified** - a real published value, read from a secondary reference
//!   and not checked against a primary formulation. The water and air tables are
//!   in this state.
//! * **Verified** - checked against the primary source by a named person.
//!
//! The middle case is the one that gets lost in practice. A value being "not a
//! placeholder" and a value being "trustworthy" are different claims, and a
//! library that collapses them is telling its users something it does not know.

use azoth_core::{AzothError, Result, WarningCode};

/// Provenance of a data value.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerifyStatus {
    /// Placeholder for software testing. Not from any source. Not engineering
    /// data.
    EstimatedDummy,
    /// Read from a secondary public reference, not checked against the primary
    /// source.
    Unverified,
    /// Read from the primary source by a named person.
    Verified,
}

impl VerifyStatus {
    /// Parse the spelling used in the data files.
    ///
    /// # Errors
    /// An unrecognised value is an error rather than a default, because
    /// guessing a provenance is the one thing this type exists to prevent.
    pub fn parse(raw: &str) -> Result<Self> {
        match raw.trim() {
            "estimated_dummy" => Ok(Self::EstimatedDummy),
            "unverified" => Ok(Self::Unverified),
            "verified" => Ok(Self::Verified),
            other => Err(AzothError::invalid_input(
                "verify_status",
                format!(
                    "unknown verify_status `{other}`; expected estimated_dummy, \
                     unverified or verified"
                ),
            )),
        }
    }

    /// Stable machine-readable form, matching the data files and the Python side.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::EstimatedDummy => "estimated_dummy",
            Self::Unverified => "unverified",
            Self::Verified => "verified",
        }
    }

    /// True when the value is a placeholder rather than a measurement.
    #[must_use]
    pub const fn is_placeholder(self) -> bool {
        matches!(self, Self::EstimatedDummy)
    }

    /// The warning a result should carry when it used a value in this state.
    ///
    /// `None` for [`Self::Verified`], which is the only state that is silent.
    ///
    /// The mapping lives here rather than in each caller because it is the one
    /// thing every data-backed calculation has to agree on. The two callers, the
    /// fitting registry and the fluid tables, had already diverged: a fitting row's
    /// status became a warning on every result, and a fluid row's became nothing at
    /// all. Which code a status raises is a property of the status, not of the file
    /// it happened to be read from.
    ///
    /// The *message* is the caller's, because a placeholder fitting and a
    /// placeholder fluid have different things to say about what to do next.
    #[must_use]
    pub const fn warning_code(self) -> Option<WarningCode> {
        match self {
            Self::EstimatedDummy => Some(WarningCode::EstimatedData),
            Self::Unverified => Some(WarningCode::UnverifiedSource),
            Self::Verified => None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn round_trips_through_its_string_form() {
        for status in [
            VerifyStatus::EstimatedDummy,
            VerifyStatus::Unverified,
            VerifyStatus::Verified,
        ] {
            assert_eq!(VerifyStatus::parse(status.as_str()).unwrap(), status);
        }
    }

    #[test]
    fn placeholder_is_only_the_dummy_state() {
        assert!(VerifyStatus::EstimatedDummy.is_placeholder());
        assert!(!VerifyStatus::Unverified.is_placeholder());
        assert!(!VerifyStatus::Verified.is_placeholder());
    }

    #[test]
    fn unknown_status_is_rejected_rather_than_defaulted() {
        // Defaulting to `verified` would silently promote unknown data to
        // trusted data, which is the failure this type is here to stop.
        assert!(VerifyStatus::parse("probably_fine").is_err());
        assert!(VerifyStatus::parse("").is_err());
    }
}
