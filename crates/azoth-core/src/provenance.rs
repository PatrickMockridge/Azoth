//! What a number rests on, carried by the number.
//!
//! A result's values say what was computed and its [`Warning`]s say which declared checks
//! could not be evaluated. Neither says *what code produced it*, and that is the gap this
//! closes: a reader holding a number should not have to find the spec and the two kernels
//! before it can say where the number came from.
//!
//! # Static half and dynamic half
//!
//! The static half is a fact about the *calculation* and never changes between calls: its
//! id, the spec and the two implementations with a SHA-256 of each, the source its equation
//! is attributed to, and how far the answer is externally checked. `tools/gen_registry.py`
//! embeds it in [`crate::provenance_gen`], so it travels with an installed crate rather than
//! only with a checkout.
//!
//! The dynamic half is a fact about *this call*: which declared checks were skipped, and
//! which warnings were raised. [`WarningCode::RangeCheckSkipped`] is the one that matters
//! most, because it is the difference between a bound that was tested and a bound that was
//! never reached - and that difference is invisible in the value.
//!
//! # Why the status is derived rather than declared
//!
//! A spec refuses a verification-status field, and this does not add one back. The status is
//! measured by `tools/provenance.py` from the evidence the tree already holds - the tests a
//! spec ships and the external cases under `validation/` - so it cannot disagree with what it
//! summarises. A declared status is a status that can be wrong, which is the failure this
//! module exists to make impossible.

use crate::warning::{Warning, WarningCode};

/// Whether a calculation's answer rests on an external source, and how good that source is.
///
/// The spellings are the ones a `validation/*.json` case already uses, so a status here and
/// a status there read the same way rather than being two dialects for one idea.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VerificationStatus {
    /// At least one external case exists and every one of them is `verified`.
    Verified,
    /// Exercised against expectations pinned in its own spec, with no independent oracle
    /// recorded for it. **The normal case rather than a defect**, and why the counts travel
    /// beside the status: the reader is told which it is rather than being left with a rank.
    PartiallyVerified,
    /// A source was sought and not found, from either tree. This outranks good evidence on
    /// purpose - an attribution that does not resolve is not made to resolve by the other
    /// evidence being sound.
    SourceNeeded,
    /// No test is run at all, so nothing can be claimed.
    Unverified,
}

impl VerificationStatus {
    /// Stable lowercase string form, as the generated table and the Python side spell it.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::PartiallyVerified => "partially_verified",
            Self::SourceNeeded => "source_needed",
            Self::Unverified => "unverified",
        }
    }

    /// Every variant, for the parity test against the generated table.
    #[must_use]
    pub const fn all() -> &'static [Self] {
        &[
            Self::Verified,
            Self::PartiallyVerified,
            Self::SourceNeeded,
            Self::Unverified,
        ]
    }

    /// Parse a spelling back to a variant.
    ///
    /// `None` rather than a default: guessing a provenance is the one thing this type exists
    /// to prevent, so an unrecognised spelling must not resolve to something plausible.
    #[must_use]
    pub fn parse(value: &str) -> Option<Self> {
        Self::all().iter().copied().find(|s| s.as_str() == value)
    }
}

impl std::fmt::Display for VerificationStatus {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// The static half, as the generated table carries it.
///
/// Every field is a `&'static str` or a `u32`, so this is `Copy` and holds no float. That is
/// load-bearing rather than incidental: the block travels inside a JSON document, and JSON
/// has no representation for a non-finite number, so a block with no floats in it cannot be
/// refused by that rule.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ProvenanceStatic {
    /// The spec id, e.g. `hydraulics.darcy_weisbach`. Equals [`crate::result::CalcResult::CALC_ID`].
    pub calc_id: &'static str,
    /// The spec's own name for the calculation.
    pub name: &'static str,
    /// The spec file the bounds and the equation come from, and its hash.
    pub spec_path: &'static str,
    pub spec_sha256: &'static str,
    /// Both implementations, each with its hash. Neither is a wrapper for the other.
    pub python_path: &'static str,
    pub python_sha256: &'static str,
    pub rust_path: &'static str,
    pub rust_sha256: &'static str,
    /// What the spec attributes its equation to, e.g. `Crane TP-410`.
    pub source: &'static str,
    /// How far the answer is externally checked.
    pub verification: VerificationStatus,
    /// How many external cases exist for this id, and the spec's own test tally.
    pub validation_cases: u32,
    pub tests_active: u32,
    pub tests_skipped: u32,
}

/// The static half with the dynamic half filled in.
///
/// Borrows the skipped fields from the result's own warnings rather than owning them, so a
/// result that already holds its warnings pays nothing to describe them.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Provenance<'a> {
    /// The generated half.
    pub static_half: ProvenanceStatic,
    /// Declared checks this call could not evaluate, by the field they are about, sorted.
    ///
    /// A warning that names no field is absent rather than present as an empty name: this is
    /// a list of fields, and "some unnamed check did not run" is not one.
    pub skipped_checks: Vec<&'a str>,
    /// Every warning raised, by code, in [`WarningCode::all`]'s order.
    ///
    /// Ordered by the vocabulary rather than by the call, so two runs raising the same
    /// warnings in a different order describe themselves identically.
    pub warning_codes: Vec<&'static str>,
    /// True when the call carried no warnings at all.
    pub clean: bool,
}

impl<'a> Provenance<'a> {
    /// Merge the generated static half with what this call actually reported.
    ///
    /// The dynamic half is read from `warnings` rather than passed separately, because a
    /// caller able to report "which checks were skipped" apart from the warnings that said so
    /// would be a second answer to a question the result has already answered.
    #[must_use]
    pub fn of(static_half: ProvenanceStatic, warnings: &'a [Warning]) -> Self {
        let mut skipped: Vec<&'a str> = warnings
            .iter()
            .filter(|warning| warning.code == WarningCode::RangeCheckSkipped)
            .filter_map(|warning| warning.field.as_deref())
            .collect();
        skipped.sort_unstable();
        skipped.dedup();

        let raised: Vec<WarningCode> = warnings.iter().map(|warning| warning.code).collect();
        let warning_codes = WarningCode::all()
            .iter()
            .copied()
            .filter(|code| raised.contains(code))
            .map(WarningCode::as_str)
            .collect();

        Self {
            static_half,
            skipped_checks: skipped,
            warning_codes,
            clean: warnings.is_empty(),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    const STATIC: ProvenanceStatic = ProvenanceStatic {
        calc_id: "hydraulics.darcy_weisbach",
        name: "Darcy-Weisbach pressure drop",
        spec_path: "specs/calcs/hydraulics/darcy_weisbach.toml",
        spec_sha256: "aa",
        python_path: "python/src/azoth/hydraulics/reference/darcy_weisbach.py",
        python_sha256: "bb",
        rust_path: "crates/azoth-hydraulics/src/darcy_weisbach.rs",
        rust_sha256: "cc",
        source: "Crane TP-410",
        verification: VerificationStatus::SourceNeeded,
        validation_cases: 1,
        tests_active: 2,
        tests_skipped: 2,
    };

    fn warning(code: WarningCode, field: Option<&str>) -> Warning {
        Warning {
            code,
            message: "a message".to_string(),
            field: field.map(str::to_string),
        }
    }

    #[test]
    fn a_clean_call_says_so_and_names_no_skipped_check() {
        let block = Provenance::of(STATIC, &[]);

        assert!(block.clean);
        assert!(block.skipped_checks.is_empty());
        assert!(block.warning_codes.is_empty());
    }

    #[test]
    fn only_a_skipped_check_naming_a_field_becomes_a_skipped_check() {
        // Three exclusions in one case: a warning whose code does not mean "a check did not
        // run", one that names no field, and a repeat of a field already named.
        let warnings = [
            warning(WarningCode::RangeCheckSkipped, Some("re")),
            warning(WarningCode::RangeCheckSkipped, Some("re")),
            warning(WarningCode::RangeCheckSkipped, None),
            warning(WarningCode::OutOfValidRange, Some("L")),
        ];
        let block = Provenance::of(STATIC, &warnings);

        assert_eq!(block.skipped_checks, vec!["re"]);
        assert!(!block.clean);
    }

    #[test]
    fn both_lists_are_ordered_by_something_other_than_the_call() {
        let warnings = [
            warning(WarningCode::OutOfValidRange, Some("v")),
            warning(WarningCode::RangeCheckSkipped, Some("v")),
            warning(WarningCode::RangeCheckSkipped, Some("D")),
        ];
        let block = Provenance::of(STATIC, &warnings);

        assert_eq!(block.skipped_checks, vec!["D", "v"]);
        assert_eq!(
            block.warning_codes,
            vec!["OUT_OF_VALID_RANGE", "RANGE_CHECK_SKIPPED"]
        );
    }

    #[test]
    fn the_status_spellings_round_trip_and_an_unknown_one_is_refused() {
        for status in VerificationStatus::all() {
            assert_eq!(VerificationStatus::parse(status.as_str()), Some(*status));
        }
        assert_eq!(VerificationStatus::parse("probably_fine"), None);
        assert_eq!(VerificationStatus::parse(""), None);
    }
}
