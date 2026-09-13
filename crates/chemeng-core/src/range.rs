//! Range checking, driven by the `valid_range` block of a calc spec.
//!
//! The spec declares bounds as typed numbers rather than as expressions like
//! `"> 0"`, which means both implementations can evaluate them with no parser and
//! no risk of the two languages disagreeing about what an expression means.
//!
//! The distinction this module exists to make is between a value that was
//! *checked and passed* and a value that was *never checked*. Those look
//! identical to a caller unless the second case says so, so
//! [`apply_checks`] emits a [`WarningCode::RangeCheckSkipped`] whenever a check
//! cannot be evaluated. That happens in practice: `hydraulics.darcy_weisbach`
//! takes viscosity as an optional input, and without it there is no Reynolds
//! number to range-check against.

use crate::error::{ChemEngError, Result};
use crate::warning::{Warning, WarningCode};

/// What happens when a range check is violated.
///
/// Mirrors `severity` in the spec schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// The value is outside the validated range but the calculation is still
    /// defined. The result is returned carrying a warning.
    Warning,
    /// The calculation is undefined or the input is meaningless. This raises.
    Error,
}

/// Which side of the interval counts as a violation.
///
/// Mirrors `when` in the spec schema.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Band {
    /// The interval is the *allowed* range; values beyond it violate. The
    /// default, and what almost every check means.
    Outside,
    /// The interval is the *forbidden* band; values inside it violate. Needed
    /// for transitional flow, which is a problem precisely between 2000 and 4000
    /// rather than outside it.
    Inside,
}

/// A single bound on a quantity.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RangeCheck {
    /// Name of the quantity being bounded, as it appears in the spec.
    pub quantity: &'static str,
    /// Lower bound of the interval, if any.
    pub min: Option<f64>,
    /// Whether the lower bound itself is inside the interval.
    pub min_inclusive: bool,
    /// Upper bound of the interval, if any.
    pub max: Option<f64>,
    /// Whether the upper bound itself is inside the interval.
    pub max_inclusive: bool,
    /// Which side of the interval violates.
    pub band: Band,
    /// What a violation means.
    pub severity: Severity,
    /// Which warning code a violation emits.
    ///
    /// Defaults to [`WarningCode::OutOfValidRange`], which is right for a plain
    /// bound. It exists so a bound with a more specific meaning can say so: the
    /// transitional-flow band emits [`WarningCode::TransitionalFlow`], letting a
    /// caller react to indeterminate friction without string-matching a message.
    pub code: WarningCode,
    /// Why the bound exists. Shown to the user, so it should explain the physics
    /// rather than restate the number.
    pub rationale: &'static str,
}

impl RangeCheck {
    /// A lower bound that must be exceeded. Violated when `value < min`.
    #[must_use]
    pub const fn above(
        quantity: &'static str,
        min: f64,
        inclusive: bool,
        severity: Severity,
        rationale: &'static str,
    ) -> Self {
        Self {
            quantity,
            min: Some(min),
            min_inclusive: inclusive,
            max: None,
            max_inclusive: true,
            band: Band::Outside,
            severity,
            code: WarningCode::OutOfValidRange,
            rationale,
        }
    }

    /// An upper bound. Violated when `value > max`.
    #[must_use]
    pub const fn below(
        quantity: &'static str,
        max: f64,
        inclusive: bool,
        severity: Severity,
        rationale: &'static str,
    ) -> Self {
        Self {
            quantity,
            min: None,
            min_inclusive: true,
            max: Some(max),
            max_inclusive: inclusive,
            band: Band::Outside,
            severity,
            code: WarningCode::OutOfValidRange,
            rationale,
        }
    }

    /// A closed interval outside which the value violates. Both bounds
    /// inclusive.
    #[must_use]
    pub const fn between(
        quantity: &'static str,
        min: f64,
        max: f64,
        severity: Severity,
        rationale: &'static str,
    ) -> Self {
        Self {
            quantity,
            min: Some(min),
            min_inclusive: true,
            max: Some(max),
            max_inclusive: true,
            band: Band::Outside,
            severity,
            code: WarningCode::OutOfValidRange,
            rationale,
        }
    }

    /// A band *inside* which the value violates. Both bounds inclusive.
    #[must_use]
    pub const fn inside_band(
        quantity: &'static str,
        min: f64,
        max: f64,
        severity: Severity,
        rationale: &'static str,
    ) -> Self {
        Self {
            band: Band::Inside,
            ..Self::between(quantity, min, max, severity, rationale)
        }
    }

    /// Override the warning code emitted on violation.
    #[must_use]
    pub const fn with_code(mut self, code: WarningCode) -> Self {
        self.code = code;
        self
    }

    /// Whether `value` violates this check. `NaN` violates everything.
    #[must_use]
    pub fn violated(&self, value: f64) -> bool {
        if value.is_nan() {
            return true;
        }
        let below = match self.min {
            Some(min) if self.min_inclusive => value < min,
            Some(min) => value <= min,
            None => false,
        };
        let above = match self.max {
            Some(max) if self.max_inclusive => value > max,
            Some(max) => value >= max,
            None => false,
        };
        match self.band {
            Band::Outside => below || above,
            Band::Inside => !(below || above),
        }
    }

    /// Human-readable form of the bound, for messages and generated docs. The
    /// docs generator renders the same text from the spec.
    #[must_use]
    pub fn describe(&self) -> String {
        let q = self.quantity;
        let lo = self.min.map(|m| {
            let op = if self.min_inclusive { ">=" } else { ">" };
            format!("{op} {m}")
        });
        let hi = self.max.map(|m| {
            let op = if self.max_inclusive { "<=" } else { "<" };
            format!("{op} {m}")
        });
        let bounds = match (lo, hi) {
            // Both parts name the quantity: "re >= 2000 and re <= 4000" reads
            // unambiguously, where ">= 2000 and re <= 4000" does not.
            (Some(l), Some(h)) => format!("{q} {l} and {q} {h}"),
            (Some(l), None) => format!("{q} {l}"),
            (None, Some(h)) => format!("{q} {h}"),
            (None, None) => format!("{q} (unbounded)"),
        };
        match self.band {
            Band::Outside => bounds,
            Band::Inside => format!("{bounds} (violation lies inside this band)"),
        }
    }

    /// Apply this check to a value.
    ///
    /// # Errors
    /// Returns [`ChemEngError::OutOfRange`] when violated with
    /// [`Severity::Error`].
    pub fn apply(&self, value: f64, warnings: &mut Vec<Warning>) -> Result<()> {
        if !self.violated(value) {
            return Ok(());
        }
        match self.severity {
            Severity::Error => Err(ChemEngError::out_of_range(
                self.quantity,
                value,
                format!("must satisfy {}. {}", self.describe(), self.rationale),
            )),
            Severity::Warning => {
                warnings.push(Warning::for_field(
                    self.code,
                    self.quantity,
                    format!(
                        "value {value} violates {}. {}",
                        self.describe(),
                        self.rationale
                    ),
                ));
                Ok(())
            }
        }
    }
}

/// Apply a set of checks, resolving each check's quantity through `resolve`.
///
/// `resolve` returns `None` for a quantity that cannot be computed from the
/// inputs actually supplied - typically because an optional input such as
/// viscosity was omitted. That is not a pass: it produces a
/// [`WarningCode::RangeCheckSkipped`] warning, so the caller can distinguish
/// "checked and fine" from "never checked".
///
/// Takes any iterator of checks, so a caller can apply the input-phase checks
/// and the derived-phase checks separately without collecting either into a
/// temporary `Vec`.
///
/// # Errors
/// Propagates the first [`Severity::Error`] violation.
pub fn apply_checks<'a>(
    checks: impl IntoIterator<Item = &'a RangeCheck>,
    resolve: impl Fn(&str) -> Option<f64>,
    warnings: &mut Vec<Warning>,
) -> Result<()> {
    for check in checks {
        match resolve(check.quantity) {
            Some(value) => check.apply(value, warnings)?,
            // Deliberately not `check.code`: a skipped check was never a
            // violation of anything, so it always reports as skipped regardless
            // of which code a violation would have carried.
            None => warnings.push(Warning::for_field(
                WarningCode::RangeCheckSkipped,
                check.quantity,
                format!(
                    "could not check `{}` (must satisfy {}) because an input it depends on \
                     was not supplied; this range is UNVERIFIED for this result",
                    check.quantity,
                    check.describe()
                ),
            )),
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn checks() -> Vec<RangeCheck> {
        vec![
            RangeCheck::above("re", 0.0, false, Severity::Error, "singular at zero"),
            RangeCheck::above("re", 4000.0, true, Severity::Warning, "turbulent only"),
            RangeCheck::inside_band("re", 2000.0, 4000.0, Severity::Warning, "transitional"),
        ]
    }

    #[test]
    fn error_severity_raises_and_warning_severity_does_not() {
        let mut w = Vec::new();
        // re = 0 violates the hard bound.
        let err = checks()[0].apply(0.0, &mut w);
        assert!(matches!(err, Err(ChemEngError::OutOfRange { .. })));
        assert!(w.is_empty(), "an error must not also push a warning");

        // re = 1000 violates only the soft bound (below 4000).
        let ok = checks()[1].apply(1000.0, &mut w);
        assert!(ok.is_ok());
        assert_eq!(w.len(), 1);
        assert_eq!(w[0].code, WarningCode::OutOfValidRange);
        assert_eq!(w[0].field.as_deref(), Some("re"));
    }

    #[test]
    fn inside_band_only_fires_within_the_band() {
        let band = checks()[2];
        assert!(!band.violated(1999.0), "laminar is outside the band");
        assert!(band.violated(2000.0), "lower edge is transitional");
        assert!(band.violated(3000.0));
        assert!(band.violated(4000.0), "upper edge is transitional");
        assert!(!band.violated(4001.0), "turbulent is outside the band");
    }

    #[test]
    fn inclusivity_is_respected() {
        // Exclusive lower bound: the bound itself violates.
        let exclusive = RangeCheck::above("x", 0.0, false, Severity::Error, "must be positive");
        assert!(exclusive.violated(0.0));
        assert!(!exclusive.violated(1e-300));

        // Inclusive upper bound: the bound itself is fine.
        let inclusive = RangeCheck::below("x", 10.0, true, Severity::Warning, "cap");
        assert!(!inclusive.violated(10.0));
        assert!(inclusive.violated(10.000_001));

        // Exclusive upper bound: the bound itself violates.
        let exclusive = RangeCheck::below("x", 10.0, false, Severity::Warning, "cap");
        assert!(exclusive.violated(10.0));
    }

    #[test]
    fn nan_violates_everything() {
        // NaN comparisons are all false, so without an explicit check a NaN
        // would sail through every bound. It must not.
        for check in checks() {
            assert!(
                check.violated(f64::NAN),
                "{} let NaN through",
                check.quantity
            );
        }
    }

    #[test]
    fn evaluate_order_resolves_then_checks() {
        let checks = checks();
        let mut warnings = Vec::new();
        apply_checks(&checks, |q| (q == "re").then_some(3000.0), &mut warnings).unwrap();
        // 3000 violates the "above 4000" warning check and is inside the band.
        let codes: Vec<_> = warnings.iter().map(|w| w.code).collect();
        assert!(codes.contains(&WarningCode::OutOfValidRange));
        assert_eq!(warnings.len(), 2);
    }

    #[test]
    fn unresolvable_quantity_warns_rather_than_passing() {
        // This is the optional-input case: darcy_weisbach without mu has no
        // Reynolds number. It must be distinguishable from a value that passed.
        let mut warnings = Vec::new();
        apply_checks(&checks(), |_| None, &mut warnings).unwrap();
        assert_eq!(warnings.len(), 3, "one skip per unresolvable check");
        for w in &warnings {
            assert_eq!(w.code, WarningCode::RangeCheckSkipped);
            assert!(w.message.contains("UNVERIFIED"), "{}", w.message);
        }
    }

    #[test]
    fn describe_reads_like_the_bound() {
        assert_eq!(checks()[0].describe(), "re > 0");
        assert_eq!(checks()[1].describe(), "re >= 4000");
        let band = RangeCheck::inside_band("re", 2000.0, 4000.0, Severity::Warning, "x");
        assert_eq!(
            band.describe(),
            "re >= 2000 and re <= 4000 (violation lies inside this band)"
        );
    }
}
