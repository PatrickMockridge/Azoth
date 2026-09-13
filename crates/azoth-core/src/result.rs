//! Result shapes shared by every calculation.

use crate::warning::Warning;

/// Flow regime of a pipe flow, under the Crane/Moody boundaries.
///
/// A note on terminology, because this genuinely trips people up. Crane and
/// Moody use "transition zone" to mean something *different* from the
/// 2000-4000 regime: for them it is the region on the Moody chart between the
/// smooth-flow line and complete turbulence, where the friction factor still
/// depends on Reynolds number. This enum uses the more common modern sense
/// (laminar / transitional / turbulent by Reynolds number), and the transition
/// zone of the Moody chart is not represented here at all.
///
/// The boundaries are approximate engineering guidance, not exact physical
/// transitions. Real transition depends on inlet geometry, vibration and
/// roughness, which is why [`FlowRegime::Transitional`] is a warning condition
/// everywhere it is used rather than a clean category.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum FlowRegime {
    /// `Re < 2000`.
    Laminar,
    /// `2000 <= Re <= 4000`. The friction factor is indeterminate here.
    Transitional,
    /// `Re > 4000`.
    Turbulent,
}

impl FlowRegime {
    /// Upper bound of the laminar band (exclusive).
    pub const LAMINAR_MAX: f64 = 2000.0;
    /// Lower bound of the turbulent band (exclusive).
    pub const TURBULENT_MIN: f64 = 4000.0;

    /// Classify a Reynolds number.
    ///
    /// Both band edges fall inside `Transitional`, so a value of exactly 2000 or
    /// exactly 4000 is transitional rather than laminar or turbulent. That is a
    /// convention, chosen because the boundaries are not sharp in reality and the
    /// conservative reading is the useful one.
    #[must_use]
    pub fn from_reynolds_number(re: f64) -> Self {
        if re < Self::LAMINAR_MAX {
            Self::Laminar
        } else if re <= Self::TURBULENT_MIN {
            Self::Transitional
        } else {
            Self::Turbulent
        }
    }

    /// Whether the friction factor is indeterminate in this regime.
    #[must_use]
    pub const fn is_indeterminate(self) -> bool {
        matches!(self, Self::Transitional)
    }

    /// Stable lowercase string form. A cross-language contract: the Python side
    /// produces the same strings.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Laminar => "laminar",
            Self::Transitional => "transitional",
            Self::Turbulent => "turbulent",
        }
    }

    /// Every variant, for the parity test against the Python side.
    #[must_use]
    pub const fn all() -> &'static [Self] {
        &[Self::Laminar, Self::Transitional, Self::Turbulent]
    }
}

impl std::fmt::Display for FlowRegime {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.write_str(self.as_str())
    }
}

/// Common shape of every calculation result.
///
/// `FIELDS` exists so the Rust shape can be asserted against the Python
/// dataclass field-for-field. Without a machine-checkable contract the two
/// languages drift, and the drift is invisible until someone accesses an
/// attribute that only exists on one side.
pub trait CalcResult {
    /// The spec id, e.g. `hydraulics.darcy_weisbach`. Must match `id` in the
    /// spec file.
    const CALC_ID: &'static str;

    /// Public field names, in declaration order. Must equal the Python result
    /// dataclass's field names.
    const FIELDS: &'static [&'static str];

    /// Warnings accumulated while computing this result.
    fn warnings(&self) -> &[Warning];

    /// True when the result carries no warnings at all.
    ///
    /// Useful precisely because it is a blunt instrument: a caller who checks
    /// this is opting out of every caveat the library attaches.
    fn is_clean(&self) -> bool {
        self.warnings().is_empty()
    }

    /// True when any warning is the given code.
    fn has_warning(&self, code: crate::warning::WarningCode) -> bool {
        self.warnings().iter().any(|w| w.code == code)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn regime_boundaries_are_as_documented() {
        assert_eq!(FlowRegime::from_reynolds_number(0.0), FlowRegime::Laminar);
        assert_eq!(
            FlowRegime::from_reynolds_number(1999.999),
            FlowRegime::Laminar
        );
        // Both edges are transitional by convention.
        assert_eq!(
            FlowRegime::from_reynolds_number(2000.0),
            FlowRegime::Transitional
        );
        assert_eq!(
            FlowRegime::from_reynolds_number(4000.0),
            FlowRegime::Transitional
        );
        assert_eq!(
            FlowRegime::from_reynolds_number(4000.001),
            FlowRegime::Turbulent
        );
        // A real case from the slice: water in a 100 mm pipe at 1.5 m/s.
        assert_eq!(
            FlowRegime::from_reynolds_number(149_401.197_6),
            FlowRegime::Turbulent
        );
    }

    #[test]
    fn only_transitional_is_indeterminate() {
        assert!(!FlowRegime::Laminar.is_indeterminate());
        assert!(FlowRegime::Transitional.is_indeterminate());
        assert!(!FlowRegime::Turbulent.is_indeterminate());
    }

    #[test]
    fn regime_strings_are_distinct() {
        let mut seen = std::collections::HashSet::new();
        for regime in FlowRegime::all() {
            assert!(seen.insert(regime.as_str()));
        }
        assert_eq!(FlowRegime::Turbulent.to_string(), "turbulent");
    }
}
