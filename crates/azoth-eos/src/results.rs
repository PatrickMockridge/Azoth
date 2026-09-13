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

/// Result of `eos.pr_alpha_ab`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrAlphaAbResult {
    /// The alpha function, where Peng-Robinson's temperature dependence lives.
    /// A square, so never negative, and exactly 1 at `Tr = 1` whatever `kappa` is.
    pub alpha: f64,
    /// `A = a*alpha*P/(R**2*T**2)`, the dimensionless attraction parameter.
    pub a_reduced: f64,
    /// `B = b*P/(R*T)`, the dimensionless repulsion parameter.
    pub b_reduced: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PrAlphaAbResult {
    const CALC_ID: &'static str = "eos.pr_alpha_ab";
    const FIELDS: &'static [&'static str] = &["alpha", "a_reduced", "b_reduced", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// How many admissible real roots the cubic had.
///
/// Reported so a caller can tell a single-root state from one where `z_min` and
/// `z_max` are genuinely two different roots, without comparing floats.
///
/// **There is no `Two`, and that is a theorem.** The polynomial at `z = b_reduced`
/// is exactly `-2*b_reduced**2` - see the spec's `root_structure` description for
/// the algebra - so `B` lies either below all three roots or between the middle
/// and the largest one. The admissible count is therefore 1 or 3 and never 2. A
/// variant that cannot occur would be a value a caller branches on and never sees,
/// which is worse than an absent one, so it is absent.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RootStructure {
    /// One admissible root: `z_min` and `z_max` are the same number.
    One,
    /// Three admissible roots, the usual subcritical case.
    Three,
}

impl RootStructure {
    /// The spec's spelling.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::One => "one_root",
            Self::Three => "three_roots",
        }
    }
}

/// Result of `eos.pr_z_factor`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrZFactorResult {
    /// The smallest admissible root.
    pub z_min: f64,
    /// The largest admissible root. Equal to `z_min` when only one is admissible.
    pub z_max: f64,
    /// How many admissible roots there were.
    pub root_structure: RootStructure,
    /// Newton steps the polish took, summed over the roots. Carried because the
    /// answer alone does not say whether the solver did any work, and because the
    /// cross-language agreement test compares iteration counts as the sharpest
    /// cheap check that both implementations ran the same scheme.
    pub iterations: u32,
    /// Whether the polish met its stopping rule.
    pub converged: bool,
    /// The largest `|x_k - x_{k-1}|` at the final polish step.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Result of `eos.prsv_kappa`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrsvKappaResult {
    /// The PRSV alpha-function coefficient. Dimensionless, and unlike
    /// Peng-Robinson's it varies with temperature.
    pub kappa: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PrsvKappaResult {
    const CALC_ID: &'static str = "eos.prsv_kappa";
    const FIELDS: &'static [&'static str] = &["kappa", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

impl CalcResult for PrZFactorResult {
    const CALC_ID: &'static str = "eos.pr_z_factor";
    const FIELDS: &'static [&'static str] = &[
        "z_min",
        "z_max",
        "root_structure",
        "iterations",
        "converged",
        "residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
