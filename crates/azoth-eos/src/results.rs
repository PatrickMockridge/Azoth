//! Result types for the equations-of-state calculations.
//!
//! The same contract as every other namespace's results: one struct per
//! calculation, field names identical to the Python result dataclass and listed in
//! [`CalcResult::FIELDS`], with a test asserting the three agree. See
//! `crates/azoth-core/src/result.rs` for why the duplication is deliberate.

use azoth_core::units::{MassDensity, MolarVolume, Pressure};
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

/// Result of `eos.pr_departure`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrDepartureResult {
    /// The logarithm of the fugacity coefficient.
    pub ln_phi: f64,
    /// The departure enthalpy over `R*T`.
    pub h_dep_rt: f64,
    /// The departure entropy over `R`.
    pub s_dep_r: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Result of `eos.vdw1f_mix_binary`.
#[derive(Debug, Clone, PartialEq)]
pub struct Vdw1fMixBinaryResult {
    /// The mixture's attraction parameter.
    pub a_mix: f64,
    /// The mixture's repulsion parameter.
    pub b_mix: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Result of `eos.pr_molar_volume`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrMolarVolumeResult {
    /// Molar volume. The namespace's only dimensioned output.
    pub v: MolarVolume,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Result of `eos.pure_saturation`.
#[derive(Debug, Clone, PartialEq)]
pub struct PureSaturationResult {
    /// The saturation pressure.
    pub p_sat: Pressure,
    /// The common value of `ln phi_L` and `ln phi_V` at the converged pressure.
    pub ln_phi: f64,
    /// Bisection steps taken.
    pub iterations: u32,
    /// The dimensionless half-width of the final bracket - the relative uncertainty
    /// in the reduced pressure, not the fugacity residual.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PureSaturationResult {
    const CALC_ID: &'static str = "eos.pure_saturation";
    const FIELDS: &'static [&'static str] =
        &["p_sat", "ln_phi", "iterations", "residual", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// What a converged flash turned out to be.
///
/// A separate type from [`RootStructure`], which counts the roots of a *pure*
/// component's cubic and says nothing about phases. The two are easy to confuse and
/// mean different things: three roots is a mathematical fact about a polynomial,
/// `two_phase` is a physical claim about a mixture.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Phase {
    /// A genuine split: `beta` in `[0, 1]` and the two compositions differ.
    TwoPhase,
    /// The feed is subcooled liquid.
    ///
    /// Reached two ways, and `beta` distinguishes them: either the converged
    /// Rachford-Rice root is negative, in which case `beta` is present as the
    /// negative-flash value and the result carries `OutOfValidRange`; or every
    /// K-value is below one, in which case **no root exists at all** and `beta` is
    /// absent.
    AllLiquid,
    /// The feed is superheated vapour. The same two routes as [`Self::AllLiquid`].
    AllVapour,
    /// The iteration converged to `x = y = z`.
    ///
    /// The feed is single phase, and **this does not say which one** - the
    /// K-values straddled one throughout, so nothing in the model ever proved which
    /// phase the feed is. That is what a tangent-plane stability analysis decides,
    /// and this model has none; see the model spec's assumptions. `beta` is absent,
    /// because at the trivial solution it is indeterminate rather than out of range.
    Trivial,
}

impl Phase {
    /// The spec's spelling.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::TwoPhase => "two_phase",
            Self::AllLiquid => "all_liquid",
            Self::AllVapour => "all_vapour",
            Self::Trivial => "trivial",
        }
    }
}

/// Result of `eos.pt_flash`.
///
/// # Why `beta` is optional here and nowhere else in this crate
///
/// At a trivial solution every `K_i` is 1, the Rachford-Rice function is identically
/// zero, and the vapour fraction is **indeterminate** rather than merely outside
/// `[0, 1]`. Successive substitution approaches the point through geometrically
/// growing `beta`, and where a bisection stops on an identically-zero function is a
/// ratio of round-off - measured at `-7.7e10` for one feed and `-2.2e11` for
/// another, neither reproducible across implementations.
///
/// Reporting a number there would be reporting a fabrication that looks exactly like
/// a real vapour fraction. `None` is the honest answer, and it is type-level: a
/// caller cannot read it without noticing.
#[derive(Debug, Clone, PartialEq)]
pub struct PtFlashResult {
    /// The vapour fraction, or `None` when the solution is trivial.
    pub beta: Option<f64>,
    /// Liquid-phase mole fractions.
    pub x: Vec<f64>,
    /// Vapour-phase mole fractions.
    pub y: Vec<f64>,
    /// `K_i = y_i / x_i`, the iterate the loop converges on.
    pub k: Vec<f64>,
    /// `ln phi_i` in the liquid phase.
    pub ln_phi_liquid: Vec<f64>,
    /// `ln phi_i` in the vapour phase.
    pub ln_phi_vapour: Vec<f64>,
    /// The liquid root of the cubic, the smallest admissible one.
    pub z_liquid: f64,
    /// The vapour root, the largest admissible one.
    pub z_vapour: f64,
    /// The smallest `T / Tc_i` over the components.
    pub min_t_over_tc: f64,
    /// What the converged state is.
    pub phase: Phase,
    /// Successive-substitution steps taken.
    pub iterations: u32,
    /// `rms_i |ln K_i - ln K_i_previous|` at the step that met the tolerance.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PtFlashResult {
    const CALC_ID: &'static str = "eos.pt_flash";
    const FIELDS: &'static [&'static str] = &[
        "beta",
        "x",
        "y",
        "k",
        "ln_phi_liquid",
        "ln_phi_vapour",
        "z_liquid",
        "z_vapour",
        "min_t_over_tc",
        "phase",
        "iterations",
        "residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

impl CalcResult for PrMolarVolumeResult {
    const CALC_ID: &'static str = "eos.pr_molar_volume";
    const FIELDS: &'static [&'static str] = &["v", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.pr_mass_density`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrMassDensityResult {
    /// Mass density.
    pub rho: MassDensity,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PrMassDensityResult {
    const CALC_ID: &'static str = "eos.pr_mass_density";
    const FIELDS: &'static [&'static str] = &["rho", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

impl CalcResult for Vdw1fMixBinaryResult {
    const CALC_ID: &'static str = "eos.vdw1f_mix_binary";
    const FIELDS: &'static [&'static str] = &["a_mix", "b_mix", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.rachford_rice_binary`.
#[derive(Debug, Clone, PartialEq)]
pub struct RachfordRiceBinaryResult {
    /// The vapour fraction that solves the Rachford-Rice equation. Outside `[0, 1]`
    /// the feed is single phase and this is the tangent-plane value rather than a
    /// phase split; the result carries `OutOfValidRange` when so.
    pub beta: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for RachfordRiceBinaryResult {
    const CALC_ID: &'static str = "eos.rachford_rice_binary";
    const FIELDS: &'static [&'static str] = &["beta", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

impl CalcResult for PrDepartureResult {
    const CALC_ID: &'static str = "eos.pr_departure";
    const FIELDS: &'static [&'static str] = &["ln_phi", "h_dep_rt", "s_dep_r", "warnings"];

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
