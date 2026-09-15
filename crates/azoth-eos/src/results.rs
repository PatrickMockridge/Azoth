//! Result types for the equations-of-state calculations.
//!
//! The same contract as every other namespace's results: one struct per
//! calculation, field names identical to the Python result dataclass and listed in
//! [`CalcResult::FIELDS`], with a test asserting the three agree. See
//! `crates/azoth-core/src/result.rs` for why the duplication is deliberate.

use azoth_core::units::{
    MassDensity, MolarEnergy, MolarHeatCapacity, MolarVolume, Pressure, ThermodynamicTemperature,
};
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
    /// Newton steps the polish took, summed over the roots. Reported because the
    /// answer alone does not say whether the solver did any work.
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
    /// The departure heat capacity over `R`.
    pub cp_dep_r: f64,
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
/// `[0, 1]`: where a bisection stops on an identically-zero function is a ratio of
/// round-off. Reporting a number there would be reporting a fabrication that looks
/// exactly like a real vapour fraction. `None` is the honest answer, and it is
/// type-level: a caller cannot read it without noticing.
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

/// Result of `eos.ph_flash`.
///
/// The state a mixture reaches when a duty is applied at a fixed pressure: the
/// temperature that satisfies the energy balance, and the phase split at it. The split
/// is reported in as much detail as [`PtFlashResult`] because it *is* one - the flash
/// evaluated at the answer - and a caller who needs the compositions should not have to
/// run it again to get them.
#[derive(Debug, Clone)]
pub struct PhFlashResult {
    /// The temperature that satisfies the enthalpy. This is the model's answer.
    pub temperature: ThermodynamicTemperature,
    /// The vapour fraction at that temperature, or `None` for a single-phase feed.
    ///
    /// `None` rather than a number outside `[0, 1]`: the flash extrapolates a split
    /// that does not exist, and reporting it would invite a caller to use it.
    pub beta: Option<f64>,
    /// Liquid-phase mole fractions at the answer.
    pub x: Vec<f64>,
    /// Vapour-phase mole fractions at the answer.
    pub y: Vec<f64>,
    /// `K_i = y_i / x_i` at the answer.
    pub k: Vec<f64>,
    /// Which phase the feed is in at the answer.
    pub phase: Phase,
    /// The liquid root of the cubic at the answer.
    pub z_liquid: f64,
    /// The vapour root.
    pub z_vapour: f64,
    /// Bisection steps taken.
    pub iterations: u32,
    /// `|H(T) - H_target| / max(|H_target|, 1)` at the answer.
    pub residual: f64,
    /// Caveats, deduplicated - the search evaluates the flash thousands of times.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PhFlashResult {
    const CALC_ID: &'static str = "eos.ph_flash";
    const FIELDS: &'static [&'static str] = &[
        "T",
        "beta",
        "x",
        "y",
        "k",
        "phase",
        "z_liquid",
        "z_vapour",
        "iterations",
        "residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.ps_flash`.
///
/// The state a stream reaches when it is expanded or compressed **isentropically** at a
/// fixed pressure. The same shape as [`PhFlashResult`], deliberately: the two models
/// differ in which property they invert and agree on everything else, so a caller
/// reading one already knows how to read the other.
#[derive(Debug, Clone)]
pub struct PsFlashResult {
    /// The temperature that satisfies the entropy. This is the model's answer.
    pub temperature: ThermodynamicTemperature,
    /// The vapour fraction at that temperature, or `None` for a single-phase feed.
    pub beta: Option<f64>,
    /// Liquid-phase mole fractions at the answer.
    pub x: Vec<f64>,
    /// Vapour-phase mole fractions at the answer.
    pub y: Vec<f64>,
    /// `K_i = y_i / x_i` at the answer.
    pub k: Vec<f64>,
    /// Which phase the feed is in at the answer.
    pub phase: Phase,
    /// The liquid root of the cubic at the answer.
    pub z_liquid: f64,
    /// The vapour root.
    pub z_vapour: f64,
    /// Bisection steps taken.
    pub iterations: u32,
    /// `|S(T) - S_target| / max(|S_target|, 1)` at the answer.
    pub residual: f64,
    /// Caveats, deduplicated.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PsFlashResult {
    const CALC_ID: &'static str = "eos.ps_flash";
    const FIELDS: &'static [&'static str] = &[
        "T",
        "beta",
        "x",
        "y",
        "k",
        "phase",
        "z_liquid",
        "z_vapour",
        "iterations",
        "residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Whether a feed is stable as a single phase.
///
/// Two values rather than a boolean because the *asymmetry* between them is the
/// point: [`Self::Unstable`] is a proof - a trial reached a stationary point below
/// the tangent plane, so a single phase is not the Gibbs minimum - while
/// [`Self::Stable`] is the absence of one, from two trials that were placed by a
/// gas-liquid correlation. A caller who reads `stable` as "no split exists" has read
/// it wrong, and a bare `true` invites exactly that.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum StabilityVerdict {
    /// Neither trial found a stationary point below the tangent plane.
    Stable,
    /// At least one trial did. A single phase is not the Gibbs minimum here.
    Unstable,
}

impl StabilityVerdict {
    /// The spec's spelling.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Stable => "stable",
            Self::Unstable => "unstable",
        }
    }
}

/// Result of `eos.stability_test`.
///
/// # Why `tm` and `w` are fixed-length vectors rather than "one per phase found"
///
/// Two trials are run, always, so both are length two in a fixed order - the
/// vapour-like trial first. A trial that converges to the feed itself still has a
/// tangent-plane distance, and it is that near-zero number which is the evidence it
/// was trivial; reporting only the trials that found something would make the
/// vector's *index* mean something different on every state.
#[derive(Debug, Clone, PartialEq)]
pub struct StabilityTestResult {
    /// Whether the feed is stable as a single phase.
    pub verdict: StabilityVerdict,
    /// The tangent-plane distance at each trial's stationary point, vapour-like
    /// trial first. Negative means that trial lies below the tangent plane.
    pub tm: Vec<f64>,
    /// The stationary-point composition of each trial, in the same order.
    pub w: Vec<Vec<f64>>,
    /// Iterations each trial took, in the same order.
    pub iterations: Vec<u32>,
    /// The smallest `T / Tc_i` over the components.
    pub min_t_over_tc: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for StabilityTestResult {
    const CALC_ID: &'static str = "eos.stability_test";
    const FIELDS: &'static [&'static str] = &[
        "verdict",
        "tm",
        "w",
        "iterations",
        "min_t_over_tc",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.bubble_pressure`.
///
/// Deliberately sharing a shape with [`DewPressureResult`] rather than one type
/// with a switch: the two differ in *which* composition is the input, and a shared
/// field would have to be named after neither.
#[derive(Debug, Clone, PartialEq)]
pub struct BubblePressureResult {
    /// The bubble-point pressure.
    pub pressure: Pressure,
    /// The composition of the vapour that first appears.
    pub incipient: Vec<f64>,
    /// K-values at the converged pressure.
    pub k: Vec<f64>,
    /// The liquid root of the cubic at the converged state.
    pub z_liquid: f64,
    /// The vapour root.
    pub z_vapour: f64,
    /// The smallest `T / Tc_i` over the components.
    pub min_t_over_tc: f64,
    /// Pressure updates taken.
    pub iterations: u32,
    /// `|sum_i x_i K_i - 1|` at the last completed step.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for BubblePressureResult {
    const CALC_ID: &'static str = "eos.bubble_pressure";
    const FIELDS: &'static [&'static str] = &[
        "pressure",
        "incipient",
        "k",
        "z_liquid",
        "z_vapour",
        "min_t_over_tc",
        "iterations",
        "residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.dew_pressure`.
#[derive(Debug, Clone, PartialEq)]
pub struct DewPressureResult {
    /// The dew-point pressure.
    pub pressure: Pressure,
    /// The composition of the liquid that first appears.
    pub incipient: Vec<f64>,
    /// K-values at the converged pressure.
    pub k: Vec<f64>,
    /// The liquid root of the cubic at the converged state.
    pub z_liquid: f64,
    /// The vapour root.
    pub z_vapour: f64,
    /// The smallest `T / Tc_i` over the components.
    pub min_t_over_tc: f64,
    /// Pressure updates taken.
    pub iterations: u32,
    /// `|sum_i y_i / K_i - 1|` at the last completed step.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for DewPressureResult {
    const CALC_ID: &'static str = "eos.dew_pressure";
    const FIELDS: &'static [&'static str] = &[
        "pressure",
        "incipient",
        "k",
        "z_liquid",
        "z_vapour",
        "min_t_over_tc",
        "iterations",
        "residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.ideal_gas_cp`.
///
/// Reports the polynomial's own dimensionless value as well as the dimensioned heat
/// capacity, because the dimensionless form is what a reader checking the arithmetic
/// by hand computes first - and because the multiplication by `R` is then visible as
/// the single step it is, rather than folded into the answer.
#[derive(Debug, Clone, PartialEq)]
pub struct IdealGasCpResult {
    /// The ideal-gas heat capacity. Carries `OutOfValidRange` when it is not
    /// positive, which means the polynomial has been evaluated outside its range.
    pub cp: MolarHeatCapacity,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Result of `eos.molar_enthalpy_entropy`.
///
/// Reports the ideal-gas and departure parts separately as well as their sum, because
/// the split is what a caller checking the answer needs: `h_ideal` carries the datum,
/// `h_departure` carries the equation of state, and a single total hides which of the
/// two a disagreement came from.
#[derive(Debug, Clone, PartialEq)]
pub struct MolarEnthalpyEntropyResult {
    /// The molar enthalpy, `h_ideal + h_departure`.
    pub h: MolarEnergy,
    /// The molar entropy, `s_ideal + s_departure`.
    pub s: MolarHeatCapacity,
    /// The ideal-gas part of the enthalpy - the reference values and the integrals.
    pub h_ideal: MolarEnergy,
    /// The ideal-gas part of the entropy.
    pub s_ideal: MolarHeatCapacity,
    /// The residual enthalpy, `R*T*h_dep_rt`.
    pub h_departure: MolarEnergy,
    /// The residual entropy, `R*s_dep_r`. Does **not** include the entropy of mixing.
    pub s_departure: MolarHeatCapacity,
    /// The composition-weighted average of the components' `psi`.
    pub psi_bar: f64,
    /// The molar heat capacity at constant pressure, `cp_ideal + cp_departure`.
    ///
    /// The derivative the isentropic and isenthalpic flashes step on: theirs is
    /// `dS/dT = cp/T`, in the entropy's case, and `dH/d(1/T) = -T**2*cp` in the
    /// enthalpy's.
    pub cp: MolarHeatCapacity,
    /// The ideal-gas part of the heat capacity - the polynomial, evaluated at `T`.
    pub cp_ideal: MolarHeatCapacity,
    /// The residual heat capacity, `R*cp_dep_r`.
    pub cp_departure: MolarHeatCapacity,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for MolarEnthalpyEntropyResult {
    const CALC_ID: &'static str = "eos.molar_enthalpy_entropy";
    const FIELDS: &'static [&'static str] = &[
        "h",
        "s",
        "h_ideal",
        "s_ideal",
        "h_departure",
        "s_departure",
        "psi_bar",
        "cp",
        "cp_ideal",
        "cp_departure",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

impl CalcResult for IdealGasCpResult {
    const CALC_ID: &'static str = "eos.ideal_gas_cp";
    const FIELDS: &'static [&'static str] = &["cp", "warnings"];

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
    const FIELDS: &'static [&'static str] =
        &["ln_phi", "h_dep_rt", "s_dep_r", "cp_dep_r", "warnings"];

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

/// Result of `eos.srk_kappa`.
#[derive(Debug, Clone, PartialEq)]
pub struct SrkKappaResult {
    /// The Soave-Redlich-Kwong alpha-function coefficient. Dimensionless.
    pub kappa: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SrkKappaResult {
    const CALC_ID: &'static str = "eos.srk_kappa";
    const FIELDS: &'static [&'static str] = &["kappa", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.pr_peneloux_shift`.
#[derive(Debug, Clone, PartialEq)]
pub struct PrPenelouxShiftResult {
    /// The Peneloux volume-translation parameter.
    pub c: MolarVolume,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for PrPenelouxShiftResult {
    const CALC_ID: &'static str = "eos.pr_peneloux_shift";
    const FIELDS: &'static [&'static str] = &["c", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.srk_peneloux_shift`.
#[derive(Debug, Clone, PartialEq)]
pub struct SrkPenelouxShiftResult {
    /// The Peneloux volume-translation parameter.
    pub c: MolarVolume,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SrkPenelouxShiftResult {
    const CALC_ID: &'static str = "eos.srk_peneloux_shift";
    const FIELDS: &'static [&'static str] = &["c", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.heat_of_vaporization`.
#[derive(Debug, Clone, PartialEq)]
pub struct HeatOfVaporizationResult {
    /// The pure-component heat of vaporisation.
    pub hov: MolarEnergy,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for HeatOfVaporizationResult {
    const CALC_ID: &'static str = "eos.heat_of_vaporization";
    const FIELDS: &'static [&'static str] = &["hov", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.liquid_heat_capacity`.
#[derive(Debug, Clone, PartialEq)]
pub struct LiquidHeatCapacityResult {
    /// The pure-component liquid heat capacity.
    pub cp: MolarHeatCapacity,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for LiquidHeatCapacityResult {
    const CALC_ID: &'static str = "eos.liquid_heat_capacity";
    const FIELDS: &'static [&'static str] = &["cp", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.antoine_vapor_pressure`.
#[derive(Debug, Clone, PartialEq)]
pub struct AntoineVaporPressureResult {
    /// The pure-component vapour pressure.
    pub p_sat: Pressure,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for AntoineVaporPressureResult {
    const CALC_ID: &'static str = "eos.antoine_vapor_pressure";
    const FIELDS: &'static [&'static str] = &["p_sat", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.rackett_molar_volume`.
#[derive(Debug, Clone, PartialEq)]
pub struct RackettMolarVolumeResult {
    /// The saturated liquid molar volume.
    pub v: MolarVolume,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for RackettMolarVolumeResult {
    const CALC_ID: &'static str = "eos.rackett_molar_volume";
    const FIELDS: &'static [&'static str] = &["v", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.costald_molar_volume`.
#[derive(Debug, Clone, PartialEq)]
pub struct CostaldMolarVolumeResult {
    /// The saturated liquid molar volume.
    pub v: MolarVolume,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for CostaldMolarVolumeResult {
    const CALC_ID: &'static str = "eos.costald_molar_volume";
    const FIELDS: &'static [&'static str] = &["v", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.srk_alpha_ab`.
#[derive(Debug, Clone, PartialEq)]
pub struct SrkAlphaAbResult {
    /// The Soave alpha function.
    pub alpha: f64,
    /// `A = a*alpha*P/(R**2*T**2)`, the dimensionless attraction parameter.
    pub a_reduced: f64,
    /// `B = b*P/(R*T)`, the dimensionless repulsion parameter.
    pub b_reduced: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SrkAlphaAbResult {
    const CALC_ID: &'static str = "eos.srk_alpha_ab";
    const FIELDS: &'static [&'static str] = &["alpha", "a_reduced", "b_reduced", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.srk_z_factor`.
#[derive(Debug, Clone, PartialEq)]
pub struct SrkZFactorResult {
    /// The smallest admissible root.
    pub z_min: f64,
    /// The largest admissible root.
    pub z_max: f64,
    /// How many admissible roots there were.
    pub root_structure: RootStructure,
    /// Newton steps the polish took.
    pub iterations: u32,
    /// Whether the polish met its stopping rule.
    pub converged: bool,
    /// The largest `|x_k - x_{k-1}|` at the final polish step.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SrkZFactorResult {
    const CALC_ID: &'static str = "eos.srk_z_factor";
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

/// Result of `eos.srk_departure`.
#[derive(Debug, Clone, PartialEq)]
pub struct SrkDepartureResult {
    /// The logarithm of the fugacity coefficient.
    pub ln_phi: f64,
    /// The departure enthalpy over `R*T`.
    pub h_dep_rt: f64,
    /// The departure entropy over `R`.
    pub s_dep_r: f64,
    /// The departure heat capacity over `R`.
    pub cp_dep_r: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for SrkDepartureResult {
    const CALC_ID: &'static str = "eos.srk_departure";
    const FIELDS: &'static [&'static str] =
        &["ln_phi", "h_dep_rt", "s_dep_r", "cp_dep_r", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.rk_alpha_ab`.
#[derive(Debug, Clone, PartialEq)]
pub struct RkAlphaAbResult {
    /// The Redlich-Kwong alpha function, `1/sqrt(Tr)`.
    pub alpha: f64,
    /// `A = a*alpha*P/(R**2*T**2)`, the dimensionless attraction parameter.
    pub a_reduced: f64,
    /// `B = b*P/(R*T)`, the dimensionless repulsion parameter.
    pub b_reduced: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for RkAlphaAbResult {
    const CALC_ID: &'static str = "eos.rk_alpha_ab";
    const FIELDS: &'static [&'static str] = &["alpha", "a_reduced", "b_reduced", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.rk_departure`.
#[derive(Debug, Clone, PartialEq)]
pub struct RkDepartureResult {
    /// The logarithm of the fugacity coefficient.
    pub ln_phi: f64,
    /// The departure enthalpy over `R*T`.
    pub h_dep_rt: f64,
    /// The departure entropy over `R`.
    pub s_dep_r: f64,
    /// The departure heat capacity over `R`.
    pub cp_dep_r: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for RkDepartureResult {
    const CALC_ID: &'static str = "eos.rk_departure";
    const FIELDS: &'static [&'static str] =
        &["ln_phi", "h_dep_rt", "s_dep_r", "cp_dep_r", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.pr78_kappa`.
#[derive(Debug, Clone, PartialEq)]
pub struct Pr78KappaResult {
    /// The 1978 Peng-Robinson alpha-function coefficient. Dimensionless.
    pub kappa: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for Pr78KappaResult {
    const CALC_ID: &'static str = "eos.pr78_kappa";
    const FIELDS: &'static [&'static str] = &["kappa", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.twu_kappa`.
#[derive(Debug, Clone, PartialEq)]
pub struct TwuKappaResult {
    /// Twu's alpha-function coefficient. Dimensionless.
    pub kappa: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TwuKappaResult {
    const CALC_ID: &'static str = "eos.twu_kappa";
    const FIELDS: &'static [&'static str] = &["kappa", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `eos.critical_point`.
///
/// The four state variables of a mixture critical point. `z_c` is here because it is the
/// quantity that distinguishes this model from the mechanical conditions: a pure
/// component's is `(1 - omega_b)/3`, and a mixture's varies with composition.
#[derive(Debug, Clone, PartialEq)]
pub struct CriticalPointResult {
    /// The critical temperature.
    pub tc: ThermodynamicTemperature,
    /// The critical pressure.
    pub pc: Pressure,
    /// The critical molar volume.
    pub vc: MolarVolume,
    /// `Pc Vc/(R Tc)`.
    pub z_c: f64,
    /// Outer iterations taken.
    pub iterations: u32,
    /// `max(|smallest eigenvalue|, |cubic form|)` at the returned state.
    pub residual: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for CriticalPointResult {
    const CALC_ID: &'static str = "eos.critical_point";
    const FIELDS: &'static [&'static str] = &[
        "tc",
        "pc",
        "vc",
        "z_c",
        "iterations",
        "residual",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
