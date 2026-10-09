//! `characterization.tbp_closure` - a cut's molar mass from its boiling point and gravity.
//!
//! Spec: `specs/models/characterization/tbp_closure.toml`. Oracle:
//! `validation/neqsim/captures/tbp_closure_probe.tsv`.
//!
//! # The bisection is NeqSim's
//!
//! Two of the four members answer in closed form and two ask "which molar mass, fed forward,
//! reproduces this boiling point?" - and that question is answered by halving the fixed bracket
//! `[0.010, 0.800]` kg/mol until the bracket is narrower than `1e-9`, at most 200 times. The
//! value returned is the **final bracket's midpoint**, not the last midpoint evaluated: the loop
//! returns `0.5 * (lower + upper)` after it stops, so the answer is a property of the bracket it
//! left behind. Reproduced exactly, including that, because a port that returned the last
//! midpoint would agree to twelve places and differ in the thirteenth.
//!
//! Where the endpoints share a sign the pair is not attainable and NeqSim refuses rather than
//! narrowing - the message names the range the bracket can reach, which is what makes the
//! refusal actionable.
//!
//! # The other direction is `characterization.tbp_density`
//!
//! Only the 1980 Riazi-Daubert pair inverts for specific gravity, which is a separate id for
//! that reason rather than an input here.

use azoth_core::units::{MassDensity, ThermodynamicTemperature, kilograms_per_mole};
use azoth_core::{AzothError, Result, apply_checks};
use std::str::FromStr;

use crate::model_gen;
use crate::results::TbpClosureResult;
use crate::tbp_cut_properties::{TbpModel, model_boiling_point};

/// The search bracket, in kg/mol, and its two stopping rules. NeqSim's own constants.
const SEARCH_LOWER: f64 = 0.010;
const SEARCH_UPPER: f64 = 0.800;
const TOLERANCE: f64 = 1.0e-9;
const MAX_ITERATIONS: u32 = 200;

/// Which of `TbpClosure`'s four members answers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum TbpClosureKind {
    /// The 1980 Riazi-Daubert pair, in closed form. The default, and the only member that also
    /// supports the density direction.
    #[default]
    RiaziDaubert1980,
    /// The 1987 Riazi-Daubert pair, in closed form.
    RiaziDaubert1987,
    /// Soreide's correlation, by bisection.
    Soreide,
    /// A TBP model's own boiling point, by bisection. Needs the `model` input.
    TbpModel,
}

impl FromStr for TbpClosureKind {
    type Err = std::convert::Infallible;

    /// An unrecognised name is the 1980 pair, which is what `getModel`-style defaults do
    /// elsewhere in this crate; the spec's own vocabulary check is what should refuse.
    fn from_str(text: &str) -> std::result::Result<Self, Self::Err> {
        Ok(match text.trim() {
            "riazi_daubert_1987" => Self::RiaziDaubert1987,
            "soreide" => Self::Soreide,
            "tbp_model" => Self::TbpModel,
            _ => Self::RiaziDaubert1980,
        })
    }
}

/// `TbpClosure.calcBoilingPointSoreide`, `[K]`.
///
/// The forward form the `SOREIDE` member inverts. `molarMass` is kg/mol and the class scales it
/// to g/mol to match the correlation's own scale, which is the same trap the ten models carry.
#[must_use]
pub fn soreide_boiling_point(molar_mass: f64, density: f64) -> f64 {
    let molar_mass_gmol = molar_mass * 1000.0;
    let rankine = 1928.3
        - 1.695e5
            * molar_mass_gmol.powf(-0.03522)
            * density.powf(3.266)
            * (-4.922e-3 * molar_mass_gmol - 4.7685 * density
                + 3.462e-3 * molar_mass_gmol * density)
                .exp();
    rankine / 1.8
}

/// The forward boiling point a member's bisection is inverting, `[K]`.
fn forward_boiling_point(
    closure: TbpClosureKind,
    molar_mass: f64,
    density: f64,
    model: Option<TbpModel>,
) -> Result<f64> {
    match closure {
        TbpClosureKind::TbpModel => {
            let model = model.ok_or_else(|| {
                AzothError::invalid_input(
                    "model",
                    "must be supplied for the `tbp_model` closure, which bisects that model's own \
                     boiling point",
                )
            })?;
            Ok(model_boiling_point(model, molar_mass * 1000.0, density))
        }
        _ => Ok(soreide_boiling_point(molar_mass, density)),
    }
}

/// The 1980 Riazi-Daubert pair, `[kg/mol]`. Specific gravity in g/cm3.
fn riazi_daubert_1980_boiling_point_to_molar_mass(boiling_point: f64, density: f64) -> f64 {
    4.5673e-5 * (boiling_point * 1.8).powf(2.1962) * density.powf(-1.0164) / 1000.0
}

/// The 1987 Riazi-Daubert pair, `[kg/mol]`.
fn riazi_daubert_1987_boiling_point_to_molar_mass(boiling_point: f64, density: f64) -> f64 {
    let molar_mass_gmol = 42.965
        * (2.097e-4 * boiling_point - 7.78712 * density + 2.08476e-3 * boiling_point * density)
            .exp()
        * boiling_point.powf(1.26007)
        * density.powf(4.98308);
    molar_mass_gmol / 1000.0
}

/// The bisection, reproducing `TbpClosure.solveMolarMass`.
fn solve_molar_mass(
    closure: TbpClosureKind,
    boiling_point: f64,
    density: f64,
    model: Option<TbpModel>,
) -> Result<f64> {
    let mut lower = SEARCH_LOWER;
    let mut upper = SEARCH_UPPER;
    let mut f_lower = forward_boiling_point(closure, lower, density, model)? - boiling_point;
    let f_upper = forward_boiling_point(closure, upper, density, model)? - boiling_point;
    if f_lower * f_upper > 0.0 {
        return Err(AzothError::invalid_input(
            "boiling_point",
            format!(
                "{boiling_point} K is not attainable with this closure at a specific gravity of \
                 {density}. Its range over the molar-mass bracket is {} to {} K, for {} to {} g/mol",
                f_lower + boiling_point,
                f_upper + boiling_point,
                SEARCH_LOWER * 1000.0,
                SEARCH_UPPER * 1000.0
            ),
        ));
    }
    let mut iterations = 0_u32;
    while iterations < MAX_ITERATIONS && (upper - lower) > TOLERANCE {
        let molar_mass = 0.5 * (lower + upper);
        let f_mid = forward_boiling_point(closure, molar_mass, density, model)? - boiling_point;
        if f_mid == 0.0 {
            return Ok(molar_mass);
        }
        if f_lower * f_mid < 0.0 {
            upper = molar_mass;
        } else {
            lower = molar_mass;
            f_lower = f_mid;
        }
        iterations += 1;
    }
    // The **bracket's** midpoint, not the last one evaluated - see the module documentation.
    Ok(0.5 * (lower + upper))
}

/// A cut's molar mass, from its normal boiling point and its specific gravity.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `closure` is `tbp_model` and no `model` is supplied, or if
///   the boiling point is not attainable over the search bracket.
/// * [`AzothError::OutOfRange`] if `boiling_point` or `density` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kelvins};
/// use azoth_characterization::{TbpClosureKind, tbp_closure};
///
/// let r = tbp_closure(
///     TbpClosureKind::RiaziDaubert1980,
///     kelvins(300.0),
///     kilograms_per_cubic_meter(650.0),
///     None,
/// )?;
/// assert!((r.molar_mass.value - 0.0709083824627520).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tbp_closure(
    closure: TbpClosureKind,
    boiling_point: ThermodynamicTemperature,
    density: MassDensity,
    model: Option<TbpModel>,
) -> Result<TbpClosureResult> {
    let spec = &model_gen::TBP_CLOSURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "boiling_point" => Some(boiling_point.value),
            "density" => Some(density.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The boundary: the correlations take g/cm3, which is the same number a thousand smaller.
    let specific_gravity = density.value / 1000.0;
    let molar_mass = match closure {
        TbpClosureKind::RiaziDaubert1980 => {
            riazi_daubert_1980_boiling_point_to_molar_mass(boiling_point.value, specific_gravity)
        }
        TbpClosureKind::RiaziDaubert1987 => {
            riazi_daubert_1987_boiling_point_to_molar_mass(boiling_point.value, specific_gravity)
        }
        TbpClosureKind::Soreide | TbpClosureKind::TbpModel => {
            solve_molar_mass(closure, boiling_point.value, specific_gravity, model)?
        }
    };

    Ok(TbpClosureResult {
        molar_mass: kilograms_per_mole(molar_mass),
        warnings,
    })
}
