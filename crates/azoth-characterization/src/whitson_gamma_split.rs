//! `characterization.whitson_gamma_split` - a plus fraction split by Whitson's gamma distribution.
//!
//! Spec: `specs/models/characterization/whitson_gamma_split.toml`. Oracle:
//! `validation/neqsim/captures/plus_fraction_probe.tsv`.
//!
//! # A window in molar mass, not a carbon number
//!
//! The cuts are windows that begin at `eta` and step by the SCN increment of 14 g/mol, with the
//! last widened to 10 000 to catch the tail. `firstCarbonNumber` and `lastPlusFractionNumber` set
//! **how many** windows there are and nothing else - `WhitsonGammaModel` never reads the plus
//! component's name, so the first carbon number stays at the outer object's own value of one - and
//! the abundance is the gamma density integrated across each window.
//!
//! # Two correlations, one of which falls back to the other
//!
//! `densityUOP` is `6.0108 * M^0.17947 * Kw^-1.18241` on the plus fraction's own Watson factor;
//! `densitySoreide` is `0.2855 + C_f*(M - 66)^0.13` clamped to `[0.6, 1.2]`, and calls the Watson
//! form when its own `C_f` comes out non-positive. `UOP` is the class's default.
//!
//! # The gamma function is an eight-term fit
//!
//! `gamma` reduces its argument to `[0, 1)`, evaluates the classical polynomial for `Gamma(1+z)`,
//! and reflects below one by dividing by `X` - which is `Gamma(X+1)/X`, and so correct. It is
//! reproduced as written rather than replaced with a library call, because the polynomial's error
//! is part of every `P0` and `P1` the split is made of.

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::WhitsonGammaSplitResult;

/// The SCN increment the windows step by, g/mol.
const SCN_INCREMENT: f64 = 14.0;

/// The upper edge the last window is widened to, g/mol - the class's own tail catcher.
const LAST_WINDOW_UPPER: f64 = 10_000.0;

/// The floor the integrated window abundance is held above, before renormalisation.
const WINDOW_FLOOR: f64 = 1.0e-15;

/// The denominator below which a window's first moment is replaced by its midpoint.
const MOMENT_FLOOR: f64 = 1.0e-15;

/// The term the incomplete-gamma series stops at.
const SERIES_TOLERANCE: f64 = 1.0e-8;

/// The incomplete-gamma series' own cap.
const SERIES_MAX_TERMS: u32 = 10_000;

/// `WhitsonGammaModel`'s constructor values for the two shape parameters.
const DEFAULT_SHAPE: f64 = 1.0;
const DEFAULT_ETA_GMOL: f64 = 90.0;

/// The Søreide correlation's own floor on `M - 66`.
const SOREIDE_MOLAR_MASS_FLOOR: f64 = 66.0;
/// The Søreide gravity's clamp, which the class applies before storing it.
const SOREIDE_LOWER: f64 = 0.6;
const SOREIDE_UPPER: f64 = 1.2;
/// The Søreide intercept, which is also the gravity of a cut at the molar-mass floor.
const SOREIDE_INTERCEPT: f64 = 0.2855;

/// Which of the class's two gravity correlations answers.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum WhitsonDensityModel {
    /// The Watson factor form. `densityModel`'s default, and the other's fallback.
    #[default]
    Uop,
    /// The Søreide correlation, which falls back to [`Self::Uop`] on a non-positive factor.
    Soreide,
}

impl std::str::FromStr for WhitsonDensityModel {
    type Err = std::convert::Infallible;

    /// An unrecognised name is `UOP`, which is what `calculateDensities`' own `equalsIgnoreCase`
    /// test does with one; the spec's vocabulary check is what should refuse.
    fn from_str(text: &str) -> std::result::Result<Self, Self::Err> {
        Ok(match text.trim() {
            "soreide" => Self::Soreide,
            _ => Self::Uop,
        })
    }
}

/// `WhitsonGammaModel.gamma`, the eight-term fit for `Gamma(1+z)` on `[0, 1)` with the reflection.
#[must_use]
pub fn gamma(x: f64) -> f64 {
    const COEFFICIENTS: [f64; 8] = [
        -0.577_191_652,
        0.988_205_891,
        -0.897_056_937,
        0.918_206_857,
        -0.756_704_078,
        0.482_199_394,
        -0.193_527_818,
        0.035_868_343,
    ];
    let mut constant = 1.0;
    let mut reduced = if x < 1.0 { x + 1.0 } else { x };
    while reduced >= 2.0 {
        constant *= reduced - 1.0;
        reduced -= 1.0;
    }
    reduced -= 1.0;
    let mut polynomial = 1.0;
    let mut power = reduced;
    for coefficient in COEFFICIENTS {
        polynomial += coefficient * power;
        power *= reduced;
    }
    let value = constant * polynomial;
    if x < 1.0 { value / x } else { value }
}

/// `WhitsonGammaModel.P0P1`: the incomplete-gamma pair at one molar mass, g/mol.
///
/// A window exactly at `eta` returns zero without evaluating, which is what makes the first
/// window's lower edge free rather than a limit.
#[must_use]
pub fn p0_p1(molar_mass: f64, eta: f64, shape: f64, scale: f64) -> [f64; 2] {
    if molar_mass == eta {
        return [0.0, 0.0];
    }
    let y = (molar_mass - eta) / scale;
    let q = (-y).exp() * y.powf(shape) / gamma(shape);
    let mut term = 1.0 / shape;
    let mut sum = term;
    for j in 1..=SERIES_MAX_TERMS {
        term *= y / (shape + f64::from(j));
        sum += term;
        if term.abs() <= SERIES_TOLERANCE {
            return [q * sum, q * (sum - 1.0 / shape)];
        }
    }
    // The class leaves `P0` and `P1` at zero when the series has not broken by the cap.
    [0.0, 0.0]
}

/// The plus fraction's Watson factor, `getWatsonKFactor`. `dens_plus` is a specific gravity.
#[must_use]
pub fn watson_factor(m_plus: f64, dens_plus: f64) -> f64 {
    4.5579 * (m_plus * 1000.0).powf(0.151_78) * dens_plus.powf(-1.182_41)
}

/// `WhitsonGammaModel.estimateAlpha`, the four bands on the Watson factor.
///
/// The class's own signature calls its second argument a density in kg/m3 and its caller passes
/// the plus fraction's specific gravity in g/cm3. The gravity is what is passed here, faithfully.
#[must_use]
pub fn estimate_shape(m_plus: f64, dens_plus: f64) -> f64 {
    let kw = watson_factor(m_plus, dens_plus);
    if kw >= 12.5 {
        0.5 + 0.1 * (kw - 12.5)
    } else if kw >= 11.5 {
        1.0 + 0.5 * (kw - 11.5)
    } else if kw >= 10.5 {
        1.5 + 0.5 * (kw - 10.5)
    } else {
        2.0 + 0.5 * (10.5 - kw)
    }
}

/// `densityUOP`: every cut's gravity from the plus fraction's Watson factor, in g/cm3.
fn density_uop(molar_masses: &[f64], kw: f64) -> Vec<f64> {
    molar_masses
        .iter()
        .map(|m| if *m > 0.0 { 6.0108 * (m * 1000.0).powf(0.179_47) * kw.powf(-1.182_41) } else { 0.0 })
        .collect()
}

/// `densitySoreide`, which falls back to [`density_uop`] on a non-positive factor.
fn density_soreide(molar_masses: &[f64], m_plus: f64, dens_plus: f64, kw: f64) -> Vec<f64> {
    let exponent = (m_plus * 1000.0 - SOREIDE_MOLAR_MASS_FLOOR).max(1.0);
    let factor = (dens_plus - SOREIDE_INTERCEPT) / exponent.powf(0.13);
    if factor <= 0.0 || factor.is_nan() {
        return density_uop(molar_masses, kw);
    }
    molar_masses
        .iter()
        .map(|m| {
            if *m <= 0.0 {
                return 0.0;
            }
            let argument = (m * 1000.0 - SOREIDE_MOLAR_MASS_FLOOR).max(1.0);
            (SOREIDE_INTERCEPT + factor * argument.powf(0.13)).clamp(SOREIDE_LOWER, SOREIDE_UPPER)
        })
        .collect()
}

/// A plus fraction split into carbon-number cuts by Whitson's three-parameter gamma.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if an input is outside its declared range.
/// * [`AzothError::InvalidInput`] if the plus fraction's molar mass is at or below `eta`, so the
///   derived scale is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::{WhitsonDensityModel, whitson_gamma_split};
///
/// let split = whitson_gamma_split(
///     kilograms_per_mole(0.4),
///     kilograms_per_cubic_meter(850.0),
///     0.1,
///     1,
///     80,
///     None,
///     None,
///     WhitsonDensityModel::Uop,
///     false,
/// )?;
/// assert_eq!(split.cut_z.len(), 79);
/// assert!((split.shape - 1.0).abs() < 1e-15);
/// assert!((split.cut_molar_mass[0].value - 0.0969473136183643).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)]
pub fn whitson_gamma_split(
    molar_mass: MolarMass,
    density: MassDensity,
    mole_fraction: f64,
    first_carbon_number: usize,
    last_carbon_number: usize,
    alpha: Option<f64>,
    eta: Option<MolarMass>,
    density_model: WhitsonDensityModel,
    auto_estimate_shape: bool,
) -> Result<WhitsonGammaSplitResult> {
    let spec = &model_gen::WHITSON_GAMMA_SPLIT_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "molar_mass" => Some(molar_mass.value),
            "density" => Some(density.value),
            "mole_fraction" => Some(mole_fraction),
            "first_carbon_number" => Some(first_carbon_number as f64),
            "last_carbon_number" => Some(last_carbon_number as f64),
            "alpha" => alpha,
            "eta" => eta.map(|q| q.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if last_carbon_number <= first_carbon_number {
        return Err(AzothError::invalid_input(
            "last_carbon_number",
            format!(
                "{last_carbon_number} leaves no cut above {first_carbon_number}: the range is \
                 half-open, so at least one cut has to lie between them"
            ),
        ));
    }

    // The class works in g/mol and g/cm3 from here on, the conversion made once.
    let m_plus = molar_mass.value * 1000.0;
    let dens_plus = density.value / 1000.0;
    let eta_gmol = eta.map_or(DEFAULT_ETA_GMOL, |q| q.value * 1000.0);
    // The estimate runs before the scale is derived, so it moves both parameters.
    let shape = if auto_estimate_shape {
        estimate_shape(molar_mass.value, dens_plus)
    } else {
        alpha.unwrap_or(DEFAULT_SHAPE)
    };
    let scale = (m_plus - eta_gmol) / shape;
    if scale <= 0.0 {
        return Err(AzothError::invalid_input(
            "eta",
            format!(
                "the derived gamma scale is {scale}, because the plus fraction's {m_plus} g/mol is \
                 not above eta's {eta_gmol} g/mol"
            ),
        ));
    }

    let count = last_carbon_number - first_carbon_number;
    let mut cut_z = Vec::with_capacity(count);
    let mut cut_molar_mass = Vec::with_capacity(count);
    // The first window opens at `eta`; each later one opens where the one before it closed.
    let mut window_lower = eta_gmol;
    for i in first_carbon_number..last_carbon_number {
        let mut window_upper = window_lower + SCN_INCREMENT;
        if i == last_carbon_number - 1 {
            window_upper = LAST_WINDOW_UPPER;
        }

        let low = p0_p1(window_lower, eta_gmol, shape, scale);
        let high = p0_p1(window_upper, eta_gmol, shape, scale);
        let mut window = high[0] - low[0];
        if window < WINDOW_FLOOR {
            window = WINDOW_FLOOR;
        }
        cut_z.push(window * mole_fraction);

        // The moment's denominator is the **unfloored** difference, so a window the floor
        // rescued above takes the midpoint.
        let denominator = high[0] - low[0];
        let mean = if denominator.abs() < MOMENT_FLOOR {
            0.5 * (window_lower + window_upper)
        } else {
            eta_gmol + shape * scale * (high[1] - low[1]) / denominator
        };
        cut_molar_mass.push(mean / 1000.0);
        window_lower = window_upper;
    }

    // The row is normalised to the plus fraction's own mole fraction, which is what makes the
    // `window * z_plus` above a weight rather than an answer.
    let total: f64 = cut_z.iter().sum();
    if total > 0.0 {
        for value in &mut cut_z {
            *value *= mole_fraction / total;
        }
    }

    let kw = watson_factor(molar_mass.value, dens_plus);
    let densities = match density_model {
        WhitsonDensityModel::Uop => density_uop(&cut_molar_mass, kw),
        WhitsonDensityModel::Soreide => density_soreide(&cut_molar_mass, molar_mass.value, dens_plus, kw),
    };

    Ok(WhitsonGammaSplitResult {
        cut_z,
        cut_molar_mass: cut_molar_mass
            .into_iter()
            .map(kilograms_per_mole)
            .collect(),
        cut_density: densities
            .into_iter()
            .map(|d| kilograms_per_cubic_meter(d * 1000.0))
            .collect(),
        shape,
        minimum_molar_mass: kilograms_per_mole(eta_gmol / 1000.0),
        scale: kilograms_per_mole(scale / 1000.0),
        warnings,
    })
}
