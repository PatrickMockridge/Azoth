//! `eos.iapws_henry_law` - the Henry constant of a gas in water.
//!
//! Spec: `specs/calcs/eos/iapws_henry_law.toml`, which records the equation, the 14 rows
//! and why a temperature outside a row's fitted window is a returned number rather than a
//! refusal.
//!
//! Ported from `thermo/component/IapwsHenryLaw`. It is the second of the two Henry arms
//! this crate carries: [`crate::henry`] is the database correlation out of `COMP.csv`, and
//! this is the guideline's. Which one a phase takes is the phase's decision - see
//! `eos.pitzer_phase`, where both are reachable and the gate between them is measured.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::{HenryStatus, IapwsHenryLawResult};
use crate::spec_gen;

/// The temperature the guideline's reduced temperature is taken against.
const WATER_CRITICAL_TEMPERATURE: f64 = 647.096;

/// Water's critical pressure in MPa, the guideline's reference pressure.
const WATER_CRITICAL_PRESSURE_MPA: f64 = 22.064;

/// The guideline's domain: liquid water.
const CORRELATION_MINIMUM_TEMPERATURE: f64 = 273.15;

/// The six coefficients of the Wagner saturation-pressure series.
const VAPOR_PRESSURE_A: [f64; 6] = [
    -7.859_517_83,
    1.844_082_59,
    -11.786_649_7,
    22.680_741_1,
    -15.961_871_9,
    1.801_225_02,
];

/// The exponents the six coefficients are raised to.
const VAPOR_PRESSURE_B: [f64; 6] = [1.0, 1.5, 3.0, 3.5, 4.0, 7.5];

/// Megapascals per bar. The guideline evaluates in MPa and NeqSim returns bar.
const MPA_TO_BAR: f64 = 10.0;

/// Pascals per bar, so that this model returns the unit the rest of the library works in.
const BAR_TO_PA: f64 = 1.0e5;

/// Water's molar mass, in kg/mol.
///
/// `IapwsHenryLaw.WATER_MOLAR_MASS_KG_PER_MOL`, and **the factor that is not applied
/// here**: `ComponentGePitzer` multiplies the constant by it to move the standard state
/// onto the molality scale, which is a phase's decision and not this correlation's.
pub const WATER_MOLAR_MASS_KG_PER_MOL: f64 = 0.018_015_28;

/// The dataset this table is, as NeqSim names it.
pub const DATASET_ID: &str = "iapws-g7-04-water-gases-v1";

/// One row of the guideline's gas table.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GasRow {
    /// The formula the row is filed under, which is the spec's spelling too.
    pub name: &'static str,
    /// The row's first fitted constant.
    pub a: f64,
    /// The row's second fitted constant.
    pub b: f64,
    /// The row's third fitted constant.
    pub c: f64,
    /// The low end of the range the row was fitted over.
    pub minimum_temperature: f64,
    /// The high end of the range the row was fitted over.
    pub maximum_temperature: f64,
    /// The row's reported root-mean-square residual in `ln kH`.
    pub rms_log_residual: f64,
}

/// The 14 gases, and the one place their constants live.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Gas {
    /// Helium.
    He,
    /// Neon.
    Ne,
    /// Argon.
    Ar,
    /// Krypton.
    Kr,
    /// Xenon.
    Xe,
    /// Hydrogen.
    H2,
    /// Nitrogen.
    N2,
    /// Oxygen.
    O2,
    /// Carbon monoxide.
    Co,
    /// Carbon dioxide.
    Co2,
    /// Hydrogen sulphide.
    H2s,
    /// Methane.
    Ch4,
    /// Ethane.
    C2h6,
    /// Sulfur hexafluoride.
    Sf6,
}

impl Gas {
    /// The row's fitted constants and range.
    #[must_use]
    pub fn row(self) -> GasRow {
        match self {
            Self::He => GasRow {
                name: "he",
                a: -3.528_39,
                b: 7.129_83,
                c: 4.477_70,
                minimum_temperature: 273.21,
                maximum_temperature: 553.18,
                rms_log_residual: 0.034_1,
            },
            Self::Ne => GasRow {
                name: "ne",
                a: -3.183_01,
                b: 5.314_48,
                c: 5.437_74,
                minimum_temperature: 273.20,
                maximum_temperature: 543.36,
                rms_log_residual: 0.057_7,
            },
            Self::Ar => GasRow {
                name: "ar",
                a: -8.409_54,
                b: 4.295_87,
                c: 10.527_79,
                minimum_temperature: 273.19,
                maximum_temperature: 568.36,
                rms_log_residual: 0.044_3,
            },
            Self::Kr => GasRow {
                name: "kr",
                a: -8.973_58,
                b: 3.615_08,
                c: 11.299_63,
                minimum_temperature: 273.19,
                maximum_temperature: 525.56,
                rms_log_residual: 0.043_4,
            },
            Self::Xe => GasRow {
                name: "xe",
                a: -14.216_35,
                b: 4.000_41,
                c: 15.609_99,
                minimum_temperature: 273.22,
                maximum_temperature: 574.85,
                rms_log_residual: 0.036_3,
            },
            Self::H2 => GasRow {
                name: "h2",
                a: -4.732_84,
                b: 6.089_54,
                c: 6.060_66,
                minimum_temperature: 273.15,
                maximum_temperature: 636.09,
                rms_log_residual: 0.051_7,
            },
            Self::N2 => GasRow {
                name: "n2",
                a: -9.675_78,
                b: 4.721_62,
                c: 11.705_85,
                minimum_temperature: 278.12,
                maximum_temperature: 636.46,
                rms_log_residual: 0.037_2,
            },
            Self::O2 => GasRow {
                name: "o2",
                a: -9.448_33,
                b: 4.438_22,
                c: 11.420_05,
                minimum_temperature: 274.15,
                maximum_temperature: 616.52,
                rms_log_residual: 0.037_7,
            },
            Self::Co => GasRow {
                name: "co",
                a: -10.528_62,
                b: 5.132_59,
                c: 12.014_21,
                minimum_temperature: 278.15,
                maximum_temperature: 588.67,
                rms_log_residual: 0.003_9,
            },
            Self::Co2 => GasRow {
                name: "co2",
                a: -8.554_45,
                b: 4.011_95,
                c: 9.523_45,
                minimum_temperature: 274.19,
                maximum_temperature: 642.66,
                rms_log_residual: 0.052_8,
            },
            Self::H2s => GasRow {
                name: "h2s",
                a: -4.514_99,
                b: 5.235_38,
                c: 4.421_26,
                minimum_temperature: 273.15,
                maximum_temperature: 533.09,
                rms_log_residual: 0.040_8,
            },
            Self::Ch4 => GasRow {
                name: "ch4",
                a: -10.447_08,
                b: 4.664_91,
                c: 12.129_86,
                minimum_temperature: 275.46,
                maximum_temperature: 633.11,
                rms_log_residual: 0.038_6,
            },
            Self::C2h6 => GasRow {
                name: "c2h6",
                a: -19.675_63,
                b: 4.512_22,
                c: 20.625_67,
                minimum_temperature: 275.44,
                maximum_temperature: 473.46,
                rms_log_residual: 0.025_9,
            },
            Self::Sf6 => GasRow {
                name: "sf6",
                a: -16.561_18,
                b: 2.152_89,
                c: 20.354_40,
                minimum_temperature: 283.14,
                maximum_temperature: 505.55,
                rms_log_residual: 0.050_5,
            },
        }
    }

    /// The spec's spelling.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        self.row().name
    }
}

impl std::str::FromStr for Gas {
    type Err = String;

    fn from_str(value: &str) -> std::result::Result<Self, Self::Err> {
        match value {
            "he" => Ok(Self::He),
            "ne" => Ok(Self::Ne),
            "ar" => Ok(Self::Ar),
            "kr" => Ok(Self::Kr),
            "xe" => Ok(Self::Xe),
            "h2" => Ok(Self::H2),
            "n2" => Ok(Self::N2),
            "o2" => Ok(Self::O2),
            "co" => Ok(Self::Co),
            "co2" => Ok(Self::Co2),
            "h2s" => Ok(Self::H2s),
            "ch4" => Ok(Self::Ch4),
            "c2h6" => Ok(Self::C2h6),
            "sf6" => Ok(Self::Sf6),
            other => Err(format!(
                "unknown gas `{other}`; the guideline's table carries `he`, `ne`, `ar`, \
                 `kr`, `xe`, `h2`, `n2`, `o2`, `co`, `co2`, `h2s`, `ch4`, `c2h6` and `sf6`"
            )),
        }
    }
}

/// The row a component name means, which is `IapwsHenryLaw.findGas`.
///
/// A row answers to its formula and to its name, case-insensitively, and to nothing else:
/// `krypton` is a row here and not a component this crate carries data for, so the
/// resolution is by string and cannot be by lookup. This is the boundary a caller with a
/// databank name crosses, and it is public because that caller is a phase.
#[must_use]
pub fn gas_from_name(name: &str) -> Option<Gas> {
    match name.trim().to_ascii_lowercase().as_str() {
        "he" | "helium" => Some(Gas::He),
        "ne" | "neon" => Some(Gas::Ne),
        "ar" | "argon" => Some(Gas::Ar),
        "kr" | "krypton" => Some(Gas::Kr),
        "xe" | "xenon" => Some(Gas::Xe),
        "h2" | "hydrogen" => Some(Gas::H2),
        "n2" | "nitrogen" => Some(Gas::N2),
        "o2" | "oxygen" => Some(Gas::O2),
        "co" | "carbon monoxide" => Some(Gas::Co),
        "co2" | "carbon dioxide" => Some(Gas::Co2),
        "h2s" | "hydrogen sulfide" | "hydrogen sulphide" => Some(Gas::H2s),
        "ch4" | "methane" => Some(Gas::Ch4),
        "c2h6" | "ethane" => Some(Gas::C2h6),
        "sf6" | "sulfur hexafluoride" | "sulphur hexafluoride" => Some(Gas::Sf6),
        _ => None,
    }
}

/// The row a component name means, as a refusal rather than an `Option`.
///
/// The same lookup as [`gas_from_name`], and the one the boundary should use: a caller that
/// has a name and no row has made an input error, and both kernels raise the same one with
/// the same message.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if no row answers to the name.
pub fn gas_by_name(name: &str) -> Result<Gas> {
    gas_from_name(name).ok_or_else(|| {
        AzothError::invalid_input(
            "gas",
            format!(
                "unknown gas `{name}`; the guideline's table carries `he`, `ne`, `ar`, `kr`, \
                 `xe`, `h2`, `n2`, `o2`, `co`, `co2`, `h2s`, `ch4`, `c2h6` and `sf6`"
            ),
        )
    })
}

/// The Wagner saturation-pressure series, `S(tau)`.
fn vapor_pressure_series(tau: f64) -> f64 {
    let mut sum = 0.0;
    for (a, b) in VAPOR_PRESSURE_A.iter().zip(VAPOR_PRESSURE_B.iter()) {
        sum += a * tau.powf(*b);
    }
    sum
}

/// `kH` in bar, the guideline's own unit, and `d(ln kH)/dT` in `1/K`.
fn evaluate(gas: Gas, temperature: f64) -> (f64, f64) {
    let GasRow { a, b, c, .. } = gas.row();
    let tr = temperature / WATER_CRITICAL_TEMPERATURE;
    let tau = 1.0 - tr;
    let series = vapor_pressure_series(tau);
    let log_pressure_mpa = WATER_CRITICAL_PRESSURE_MPA.ln() + series / tr;
    let log_henry_mpa =
        log_pressure_mpa + a / tr + b * tau.powf(0.355) / tr + c * tr.powf(-0.41) * tau.exp();
    let henry_bar = log_henry_mpa.exp() * MPA_TO_BAR;

    // The logarithmic derivative, `getLnHenryCoefficientTemperatureDerivative`. It is
    // written against `Tr` and divided once at the end, which is how NeqSim writes it.
    let mut series_derivative = 0.0;
    for (a, b) in VAPOR_PRESSURE_A.iter().zip(VAPOR_PRESSURE_B.iter()) {
        series_derivative -= a * b * tau.powf(b - 1.0);
    }
    let mut derivative_by_tr = series_derivative / tr - series / (tr * tr) - a / (tr * tr)
        + b * (-0.355 * tau.powf(-0.645) / tr - tau.powf(0.355) / (tr * tr));
    let final_term = c * tr.powf(-0.41) * tau.exp();
    derivative_by_tr += final_term * (-0.41 / tr - 1.0);

    (henry_bar, derivative_by_tr / WATER_CRITICAL_TEMPERATURE)
}

/// The Henry constant of a gas in water, and the range it was fitted over.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` is outside liquid water, `[273.15, 647.096)` K.
///
/// A temperature inside liquid water but outside the row's own fitted window is
/// **returned**, with `status` = `guideline_extrapolation`: the equation is defined there
/// and the guideline's own extrapolating entry point exists for it. What an extrapolation
/// means is the consumer's decision, and NeqSim's consumers turn it into the insoluble
/// limit.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::{Gas, iapws_henry_law};
///
/// let r = iapws_henry_law(Gas::Ch4, kelvins(298.15))?;
/// assert!((r.henry.value - 3.947_965_646_000_571e9).abs() < 1.0);
/// assert!((r.d_ln_henry_d_t - 0.016_661_113_277_496_545).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn iapws_henry_law(gas: Gas, t: ThermodynamicTemperature) -> Result<IapwsHenryLawResult> {
    let spec = &spec_gen::IAPWS_HENRY_LAW_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The guideline's domain is liquid water, and it is checked here rather than left to
    // the range table: the table's two bounds are this comparison, and a caller reading
    // the code should see the domain stated once.
    if !t.value.is_finite()
        || t.value < CORRELATION_MINIMUM_TEMPERATURE
        || t.value >= WATER_CRITICAL_TEMPERATURE
    {
        return Err(AzothError::out_of_range(
            "T",
            t.value,
            format!(
                "the guideline is stated for liquid water, [{CORRELATION_MINIMUM_TEMPERATURE}, \
                 {WATER_CRITICAL_TEMPERATURE}) K, and the reference state ceases to exist above \
                 it"
            ),
        ));
    }

    let row = gas.row();
    let (henry_bar, d_ln_henry_d_t) = evaluate(gas, t.value);
    let status = if t.value >= row.minimum_temperature && t.value <= row.maximum_temperature {
        HenryStatus::WithinFittedRange
    } else {
        HenryStatus::GuidelineExtrapolation
    };

    // `ln_henry` is the logarithm of the value this model returns and not of NeqSim's bar
    // figure: the two differ by `ln 1e5`, and a result whose parts are in two scales is a
    // trap for whoever reads one against the other. The derivative is scale-invariant, so
    // it is unaffected.
    let henry_pa = henry_bar * BAR_TO_PA;

    Ok(IapwsHenryLawResult {
        henry: pascals(henry_pa),
        ln_henry: henry_pa.ln(),
        d_ln_henry_d_t,
        status,
        rms_log_residual: row.rms_log_residual,
        warnings,
    })
}
