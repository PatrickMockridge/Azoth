//! `eos.antoine_vapor_pressure` - NeqSim's pure-component vapour-pressure correlation.
//!
//! NeqSim's `getAntoineVaporPressure` dispatches on a string label and evaluates one of
//! four correlations. The label vocabulary in the vendored data is messy - `exp` and
//! `log` are one formula under two names, and `loglog`/`log10` have no branch and fall
//! through to Wagner - so the label is cleaned onto [`AntoineForm`] by
//! [`form_from_type`] before the correlation is evaluated.
//!
//! Spec: `specs/calcs/eos/antoine_vapor_pressure.toml`, which records the four formulas,
//! the cleaning rule, and the dead DIPPR-101 branch.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};

use crate::results::AntoineVaporPressureResult;
use crate::spec_gen;

/// The four vapour-pressure correlations NeqSim's `getAntoineVaporPressure` evaluates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AntoineForm {
    /// `1e5 * 10**(A - B/(T + C - 273.15))`.
    Pow10,
    /// `10**(A - B/(T + C))`.
    Pow10Kpa,
    /// `1e5 * exp(A - B/(T + C))` - NeqSim's `exp` and `log` labels share it.
    Exp,
    /// `exp((A x + B x**1.5 + C x**3 + D x**6)/(1-x)) * Pc`, `x = 1 - T/Tc`.
    Wagner,
}

impl AntoineForm {
    /// The short name that crosses the Python boundary.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            AntoineForm::Pow10 => "pow10",
            AntoineForm::Pow10Kpa => "pow10kpa",
            AntoineForm::Exp => "exp",
            AntoineForm::Wagner => "wagner",
        }
    }
}

impl std::str::FromStr for AntoineForm {
    type Err = String;

    fn from_str(s: &str) -> std::result::Result<Self, Self::Err> {
        match s {
            "pow10" => Ok(AntoineForm::Pow10),
            "pow10kpa" => Ok(AntoineForm::Pow10Kpa),
            "exp" => Ok(AntoineForm::Exp),
            "wagner" => Ok(AntoineForm::Wagner),
            other => Err(format!(
                "unknown Antoine form `{other}`; expected `pow10`, `pow10kpa`, `exp` or `wagner`"
            )),
        }
    }
}

/// Map NeqSim's raw `AntoineVapPresLiqType` label onto [`AntoineForm`].
///
/// `exp` and `log` are one formula under two names, and `loglog`/`log10` have no
/// branch in NeqSim's dispatch, so they fall through to Wagner - a defect this
/// reproduces rather than silently repairs.
#[must_use]
pub fn form_from_type(label: &str) -> AntoineForm {
    match label {
        "pow10" => AntoineForm::Pow10,
        "pow10KPa" => AntoineForm::Pow10Kpa,
        "exp" | "log" => AntoineForm::Exp,
        _ => AntoineForm::Wagner,
    }
}

/// The pure-component vapour pressure at a temperature, from NeqSim's correlation.
///
/// `A`-`E` are the raw `ANTOINEA`-`ANTOINEE` and `Tc`/`Pc` the critical constants, all
/// caller-supplied. `E` is the dead DIPPR-101 exponent and is unused; `Tc` and `Pc`
/// are used only by the Wagner form.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `T`, `Tc` or `Pc` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::antoine_vapor_pressure;
/// use azoth_eos::AntoineForm;
///
/// let r = antoine_vapor_pressure(5.23243, 891.0098, 332.0975, 0.0, 0.0,
///     AntoineForm::Pow10, kelvins(190.56), pascals(4_599_000.0), kelvins(300.0))?;
/// assert!((r.p_sat.value - 56_252_981.19539252).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn antoine_vapor_pressure(
    A: f64,
    B: f64,
    C: f64,
    D: f64,
    E: f64,
    form: AntoineForm,
    Tc: ThermodynamicTemperature,
    Pc: Pressure,
    T: ThermodynamicTemperature,
) -> Result<AntoineVaporPressureResult> {
    let spec = &spec_gen::ANTOINE_VAPOR_PRESSURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "A" => Some(A),
            "B" => Some(B),
            "C" => Some(C),
            "D" => Some(D),
            "E" => Some(E),
            "Tc" => Some(Tc.value),
            "Pc" => Some(Pc.value),
            "T" => Some(T.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let t = T.value;
    let p_sat = match form {
        AntoineForm::Pow10 => 1e5 * (10.0_f64).powf(A - B / (t + C - 273.15)),
        AntoineForm::Pow10Kpa => (10.0_f64).powf(A - B / (t + C)),
        AntoineForm::Exp => 1e5 * (A - B / (t + C)).exp(),
        AntoineForm::Wagner => {
            let x = 1.0 - t / Tc.value;
            ((A * x + B * x.powf(1.5) + C * x.powi(3) + D * x.powi(6)) / (1.0 - x)).exp() * Pc.value
        }
    };

    apply_checks(
        spec.derived_checks(),
        |name| (name == "p_sat").then_some(p_sat),
        &mut warnings,
    )?;

    Ok(AntoineVaporPressureResult {
        p_sat: pascals(p_sat),
        warnings,
    })
}
