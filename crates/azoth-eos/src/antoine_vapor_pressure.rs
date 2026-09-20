//! `eos.antoine_vapor_pressure` - NeqSim's pure-component vapour-pressure correlation.
//!
//! NeqSim's `getAntoineVaporPressure` dispatches on a string label and evaluates one of
//! five correlations. The label vocabulary in the vendored data is messy - `exp` and
//! `log` are one formula under two names, `loglog`/`log10` have no branch and fall
//! through to Wagner - so the label is cleaned onto [`AntoineForm`] by
//! [`form_from_type`] before the correlation is evaluated. That function also takes the
//! fifth coefficient, because for twenty components the label names the wrong form.
//!
//! Spec: `specs/calcs/eos/antoine_vapor_pressure.toml`, which records the four formulas,
//! the cleaning rule, and the dead DIPPR-101 branch.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::AntoineVaporPressureResult;
use crate::spec_gen;

/// The five vapour-pressure correlations NeqSim's `getAntoineVaporPressure` evaluates.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum AntoineForm {
    /// `exp(A + B/T + C ln T + D T**E)` in pascals - the DIPPR-101 form.
    ///
    /// Selected by a non-zero `E` rather than by the label: the rows carrying one are
    /// labelled `log`, which names a different correlation. See [`form_from_type`].
    Dippr101,
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
            AntoineForm::Dippr101 => "dippr101",
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
            "dippr101" => Ok(AntoineForm::Dippr101),
            "pow10" => Ok(AntoineForm::Pow10),
            "pow10kpa" => Ok(AntoineForm::Pow10Kpa),
            "exp" => Ok(AntoineForm::Exp),
            "wagner" => Ok(AntoineForm::Wagner),
            other => Err(format!(
                "unknown Antoine form `{other}`; expected `dippr101`, `pow10`, `pow10kpa`, \
                 `exp` or `wagner`"
            )),
        }
    }
}

/// Map NeqSim's raw `AntoineVapPresLiqType` label onto [`AntoineForm`].
///
/// **`E` is part of the question, not only the label.** Twenty rows in NeqSim's
/// `COMP.csv` carry DIPPR-101 coefficients - `exp(A + B/T + C ln T + D T**E)` - under
/// the label `log`, which names the two-term exponential instead. Reading the label
/// alone therefore picks the wrong correlation for them, by thirty to ninety orders of
/// magnitude; `i-pentane` came out at 8.3e38 bar.
///
/// The rule is NeqSim's own, `Component.usesDipprVaporPressureCorrelation`: a non-zero
/// `E` decides, except that `pow10` and `pow10KPa` keep precedence because those
/// coefficients are log10-based and would not survive the exponential form.
///
/// `exp` and `log` are otherwise one formula under two names, and `loglog`/`log10`
/// have no branch in NeqSim's dispatch, so they fall through to Wagner - a defect this
/// reproduces rather than silently repairs.
///
/// **`none` is not a form, and it is `None` here rather than a fall-through.** It is the
/// marker upstream added in `83b64e5` (PR #3775) for a row whose correlation is
/// *unavailable*: 313 of its 389 rows now carry it, with all five coefficients zero. Sent
/// to Wagner those coefficients give `exp(0) * Pc = Pc`, so an unavailable correlation
/// would come back as the component's **critical pressure** - measured, `1.82e6 Pa` for
/// `nc12`, whose vapour pressure there is about `42 Pa`. A plausible number wrong by four
/// orders of magnitude is the failure this library exists to make impossible, so the
/// marker is returned as an absence and the caller refuses.
#[must_use]
pub fn form_from_type(label: &str, e: f64) -> Option<AntoineForm> {
    if label == "none" {
        return None;
    }
    if e.abs() > 1e-12 && label != "pow10" && label != "pow10KPa" {
        return Some(AntoineForm::Dippr101);
    }
    Some(match label {
        "pow10" => AntoineForm::Pow10,
        "pow10KPa" => AntoineForm::Pow10Kpa,
        "exp" | "log" => AntoineForm::Exp,
        _ => AntoineForm::Wagner,
    })
}

/// The pure-component vapour pressure at a temperature, from NeqSim's correlation.
///
/// `A`-`E` are the raw `ANTOINEA`-`ANTOINEE` and `Tc`/`Pc` the critical constants, all
/// caller-supplied. `E` is the DIPPR-101 exponent, used by [`AntoineForm::Dippr101`]
/// and ignored by the others; `Tc` and `Pc` are used only by the Wagner form.
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
        // NeqSim returns this one in pascals already - `exp(...) / 100000` bar - so
        // unlike the other bar-returning forms there is no `1e5` factor here.
        AntoineForm::Dippr101 => (A + B / t + C * t.ln() + D * t.powf(E)).exp(),
        AntoineForm::Pow10 => 1e5 * (10.0_f64).powf(A - B / (t + C - 273.15)),
        AntoineForm::Pow10Kpa => (10.0_f64).powf(A - B / (t + C)),
        AntoineForm::Exp => 1e5 * (A - B / (t + C)).exp(),
        AntoineForm::Wagner => {
            // **The form is defined on `0 < T <= Tc`, and the kernel says so.** `x = 1 -
            // T/Tc` is negative above the critical temperature, and `x.powf(1.5)` is then
            // `NaN` - which flows through the `exp` into `p_sat` and is returned as a
            // *result*. A caller checking `p_sat > 0` sees `false` and has nothing to
            // attribute it to, so the out-of-domain state is refused rather than computed.
            // Above `Tc` there is no saturation pressure to report, which is why this is a
            // refusal and not a clamp.
            let x = 1.0 - t / Tc.value;
            if x < 0.0 {
                return Err(AzothError::out_of_range(
                    "T",
                    t,
                    format!(
                        "the Wagner form is defined up to the critical temperature, and this                          state is above it: `1 - T/Tc` is {x}, whose 1.5 power is not a real                          number. There is no saturation pressure above `Tc` to return, so the                          state is refused rather than reported as `NaN`"
                    ),
                ));
            }
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
