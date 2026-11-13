//! `process.probit` - the dose-response curve as a registered id.
//!
//! Spec: `specs/models/process/probit.toml`. The arithmetic is
//! `process/safety/dispersion/ProbitModel.java` at the pin: the probit, the normal CDF over
//! the class's Abramowitz-Stegun `erf`, and the five factory triples.

use azoth_core::{AzothError, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::warnings as wire_warnings;
use crate::model_gen;

/// Result of `process.probit`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct ProbitResult {
    /// `Y = a + b n ln(dose)`; absent when no dose was stated.
    pub probit: Option<f64>,
    /// The probability at the stated or computed probit.
    pub probability: f64,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

/// Abramowitz & Stegun 7.1.26, the class's approximation.
fn erf(x: f64) -> f64 {
    let sign = if x < 0.0 { -1.0 } else { 1.0 };
    let x = x.abs();
    let t = 1.0 / (1.0 + 0.327_591_1 * x);
    let y = 1.0
        - (((((1.061_405_429 * t - 1.453_152_027) * t) + 1.421_413_741) * t - 0.284_496_736) * t
            + 0.254_829_592)
            * t
            * (-x * x).exp();
    sign * y
}

/// The probit `Y` of a positive dose.
fn probit_of(a: f64, b: f64, n: f64, dose: f64) -> f64 {
    a + b * n * dose.ln()
}

/// The standard normal CDF over [`erf`].
fn probability_of(probit: f64) -> f64 {
    0.5 * (1.0 + erf((probit - 5.0) / 2.0f64.sqrt()))
}

/// The probit of a hazard dose, and the probability a probit value implies.
///
/// # Errors
/// `InvalidInput` when neither `dose` nor `probit_value` is stated, or when the stated dose
/// is not positive.
pub fn probit(
    a: f64,
    b: f64,
    n: f64,
    dose: Option<f64>,
    probit_value: Option<f64>,
) -> Result<ProbitResult> {
    let spec = &model_gen::PROBIT_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "dose" => dose,
            _ => None,
        },
        &mut warnings,
    )?;
    if dose.is_none() && probit_value.is_none() {
        return Err(AzothError::invalid_input(
            "dose",
            "a probit needs a dose or a probit value, and neither was stated",
        ));
    }
    let probit = dose.map(|value| probit_of(a, b, n, value));
    let evaluated = probit_value.or(probit).unwrap_or(f64::NAN);
    Ok(ProbitResult {
        probit,
        probability: probability_of(evaluated),
        warnings,
    })
}
