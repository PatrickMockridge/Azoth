//! `eos.twucoon_alpha` - the Twu-Coon alpha function.
//!
//! ```text
//! alpha = Tr**a*exp(b*(1 - Tr**c)) + omega*(Tr**d*exp(e*(1 - Tr**f)) - Tr**a*exp(b*(1 - Tr**c)))
//! ```
//!
//! Spec: `specs/calcs/eos/twucoon_alpha.toml`, which carries the six fitted constants
//! and the provenance.

use azoth_core::{Result, apply_checks};

use crate::results::TwucoonAlphaResult;
use crate::spec_gen;

/// The six Twu-Coon constants, in the order the equation reads them.
const A: f64 = -0.201_158;
const B: f64 = 0.141_599;
const C: f64 = 2.295_28;
const D: f64 = -0.660_145;
const E: f64 = 0.500_315;
const F: f64 = 2.631_65;

/// The Twu-Coon alpha function for a pure component.
///
/// `omega` is the acentric factor; `Tr` the reduced temperature. The alpha is
/// cubic-agnostic - the `a` and `b` it scales are the cubic's own.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes fractional
///   powers of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::twucoon_alpha;
///
/// let r = twucoon_alpha(0.152, 0.7)?;
/// assert!((r.alpha - 1.2469731129671788).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn twucoon_alpha(omega: f64, Tr: f64) -> Result<TwucoonAlphaResult> {
    let spec = &spec_gen::TWUCOON_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let low = Tr.powf(A) * (B * (1.0 - Tr.powf(C))).exp();
    let high = Tr.powf(D) * (E * (1.0 - Tr.powf(F))).exp();
    let alpha = low + omega * (high - low);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(TwucoonAlphaResult { alpha, warnings })
}
