//! `eos.twucoon_param_alpha` - the Twu-Coon parameter alpha function.
//!
//! ```text
//! alpha = Tr**(c*(b-1)) * exp(a*(1 - Tr**(b*c)))
//! ```
//!
//! Spec: `specs/calcs/eos/twucoon_param_alpha.toml`, which carries the provenance and why
//! the three parameters are the caller's rather than derived.

use azoth_core::{Result, apply_checks};

use crate::results::TwucoonParamAlphaResult;
use crate::spec_gen;

/// The Twu-Coon parameter alpha function for a pure component.
///
/// `a`, `b` and `c` are the fitted Twu-Coon parameters; `Tr` the reduced temperature. The
/// alpha is cubic-agnostic - the `a` and `b` it scales are the cubic's own.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes fractional
///   powers of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::twucoon_param_alpha;
///
/// let r = twucoon_param_alpha(0.1, 0.5, 2.0, 0.7)?;
/// assert!((r.alpha - 1.4720779056478814).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn twucoon_param_alpha(a: f64, b: f64, c: f64, Tr: f64) -> Result<TwucoonParamAlphaResult> {
    let spec = &spec_gen::TWUCOON_PARAM_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "a" => Some(a),
            "b" => Some(b),
            "c" => Some(c),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let alpha = Tr.powf(c * (b - 1.0)) * (a * (1.0 - Tr.powf(b * c))).exp();

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(TwucoonParamAlphaResult { alpha, warnings })
}
