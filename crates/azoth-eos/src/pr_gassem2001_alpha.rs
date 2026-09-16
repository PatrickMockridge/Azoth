//! `eos.pr_gassem2001_alpha` - the Gassem (2001) alpha function.
//!
//! ```text
//! alpha = exp((A + B*Tr)*(1 - Tr**(C + D*omega + E*omega**2)))
//! ```
//!
//! Spec: `specs/calcs/eos/pr_gassem2001_alpha.toml`, which carries the five fitted
//! constants and the provenance.

use azoth_core::{Result, apply_checks};

use crate::results::PrGassem2001AlphaResult;
use crate::spec_gen;

/// The five Gassem (2001) constants.
const A: f64 = 2.0;
const B: f64 = 0.836;
const C: f64 = 0.134;
const D: f64 = 0.508;
const E: f64 = -0.0467;

/// The Gassem (2001) alpha function for a pure component.
///
/// `omega` is the acentric factor; `Tr` the reduced temperature. The alpha is
/// cubic-agnostic - the `a` and `b` it scales are the cubic's own.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a fractional
///   power of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::pr_gassem2001_alpha;
///
/// let r = pr_gassem2001_alpha(0.152, 0.7)?;
/// assert!((r.alpha - 1.2052404595821262).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_gassem2001_alpha(omega: f64, Tr: f64) -> Result<PrGassem2001AlphaResult> {
    let spec = &spec_gen::PR_GASSEM2001_ALPHA_SPEC;
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

    let exponent = C + D * omega + E * omega * omega;
    let alpha = ((A + B * Tr) * (1.0 - Tr.powf(exponent))).exp();

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(PrGassem2001AlphaResult { alpha, warnings })
}
