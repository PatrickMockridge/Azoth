//! `eos.schwartzentruber_alpha` - the Schwartzentruber-Renon alpha function.
//!
//! ```text
//! alpha = (1 + m*(1 - sqrt(Tr)) - p1*(1 - Tr)*(1 + p2*Tr + p3*Tr**2))**2
//! m     = 0.48508 + 1.55191*omega - 0.15613*omega**2
//! ```
//!
//! Spec: `specs/calcs/eos/schwartzentruber_alpha.toml`, which carries the provenance and
//! why the three parameters are the caller's rather than derived.

use azoth_core::{Result, apply_checks};

use crate::results::SchwartzentruberAlphaResult;
use crate::spec_gen;

/// The Schwartzentruber-Renon alpha function for a pure component.
///
/// `p1`, `p2` and `p3` are the fitted Schwartzentruber-Renon parameters; `omega` the
/// acentric factor; `Tr` the reduced temperature.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::schwartzentruber_alpha;
///
/// let r = schwartzentruber_alpha(0.1, 0.3, 0.2, 0.1, 0.7)?;
/// assert!((r.alpha - 0.9946408503265206).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn schwartzentruber_alpha(
    omega: f64,
    p1: f64,
    p2: f64,
    p3: f64,
    Tr: f64,
) -> Result<SchwartzentruberAlphaResult> {
    let spec = &spec_gen::SCHWARTZENTRUBER_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "p1" => Some(p1),
            "p2" => Some(p2),
            "p3" => Some(p3),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let m = 0.48508 + 1.55191 * omega - 0.15613 * omega * omega;
    let inner = 1.0 + m * (1.0 - Tr.sqrt()) - p1 * (1.0 - Tr) * (1.0 + p2 * Tr + p3 * Tr * Tr);
    let alpha = inner * inner;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(SchwartzentruberAlphaResult { alpha, warnings })
}
