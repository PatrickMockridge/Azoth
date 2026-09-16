//! `eos.mollerup_alpha` - the Mollerup alpha function.
//!
//! ```text
//! alpha = 1 + p1*(1/Tr - 1) + p2*Tr*ln(Tr) + p3*(Tr - 1)
//! ```
//!
//! Spec: `specs/calcs/eos/mollerup_alpha.toml`, which carries the provenance and why the
//! three parameters are the caller's rather than derived.

use azoth_core::{Result, apply_checks};

use crate::results::MollerupAlphaResult;
use crate::spec_gen;

/// The Mollerup alpha function for a pure component.
///
/// `p1`, `p2` and `p3` are the fitted Mollerup parameters; `Tr` the reduced
/// temperature. The alpha is cubic-agnostic - the `a` and `b` it scales are the cubic's
/// own.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes `1/Tr` and
///   `ln(Tr)`.
///
/// # Example
/// ```
/// use azoth_eos::mollerup_alpha;
///
/// let r = mollerup_alpha(0.1, 0.05, 0.02, 0.7)?;
/// assert!((r.alpha - 1.0243735198192874).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn mollerup_alpha(p1: f64, p2: f64, p3: f64, Tr: f64) -> Result<MollerupAlphaResult> {
    let spec = &spec_gen::MOLLERUP_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "p1" => Some(p1),
            "p2" => Some(p2),
            "p3" => Some(p3),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let alpha = 1.0 + p1 * (1.0 / Tr - 1.0) + p2 * Tr * Tr.ln() + p3 * (Tr - 1.0);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(MollerupAlphaResult { alpha, warnings })
}
