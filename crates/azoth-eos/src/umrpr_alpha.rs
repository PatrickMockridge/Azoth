//! `eos.umrpr_alpha` - the UMR-PR alpha function.
//!
//! ```text
//! alpha = (1 + m*(1 - sqrt(Tr)))**2
//! m     = 0.384401 + 1.52276*omega - 0.213808*omega**2 + 0.034616*omega**3
//!         - 0.001976*omega**4
//! ```
//!
//! Spec: `specs/calcs/eos/umrpr_alpha.toml`, which carries the provenance and why the
//! m-factor is the UMR-PR kappa.

use azoth_core::{Result, apply_checks};

use crate::results::UmrprAlphaResult;
use crate::spec_gen;

/// The UMR-PR alpha function for a pure component.
///
/// `omega` is the acentric factor; `Tr` the reduced temperature. The alpha is
/// cubic-agnostic - the `a` and `b` it scales are the cubic's own.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::umrpr_alpha;
///
/// let r = umrpr_alpha(0.1, 0.7)?;
/// assert!((r.alpha - 1.1822586823466172).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn umrpr_alpha(omega: f64, Tr: f64) -> Result<UmrprAlphaResult> {
    let spec = &spec_gen::UMRPR_ALPHA_SPEC;
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

    let m = 0.384_401 + 1.522_76 * omega - 0.213_808 * omega * omega
        + 0.034_616 * omega * omega * omega
        - 0.001_976 * omega * omega * omega * omega;
    let t = 1.0 + m * (1.0 - Tr.sqrt());
    let alpha = t * t;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(UmrprAlphaResult { alpha, warnings })
}
