//! `eos.matcop_prumr_alpha` - the Mathias-Copeman alpha function with the UMR-PR
//! fallback.
//!
//! ```text
//! alpha = (1 + mc1*(1 - sqrt(Tr)) + mc2*(1 - sqrt(Tr))**2 + mc3*(1 - sqrt(Tr))**3)**2
//!         falling back to (1 + m*(1 - sqrt(Tr)))**2 when mc1 < 1e-20
//! m     = 0.384401 + 1.52276*omega - 0.213808*omega**2 + 0.034616*omega**3
//!         - 0.001976*omega**4
//! ```
//!
//! Spec: `specs/calcs/eos/matcop_prumr_alpha.toml`, which carries the provenance and why
//! the three coefficients are the caller's rather than derived.

use azoth_core::{Result, apply_checks};

use crate::results::MatcopPrumrAlphaResult;
use crate::spec_gen;

/// The Mathias-Copeman alpha function for a pure component, with the standard
/// Peng-Robinson alpha when the first coefficient is essentially zero.
///
/// `mc1`, `mc2` and `mc3` are the fitted Mathias-Copeman coefficients; `omega` the
/// acentric factor, used only for the fallback; `Tr` the reduced temperature.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::matcop_prumr_alpha;
///
/// let r = matcop_prumr_alpha(0.1, 0.5, 0.2, -0.1, 0.7)?;
/// assert!((r.alpha - 1.180634768967036).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn matcop_prumr_alpha(
    omega: f64,
    mc1: f64,
    mc2: f64,
    mc3: f64,
    Tr: f64,
) -> Result<MatcopPrumrAlphaResult> {
    let spec = &spec_gen::MATCOP_PRUMR_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "mc1" => Some(mc1),
            "mc2" => Some(mc2),
            "mc3" => Some(mc3),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let alpha = if mc1 < 1e-20 {
        let m = 0.384_401 + 1.522_76 * omega - 0.213_808 * omega * omega
            + 0.034_616 * omega * omega * omega
            - 0.001_976 * omega * omega * omega * omega;
        let t = 1.0 + m * (1.0 - Tr.sqrt());
        t * t
    } else {
        let root = 1.0 - Tr.sqrt();
        let t = 1.0 + mc1 * root + mc2 * root * root + mc3 * root * root * root;
        t * t
    };

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(MatcopPrumrAlphaResult { alpha, warnings })
}
