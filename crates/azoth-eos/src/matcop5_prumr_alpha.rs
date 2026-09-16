//! `eos.matcop5_prumr_alpha` - the five-parameter Mathias-Copeman alpha function.
//!
//! ```text
//! alpha = (1 + mc1*u + mc2*u**2 + mc3*u**3 + mc4*u**4 + mc5*u**5)**2,  u = 1 - sqrt(Tr)
//!         falling back to (1 + m*(1 - sqrt(Tr)))**2 when every coefficient is zero
//! m     = 0.37464 + 1.54226*omega - 0.26992*omega**2
//! ```
//!
//! Spec: `specs/calcs/eos/matcop5_prumr_alpha.toml`, which carries the provenance and why
//! the five coefficients are the caller's rather than derived.

use azoth_core::{Result, apply_checks};

use crate::results::Matcop5PrumrAlphaResult;
use crate::spec_gen;

/// The five-parameter Mathias-Copeman alpha function for a pure component.
///
/// `mc1` .. `mc5` are the fitted Mathias-Copeman coefficients; `omega` the acentric
/// factor, used only for the fallback; `Tr` the reduced temperature.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::matcop5_prumr_alpha;
///
/// let r = matcop5_prumr_alpha(0.1, 0.5, 0.2, -0.1, 0.05, 0.01, 0.7)?;
/// assert!((r.alpha - 1.1807146411895904).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn matcop5_prumr_alpha(
    omega: f64,
    mc1: f64,
    mc2: f64,
    mc3: f64,
    mc4: f64,
    mc5: f64,
    Tr: f64,
) -> Result<Matcop5PrumrAlphaResult> {
    let spec = &spec_gen::MATCOP5_PRUMR_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "mc1" => Some(mc1),
            "mc2" => Some(mc2),
            "mc3" => Some(mc3),
            "mc4" => Some(mc4),
            "mc5" => Some(mc5),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let coeffs = [mc1, mc2, mc3, mc4, mc5];
    let all_zero = coeffs.iter().all(|c| c.abs() < 1e-20);

    let alpha = if all_zero {
        let m = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;
        let t = 1.0 + m * (1.0 - Tr.sqrt());
        t * t
    } else {
        let u = 1.0 - Tr.sqrt();
        let s = 1.0
            + mc1 * u
            + mc2 * u * u
            + mc3 * u * u * u
            + mc4 * u * u * u * u
            + mc5 * u * u * u * u * u;
        s * s
    };

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(Matcop5PrumrAlphaResult { alpha, warnings })
}
