//! `eos.prsv_kappa` - the Peng-Robinson-Stryjek-Vera alpha-function coefficient.
//!
//! ```text
//! kappa = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3
//!       + kappa1*(1 + Tr**0.5)*(0.7 - Tr)
//! ```
//!
//! Spec: `specs/calcs/eos/prsv_kappa.toml`, which carries the provenance, the
//! parameter `kappa1` this library does not ship, and why no bound asserts the fitted
//! range.
//!
//! PRSV changes the temperature dependence of the attraction coefficient and nothing
//! else, so [`crate::pr_alpha_ab`] serves it unchanged. Unlike Peng-Robinson's
//! coefficient this one varies with `Tr`, so it has to be recomputed at each
//! temperature rather than passed as a fixed number.

use azoth_core::{Result, apply_checks};

use crate::results::PrsvKappaResult;
use crate::spec_gen;

/// The PRSV alpha-function coefficient for a pure component.
///
/// `omega` is the acentric factor and `kappa1` the pure-compound parameter PRSV
/// adds; both are the caller's, because this library ships no component data.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the temperature factor
///   takes a square root, and a reduced temperature at or below zero is not a state.
///
/// # Example
/// ```
/// use azoth_eos::prsv_kappa;
///
/// // Tr = 0.7 is the anchor: the kappa1 term vanishes there whatever kappa1 is.
/// let at_anchor = prsv_kappa(0.152, 0.7, 0.05)?;
/// assert!((at_anchor.kappa - 0.6014406094290431).abs() < 1e-15);
///
/// let elsewhere = prsv_kappa(0.152, 0.8, 0.05)?;
/// assert!((elsewhere.kappa - 0.5919684734740435).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn prsv_kappa(omega: f64, Tr: f64, kappa1: f64) -> Result<PrsvKappaResult> {
    let spec = &spec_gen::PRSV_KAPPA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            "Tr" => Some(Tr),
            "kappa1" => Some(kappa1),
            _ => None,
        },
        &mut warnings,
    )?;

    // Written in the equation's own order, with `omega * omega` for the square
    // rather than `powi(2)`. The two agree bit-for-bit for these inputs - measured
    // - but the explicit form is the one the Python side can match without relying
    // on how each language lowers `**`.
    let acentric_only = 0.378893 + 1.4897153 * omega - 0.17131848 * (omega * omega)
        + 0.0196554 * (omega * omega * omega);
    // `0.7 - Tr` is zero at Tr = 0.7 whatever `kappa1` is, which is the reduced
    // temperature the acentric factor is defined at - the correlation's anchor.
    let kappa = acentric_only + kappa1 * (1.0 + Tr.sqrt()) * (0.7 - Tr);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "kappa").then_some(kappa),
        &mut warnings,
    )?;

    Ok(PrsvKappaResult { kappa, warnings })
}
