//! `eos.srk_kappa` - the Soave-Redlich-Kwong attraction-parameter coefficient.
//!
//! ```text
//! kappa = 0.48 + 1.574*omega - 0.176*omega**2
//! ```
//!
//! Spec: `specs/calcs/eos/srk_kappa.toml`, which carries the citation and why a
//! negative `kappa` is returned with a warning rather than refused.
//!
//! Soave's coefficient is the whole temperature dependence of the SRK attraction term,
//! in one number, and it is the same *form* as Peng-Robinson's with three different
//! constants - which is why the alpha functions themselves share [`crate::alpha_term`].

use azoth_core::{Result, apply_checks};

use crate::results::SrkKappaResult;
use crate::spec_gen;

/// The Soave-Redlich-Kwong alpha-function coefficient for a pure component.
///
/// `omega` is the Pitzer acentric factor, supplied by the caller. The returned
/// [`SrkKappaResult`] carries warnings but no failure path beyond a malformed spec: the
/// polynomial is defined for every real `omega`.
///
/// # Example
/// ```
/// use azoth_eos::srk_kappa;
///
/// // Propane-like.
/// let r = srk_kappa(0.152)?;
/// assert!((r.kappa - 0.715181696).abs() < 1e-12);
/// assert!(r.warnings.is_empty());
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn srk_kappa(omega: f64) -> Result<SrkKappaResult> {
    let spec = &spec_gen::SRK_KAPPA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            _ => None,
        },
        &mut warnings,
    )?;

    let kappa = 0.48 + 1.574 * omega - 0.176 * omega * omega;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "kappa").then_some(kappa),
        &mut warnings,
    )?;

    Ok(SrkKappaResult { kappa, warnings })
}
