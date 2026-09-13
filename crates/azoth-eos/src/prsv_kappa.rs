//! `eos.prsv_kappa` - the Peng-Robinson-Stryjek-Vera alpha-function coefficient.
//!
//! ```text
//! kappa = 0.378893 + 1.4897153*omega - 0.17131848*omega**2 + 0.0196554*omega**3
//!       + kappa1*(1 + Tr**0.5)*(0.7 - Tr)
//! ```
//!
//! Stryjek, R.; Vera, J. H. (1986). "PRSV: An improved Peng-Robinson equation of
//! state for pure compounds and mixtures." Can. J. Chem. Eng. 64(2), 323-333.
//! DOI 10.1002/cjce.5450640224
//!
//! Spec: `specs/calcs/eos/prsv_kappa.yaml`
//!
//! # What PRSV changes, and what it does not
//!
//! It changes the temperature dependence of the attraction coefficient and nothing
//! else. The alpha function it feeds is Peng-Robinson's, unchanged, which is why
//! [`crate::pr_alpha_ab`] serves both and why this calc is a coefficient rather
//! than a second equation of state.
//!
//! Peng-Robinson's coefficient is a constant per substance. This one is not: the
//! `(1 + sqrt(Tr))*(0.7 - Tr)` term varies with temperature, and `kappa1` is fitted
//! per component so that the variation matches that component's vapour pressure. It
//! therefore cannot be passed to [`crate::pr_alpha_ab`] as a fixed number the way
//! Peng-Robinson's can - it has to be recomputed at each temperature.
//!
//! # The anchor at Tr = 0.7
//!
//! The temperature factor is exactly zero at `Tr = 0.7`, for any `kappa1`, because
//! `0.7 - Tr` is. So the coefficient collapses to its acentric-only part there -
//! which is the temperature at which the acentric factor is defined. PRSV is pinned
//! to the acentric-only fit at Tr = 0.7 and departs from it elsewhere by whatever
//! `kappa1` says. The spec's `at_tr_0_7_the_kappa1_term_vanishes` case asserts it
//! with a deliberately non-zero `kappa1`, so it cannot pass by the term being
//! absent.
//!
//! # The parameter this crate does not ship
//!
//! `kappa1` is fitted per substance, and a table of fitted parameters is the
//! databank this library deliberately has none of. It is an input. A caller who
//! expected PRSV's published accuracy for free will not get it, and the spec says
//! so rather than leaving it to be discovered.

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
    let kappa = acentric_only + kappa1 * (1.0 + Tr.sqrt()) * (0.7 - Tr);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "kappa").then_some(kappa),
        &mut warnings,
    )?;

    Ok(PrsvKappaResult { kappa, warnings })
}
