//! `eos.pr_kappa` - the Peng-Robinson attraction-parameter coefficient.
//!
//! ```text
//! kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2
//! ```
//!
//! Peng, D. Y.; Robinson, D. B. (1976). "A New Two-Constant Equation of State."
//! Ind. Eng. Chem. Fundam. 15(1), 59-64. DOI 10.1021/i160057a011
//!
//! Spec: `specs/calcs/eos/pr_kappa.yaml`
//!
//! # Where this sits in the equation of state
//!
//! `kappa` is the whole temperature dependence of the Peng-Robinson attraction
//! term, in one number:
//!
//! ```text
//! alpha(T) = (1 + kappa*(1 - sqrt(Tr)))**2
//! a(T)     = 0.45724 * R**2 * Tc**2 / Pc * alpha(T)
//! ```
//!
//! It is a property of the substance alone - no temperature, no pressure - which is
//! what makes it worth a calculation of its own rather than a line inside one. A
//! transposed digit in 0.26992 is not caught by anything downstream: the Z factor,
//! the fugacity coefficient and the phase split it feeds all still converge and all
//! still pass their own consistency checks, and all are slightly wrong.
//!
//! # The sign of kappa
//!
//! `kappa` is negative for `omega < -0.23338349942403008`, which the polynomial's
//! quadratic term makes reachable - helium is at -0.385. A negative coefficient
//! makes `alpha` *grow* with temperature, which is not a statement about any
//! fluid. It is returned rather than refused, carrying a warning: the arithmetic
//! is well defined and a caller inspecting the limit deliberately should not be
//! stopped. See the spec's `valid_range` rationale.

use azoth_core::{Result, apply_checks};

use crate::results::PrKappaResult;
use crate::spec_gen;

/// The Peng-Robinson alpha-function coefficient for a pure component.
///
/// `omega` is the Pitzer acentric factor - dimensionless, and supplied by the
/// caller, because this library ships no component databank.
///
/// The returned [`PrKappaResult`] carries warnings but no failure path beyond a
/// malformed spec: the polynomial is defined for every real `omega`, so there is
/// no input for which a refusal would be the right answer. A `kappa` below zero is
/// out of the range the coefficient is physical over and says so, rather than
/// being an error.
///
/// # Errors
/// In practice none. [`azoth_core::AzothError::OutOfRange`] would require a spec
/// bound with `error` severity, and the spec deliberately declares none - see its
/// `notes` for why an error bound would be an accuracy claim in
/// disguise.
///
/// # Example
/// ```
/// use azoth_eos::pr_kappa;
///
/// // Propane-like.
/// let r = pr_kappa(0.152)?;
/// assert!((r.kappa - 0.60282728832).abs() < 1e-15);
/// assert!(r.warnings.is_empty());
///
/// // Helium's acentric factor gives a negative coefficient, which is warned about
/// // rather than refused.
/// let helium = pr_kappa(-0.385)?;
/// assert!((helium.kappa + 0.259138992).abs() < 1e-12);
/// assert!(!helium.warnings.is_empty());
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pr_kappa(omega: f64) -> Result<PrKappaResult> {
    let spec = &spec_gen::PR_KAPPA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            _ => None,
        },
        &mut warnings,
    )?;

    let kappa = 0.37464 + 1.54226 * omega - 0.26992 * omega * omega;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "kappa").then_some(kappa),
        &mut warnings,
    )?;

    Ok(PrKappaResult { kappa, warnings })
}
