//! `eos.pr78_kappa` - the Peng-Robinson (1978) attraction-parameter coefficient.
//!
//! ```text
//! kappa = 0.379642 + 1.48503*omega - 0.164423*omega**2 + 0.01666*omega**3   if omega > 0.49
//! kappa = 0.37464 + 1.54226*omega - 0.26992*omega**2                        otherwise
//! ```
//!
//! Spec: `specs/calcs/eos/pr78_kappa.toml`. The 1978 revision of the Peng-Robinson
//! coefficient: the 1976 form for light components, a heavier-acentric branch for
//! `omega` above 0.49. Both feed the same Soave alpha form.

use azoth_core::{Result, apply_checks};

use crate::results::Pr78KappaResult;
use crate::spec_gen;

/// The 1978 Peng-Robinson alpha-function coefficient for a pure component.
///
/// # Example
/// ```
/// use azoth_eos::pr78_kappa;
///
/// let heavy = pr78_kappa(0.6)?;
/// assert!((heavy.kappa - 1.2150662799999998).abs() < 1e-12);
/// let light = pr78_kappa(0.152)?;
/// assert!((light.kappa - 0.60282728832).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn pr78_kappa(omega: f64) -> Result<Pr78KappaResult> {
    let spec = &spec_gen::PR78_KAPPA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            _ => None,
        },
        &mut warnings,
    )?;

    let kappa = if omega > 0.49 {
        0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.01666 * omega * omega * omega
    } else {
        0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    };

    apply_checks(
        spec.derived_checks(),
        |name| (name == "kappa").then_some(kappa),
        &mut warnings,
    )?;

    Ok(Pr78KappaResult { kappa, warnings })
}
