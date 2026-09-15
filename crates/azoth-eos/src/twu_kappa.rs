//! `eos.twu_kappa` - the Twu attraction-parameter coefficient.
//!
//! ```text
//! kappa = 0.48 + 1.574*omega - 0.175*omega**2
//! ```
//!
//! Spec: `specs/calcs/eos/twu_kappa.toml`. Twu's coefficient, which differs from
//! Soave's only in the last constant - 0.175 against 0.176 - and feeds the same Soave
//! alpha form.

use azoth_core::{Result, apply_checks};

use crate::results::TwuKappaResult;
use crate::spec_gen;

/// Twu's alpha-function coefficient for a pure component.
///
/// # Example
/// ```
/// use azoth_eos::twu_kappa;
///
/// let r = twu_kappa(0.152)?;
/// assert!((r.kappa - 0.7152048).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn twu_kappa(omega: f64) -> Result<TwuKappaResult> {
    let spec = &spec_gen::TWU_KAPPA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "omega" => Some(omega),
            _ => None,
        },
        &mut warnings,
    )?;

    let kappa = 0.48 + 1.574 * omega - 0.175 * omega * omega;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "kappa").then_some(kappa),
        &mut warnings,
    )?;

    Ok(TwuKappaResult { kappa, warnings })
}
