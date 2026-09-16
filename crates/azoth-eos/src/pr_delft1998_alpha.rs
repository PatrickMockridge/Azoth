//! `eos.pr_delft1998_alpha` - the Peng-Robinson alpha function, Delft (1998).
//!
//! ```text
//! alpha = (1 + m*(1 - sqrt(Tr)))**2
//! m     = 0.379642 + 1.48503*omega - 0.164423*omega**2 + 0.01666*omega**3  if omega > 0.49
//! m     = 0.37464 + 1.54226*omega - 0.26992*omega**2                        otherwise
//! ```
//!
//! Spec: `specs/calcs/eos/pr_delft1998_alpha.toml`, which carries the provenance and the
//! methane-specific branch this general form does not express.

use azoth_core::{Result, apply_checks};

use crate::results::PrDelft1998AlphaResult;
use crate::spec_gen;

/// The Peng-Robinson alpha function, Delft (1998), for a pure component.
///
/// `omega` is the acentric factor; `Tr` the reduced temperature. The `m` is the 1978
/// Peng-Robinson coefficient, the same one `pr78_kappa` computes.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` - the formula takes a square
///   root of `Tr`.
///
/// # Example
/// ```
/// use azoth_eos::pr_delft1998_alpha;
///
/// let r = pr_delft1998_alpha(0.1, 0.7)?;
/// assert!((r.alpha - 1.1792745256672486).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_delft1998_alpha(omega: f64, Tr: f64) -> Result<PrDelft1998AlphaResult> {
    let spec = &spec_gen::PR_DELFT1998_ALPHA_SPEC;
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

    let m = if omega > 0.49 {
        0.379642 + 1.48503 * omega - 0.164423 * omega * omega + 0.01666 * omega * omega * omega
    } else {
        0.37464 + 1.54226 * omega - 0.26992 * omega * omega
    };
    let t = 1.0 + m * (1.0 - Tr.sqrt());
    let alpha = t * t;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(PrDelft1998AlphaResult { alpha, warnings })
}
