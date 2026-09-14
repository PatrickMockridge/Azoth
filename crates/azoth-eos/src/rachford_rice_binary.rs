//! `eos.rachford_rice_binary` - the Rachford-Rice vapour fraction, for a binary.
//!
//! ```text
//! beta = -(z1*(K1 - 1) + z2*(K2 - 1)) / ((K1 - 1)*(K2 - 1))     z2 = 1 - z1
//! ```
//!
//! Spec: `specs/calcs/eos/rachford_rice_binary.toml`, which carries the provenance, the
//! derivation of the closed form, and why a `beta` outside `[0, 1]` is a warning rather
//! than an error.

use azoth_core::{Result, apply_checks};

use crate::results::RachfordRiceBinaryResult;
use crate::spec_gen;

/// The vapour fraction that solves the Rachford-Rice equation for two components.
///
/// `z1` is the overall mole fraction of component 1; component 2 is `1 - z1`.
/// `K1` and `K2` are the K-values, `y/x`, for the same temperature and pressure.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `z1` is outside `[0, 1]`, if either
///   K-value is not positive, or if either equals 1 - which makes `K - 1` a divisor
///   and the component degenerate, dropping out of the sum entirely.
///
/// A `beta` outside `[0, 1]` is *not* an error: it is returned with an
/// `OUT_OF_VALID_RANGE` warning.
///
/// # Example
/// ```
/// use azoth_eos::rachford_rice_binary;
///
/// let r = rachford_rice_binary(0.6, 4.0, 0.25)?;
/// assert!((r.beta - 0.6666666666666665).abs() < 1e-15);
/// assert!(r.warnings.is_empty());
///
/// // A feed that does not split: both K-values exceed 1.
/// let single = rachford_rice_binary(0.5, 2.0, 1.5)?;
/// assert!((single.beta + 1.5).abs() < 1e-15);
/// assert!(!single.warnings.is_empty());
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `K1` and `K2` are the symbols in the published equation
pub fn rachford_rice_binary(z1: f64, K1: f64, K2: f64) -> Result<RachfordRiceBinaryResult> {
    let spec = &spec_gen::RACHFORD_RICE_BINARY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "z1" => Some(z1),
            "K1" => Some(K1),
            "K2" => Some(K2),
            _ => None,
        },
        &mut warnings,
    )?;

    // Hoisted but one operation each, so this is the same double as writing
    // `K1 - 1.0` inline four times - which is what the spec's equation shows.
    let a = K1 - 1.0;
    let b = K2 - 1.0;
    // Guarded by the `equals: 1` bounds above, so neither divisor is zero.
    let beta = -(z1 * a + (1.0 - z1) * b) / (a * b);

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "beta" => Some(beta),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(RachfordRiceBinaryResult { beta, warnings })
}
