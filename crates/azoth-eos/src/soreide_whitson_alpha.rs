//! `eos.soreide_whitson_alpha` - the Soreide-Whitson alpha function for water.
//!
//! ```text
//! alpha = (1 + 0.453*(1 - Tr*(1 - 0.0103*salinity**1.1)) + 0.0034*((1/Tr)**3 - 1))**2
//! ```
//!
//! Spec: `specs/calcs/eos/soreide_whitson_alpha.toml`, which carries the provenance and
//! why the salinity is molality and the water-only scope.

use azoth_core::{Result, apply_checks};

use crate::results::SoreideWhitsonAlphaResult;
use crate::spec_gen;

/// The Soreide-Whitson alpha function for water.
///
/// `salinity` is the molality (mol NaCl / kg H2O); `Tr` the reduced temperature. For a
/// non-water component NeqSim's class delegates to the 1978 Peng-Robinson alpha, which
/// is ported separately as `pr78_kappa` and the standard alpha form.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tr <= 0` or `salinity < 0`.
///
/// # Example
/// ```
/// use azoth_eos::soreide_whitson_alpha;
///
/// let r = soreide_whitson_alpha(1.0, 0.7)?;
/// assert!((r.alpha - 1.3125796067429514).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn soreide_whitson_alpha(salinity: f64, Tr: f64) -> Result<SoreideWhitsonAlphaResult> {
    let spec = &spec_gen::SOREIDE_WHITSON_ALPHA_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "salinity" => Some(salinity),
            "Tr" => Some(Tr),
            _ => None,
        },
        &mut warnings,
    )?;

    let a = 1.0
        + 0.453 * (1.0 - Tr * (1.0 - 0.0103 * salinity.powf(1.1)))
        + 0.0034 * ((1.0 / Tr).powf(3.0) - 1.0);
    let alpha = a * a;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "alpha").then_some(alpha),
        &mut warnings,
    )?;

    Ok(SoreideWhitsonAlphaResult { alpha, warnings })
}
