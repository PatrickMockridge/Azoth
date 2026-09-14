//! `hydraulics.friction_factor_haaland` - the explicit Haaland friction factor.
//!
//! ```text
//! f = (1.0 / (-1.8 * log10((relative_roughness / 3.7)**1.11 + 6.9 / Re)))**2
//! ```
//!
//! Spec: `specs/calcs/hydraulics/friction_factor_haaland.toml`, which carries the
//! citation.
//!
//! This is an *approximation to* Colebrook, not a more correct alternative to it,
//! and it is the second such approximation here alongside Swamee-Jain. They were
//! fitted differently and are accurate to about 1% and 2% respectively, so a
//! caller comparing them at the same inputs gets a cheap sense of how much the
//! choice of explicit form matters.
//!
//! The attribution is unconfirmed: the equation is not in doubt, but nobody has
//! opened the paper and checked this against it. See the spec's `notes`, which
//! also record that the valid range below is the framework's rather than a range
//! read from the source.

use crate::results::HaalandResult;
use crate::spec_gen;
use azoth_core::{Result, apply_checks};

/// Haaland explicit friction factor. `relative_roughness` is `epsilon / D`.
///
/// Returns immediately: no iteration, and so no convergence report.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Re <= 0` (the `6.9/Re` term is
///   singular at the origin) or if `relative_roughness < 0`.
///
/// Outside the range in which the correlation describes turbulent pipe flow the
/// value is still returned, carrying an `OutOfValidRange` warning: the equation
/// is well defined there, it simply is not covered by the accuracy claim.
///
/// # Example
/// ```
/// use azoth_hydraulics::friction_factor_haaland;
///
/// let r = friction_factor_haaland(100_000.0, 4.6e-4)?;
/// assert!((r.f - 0.019898058).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn friction_factor_haaland(re: f64, relative_roughness: f64) -> Result<HaalandResult> {
    let spec = &spec_gen::FRICTION_FACTOR_HAALAND_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "re" => Some(re),
            "relative_roughness" => Some(relative_roughness),
            _ => None,
        },
        &mut warnings,
    )?;

    // `1/sqrt(f) = -1.8 * log10(...)`, so `f` is the reciprocal of that, squared.
    // The square is `x * x` rather than `x.powi(2)`, matching Swamee-Jain and the
    // Python side.
    let inner = (relative_roughness / 3.7).powf(1.11) + 6.9 / re;
    let inv_sqrt_f = -1.8 * inner.log10();
    let f = 1.0 / (inv_sqrt_f * inv_sqrt_f);

    apply_checks(
        spec.derived_checks(),
        |q| (q == "f").then_some(f),
        &mut warnings,
    )?;

    Ok(HaalandResult { f, warnings })
}
