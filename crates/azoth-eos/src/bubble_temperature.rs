//! `eos.bubble_temperature` - the temperature at which a liquid first gives off vapour.
//!
//! Spec: `specs/models/eos/bubble_temperature.toml`
//!
//! The iteration lives in [`crate::saturation_temperature`], shared with
//! `eos.dew_temperature`.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{Result, apply_checks};

use crate::mixture::Mixture;
use crate::phase_boundary::Incipient;
use crate::results::BubbleTemperatureResult;
use crate::saturation_temperature::phase_boundary_temperature;
use crate::{model_gen, phase_boundary};

/// The temperature at which a liquid of composition `x` first gives off vapour at a
/// pressure.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `P` is not positive, or if the mixture has no bubble
///   point at this pressure.
/// * [`AzothError::InvalidInput`] if the mixture has one component, or `x` is not a
///   composition.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{pascals, kelvins};
/// use azoth_eos::{bubble_pressure, bubble_temperature, databank};
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves")
///     .0;
/// // The round trip: the bubble pressure at 300 K, held at that pressure.
/// let bp = bubble_pressure(&mixture, kelvins(300.0), &[0.2, 0.8])?;
/// let r = bubble_temperature(&mixture, pascals(bp.pressure.value), &[0.2, 0.8])?;
/// assert!((r.temperature.value - 300.0).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn bubble_temperature(
    mixture: &Mixture,
    p: Pressure,
    x: &[f64],
) -> Result<BubbleTemperatureResult> {
    let spec = &model_gen::BUBBLE_TEMPERATURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "P").then_some(p.value),
        &mut warnings,
    )?;

    let boundary = phase_boundary_temperature(mixture, p, x, Incipient::Vapour, None)?;
    warnings.extend(boundary.warnings);

    let min_t_over_tc = mixture
        .components()
        .iter()
        .map(|c| boundary.temperature / c.tc.value)
        .fold(f64::INFINITY, f64::min);
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "min_t_over_tc").then_some(min_t_over_tc),
        &mut warnings,
    )?;

    Ok(BubbleTemperatureResult {
        temperature: kelvins(boundary.temperature),
        incipient: boundary.incipient,
        k: boundary.k,
        z_liquid: boundary.z_held,
        z_vapour: boundary.z_incipient,
        min_t_over_tc,
        iterations: boundary.iterations,
        residual: boundary.residual,
        warnings,
    })
}

/// The trivial-solution threshold, re-exported.
pub const TRIVIAL_TOLERANCE: f64 = phase_boundary::TRIVIAL_TOLERANCE;
