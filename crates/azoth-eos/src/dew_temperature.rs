//! `eos.dew_temperature` - the temperature at which a vapour first gives off liquid.
//!
//! Spec: `specs/models/eos/dew_temperature.toml`
//!
//! The iteration lives in [`crate::saturation_temperature`], shared with
//! `eos.bubble_temperature`.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{Result, apply_checks};

use crate::mixture::Mixture;
use crate::phase_boundary::Incipient;
use crate::results::DewTemperatureResult;
use crate::saturation_temperature::phase_boundary_temperature;
use crate::{model_gen, phase_boundary};

/// The temperature at which a vapour of composition `y` first gives off liquid at a
/// pressure.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `P` is not positive, or if the mixture has no dew
///   point at this pressure.
/// * [`AzothError::InvalidInput`] if the mixture has one component, or `y` is not a
///   composition.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{pascals, kelvins};
/// use azoth_eos::{dew_pressure, dew_temperature, databank};
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves")
///     .0;
/// // The round trip: the dew pressure at 300 K, held at that pressure.
/// let dp = dew_pressure(&mixture, kelvins(300.0), &[0.8, 0.2])?;
/// let r = dew_temperature(&mixture, pascals(dp.pressure.value), &[0.8, 0.2])?;
/// assert!((r.temperature.value - 300.0).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn dew_temperature(mixture: &Mixture, p: Pressure, y: &[f64]) -> Result<DewTemperatureResult> {
    let spec = &model_gen::DEW_TEMPERATURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "P").then_some(p.value),
        &mut warnings,
    )?;

    let boundary = phase_boundary_temperature(mixture, p, y, Incipient::Liquid)?;
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

    Ok(DewTemperatureResult {
        temperature: kelvins(boundary.temperature),
        incipient: boundary.incipient,
        k: boundary.k,
        z_liquid: boundary.z_incipient,
        z_vapour: boundary.z_held,
        min_t_over_tc,
        iterations: boundary.iterations,
        residual: boundary.residual,
        warnings,
    })
}

/// The trivial-solution threshold, re-exported.
pub const TRIVIAL_TOLERANCE: f64 = phase_boundary::TRIVIAL_TOLERANCE;
