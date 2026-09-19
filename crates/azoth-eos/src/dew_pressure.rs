//! `eos.dew_pressure` - the pressure at which a vapour first condenses.
//!
//! Spec: `specs/models/eos/dew_pressure.toml`
//!
//! The companion of [`crate::bubble_pressure`], and the same iteration with the
//! phases exchanged. The loop, the initialisation and the guard against the trivial
//! solution live in [`crate::phase_boundary`] so that they are written once - the
//! two models differ in which composition is the input and which is solved for, and
//! in one line of the pressure update, not in their procedure.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};

use crate::mixture::Mixture;
use crate::model_gen;
use crate::phase_boundary::{Incipient, phase_boundary_pressure};
use crate::results::DewPressureResult;

/// The pressure at which a vapour of composition `y` first condenses.
///
/// `y` is the vapour's composition and is taken as given: this model does not ask
/// whether that vapour is stable, only where its dew point is.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` is not positive, or if the mixture has no dew
///   point at this temperature - reported on `min_t_over_tc`, meaning the mixture is
///   at or above its critical condition.
/// * [`AzothError::InvalidInput`] if the mixture has one component, or if `y` is the
///   wrong length, has a negative entry, or does not sum to one.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
///
/// # Example
/// ```
/// use azoth_eos::Cubic;
/// use azoth_core::units::kelvins;
/// use azoth_eos::{databank, dew_pressure};
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
///     .expect("the pair resolves")
///     .0;
/// let r = dew_pressure(&mixture, kelvins(300.0), &[0.8, 0.2])?;
/// assert!((r.pressure.value / 1e6 - 1.5673473649394323).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn dew_pressure(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    y: &[f64],
) -> Result<DewPressureResult> {
    let spec = &model_gen::DEW_PRESSURE_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| (quantity == "T").then_some(t.value),
        &mut warnings,
    )?;

    let min_t_over_tc = mixture
        .components()
        .iter()
        .map(|c| t.value / c.tc.value)
        .fold(f64::INFINITY, f64::min);
    apply_checks(
        spec.derived_checks(),
        |quantity| (quantity == "min_t_over_tc").then_some(min_t_over_tc),
        &mut warnings,
    )?;

    let boundary = phase_boundary_pressure(mixture, t, y, Incipient::Liquid)?;
    warnings.extend(boundary.warnings);

    Ok(DewPressureResult {
        pressure: pascals(boundary.pressure),
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
