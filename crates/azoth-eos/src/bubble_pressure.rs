//! `eos.bubble_pressure` - the pressure at which a liquid first gives off vapour.
//!
//! Spec: `specs/models/eos/bubble_pressure.yaml`
//!
//! The iteration itself lives in [`crate::phase_boundary`], because `eos.dew_pressure`
//! runs the same loop with the phases exchanged and the guard against the trivial
//! solution has to be written once. What is here is the spec lookup, the range
//! checks and the result.

use azoth_core::units::{ThermodynamicTemperature, pascals};
use azoth_core::{Result, apply_checks};

use crate::mixture::Mixture;
use crate::phase_boundary::{Incipient, phase_boundary_pressure};
use crate::results::BubblePressureResult;
use crate::{model_gen, phase_boundary};

/// The pressure at which a liquid of composition `x` first gives off vapour.
///
/// `x` is the liquid's composition and is taken as given: this model does not ask
/// whether that liquid is stable, only where its bubble point is.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if `T` is not positive, or if the mixture has no
///   bubble point at this temperature - which is reported on `min_t_over_tc` and
///   means the mixture is at or above its critical condition.
/// * [`AzothError::InvalidInput`] if the mixture has one component, or if `x` is the
///   wrong length, has a negative entry, or does not sum to one.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::{bubble_pressure, databank};
///
/// let mixture = databank::mixture_of(&["methane", "n-butane"], None)
///     .expect("the pair resolves")
///     .0;
/// let r = bubble_pressure(&mixture, kelvins(300.0), &[0.2, 0.8])?;
/// assert!((r.pressure.value / 1e6 - 3.9509604937437346).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn bubble_pressure(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    x: &[f64],
) -> Result<BubblePressureResult> {
    let spec = &model_gen::BUBBLE_PRESSURE_SPEC;
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

    let boundary = phase_boundary_pressure(mixture, t, x, Incipient::Vapour)?;
    warnings.extend(boundary.warnings);

    Ok(BubblePressureResult {
        pressure: pascals(boundary.pressure),
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

/// The trivial-solution threshold, re-exported so the spec page and the tests can
/// name the same number.
///
/// Public for the same reason [`crate::pr_alpha_ab::OMEGA_A`] is: a reader checking
/// the model's behaviour against its documentation should be able to see that the
/// constant the documentation quotes is the constant the code uses.
pub const TRIVIAL_TOLERANCE: f64 = phase_boundary::TRIVIAL_TOLERANCE;
