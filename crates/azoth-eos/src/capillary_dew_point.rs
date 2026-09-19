//! `eos.capillary_dew_point` - the temperature at which a vapour condenses inside a pore.
//!
//! Spec: `specs/models/eos/capillary_dew_point.toml`. NeqSim's `CapillaryDewPointFlash`, reached
//! through `ThermodynamicOperations.capillaryDewPointTemperatureFlash(r[, theta])`.
//!
//! The same boundary [`crate::dew_temperature`] finds, drawn on a curved interface instead of a
//! flat one. Young-Laplace puts the liquid inside a pore at a pressure above the vapour outside
//! it by `2 sigma cos(theta) / r`, and the Kelvin equation turns that into a shift on every
//! K-value:
//!
//! ```text
//! K_cap,i = K_i exp(-Vm_L dP_cap / (R T))
//! ```
//!
//! so the incipient liquid is *stabilised* - the dew point moves **up**, not down, which is the
//! direction that matters for condensation in tight rock. Measured on methane/n-butane 50/50 at
//! 300 K and 20 bar with `sigma = 0.005 N/m`: the bulk dew point is `345.1416 K`, and the
//! capillary one is `347.6659 K` at 10 nm, `345.4016 K` at 100 nm and `345.1677 K` at 1 µm. The
//! shift scales as `1/r` and as `cos(theta)`, and both were checked.
//!
//! **The surface tension is the caller's, not the mixture's.** NeqSim reads it off the phase's
//! interphase properties and falls back to a constant `0.005 N/m` when that returns nothing;
//! this takes it as an argument, because a surface tension is a fitted quantity with its own
//! provenance - `eos.parachor_surface_tension` computes one - and reading it from a correlation
//! inside a model whose other inputs are all stated would hide which one was used.

use azoth_core::units::{Pressure, kelvins, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::Mixture;
use crate::model_gen;
use crate::phase_boundary::Incipient;
use crate::results::CapillaryDewPointResult;
use crate::saturation_temperature::{Curvature, phase_boundary_temperature};

/// The dew-point temperature of a vapour held in a pore of a stated radius.
///
/// `y` is the vapour's composition and is taken as given: this model does not ask whether that
/// vapour is stable, only where its dew point is. `pore_radius` is in metres, `contact_angle` in
/// radians - zero is a perfectly wetting liquid - and `surface_tension` in newtons per metre.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `pore_radius` is not positive, `surface_tension` is
///   negative, `y` is the wrong length, has a negative entry, or does not sum to one.
/// * [`AzothError::OutOfRange`] if `P` is not positive, or if the mixture has no dew point at
///   this pressure.
/// * [`AzothError::SolverNotConverged`] if the iteration hits its cap.
///
/// # Example
/// ```
/// use azoth_core::units::{pascals};
/// use azoth_eos::Cubic;
/// use azoth_eos::databank::mixture_of;
/// use azoth_eos::capillary_dew_point;
///
/// let (mixture, _) = mixture_of(&["methane", "n-butane"], Cubic::Pr, None)?;
/// // A wide pore is the flat interface, and the shift goes to nothing with it.
/// let r = capillary_dew_point(&mixture, pascals(2.0e6), &[0.5, 0.5], 1.0e-3, 0.0, 0.005)?;
/// assert!((r.temperature.value - 345.1416).abs() < 0.01);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `P` is the symbol in the chemistry
pub fn capillary_dew_point(
    mixture: &Mixture,
    P: Pressure,
    y: &[f64],
    pore_radius: f64,
    contact_angle: f64,
    surface_tension: f64,
) -> Result<CapillaryDewPointResult> {
    let spec = &model_gen::CAPILLARY_DEW_POINT_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(P.value),
            "pore_radius" => Some(pore_radius),
            "surface_tension" => Some(surface_tension),
            _ => None,
        },
        &mut warnings,
    )?;

    // The radius and the tension are range-checked by the spec's `valid_range` blocks, through
    // `apply_checks` above - a non-positive radius divides the Young-Laplace pressure and a
    // negative tension moves the dew point the wrong way, and both are declared there with
    // their reasons rather than repeated here.
    //
    // The contact angle has no range: every real angle is admissible, `0` and `pi/2` included,
    // and the two outside the first quadrant are the same interface by symmetry. Only a
    // non-finite one has no meaning.
    if !contact_angle.is_finite() {
        return Err(AzothError::invalid_input(
            "contact_angle",
            "the contact angle is not a finite number",
        ));
    }

    let curvature = Curvature {
        pore_radius,
        contact_angle,
        surface_tension,
    };
    let boundary = phase_boundary_temperature(mixture, P, y, Incipient::Liquid, Some(curvature))?;
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

    Ok(CapillaryDewPointResult {
        temperature: kelvins(boundary.temperature),
        incipient: boundary.incipient,
        k: boundary.k,
        z_liquid: boundary.z_incipient,
        z_vapour: boundary.z_held,
        capillary_pressure: pascals(curvature.capillary_pressure()),
        min_t_over_tc,
        iterations: boundary.iterations,
        residual: boundary.residual,
        warnings,
    })
}
