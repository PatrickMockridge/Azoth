//! `hydraulics.reynolds_number` - Reynolds number and flow regime.
//!
//! ```text
//! Re = rho * v * D / mu
//! ```
//!
//! Spec: `specs/calcs/hydraulics/reynolds_number.yaml`
//!
//! Every other calc in this crate consumes the Reynolds number, and two of them
//! consume the regime as well, so this is the entry point to the slice. The
//! regime boundaries are the Crane/Moody convention and are approximate: real
//! transition depends on inlet geometry, vibration and roughness. That is why
//! the transitional band produces a warning rather than a clean label.

use crate::results::ReynoldsNumberResult;
use crate::spec_gen;
use chemeng_core::units::{DynamicViscosity, Length, MassDensity, Velocity};
use chemeng_core::{FlowRegime, Result, Warning, apply_checks};

/// Reynolds number for flow in a circular pipe.
///
/// Returns the Reynolds number and its [`FlowRegime`]. When the flow is
/// transitional the result carries a
/// [`chemeng_core::WarningCode::TransitionalFlow`] warning, because the friction
/// factor - and therefore any pressure drop derived from it - is indeterminate
/// in that band rather than merely uncertain.
///
/// # Errors
/// Returns [`chemeng_core::ChemEngError::OutOfRange`] if density, diameter or
/// viscosity is not positive, or if velocity is negative. Those make the
/// Reynolds number undefined or meaningless, as opposed to merely out of range.
///
/// # Example
/// ```
/// use chemeng_core::units::{kilograms_per_cubic_meter, meters, meters_per_second, pascal_seconds};
/// use chemeng_hydraulics::reynolds_number;
///
/// let r = reynolds_number(
///     kilograms_per_cubic_meter(998.0),
///     meters_per_second(1.5),
///     meters(0.1),
///     pascal_seconds(1.002e-3),
/// )?;
/// assert!((r.re - 149_401.2).abs() < 0.5);
/// assert_eq!(r.regime.to_string(), "turbulent");
/// # Ok::<(), chemeng_core::ChemEngError>(())
/// ```
#[allow(non_snake_case)] // `D` is the symbol in the published equation
pub fn reynolds_number(
    rho: MassDensity,
    v: Velocity,
    D: Length,
    mu: DynamicViscosity,
) -> Result<ReynoldsNumberResult> {
    let spec = &spec_gen::REYNOLDS_NUMBER_SPEC;
    let mut warnings = Vec::new();

    let (rho_v, v_v, d_v, mu_v) = (rho.value, v.value, D.value, mu.value);

    // Input checks first: they guard the arithmetic below, and their severity is
    // `error`, so a violation here means we must not divide at all.
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "rho" => Some(rho_v),
            "v" => Some(v_v),
            "D" => Some(d_v),
            "mu" => Some(mu_v),
            _ => None,
        },
        &mut warnings,
    )?;

    let re = rho_v * v_v * d_v / mu_v;

    // Output and derived checks. `re` is an output here, so this is where the
    // transitional band is detected - and the spec names that check's warning
    // code as TRANSITIONAL_FLOW, so the caller gets the specific code rather
    // than a generic out-of-range.
    apply_checks(
        spec.derived_checks(),
        |q| (q == "re").then_some(re),
        &mut warnings,
    )?;

    Ok(ReynoldsNumberResult {
        re,
        regime: FlowRegime::from_reynolds_number(re),
        warnings,
    })
}

/// Classify a Reynolds number without computing it.
///
/// Provided so a caller that already has `Re` - from a measured flow, or from
/// another tool - can get the regime on the same convention this crate uses,
/// rather than reimplementing the boundaries and drifting from them.
#[must_use]
pub fn regime_for(re: f64) -> FlowRegime {
    FlowRegime::from_reynolds_number(re)
}

/// Warnings a result would carry purely because of the flow regime.
///
/// Exposed for the CLI, which composes several calcs and needs to explain the
/// regime once rather than repeat it per stage.
#[must_use]
pub fn regime_warning(regime: FlowRegime) -> Option<Warning> {
    regime.is_indeterminate().then(|| {
        Warning::new(
            chemeng_core::WarningCode::TransitionalFlow,
            format!(
                "flow is transitional (Re between {} and {}): the friction factor is \
                 indeterminate here, so pressure drop derived from it is not dependable",
                FlowRegime::LAMINAR_MAX,
                FlowRegime::TURBULENT_MIN
            ),
        )
    })
}
