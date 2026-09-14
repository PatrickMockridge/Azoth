//! `hydraulics.orifice_flow` - flow through an orifice from its pressure difference.
//!
//! ```text
//! q = Cd * (pi * d**2 / 4) * sqrt(2 * dP / rho)
//! ```
//!
//! Spec: `specs/calcs/hydraulics/orifice_flow.toml`
//!
//! # The discharge coefficient is an input, deliberately
//!
//! `Cd` is supplied by the caller; this calc does not compute it. ISO 5167's
//! discharge-coefficient equation is a long fitted expression whose constants come
//! from a table of experimental results, and reproducing that table is what this
//! project's copyright rule forbids - "just implementing the equation" being exactly
//! how it would happen by accident. Taking `Cd` as an input keeps this calc to the
//! relation that is genuinely public.
//!
//! # Which `Cd`
//!
//! This calc applies the equation exactly as written and adds no velocity-of-approach
//! factor. ISO 5167's coefficient already includes `1/sqrt(1 - beta**4)`; older texts
//! write the two separately. The caller must supply the coefficient for *this* form.

use azoth_core::units::{Length, MassDensity, Pressure, cubic_meters_per_second};
use azoth_core::{Result, apply_checks};

use crate::results::OrificeFlowResult;
use crate::spec_gen;

/// Volumetric flow through an orifice. `d` is the bore, `Cd` is dimensionless.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `d` or `rho` is not positive, if
///   `dP` is negative, or if `Cd` is outside `(0, 1]`.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, millimeters, pascals};
/// use azoth_hydraulics::orifice_flow;
///
/// let r = orifice_flow(
///     millimeters(50.0),
///     pascals(25000.0),
///     kilograms_per_cubic_meter(998.0),
///     0.62,
/// )?;
/// assert!((r.q.value - 0.008616706712060997).abs() < 1e-15);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Cd` is the symbol in the published equation
pub fn orifice_flow(
    d: Length,
    dP: Pressure,
    rho: MassDensity,
    Cd: f64,
) -> Result<OrificeFlowResult> {
    let spec = &spec_gen::ORIFICE_FLOW_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            // `.value` is the SI base magnitude, so a 50 mm bore is 0.05 here -
            // the same number the Python side converts to.
            "d" => Some(d.value),
            "dP" => Some(dP.value),
            "rho" => Some(rho.value),
            // Dimensionless, so it arrives as a plain f64 with no unit to extract.
            "Cd" => Some(Cd),
            _ => None,
        },
        &mut warnings,
    )?;

    // Guarded by the checks above, so d and rho are positive and dP is not negative.
    let area = std::f64::consts::PI * d.value * d.value / 4.0;
    let q = Cd * area * (2.0 * dP.value / rho.value).sqrt();

    apply_checks(
        spec.derived_checks(),
        |name| (name == "q").then_some(q),
        &mut warnings,
    )?;

    Ok(OrificeFlowResult {
        q: cubic_meters_per_second(q),
        warnings,
    })
}
