//! `eos.heat_of_vaporization` - NeqSim's pure-component heat-of-vaporisation
//! correlation.
//!
//! ```text
//! hov = 1e-3*c0*(1 - Tr)**(c1 + c2*Tr + c3*Tr**2)
//! ```
//!
//! Spec: `specs/calcs/eos/heat_of_vaporization.toml`, which records the internal-unit
//! coefficients and the `1e-3` that recovers `J/mol`.

use azoth_core::units::{ThermodynamicTemperature, joules_per_mole};
use azoth_core::{Result, apply_checks};

use crate::results::HeatOfVaporizationResult;
use crate::spec_gen;

/// The pure-component heat of vaporisation at a temperature, from NeqSim's
/// correlation.
///
/// `c0..c3` are the raw `heatofvaporizationcoefs1..4` and `Tc` the critical
/// temperature, all caller-supplied.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `Tc` or `T` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::heat_of_vaporization;
///
/// let r = heat_of_vaporization(52_100_000.0, 0.32, -0.212, 0.258, kelvins(425.12), kelvins(300.0))?;
/// assert!((r.hov.value - 36147.57779141032).abs() < 1e-8);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `Tc` is the symbol in the published equation
pub fn heat_of_vaporization(
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    Tc: ThermodynamicTemperature,
    T: ThermodynamicTemperature,
) -> Result<HeatOfVaporizationResult> {
    let spec = &spec_gen::HEAT_OF_VAPORIZATION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "c0" => Some(c0),
            "c1" => Some(c1),
            "c2" => Some(c2),
            "c3" => Some(c3),
            "Tc" => Some(Tc.value),
            "T" => Some(T.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let tr = T.value / Tc.value;
    let exponent = c1 + c2 * tr + c3 * tr * tr;
    let hov = 1e-3 * c0 * (1.0 - tr).powf(exponent);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "hov").then_some(hov),
        &mut warnings,
    )?;

    Ok(HeatOfVaporizationResult {
        hov: joules_per_mole(hov),
        warnings,
    })
}
