//! `hydraulics.pump_power` - shaft power from flow, head and efficiency.
//!
//! ```text
//! power = rho * g * q * H / eta
//! ```
//!
//! Spec: `specs/calcs/hydraulics/pump_power.yaml`
//!
//! # `g` is a constant here, not an input
//!
//! Standard gravity is a *defined* value - the CGPM fixed it exactly in 1901 - so
//! it belongs to the equation in the way pi does, rather than being a measured
//! property assumed behind the caller's back. The spec records the alternative of
//! taking it as an input and why that was rejected.
//!
//! # What `H` is, and what it is not
//!
//! `H` is the head the pump *delivers*, in metres of the pumped fluid. It is not
//! the head the impeller generates before the pump's own internal losses: those are
//! what `eta` accounts for, and counting them in `H` as well would apply them twice
//! and report a power that is too high.

use azoth_core::units::{Length, MassDensity, VolumeRate, watts};
use azoth_core::{Result, apply_checks};

use crate::results::PumpPowerResult;
use crate::spec_gen;

/// Standard gravity, in metres per second squared.
///
/// Exact by definition: the CGPM fixed it in 1901, which is what makes it a
/// constant of the equation rather than an input. Named rather than written inline
/// so it is greppable and so the Python side can be compared against it.
pub const STANDARD_GRAVITY_M_S2: f64 = 9.80665;

/// Shaft power a pump must be supplied with. `eta` is dimensionless.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `rho` is not positive, if `q` or `H`
///   is negative, or if `eta` is outside `(0, 1]`.
///
/// # Example
/// ```
/// use azoth_core::units::{cubic_meters_per_second, kilograms_per_cubic_meter, meters};
/// use azoth_hydraulics::pump_power;
///
/// let r = pump_power(
///     kilograms_per_cubic_meter(998.0),
///     cubic_meters_per_second(0.01),
///     meters(30.0),
///     0.75,
/// )?;
/// assert!((r.power.value - 3914.81468).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `H` is the symbol in the published equation
pub fn pump_power(rho: MassDensity, q: VolumeRate, H: Length, eta: f64) -> Result<PumpPowerResult> {
    let spec = &spec_gen::PUMP_POWER_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "rho" => Some(rho.value),
            "q" => Some(q.value),
            "H" => Some(H.value),
            // Dimensionless, so it arrives as a plain f64 with no unit to extract -
            // the same treatment colebrook's `re` gets.
            "eta" => Some(eta),
            _ => None,
        },
        &mut warnings,
    )?;

    let power = rho.value * STANDARD_GRAVITY_M_S2 * q.value * H.value / eta;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "power").then_some(power),
        &mut warnings,
    )?;

    Ok(PumpPowerResult {
        power: watts(power),
        warnings,
    })
}
