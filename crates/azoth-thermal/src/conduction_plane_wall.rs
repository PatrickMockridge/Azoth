//! `thermal.conduction_plane_wall` - steady conduction through a plane wall.
//!
//! ```text
//! q = k * A * dT / L
//! ```
//!
//! Fourier, J. (1822). "Théorie analytique de la chaleur." Paris: Firmin Didot.
//!
//! Spec: `specs/calcs/thermal/conduction_plane_wall.yaml`
//!
//! # The sign convention
//!
//! `dT` may be negative and `q` follows its sign. This calc models a temperature
//! difference across a slab, not a named hot face and cold face, so it has no way
//! to say which side of the wall is which - and inventing a convention that the
//! inputs do not carry would be worse than returning the signed answer and letting
//! the caller say what it means. The spec argues the alternative and why it was
//! not taken.
//!
//! # Why `dT` is a `TemperatureInterval`
//!
//! It is a difference, not an absolute temperature, and `uom` distinguishes them
//! for a reason: a 30 K interval is a 30 degC interval, while an absolute 30 K is
//! -243.15 degC. Taking the interval type means a caller cannot hand this function
//! an absolute temperature and be silently offset by 273.15. See
//! `azoth_core::units::kelvin_intervals`.

use azoth_core::units::{Area, Length, TemperatureInterval, ThermalConductivity, watts};
use azoth_core::{Result, apply_checks};

use crate::results::ConductionPlaneWallResult;
use crate::spec_gen;

/// Steady heat flow through a plane wall.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if `k`, `A` or `L` is not positive.
///   `dT` is deliberately unbounded - see the module documentation.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kelvin_intervals, square_meters, meters, watts_per_meter_kelvin,
/// };
/// use azoth_thermal::conduction_plane_wall;
///
/// let r = conduction_plane_wall(
///     watts_per_meter_kelvin(45.0),
///     square_meters(2.0),
///     kelvin_intervals(30.0),
///     meters(0.05),
/// )?;
/// assert!((r.q.value - 54000.0).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `A`, `dT` and `L` are the symbols in the published equation
pub fn conduction_plane_wall(
    k: ThermalConductivity,
    A: Area,
    dT: TemperatureInterval,
    L: Length,
) -> Result<ConductionPlaneWallResult> {
    let spec = &spec_gen::CONDUCTION_PLANE_WALL_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "k" => Some(k.value),
            "A" => Some(A.value),
            "dT" => Some(dT.value),
            "L" => Some(L.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // Guarded by the checks above, so `L` is non-zero here.
    let q = k.value * A.value * dT.value / L.value;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "q").then_some(q),
        &mut warnings,
    )?;

    Ok(ConductionPlaneWallResult {
        q: watts(q),
        warnings,
    })
}
