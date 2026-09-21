//! `eos.hydrate_equilibrium_line` - a hydrate curve over a pressure grid.
//!
//! NeqSim's `HydrateEquilibriumLine`. It is a loop and not a model: ten equally spaced
//! pressures from a minimum to a maximum, and `eos.hydrate_formation_temperature` solved at
//! each. The equality, the cages and the structure are that model's; what this adds is the
//! grid and the vectors.
//!
//! # The seeding, and why it is not here
//!
//! NeqSim sets the system's temperature to the previous point's answer before each solve,
//! with the comment "hydrate T increases with P". Measured against the same ten pressures
//! solved independently on a fresh system, the two vectors agree to every printed digit
//! (`validation/neqsim/captures/hydrate_equilibrium_probe.tsv`, `worst_relative_gap = 0`), so
//! the seeding is a starting guess and not part of the answer - and this does without it. A
//! port that reproduced the seeding would be reproducing a speed-up.
//!
//! # The count is NeqSim's, and it is fixed
//!
//! `HydrateEquilibriumLine.numberOfPoints` is a field initialised to ten that no constructor
//! or setter ever writes, so the grid is ten points whatever a caller asks for. That is
//! NeqSim's surface rather than this one's: a caller wanting a different grid wants
//! [`crate::hydrate_formation_temperature`] over pressures of their own choosing, which is
//! the same arithmetic without the wrap.

use azoth_core::units::{Pressure, pascals};
use azoth_core::{AzothError, Result};

use crate::hydrate_formation_temperature::hydrate_formation_temperature;
use crate::mixture::Mixture;
use crate::results::HydrateEquilibriumLineResult;

/// The points NeqSim's line always solves, whatever the bounds.
const POINTS: usize = 10;

/// The formation temperature at ten pressures between a minimum and a maximum.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if either bound is not positive or the maximum is not above
///   the minimum.
/// * [`AzothError`] from any point's own solve. NeqSim's loop catches a failed point and then
///   reports the system's temperature anyway, which repeats the previous point's answer; this
///   refuses instead, so a grid that could not be solved is not a grid.
pub fn hydrate_equilibrium_line(
    mixture: &Mixture,
    p_min: Pressure,
    p_max: Pressure,
    z: &[f64],
) -> Result<HydrateEquilibriumLineResult> {
    let spec = &crate::model_gen::HYDRATE_EQUILIBRIUM_LINE_SPEC;
    let mut warnings = Vec::new();
    azoth_core::apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P_min" => Some(p_min.value),
            "P_max" => Some(p_max.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if p_min.value <= 0.0 || p_max.value <= 0.0 || p_max.value <= p_min.value {
        return Err(AzothError::invalid_input(
            "P_max",
            format!(
                "the grid runs from a positive minimum up to a larger maximum, and this was \
                 given {} and {} Pa",
                p_min.value, p_max.value
            ),
        ));
    }

    // `min + dp i` over `0..points`, which lands on the maximum at the last point - NeqSim's
    // own stepping, and why it never sets the maximum separately.
    let dp = (p_max.value - p_min.value) / (POINTS - 1) as f64;
    let mut temperature = Vec::with_capacity(POINTS);
    let mut pressure = Vec::with_capacity(POINTS);

    for point in 0..POINTS {
        let p = p_min.value + dp * point as f64;
        let solved = hydrate_formation_temperature(mixture, pascals(p), z)?;
        temperature.push(solved.temperature.value);
        pressure.push(p);
        warnings.extend(solved.warnings.iter().cloned());
    }

    Ok(HydrateEquilibriumLineResult {
        temperature,
        pressure,
        warnings,
    })
}
