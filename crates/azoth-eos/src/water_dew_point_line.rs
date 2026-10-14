//! `eos.water_dew_point_line` - the water dew point against pressure.
//!
//! NeqSim's `WaterDewPointEquilibriumLine`, reached through
//! `ThermodynamicOperations.waterDewPointLine(minimumPressure, maximumPressure)`. Spec:
//! `specs/models/eos/water_dew_point_line.toml`. Oracle:
//! `validation/neqsim/captures/water_dew_point_probe.tsv`.
//!
//! # It is the multiphase route, ten times
//!
//! The class holds no arithmetic of its own. It clones the fluid, steps the pressure from
//! `minPressure` to `maxPressure` in `numberOfPoints = 10` equal increments, calls
//! `waterDewPointTemperatureMultiphaseFlash()` at each, and records the temperature the fluid
//! came back at. So this composes [`crate::water_dew_point_multiphase`] rather than restating
//! it, and `eos.hydrate_equilibrium_line` is the same shape one family over.
//!
//! # The catch is upstream's, and so is the silence it leaves
//!
//! Upstream wraps each point in a `try`, logs the failure and records **whatever temperature
//! the fluid is at** - which for a failed point is the previous point's. So a line with a bad
//! point is a line with a repeated temperature in it, and nothing in the answer says so. The
//! port reproduces the number and reports `failed_points` beside it, because a caller reading
//! the curve cannot otherwise tell a flat stretch from a point that did not solve.

use azoth_core::units::{Pressure, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Cubic;
use crate::model_gen;
use crate::results::WaterDewPointLineResult;
use crate::water_dew_point_multiphase::water_dew_point_multiphase;

/// How many points the class sweeps. Its own field, with no setter - so this is the class's
/// number and not a caller's choice.
const NUMBER_OF_POINTS: usize = 10;

/// The water dew point at each of ten pressures from `minimum` to `maximum`.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `minimum_pressure` is not positive, if it is above
///   `maximum_pressure`, or if the fluid names no water - each of which makes every point of
///   the line the same failure.
/// * Propagates the route's and the kernels' range checks.
pub fn water_dew_point_line(
    components: &[String],
    z: &[f64],
    minimum_pressure: Pressure,
    maximum_pressure: Pressure,
    cubic: Cubic,
) -> Result<WaterDewPointLineResult> {
    let spec = &model_gen::WATER_DEW_POINT_LINE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "minimum_pressure" => Some(minimum_pressure.value),
            "maximum_pressure" => Some(maximum_pressure.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if minimum_pressure.value <= 0.0 || maximum_pressure.value <= 0.0 {
        return Err(AzothError::invalid_input(
            "minimum_pressure",
            "a pressure is a divisor in the reduced variables and both ends of the sweep are \
             finite and positive by definition"
                .to_string(),
        ));
    }
    if minimum_pressure.value > maximum_pressure.value {
        return Err(AzothError::invalid_input(
            "minimum_pressure",
            format!(
                "the sweep starts at {} Pa and ends at {} Pa, so it runs backwards",
                minimum_pressure.value, maximum_pressure.value
            ),
        ));
    }

    let step = (maximum_pressure.value - minimum_pressure.value) / (NUMBER_OF_POINTS - 1) as f64;
    let mut temperatures = Vec::with_capacity(NUMBER_OF_POINTS);
    let mut pressures = Vec::with_capacity(NUMBER_OF_POINTS);
    let mut failed = 0_u32;
    // The class's own last-resort behaviour: a point that does not solve keeps the temperature
    // the fluid was already at, which is the previous point's.
    let mut last = f64::NAN;
    for index in 0..NUMBER_OF_POINTS {
        let p = pascals(minimum_pressure.value + step * index as f64);
        match water_dew_point_multiphase(components, z, p, cubic) {
            Ok(point) => {
                last = point.temperature.value;
                warnings.extend(point.warnings.iter().cloned());
            }
            Err(_) => failed += 1,
        }
        temperatures.push(last);
        pressures.push(p.value);
    }

    Ok(WaterDewPointLineResult {
        temperature: temperatures,
        pressure: pressures,
        failed_points: failed,
        warnings,
    })
}
