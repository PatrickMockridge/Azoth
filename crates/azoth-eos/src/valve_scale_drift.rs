//! `eos.valve_scale_drift` - NeqSim's `ValveScaleDrift`.
//!
//! A uniform radial deposit on a circular trim port: the open-area fraction, the fouling
//! fraction and the effective Kv at an elapsed time, with the time to plug the port and the
//! time for the valve opening to pin at 100 percent.

use azoth_core::Result;
use azoth_core::apply_checks;
use azoth_core::units::{Length, Time, Velocity, meters, seconds};

use crate::model_gen;
use crate::results::ValveScaleDriftResult;

/// The fouling fraction's cap, as the class's `getFoulingFraction` states it.
const FOULING_CAP: f64 = 0.999_999;

/// The trim fouling a deposit growth rate produces at an elapsed time.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
///
/// # Example
/// ```
/// use azoth_core::units::{meters, meters_per_second, seconds};
/// use azoth_eos::valve_scale_drift::valve_scale_drift;
///
/// let r = valve_scale_drift(
///     meters(0.05),
///     meters_per_second(3.168808781402895e-11), // 1 mm/yr
///     seconds(2_592_000.0),                     // 30 days
///     50.0,
///     100.0,
/// )?;
/// assert!((r.open_area_fraction - 0.9934399521016661).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn valve_scale_drift(
    port_diameter: Length,
    growth_rate: Velocity,
    elapsed_time: Time,
    clean_opening_percent: f64,
    valve_kv: f64,
) -> Result<ValveScaleDriftResult> {
    let spec = &model_gen::VALVE_SCALE_DRIFT_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "port_diameter" => Some(port_diameter.value),
            "growth_rate" => Some(growth_rate.value),
            "elapsed_time" => Some(elapsed_time.value),
            "clean_opening_percent" => Some(clean_opening_percent),
            "valve_kv" => Some(valve_kv),
            _ => None,
        },
        &mut warnings,
    )?;

    // The setters' clamps: a negative growth rate is zero, and a negative elapsed step is
    // ignored, which a total states as zero.
    let growth = growth_rate.value.max(0.0);
    let elapsed = elapsed_time.value.max(0.0);
    let port = port_diameter.value;
    let thickness = growth * elapsed;

    let open_area = if port <= 0.0 {
        1.0
    } else {
        let open_diameter = port - 2.0 * thickness;
        if open_diameter <= 0.0 {
            0.0
        } else {
            let ratio = open_diameter / port;
            ratio * ratio
        }
    };

    let fouling = (1.0 - open_area).clamp(0.0, FOULING_CAP);

    let time_to_plug = if growth <= 0.0 || port <= 0.0 {
        f64::INFINITY
    } else {
        ((port / 2.0 - thickness).max(0.0) / growth).max(0.0)
    };

    let opening_percent = if open_area <= 0.0 {
        100.0
    } else {
        (clean_opening_percent / open_area).min(100.0)
    };

    let time_to_pin = if growth <= 0.0 || port <= 0.0 {
        f64::INFINITY
    } else if clean_opening_percent >= 100.0 {
        0.0
    } else if clean_opening_percent <= 0.0 {
        f64::INFINITY
    } else {
        let pin_thickness = port / 2.0 * (1.0 - (clean_opening_percent / 100.0).sqrt());
        let remaining = pin_thickness - thickness;
        if remaining <= 0.0 {
            0.0
        } else {
            remaining / growth
        }
    };

    Ok(ValveScaleDriftResult {
        deposit_thickness: meters(thickness),
        elapsed_time: seconds(elapsed),
        open_area_fraction: open_area,
        fouling_fraction: fouling,
        effective_kv: valve_kv * open_area,
        time_to_plug: seconds(time_to_plug),
        opening_percent,
        time_to_pin: seconds(time_to_pin),
        warnings,
    })
}
