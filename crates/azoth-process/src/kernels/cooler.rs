//! `unit_ops.cooler` - `unit_ops.heater`'s arithmetic, with the sign a cooler implies.

use azoth_core::Result;
use azoth_core::units::{Power, Pressure, ThermodynamicTemperature, kelvins, pascals};
use azoth_core::AzothError;

use crate::kernels::heater::{HeaterOutcome, heater};
use crate::stream::Stream;

/// Cool a stream to a stated temperature, or by a stated duty.
///
/// **This is a delegation and the delegation is the port.** `Cooler extends Heater` and
/// overrides `runTransient`, the mechanical design and a handful of getters - **not `run`**,
/// which is the whole of the steady-state arithmetic. Its own 258 lines are a dynamic
/// valve/NTU temperature-control model that does not touch the outlet record, and the
/// caption to that claim is a measurement rather than a reading: `ProcessProbe cooler` runs
/// the same six rows as `ProcessProbe heater` through this class and the two captures are
/// **byte-identical**.
///
/// Writing the arithmetic out again here would be a second implementation of one machine's
/// physics and the two could then disagree - which is the same reason `GasScrubber`'s
/// stream-side kernel is `Separator`'s and not a copy of it.
pub fn cooler(
    feed: &Stream,
    outlet_temperature: Option<ThermodynamicTemperature>,
    duty: Option<Power>,
    pressure_drop: Option<Pressure>,
) -> Result<HeaterOutcome> {
    heater(feed, outlet_temperature, duty, pressure_drop)
}

/// `Cooler.configureDynamicTemperatureControl`'s design point and one `runTransient` step.
#[derive(Debug, Clone, Copy)]
pub struct CoolerTransientSettings {
    /// The step interval, s.
    pub time_step_s: f64,
    /// The design process mass flow, kg/h.
    pub design_mass_flow_kg_per_hour: f64,
    /// The design inlet temperature, K.
    pub design_inlet_temperature_k: f64,
    /// The design outlet temperature, K.
    pub design_outlet_temperature_k: f64,
    /// The cooling medium's supply temperature, K.
    pub cooling_medium_temperature_k: f64,
    /// The valve opening at the design point, percent.
    pub design_valve_opening_percent: f64,
    /// The utility-valve actuator's first-order constant, s; zero is instantaneous.
    pub actuator_time_constant_s: f64,
    /// The process side's first-order constant, s; zero is instantaneous.
    pub thermal_time_constant_s: f64,
    /// The carried valve opening, percent.
    pub cooling_valve_opening_percent: f64,
    /// The controller's response for this step; absent is the step-initial opening.
    pub requested_valve_opening_percent: Option<f64>,
    /// The carried outlet temperature; absent or NaN seeds to the equilibrium.
    pub dynamic_outlet_temperature_k: Option<f64>,
    /// The pressure loss across the cooler, Pa.
    pub pressure_drop_pa: Option<f64>,
}

/// What one dynamic-cooler step answered.
#[derive(Debug, Clone)]
pub struct CoolerTransientOutcome {
    /// The outlet record: the steady heater's run at the computed temperature.
    pub outlet: Stream,
    /// The duty it moved, W.
    pub duty: Power,
    /// The post-step valve opening, percent.
    pub cooling_valve_opening_percent: f64,
    /// The equilibrium outlet the thermal lag moved toward, K.
    pub equilibrium_outlet_temperature: ThermodynamicTemperature,
}

/// `Cooler.runTransient`'s configured branch, one step.
///
/// The valve moves by `dt / (tau_actuator + dt)` toward the requested opening, an NTU scaled
/// by opening and inverse mass flow gives the equilibrium outlet, and a second first-order lag
/// moves the outlet toward it; the outlet is then the steady `Heater.run` at that temperature.
/// The class's per-step UUID cache and its controller are session state, outside this step.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a non-positive side flow, a design point whose temperatures
/// are not `inlet > outlet > medium`, a non-positive design mass flow or valve opening, a
/// negative time constant, a non-finite step, or whatever the heater refuses.
pub fn cooler_transient(
    feed: &Stream,
    settings: &CoolerTransientSettings,
) -> Result<CoolerTransientOutcome> {
    if !feed.n.is_finite() || feed.n <= 0.0 {
        return Err(AzothError::invalid_input(
            "inlet_n",
            format!("a unit operation needs a positive molar flow, and {} mol/s is not", feed.n),
        ));
    }
    let dt = settings.time_step_s;
    if !dt.is_finite() || dt < 0.0 {
        return Err(AzothError::invalid_input(
            "time_step",
            format!("a forward-Euler interval is finite and non-negative, and {dt} s is not"),
        ));
    }
    let design_mass = settings.design_mass_flow_kg_per_hour;
    if !design_mass.is_finite() || design_mass <= 0.0 {
        return Err(AzothError::invalid_input(
            "design_mass_flow",
            format!("the design point needs a positive mass flow, and {design_mass} kg/h is not"),
        ));
    }
    let design_inlet = settings.design_inlet_temperature_k;
    let design_outlet = settings.design_outlet_temperature_k;
    let medium = settings.cooling_medium_temperature_k;
    if !design_inlet.is_finite() || !design_outlet.is_finite() || !medium.is_finite() {
        return Err(AzothError::invalid_input(
            "design_inlet_temperature",
            "the design point's three temperatures are finite",
        ));
    }
    if !(design_inlet > design_outlet && design_outlet > medium) {
        return Err(AzothError::invalid_input(
            "design_outlet_temperature",
            format!(
                "the design point is inlet > outlet > medium, and {design_inlet} / {design_outlet} / {medium} K is not"
            ),
        ));
    }
    let design_opening = settings.design_valve_opening_percent;
    if !design_opening.is_finite() || design_opening <= 0.0 || design_opening > 100.0 {
        return Err(AzothError::invalid_input(
            "design_valve_opening",
            format!("the design opening is in (0, 100] percent, and {design_opening} is not"),
        ));
    }
    for (field, value) in [
        ("actuator_time_constant", settings.actuator_time_constant_s),
        ("thermal_time_constant", settings.thermal_time_constant_s),
    ] {
        if !value.is_finite() || value < 0.0 {
            return Err(AzothError::invalid_input(
                field,
                format!("a time constant is finite and non-negative, and {value} s is not"),
            ));
        }
    }
    if !settings.cooling_valve_opening_percent.is_finite() {
        return Err(AzothError::invalid_input(
            "cooling_valve_opening",
            "the carried opening is finite",
        ));
    }
    if let Some(requested) = settings.requested_valve_opening_percent {
        if !requested.is_finite() {
            return Err(AzothError::invalid_input(
                "requested_valve_opening",
                "the controller's requested opening is finite",
            ));
        }
    }
    if let Some(carried) = settings.dynamic_outlet_temperature_k {
        if carried.is_infinite() {
            return Err(AzothError::invalid_input(
                "dynamic_outlet_temperature",
                format!("a carried temperature is finite or NaN, and {carried} is neither"),
            ));
        }
    }

    let design_ntu = -((design_outlet - medium) / (design_inlet - medium)).ln();
    let mut opening = settings.cooling_valve_opening_percent.clamp(0.0, 100.0);
    let requested = settings
        .requested_valve_opening_percent
        .unwrap_or(opening)
        .clamp(0.0, 100.0);
    let actuator = if settings.actuator_time_constant_s > 0.0 {
        dt / (settings.actuator_time_constant_s + dt)
    } else {
        1.0
    };
    opening += actuator * (requested - opening);

    let mass_flow = (feed.mass_flow()? * 3600.0).max(1.0e-12);
    let ntu = (design_ntu * opening / design_opening * design_mass / mass_flow).clamp(0.0, 20.0);
    let inlet_t = feed.t.value;
    let equilibrium = if inlet_t > medium {
        medium + (inlet_t - medium) * (-ntu).exp()
    } else {
        inlet_t
    };
    // **The absent state is the design outlet**, which is the state
    // `configureDynamicTemperatureControl` writes before the first step; the class's NaN
    // branch is the equilibrium only where no run established one.
    let previous = settings
        .dynamic_outlet_temperature_k
        .filter(|value| !value.is_nan())
        .unwrap_or(design_outlet);
    let thermal = if settings.thermal_time_constant_s > 0.0 {
        dt / (settings.thermal_time_constant_s + dt)
    } else {
        1.0
    };
    let outlet_t = (previous + thermal * (equilibrium - previous)).clamp(medium.min(inlet_t), inlet_t);

    let outcome = heater(
        feed,
        Some(kelvins(outlet_t)),
        None,
        settings.pressure_drop_pa.map(pascals),
    )?;
    Ok(CoolerTransientOutcome {
        outlet: outcome.outlet,
        duty: outcome.duty,
        cooling_valve_opening_percent: opening,
        equilibrium_outlet_temperature: kelvins(equilibrium),
    })
}
