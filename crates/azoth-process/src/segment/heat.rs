//! The interphase heat step: the Chilton-Colburn coefficient and the heat it moves.

use super::film::{clamp, finite_positive};

/// `calculateVolumetricHeatTransferCoefficient`: the Chilton-Colburn analogy.
///
/// **The fractional exponent is guarded, and the guard is the class's.** `(Sc/Pr)^(2/3)` is a
/// partial function at a non-positive base, so both numbers are required positive and finite
/// before it is evaluated; a state that fails either answers `0.0` rather than a `NaN` carried
/// into the segment's heat balance.
#[allow(clippy::too_many_arguments)] // The signature is the class's own argument list.
pub fn volumetric_heat_transfer_coefficient(
    mass_transfer_coefficient: f64,
    density: f64,
    heat_capacity: f64,
    viscosity: f64,
    diffusivity: f64,
    conductivity: f64,
    model_none: bool,
    correction: f64,
) -> f64 {
    let none = |value: f64| !(value > 0.0 && value.is_finite());
    if model_none
        || none(mass_transfer_coefficient)
        || none(density)
        || none(heat_capacity)
        || none(viscosity)
        || none(diffusivity)
        || none(conductivity)
    {
        return 0.0;
    }
    let prandtl = heat_capacity * viscosity / conductivity;
    let schmidt = viscosity / (density * diffusivity);
    if none(prandtl) || none(schmidt) {
        return 0.0;
    }
    let analogy = (schmidt / prandtl).powf(2.0 / 3.0);
    mass_transfer_coefficient * density * heat_capacity * analogy * correction
}

/// `combineHeatTransferCoefficients`: the two sides in series.
pub fn combine_heat_transfer_coefficients(gas: f64, liquid: f64) -> f64 {
    if !(gas > 0.0 && gas.is_finite() && liquid > 0.0 && liquid.is_finite()) {
        return 0.0;
    }
    1.0 / (1.0 / gas + 1.0 / liquid)
}

/// `calculateInterfaceTemperature`: the resistance-weighted mean of the two bulk temperatures.
pub fn interface_temperature(gas_t: f64, liquid_t: f64, gas: f64, liquid: f64) -> f64 {
    if !(gas > 0.0 && gas.is_finite() && liquid > 0.0 && liquid.is_finite()) {
        return 0.5 * (gas_t + liquid_t);
    }
    (gas * gas_t + liquid * liquid_t) / (gas + liquid)
}

/// `heatCapacityRate`: the phase's mass flow times its heat capacity, W/K.
pub fn heat_capacity_rate(mass_flow: f64, heat_capacity: f64) -> f64 {
    finite_positive(mass_flow, 0.0) * heat_capacity
}

/// What the heat step did to the two temperatures.
pub struct HeatStep {
    /// The rate moved, in W, positive from gas to liquid.
    pub rate: f64,
    /// The gas temperature after the step, K.
    pub gas_temperature: f64,
    /// The liquid temperature after it, K.
    pub liquid_temperature: f64,
}

/// `applyInterphaseHeatTransfer`.
///
/// **The `1.0 K` floor is the class's and it is reachable**: each side is moved by
/// `Q / (m cp)` and floored at one kelvin, so a cold liquid with a small capacity rate stops
/// where the class stops rather than going below a temperature this library would refuse.
#[allow(clippy::too_many_arguments)] // The signature is the class's own argument list.
pub fn apply_interphase_heat_transfer(
    gas_temperature: f64,
    liquid_temperature: f64,
    gas_heat_capacity_rate: f64,
    liquid_heat_capacity_rate: f64,
    overall_coefficient: f64,
    segment_volume: f64,
    model_none: bool,
    max_fraction: f64,
) -> HeatStep {
    let unchanged = |rate: f64| HeatStep {
        rate,
        gas_temperature,
        liquid_temperature,
    };
    if model_none || !(overall_coefficient > 0.0 && overall_coefficient.is_finite()) {
        return unchanged(0.0);
    }
    let difference = gas_temperature - liquid_temperature;
    if difference.abs() < 1.0e-12 {
        return unchanged(0.0);
    }
    if !(gas_heat_capacity_rate > 0.0
        && gas_heat_capacity_rate.is_finite()
        && liquid_heat_capacity_rate > 0.0
        && liquid_heat_capacity_rate.is_finite())
    {
        return unchanged(0.0);
    }
    let rate = overall_coefficient * segment_volume * difference;
    let maximum =
        gas_heat_capacity_rate.min(liquid_heat_capacity_rate) * difference.abs() * max_fraction;
    let rate = rate.signum() * rate.abs().min(maximum);
    if rate == 0.0 {
        return unchanged(0.0);
    }
    HeatStep {
        rate,
        gas_temperature: clamp(
            gas_temperature - rate / gas_heat_capacity_rate,
            1.0,
            f64::MAX,
        ),
        liquid_temperature: clamp(
            liquid_temperature + rate / liquid_heat_capacity_rate,
            1.0,
            f64::MAX,
        ),
    }
}
