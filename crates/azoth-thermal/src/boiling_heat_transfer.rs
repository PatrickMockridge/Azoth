//! `thermal.boiling_heat_transfer` - Chen and Gungor-Winterton flow boiling.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/BoilingHeatTransfer.java` at
//! NeqSim `f0c7436c69`, composed as `IncrementalZoneAnalysis` uses it: the average
//! Gungor-Winterton coefficient over a quality range, with the Chen local, the two factors
//! and the horizontal correction beside it.
//!
//! Spec: `specs/calcs/thermal/boiling_heat_transfer.toml`.

use azoth_core::units::{
    DynamicViscosity, HeatFluxDensity, Length, MassDensity, MassFlux, Pressure,
    SpecificHeatCapacity, SurfaceTension, TemperatureInterval, ThermalConductivity,
    watts_per_square_meter_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::BoilingHeatTransferResult;
use crate::spec_gen;

/// The class's boiling quality floor.
pub const MIN_QUALITY: f64 = 0.001;
/// The class's boiling quality ceiling.
pub const MAX_QUALITY: f64 = 0.999;
/// The class's gravity constant in the horizontal Froude number.
const GRAVITY: f64 = 9.81;

/// The composed boiling request, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct BoilingRequest {
    pub mass_flux: f64,
    pub tube_inner_diameter: f64,
    pub liquid_density: f64,
    pub vapor_density: f64,
    pub liquid_viscosity: f64,
    pub vapor_viscosity: f64,
    pub liquid_cp: f64,
    pub liquid_conductivity: f64,
    pub surface_tension: f64,
    pub heat_of_vaporization: f64,
    pub wall_superheat: f64,
    pub saturation_pressure_difference: f64,
    pub heat_flux: f64,
    pub vapor_quality: f64,
    pub quality_in: f64,
    pub quality_out: f64,
    pub intervals: Option<f64>,
    pub horizontal: Option<bool>,
}

/// The Lockhart-Martinelli parameter, answering 1 at an endpoint.
#[must_use]
pub fn martinelli_parameter(
    x: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    vapor_viscosity: f64,
) -> f64 {
    if x <= 0.0
        || x >= 1.0
        || vapor_density <= 0.0
        || liquid_density <= 0.0
        || vapor_viscosity <= 0.0
    {
        return 1.0;
    }
    ((1.0 - x) / x).powf(0.9)
        * (vapor_density / liquid_density).powf(0.5)
        * (liquid_viscosity / vapor_viscosity).powf(0.1)
}

/// The Chen enhancement factor.
#[must_use]
pub fn chen_enhancement_factor(x_tt: f64) -> f64 {
    if x_tt <= 0.0 {
        return 1.0;
    }
    let inverse = 1.0 / x_tt;
    if inverse <= 0.1 {
        return 1.0;
    }
    2.35 * (0.213 + inverse).powf(0.736)
}

/// The Chen suppression factor.
#[must_use]
pub fn chen_suppression_factor(re_two_phase: f64) -> f64 {
    if re_two_phase <= 0.0 {
        return 1.0;
    }
    1.0 / (1.0 + 2.53e-6 * re_two_phase.powf(1.17))
}

/// The Chen local coefficient.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn chen_htc(
    mass_flux: f64,
    vapor_quality: f64,
    tube_id: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    vapor_viscosity: f64,
    liquid_cp: f64,
    liquid_conductivity: f64,
    surface_tension: f64,
    heat_of_vaporization: f64,
    wall_superheat: f64,
    saturation_pressure_difference: f64,
) -> f64 {
    if liquid_viscosity <= 0.0 || liquid_conductivity <= 0.0 || tube_id <= 0.0 || mass_flux <= 0.0 {
        return 0.0;
    }
    let x = vapor_quality.clamp(MIN_QUALITY, MAX_QUALITY);
    let re_liquid = mass_flux * (1.0 - x) * tube_id / liquid_viscosity;
    let pr_liquid = liquid_cp * liquid_viscosity / liquid_conductivity;
    let h_liquid =
        0.023 * re_liquid.max(1.0).powf(0.8) * pr_liquid.powf(0.4) * liquid_conductivity / tube_id;
    let x_tt = martinelli_parameter(
        x,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        vapor_viscosity,
    );
    let enhancement = chen_enhancement_factor(x_tt);

    let mut nucleate = 0.0;
    if wall_superheat > 0.0
        && saturation_pressure_difference > 0.0
        && surface_tension > 0.0
        && heat_of_vaporization > 0.0
    {
        let numerator =
            liquid_conductivity.powf(0.79) * liquid_cp.powf(0.45) * liquid_density.powf(0.49);
        let denominator = surface_tension.powf(0.5)
            * liquid_viscosity.powf(0.29)
            * heat_of_vaporization.powf(0.24)
            * vapor_density.powf(0.24);
        if denominator > 0.0 {
            nucleate = 0.00122
                * (numerator / denominator)
                * wall_superheat.powf(0.24)
                * saturation_pressure_difference.powf(0.75);
        }
    }
    let re_two_phase = re_liquid * enhancement.powf(1.25);
    let suppression = chen_suppression_factor(re_two_phase);
    h_liquid * enhancement + nucleate * suppression
}

/// The Gungor-Winterton local coefficient.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn gungor_winterton_htc(
    mass_flux: f64,
    vapor_quality: f64,
    tube_id: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    liquid_cp: f64,
    liquid_conductivity: f64,
    heat_flux: f64,
    heat_of_vaporization: f64,
) -> f64 {
    if liquid_viscosity <= 0.0 || liquid_conductivity <= 0.0 || tube_id <= 0.0 || mass_flux <= 0.0 {
        return 0.0;
    }
    let x = vapor_quality.clamp(MIN_QUALITY, MAX_QUALITY);
    let re_liquid = mass_flux * (1.0 - x) * tube_id / liquid_viscosity;
    let pr_liquid = liquid_cp * liquid_viscosity / liquid_conductivity;
    let h_liquid =
        0.023 * re_liquid.max(1.0).powf(0.8) * pr_liquid.powf(0.4) * liquid_conductivity / tube_id;
    let boiling_number = if heat_of_vaporization > 0.0 && mass_flux > 0.0 {
        heat_flux / (mass_flux * heat_of_vaporization)
    } else {
        0.0
    };
    let convective = if x < MAX_QUALITY && vapor_density > 0.0 {
        1.12 * (x / (1.0 - x)).powf(0.75) * (liquid_density / vapor_density).powf(0.41)
    } else {
        0.0
    };
    let nucleate = 3000.0 * boiling_number.max(0.0).powf(0.86);
    h_liquid * (1.0 + nucleate + convective)
}

/// The Gungor-Winterton coefficient with the horizontal stratified correction.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn gungor_winterton_corrected_htc(
    mass_flux: f64,
    vapor_quality: f64,
    tube_id: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    liquid_cp: f64,
    liquid_conductivity: f64,
    heat_flux: f64,
    heat_of_vaporization: f64,
    horizontal: bool,
) -> f64 {
    let mut h = gungor_winterton_htc(
        mass_flux,
        vapor_quality,
        tube_id,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        liquid_cp,
        liquid_conductivity,
        heat_flux,
        heat_of_vaporization,
    );
    if horizontal && liquid_density > 0.0 {
        let x = vapor_quality.clamp(MIN_QUALITY, MAX_QUALITY);
        let liquid_flux = mass_flux * (1.0 - x);
        let liquid_velocity = liquid_flux / liquid_density;
        let froude = liquid_velocity * liquid_velocity / (GRAVITY * tube_id);
        if froude < 0.05 {
            h *= froude.powf(0.1 - 2.0 * froude);
        }
    }
    h
}

/// Simpson's-rule average of the uncorrected Gungor-Winterton coefficient.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn average_htc(
    mass_flux: f64,
    tube_id: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    liquid_cp: f64,
    liquid_conductivity: f64,
    heat_flux: f64,
    heat_of_vaporization: f64,
    quality_in: f64,
    quality_out: f64,
    intervals: i32,
) -> f64 {
    if mass_flux <= 0.0 || liquid_viscosity <= 0.0 || tube_id <= 0.0 {
        return 0.0;
    }
    let x_in = quality_in.clamp(MIN_QUALITY, MAX_QUALITY);
    let x_out = quality_out.clamp(MIN_QUALITY, MAX_QUALITY);
    if (x_in - x_out).abs() < 1e-10 {
        return gungor_winterton_htc(
            mass_flux,
            x_in,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            liquid_cp,
            liquid_conductivity,
            heat_flux,
            heat_of_vaporization,
        );
    }
    let mut n = intervals.max(4);
    if n % 2 != 0 {
        n += 1;
    }
    let dx = (x_out - x_in) / f64::from(n);
    let local = |x: f64| {
        gungor_winterton_htc(
            mass_flux,
            x,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            liquid_cp,
            liquid_conductivity,
            heat_flux,
            heat_of_vaporization,
        )
    };
    let mut sum = local(x_in) + local(x_out);
    for i in 1..n {
        let x = x_in + f64::from(i) * dx;
        let h = local(x);
        sum += if i % 2 == 0 { 2.0 * h } else { 4.0 * h };
    }
    (sum * dx / 3.0 / (x_out - x_in)).abs()
}

/// The composed boiling calculation.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a non-integral interval count, and
/// [`AzothError::OutOfRange`] for the spec's own bounds.
pub fn rate(request: &BoilingRequest) -> Result<BoilingHeatTransferResult> {
    let spec = &spec_gen::BOILING_HEAT_TRANSFER_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "mass_flux" => Some(request.mass_flux),
            "tube_inner_diameter" => Some(request.tube_inner_diameter),
            "liquid_density" => Some(request.liquid_density),
            "vapor_density" => Some(request.vapor_density),
            "liquid_viscosity" => Some(request.liquid_viscosity),
            "vapor_viscosity" => Some(request.vapor_viscosity),
            "liquid_cp" => Some(request.liquid_cp),
            "liquid_conductivity" => Some(request.liquid_conductivity),
            "surface_tension" => Some(request.surface_tension),
            "heat_of_vaporization" => Some(request.heat_of_vaporization),
            "wall_superheat" => Some(request.wall_superheat),
            "saturation_pressure_difference" => Some(request.saturation_pressure_difference),
            "heat_flux" => Some(request.heat_flux),
            _ => None,
        },
        &mut warnings,
    )?;
    let intervals = match request.intervals {
        Some(value) => {
            if !value.is_finite() || (value - value.round()).abs() > 0.0 {
                return Err(AzothError::invalid_input(
                    "intervals",
                    format!("{value} is not a whole interval count"),
                ));
            }
            value.round() as i32
        }
        None => 20,
    };
    let horizontal = request.horizontal.unwrap_or(false);
    // **The published factors are the class's public helpers at the stated quality**, while
    // the Chen local clamps the quality first; at a clamp the two differ, and the helper's
    // raw-quality answer is the one the class's own methods give a caller.
    let x_tt = martinelli_parameter(
        request.vapor_quality,
        request.liquid_density,
        request.vapor_density,
        request.liquid_viscosity,
        request.vapor_viscosity,
    );
    let x = request.vapor_quality.clamp(MIN_QUALITY, MAX_QUALITY);
    let enhancement = chen_enhancement_factor(x_tt);
    let re_liquid =
        request.mass_flux * (1.0 - x) * request.tube_inner_diameter / request.liquid_viscosity;
    Ok(BoilingHeatTransferResult {
        chen_htc: watts_per_square_meter_kelvin(chen_htc(
            request.mass_flux,
            request.vapor_quality,
            request.tube_inner_diameter,
            request.liquid_density,
            request.vapor_density,
            request.liquid_viscosity,
            request.vapor_viscosity,
            request.liquid_cp,
            request.liquid_conductivity,
            request.surface_tension,
            request.heat_of_vaporization,
            request.wall_superheat,
            request.saturation_pressure_difference,
        )),
        gungor_winterton_htc: watts_per_square_meter_kelvin(gungor_winterton_htc(
            request.mass_flux,
            request.vapor_quality,
            request.tube_inner_diameter,
            request.liquid_density,
            request.vapor_density,
            request.liquid_viscosity,
            request.liquid_cp,
            request.liquid_conductivity,
            request.heat_flux,
            request.heat_of_vaporization,
        )),
        gungor_winterton_corrected_htc: watts_per_square_meter_kelvin(
            gungor_winterton_corrected_htc(
                request.mass_flux,
                request.vapor_quality,
                request.tube_inner_diameter,
                request.liquid_density,
                request.vapor_density,
                request.liquid_viscosity,
                request.liquid_cp,
                request.liquid_conductivity,
                request.heat_flux,
                request.heat_of_vaporization,
                horizontal,
            ),
        ),
        average_htc: watts_per_square_meter_kelvin(average_htc(
            request.mass_flux,
            request.tube_inner_diameter,
            request.liquid_density,
            request.vapor_density,
            request.liquid_viscosity,
            request.liquid_cp,
            request.liquid_conductivity,
            request.heat_flux,
            request.heat_of_vaporization,
            request.quality_in,
            request.quality_out,
            intervals,
        )),
        martinelli_parameter: x_tt,
        chen_enhancement_factor: enhancement,
        chen_suppression_factor: chen_suppression_factor(re_liquid * enhancement.powf(1.25)),
        warnings,
    })
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn boiling_heat_transfer(
    mass_flux: MassFlux,
    tube_inner_diameter: Length,
    liquid_density: MassDensity,
    vapor_density: MassDensity,
    liquid_viscosity: DynamicViscosity,
    vapor_viscosity: DynamicViscosity,
    liquid_cp: SpecificHeatCapacity,
    liquid_conductivity: ThermalConductivity,
    surface_tension: SurfaceTension,
    heat_of_vaporization: f64,
    wall_superheat: TemperatureInterval,
    saturation_pressure_difference: Pressure,
    heat_flux: HeatFluxDensity,
    vapor_quality: f64,
    quality_in: f64,
    quality_out: f64,
    intervals: Option<f64>,
    horizontal: Option<bool>,
) -> Result<BoilingHeatTransferResult> {
    rate(&BoilingRequest {
        mass_flux: mass_flux.value,
        tube_inner_diameter: tube_inner_diameter.value,
        liquid_density: liquid_density.value,
        vapor_density: vapor_density.value,
        liquid_viscosity: liquid_viscosity.value,
        vapor_viscosity: vapor_viscosity.value,
        liquid_cp: liquid_cp.value,
        liquid_conductivity: liquid_conductivity.value,
        surface_tension: surface_tension.value,
        heat_of_vaporization,
        wall_superheat: wall_superheat.value,
        saturation_pressure_difference: saturation_pressure_difference.value,
        heat_flux: heat_flux.value,
        vapor_quality,
        quality_in,
        quality_out,
        intervals,
        horizontal,
    })
}
