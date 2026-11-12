//! `thermal.two_phase_pressure_drop` - Friedel and Muller-Steinhagen-Heck.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/TwoPhasePressureDrop.java` at
//! NeqSim `118d9aff`: the Friedel frictional gradient and its total and quality-averaged
//! drops, the Muller-Steinhagen-Heck gradient, and the gravitational and acceleration terms
//! of a homogeneous two-phase mixture.
//!
//! Spec: `specs/calcs/thermal/two_phase_pressure_drop.toml`.

use azoth_core::units::{DynamicViscosity, Length, MassDensity, MassFlux, SurfaceTension, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::TwoPhasePressureDropResult;
use crate::spec_gen;

/// The class's gravity constant.
const GRAVITY: f64 = 9.81;

/// The composed request, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct TwoPhasePressureDropRequest {
    pub mass_flux: f64,
    pub vapor_quality: f64,
    pub tube_inner_diameter: f64,
    pub tube_length: f64,
    pub liquid_density: f64,
    pub vapor_density: f64,
    pub liquid_viscosity: f64,
    pub vapor_viscosity: f64,
    pub surface_tension: f64,
    pub quality_in: f64,
    pub quality_out: f64,
    pub intervals: Option<f64>,
}

/// The Fanning friction factor for a smooth tube.
#[must_use]
pub fn fanning_friction(reynolds: f64) -> f64 {
    if reynolds <= 0.0 {
        return 0.0;
    }
    if reynolds < 2300.0 {
        return 16.0 / reynolds;
    }
    0.079 / reynolds.powf(0.25)
}

/// The Friedel two-phase frictional gradient.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn friedel_gradient(
    mass_flux: f64,
    vapor_quality: f64,
    tube_id: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    vapor_viscosity: f64,
    surface_tension: f64,
) -> f64 {
    if mass_flux <= 0.0
        || tube_id <= 0.0
        || liquid_density <= 0.0
        || vapor_density <= 0.0
        || liquid_viscosity <= 0.0
        || vapor_viscosity <= 0.0
    {
        return 0.0;
    }

    let x = vapor_quality.clamp(0.0, 1.0);
    let re_lo = mass_flux * tube_id / liquid_viscosity;
    let f_lo = fanning_friction(re_lo);
    let dpdz_lo = 4.0 * f_lo * mass_flux * mass_flux / (2.0 * liquid_density * tube_id);
    if x < 1e-10 {
        return dpdz_lo;
    }

    let re_vo = mass_flux * tube_id / vapor_viscosity;
    let f_vo = fanning_friction(re_vo);

    let e = (1.0 - x).powi(2) + x * x * (liquid_density * f_vo) / (vapor_density * f_lo);
    let f = x.powf(0.78) * (1.0 - x).powf(0.224);
    let h = (liquid_density / vapor_density).powf(0.91)
        * (vapor_viscosity / liquid_viscosity).powf(0.19)
        * (1.0 - vapor_viscosity / liquid_viscosity).powf(0.7);

    let rho_h = 1.0 / (x / vapor_density + (1.0 - x) / liquid_density);
    let fr = mass_flux * mass_flux / (rho_h * rho_h * GRAVITY * tube_id);
    let we = if surface_tension > 0.0 {
        mass_flux * mass_flux * tube_id / (rho_h * surface_tension)
    } else {
        1e6
    };
    let fr_we = fr.max(1e-10).powf(0.045) * we.max(1e-10).powf(0.035);
    let phi_lo_sq = e + 3.24 * f * h / fr_we;
    phi_lo_sq * dpdz_lo
}

/// The Friedel drop over a tube length at a constant quality.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn friedel_pressure_drop(
    mass_flux: f64,
    vapor_quality: f64,
    tube_id: f64,
    tube_length: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    vapor_viscosity: f64,
    surface_tension: f64,
) -> f64 {
    friedel_gradient(
        mass_flux,
        vapor_quality,
        tube_id,
        liquid_density,
        vapor_density,
        liquid_viscosity,
        vapor_viscosity,
        surface_tension,
    ) * tube_length
}

/// Simpson's-rule average of the Friedel gradient over a quality range, times the length.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn friedel_average_pressure_drop(
    mass_flux: f64,
    quality_in: f64,
    quality_out: f64,
    tube_id: f64,
    tube_length: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    vapor_viscosity: f64,
    surface_tension: f64,
    intervals: i32,
) -> f64 {
    if tube_length <= 0.0 {
        return 0.0;
    }
    let x_in = quality_in.clamp(0.0, 1.0);
    let x_out = quality_out.clamp(0.0, 1.0);
    if (x_in - x_out).abs() < 1e-10 {
        return friedel_pressure_drop(
            mass_flux,
            x_in,
            tube_id,
            tube_length,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            vapor_viscosity,
            surface_tension,
        );
    }
    let mut n = intervals.max(4);
    if n % 2 != 0 {
        n += 1;
    }
    let dx = (x_out - x_in) / f64::from(n);
    let local = |x: f64| {
        friedel_gradient(
            mass_flux,
            x,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            vapor_viscosity,
            surface_tension,
        )
    };
    let mut sum = local(x_in) + local(x_out);
    for i in 1..n {
        let x = x_in + f64::from(i) * dx;
        let gradient = local(x);
        sum += if i % 2 == 0 {
            2.0 * gradient
        } else {
            4.0 * gradient
        };
    }
    (sum * dx / 3.0 / (x_out - x_in)).abs() * tube_length
}

/// The Muller-Steinhagen-Heck two-phase frictional gradient.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn muller_steinhagen_heck_gradient(
    mass_flux: f64,
    vapor_quality: f64,
    tube_id: f64,
    liquid_density: f64,
    vapor_density: f64,
    liquid_viscosity: f64,
    vapor_viscosity: f64,
) -> f64 {
    if mass_flux <= 0.0
        || tube_id <= 0.0
        || liquid_density <= 0.0
        || vapor_density <= 0.0
        || liquid_viscosity <= 0.0
        || vapor_viscosity <= 0.0
    {
        return 0.0;
    }
    let x = vapor_quality.clamp(0.0, 1.0);
    let re_lo = mass_flux * tube_id / liquid_viscosity;
    let f_lo = fanning_friction(re_lo);
    let a = 4.0 * f_lo * mass_flux * mass_flux / (2.0 * liquid_density * tube_id);
    let re_vo = mass_flux * tube_id / vapor_viscosity;
    let f_vo = fanning_friction(re_vo);
    let b = 4.0 * f_vo * mass_flux * mass_flux / (2.0 * vapor_density * tube_id);
    let g_factor = a + 2.0 * (b - a) * x;
    g_factor * (1.0 - x).powf(1.0 / 3.0) + b * x * x * x
}

/// The homogeneous void fraction at a quality.
fn homogeneous_void_fraction(x: f64, liquid_density: f64, vapor_density: f64) -> f64 {
    if x < 1e-10 {
        0.0
    } else if x > 1.0 - 1e-10 {
        1.0
    } else {
        1.0 / (1.0 + (1.0 - x) / x * vapor_density / liquid_density)
    }
}

/// The gravitational gradient of the homogeneous mixture.
#[must_use]
pub fn gravitational_gradient(vapor_quality: f64, liquid_density: f64, vapor_density: f64) -> f64 {
    if liquid_density <= 0.0 || vapor_density <= 0.0 {
        return 0.0;
    }
    let x = vapor_quality.clamp(0.0, 1.0);
    let alpha = homogeneous_void_fraction(x, liquid_density, vapor_density);
    let rho_tp = liquid_density * (1.0 - alpha) + vapor_density * alpha;
    rho_tp * GRAVITY
}

/// The momentum-flux term at a quality.
fn momentum_flux(quality: f64, liquid_density: f64, vapor_density: f64) -> f64 {
    let x = quality.clamp(0.0, 1.0);
    if x < 1e-10 {
        return 1.0 / liquid_density;
    }
    if x > 1.0 - 1e-10 {
        return 1.0 / vapor_density;
    }
    let alpha = homogeneous_void_fraction(x, liquid_density, vapor_density);
    let liquid_term = (1.0 - x).powi(2) / ((1.0 - alpha) * liquid_density);
    let vapor_term = x * x / (alpha * vapor_density);
    liquid_term + vapor_term
}

/// The acceleration pressure drop over a quality change.
#[must_use]
pub fn acceleration_pressure_drop(
    mass_flux: f64,
    quality_in: f64,
    quality_out: f64,
    liquid_density: f64,
    vapor_density: f64,
) -> f64 {
    if mass_flux <= 0.0 || liquid_density <= 0.0 || vapor_density <= 0.0 {
        return 0.0;
    }
    let flux_out = momentum_flux(quality_out, liquid_density, vapor_density);
    let flux_in = momentum_flux(quality_in, liquid_density, vapor_density);
    mass_flux * mass_flux * (flux_out - flux_in)
}

/// The composed two-phase pressure-drop calculation.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a non-integral interval count, and
/// [`AzothError::OutOfRange`] for the spec's own bounds.
pub fn rate(request: &TwoPhasePressureDropRequest) -> Result<TwoPhasePressureDropResult> {
    let spec = &spec_gen::TWO_PHASE_PRESSURE_DROP_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "mass_flux" => Some(request.mass_flux),
            "tube_inner_diameter" => Some(request.tube_inner_diameter),
            "tube_length" => Some(request.tube_length),
            "liquid_density" => Some(request.liquid_density),
            "vapor_density" => Some(request.vapor_density),
            "liquid_viscosity" => Some(request.liquid_viscosity),
            "vapor_viscosity" => Some(request.vapor_viscosity),
            "surface_tension" => Some(request.surface_tension),
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
    Ok(TwoPhasePressureDropResult {
        friedel_gradient: pascals(friedel_gradient(
            request.mass_flux,
            request.vapor_quality,
            request.tube_inner_diameter,
            request.liquid_density,
            request.vapor_density,
            request.liquid_viscosity,
            request.vapor_viscosity,
            request.surface_tension,
        )),
        friedel_pressure_drop: pascals(friedel_pressure_drop(
            request.mass_flux,
            request.vapor_quality,
            request.tube_inner_diameter,
            request.tube_length,
            request.liquid_density,
            request.vapor_density,
            request.liquid_viscosity,
            request.vapor_viscosity,
            request.surface_tension,
        )),
        friedel_average_pressure_drop: pascals(friedel_average_pressure_drop(
            request.mass_flux,
            request.quality_in,
            request.quality_out,
            request.tube_inner_diameter,
            request.tube_length,
            request.liquid_density,
            request.vapor_density,
            request.liquid_viscosity,
            request.vapor_viscosity,
            request.surface_tension,
            intervals,
        )),
        muller_steinhagen_heck_gradient: pascals(muller_steinhagen_heck_gradient(
            request.mass_flux,
            request.vapor_quality,
            request.tube_inner_diameter,
            request.liquid_density,
            request.vapor_density,
            request.liquid_viscosity,
            request.vapor_viscosity,
        )),
        gravitational_gradient: pascals(gravitational_gradient(
            request.vapor_quality,
            request.liquid_density,
            request.vapor_density,
        )),
        acceleration_pressure_drop: pascals(acceleration_pressure_drop(
            request.mass_flux,
            request.quality_in,
            request.quality_out,
            request.liquid_density,
            request.vapor_density,
        )),
        warnings,
    })
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn two_phase_pressure_drop(
    mass_flux: MassFlux,
    vapor_quality: f64,
    tube_inner_diameter: Length,
    tube_length: Length,
    liquid_density: MassDensity,
    vapor_density: MassDensity,
    liquid_viscosity: DynamicViscosity,
    vapor_viscosity: DynamicViscosity,
    surface_tension: SurfaceTension,
    quality_in: f64,
    quality_out: f64,
    intervals: Option<f64>,
) -> Result<TwoPhasePressureDropResult> {
    rate(&TwoPhasePressureDropRequest {
        mass_flux: mass_flux.value,
        vapor_quality,
        tube_inner_diameter: tube_inner_diameter.value,
        tube_length: tube_length.value,
        liquid_density: liquid_density.value,
        vapor_density: vapor_density.value,
        liquid_viscosity: liquid_viscosity.value,
        vapor_viscosity: vapor_viscosity.value,
        surface_tension: surface_tension.value,
        quality_in,
        quality_out,
        intervals,
    })
}
