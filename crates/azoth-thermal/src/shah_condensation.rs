//! `thermal.shah_condensation` - Shah's in-tube condensation correlation.
//!
//! Transcribed from `process/mechanicaldesign/heatexchanger/ShahCondensation.java` at
//! NeqSim `118d9aff`, composed as `IncrementalZoneAnalysis` uses it: the liquid-only
//! coefficient feeds the local, average and vertical entries.
//!
//! Spec: `specs/calcs/thermal/shah_condensation.toml`.

use azoth_core::units::watts_per_square_meter_kelvin;
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::ShahCondensationResult;
use crate::spec_gen;

/// The correlation's reduced-pressure floor.
pub const MIN_REDUCED_PRESSURE: f64 = 0.001;
/// The correlation's reduced-pressure ceiling.
pub const MAX_REDUCED_PRESSURE: f64 = 0.44;
/// The class's gravity constant in the vertical `Jg`.
const GRAVITY: f64 = 9.81;

/// The composed condensation request, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct ShahRequest {
    pub mass_flux: f64,
    pub tube_inner_diameter: f64,
    pub liquid_density: f64,
    pub liquid_viscosity: f64,
    pub liquid_cp: f64,
    pub liquid_conductivity: f64,
    pub vapor_density: f64,
    pub reduced_pressure: f64,
    pub vapor_quality: f64,
    pub quality_in: f64,
    pub quality_out: f64,
    pub intervals: Option<f64>,
}

/// The all-liquid Dittus-Boelter coefficient.
#[must_use]
pub fn liquid_only_htc(
    mass_flux: f64,
    tube_id: f64,
    _liquid_density: f64,
    liquid_viscosity: f64,
    liquid_cp: f64,
    liquid_conductivity: f64,
) -> f64 {
    if liquid_viscosity <= 0.0 || liquid_conductivity <= 0.0 || tube_id <= 0.0 || mass_flux <= 0.0 {
        return 0.0;
    }
    let re_lo = mass_flux * tube_id / liquid_viscosity;
    let pr_l = liquid_cp * liquid_viscosity / liquid_conductivity;
    let nu_lo = 0.023 * re_lo.powf(0.8) * pr_l.powf(0.4);
    nu_lo * liquid_conductivity / tube_id
}

/// The local Shah coefficient at a quality and reduced pressure.
#[must_use]
pub fn local_htc(h_lo: f64, vapor_quality: f64, reduced_pressure: f64) -> f64 {
    if h_lo <= 0.0 || reduced_pressure <= 0.0 {
        return 0.0;
    }
    let x = vapor_quality.clamp(0.0, 1.0);
    let pr = reduced_pressure.clamp(MIN_REDUCED_PRESSURE, MAX_REDUCED_PRESSURE);
    let term1 = (1.0 - x).powf(0.8);
    let term2 = 3.8 * x.powf(0.76) * (1.0 - x).powf(0.04) / pr.powf(0.38);
    h_lo * (term1 + term2)
}

/// Simpson's-rule average over a quality range.
#[must_use]
pub fn average_htc(
    h_lo: f64,
    reduced_pressure: f64,
    quality_in: f64,
    quality_out: f64,
    intervals: i32,
) -> f64 {
    if h_lo <= 0.0 || reduced_pressure <= 0.0 {
        return 0.0;
    }
    let x_in = quality_in.clamp(0.0, 1.0);
    let x_out = quality_out.clamp(0.0, 1.0);
    if (x_in - x_out).abs() < 1e-10 {
        return local_htc(h_lo, x_in, reduced_pressure);
    }
    let mut n = intervals.max(4);
    if n % 2 != 0 {
        n += 1;
    }
    let dx = (x_out - x_in) / f64::from(n);
    let mut sum =
        local_htc(h_lo, x_in, reduced_pressure) + local_htc(h_lo, x_out, reduced_pressure);
    for i in 1..n {
        let x = x_in + f64::from(i) * dx;
        let h = local_htc(h_lo, x, reduced_pressure);
        sum += if i % 2 == 0 { 2.0 * h } else { 4.0 * h };
    }
    // The class's own expression is `sum*|dx|/3/|x_out-x_in|*|x_out-x_in|`; the two range
    // factors cancel, and the literal form is carried so the rounding is the same.
    sum * dx.abs() / 3.0 / (x_out - x_in).abs() * (x_out - x_in).abs()
}

/// The vertical-tube coefficient, with the class's gravity branch.
#[must_use]
pub fn vertical_htc(
    h_lo: f64,
    vapor_quality: f64,
    reduced_pressure: f64,
    mass_flux: f64,
    tube_id: f64,
    vapor_density: f64,
    liquid_density: f64,
) -> f64 {
    if h_lo <= 0.0 || reduced_pressure <= 0.0 || tube_id <= 0.0 {
        return 0.0;
    }
    let x = vapor_quality.clamp(0.0, 1.0);
    let argument = GRAVITY * tube_id * vapor_density * (liquid_density - vapor_density);
    // A non-positive argument is the class's `Double.MAX_VALUE` sentinel for `Jg` (its
    // `sqrt` answers NaN and the guard is `> 0`), which then takes the shear-dominated
    // branch; this returns it directly.
    if argument <= 0.0 {
        return local_htc(h_lo, x, reduced_pressure);
    }
    let denominator = argument.sqrt();
    let jg = x * mass_flux / denominator;
    let h_shah = local_htc(h_lo, x, reduced_pressure);
    if jg >= 0.98 {
        return h_shah;
    }
    let pr = reduced_pressure.max(MIN_REDUCED_PRESSURE);
    let h_nu = 1.32 * pr.powf(-1.0 / 3.0) * h_lo;
    let h_grav = h_nu * (1.0 - x).powf(0.8);
    h_shah.max(h_grav)
}

/// The class's validated-range predicate, over the stated values.
#[must_use]
pub fn is_in_valid_range(reduced_pressure: f64, re_liquid_only: f64, mass_flux: f64) -> bool {
    if !(MIN_REDUCED_PRESSURE..=MAX_REDUCED_PRESSURE).contains(&reduced_pressure) {
        return false;
    }
    if re_liquid_only < 350.0 {
        return false;
    }
    (4.0..=820.0).contains(&mass_flux)
}

/// The composed condensation calculation.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a non-integral interval count, and
/// [`AzothError::OutOfRange`] for the spec's own bounds.
pub fn rate(request: &ShahRequest) -> Result<ShahCondensationResult> {
    let spec = &spec_gen::SHAH_CONDENSATION_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "mass_flux" => Some(request.mass_flux),
            "tube_inner_diameter" => Some(request.tube_inner_diameter),
            "liquid_density" => Some(request.liquid_density),
            "liquid_viscosity" => Some(request.liquid_viscosity),
            "liquid_cp" => Some(request.liquid_cp),
            "liquid_conductivity" => Some(request.liquid_conductivity),
            "vapor_density" => Some(request.vapor_density),
            "reduced_pressure" => Some(request.reduced_pressure),
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
    let h_lo = liquid_only_htc(
        request.mass_flux,
        request.tube_inner_diameter,
        request.liquid_density,
        request.liquid_viscosity,
        request.liquid_cp,
        request.liquid_conductivity,
    );
    let re_lo = request.mass_flux * request.tube_inner_diameter / request.liquid_viscosity;
    Ok(ShahCondensationResult {
        liquid_only_htc: watts_per_square_meter_kelvin(h_lo),
        liquid_only_reynolds: re_lo,
        local_htc: watts_per_square_meter_kelvin(local_htc(
            h_lo,
            request.vapor_quality,
            request.reduced_pressure,
        )),
        average_htc: watts_per_square_meter_kelvin(average_htc(
            h_lo,
            request.reduced_pressure,
            request.quality_in,
            request.quality_out,
            intervals,
        )),
        vertical_htc: watts_per_square_meter_kelvin(vertical_htc(
            h_lo,
            request.vapor_quality,
            request.reduced_pressure,
            request.mass_flux,
            request.tube_inner_diameter,
            request.vapor_density,
            request.liquid_density,
        )),
        in_valid_range: is_in_valid_range(request.reduced_pressure, re_lo, request.mass_flux),
        warnings,
    })
}

/// The flat entry point the binding and the case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn shah_condensation(
    mass_flux: f64,
    tube_inner_diameter: f64,
    liquid_density: f64,
    liquid_viscosity: f64,
    liquid_cp: f64,
    liquid_conductivity: f64,
    vapor_density: f64,
    reduced_pressure: f64,
    vapor_quality: f64,
    quality_in: f64,
    quality_out: f64,
    intervals: Option<f64>,
) -> Result<ShahCondensationResult> {
    rate(&ShahRequest {
        mass_flux,
        tube_inner_diameter,
        liquid_density,
        liquid_viscosity,
        liquid_cp,
        liquid_conductivity,
        vapor_density,
        reduced_pressure,
        vapor_quality,
        quality_in,
        quality_out,
        intervals,
    })
}
