//! `unit_ops.filter` - a pressure drop, flashed at the feed's temperature.
//!
//! **It holds the temperature, which is what separates this from `throttling_valve`.**
//! `Filter.run` sets the reduced pressure and runs a `TPflash`, so the outlet is at the
//! *feed's* temperature and its enthalpy moves with the pressure; a valve is a `PHflash` and
//! holds the enthalpy instead. On the captured fluid the two differ by about `25` J/mol over
//! one bar, so a port that took the ordinary isenthalpic reading of "a pressure drop across a
//! restriction" would be wrong by a visible amount.
//!
//! **A drop larger than the inlet pressure is clamped, not refused.** `run` applies
//! `min(max(0, dP), max(0, P_in - 1e-6 bar))` and logs a warning, so the outlet lands at
//! `1e-6` bar above vacuum rather than at a negative pressure. That is the class's own
//! behaviour and this reproduces it, returning the applied drop so the caller can say so.
//!
//! **The clean drop is one of four models.** `FilterPressureDropModel`'s `FIXED` returns the
//! configured drop; `FLOW_SCALED` multiplies it by `(Q/Q_ref)^n`; `TABULATED` interpolates a
//! measured curve; `ERGUN` evaluates the packed-bed correlation. The particle capture is the
//! class's `updateParticleCapturePerformance`: a beta-ratio curve or a scalar efficiency, and
//! the inlet concentration the caller states.

use azoth_core::Result;
use azoth_core::units::{MassRate, Pressure, kilograms_per_hour, pascals};

use crate::segment::phase::{system_mass_density, system_volumetric_flow};
use crate::stream::Stream;

/// `FilterPressureDropModel`: the four clean-drop branches.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PressureDropModel {
    /// `FIXED`: the configured drop, independent of flow.
    Fixed,
    /// `FLOW_SCALED`: the configured drop scaled by `(Q/Q_ref)^n`.
    FlowScaled,
    /// `TABULATED`: a measured drop-versus-flow curve.
    Tabulated,
    /// `ERGUN`: the packed-bed correlation.
    Ergun,
}

/// One point of a `TABULATED` clean-drop curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DropPoint {
    /// The actual volumetric flow, m3/hr.
    pub flow_m3h: f64,
    /// The clean drop there, Pa.
    pub drop_pa: f64,
}

/// One point of a beta-ratio performance curve.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct BetaPoint {
    /// The particle size, m.
    pub particle_size_m: f64,
    /// The beta ratio at that size, at least one.
    pub beta_ratio: f64,
}

/// `setMediaGeometry`'s four numbers, with the class's own defaults.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Media {
    /// The bed cross-sectional area, m2.
    pub area_m2: f64,
    /// The bed depth, m.
    pub bed_depth_m: f64,
    /// The representative particle diameter, m.
    pub particle_diameter_m: f64,
    /// The bed void fraction in `(0, 1)`.
    pub void_fraction: f64,
}

impl Default for Media {
    fn default() -> Self {
        Self {
            area_m2: 1.0,
            bed_depth_m: 1.0,
            particle_diameter_m: 0.001,
            void_fraction: 0.40,
        }
    }
}

/// The settings the clean-drop model and the capture read.
#[derive(Debug, Clone)]
pub struct FilterSettings<'a> {
    /// The selected clean-drop model.
    pub model: PressureDropModel,
    /// The `FLOW_SCALED` reference flow, m3/hr; absent takes the actual flow.
    pub reference_flow_m3h: Option<f64>,
    /// The `FLOW_SCALED` exponent, positive.
    pub flow_exponent: f64,
    /// The `TABULATED` curve.
    pub drop_curve: &'a [DropPoint],
    /// The `ERGUN` media geometry.
    pub media: Media,
    /// The `ERGUN` viscosity, Pa s; absent falls back to the configured drop.
    pub viscosity_pa_s: Option<f64>,
    /// The beta-ratio curve; empty uses `nominal_efficiency`.
    pub beta_curve: &'a [BetaPoint],
    /// The particle size the beta curve is evaluated at, m.
    pub particle_size_m: f64,
    /// The scalar removal efficiency, a fraction.
    pub nominal_efficiency: f64,
}

/// What a filter did: the outlet record, the drop it applied, and the capture's answers.
#[derive(Debug, Clone)]
pub struct FilterOutcome {
    /// The outlet record.
    pub outlet: Stream,
    /// The clean drop the selected model computed, Pa.
    pub clean_drop: Pressure,
    /// The drop applied, Pa - the clean one unless it was clamped.
    pub applied_drop: Pressure,
    /// The nominal removal efficiency, a fraction.
    pub nominal_removal_efficiency: f64,
    /// The outlet particle concentration as a mass fraction, zero when the model is off.
    pub outlet_particle_concentration: f64,
    /// The captured mass rate, kg/s.
    pub captured_rate: MassRate,
}

/// The class's `FIXED` model, which is what the palette's filter runs.
///
/// # Errors
/// Whatever the flash refuses.
pub fn filter(feed: &Stream, pressure_drop: Pressure) -> Result<FilterOutcome> {
    let settings = FilterSettings {
        model: PressureDropModel::Fixed,
        reference_flow_m3h: None,
        flow_exponent: 1.0,
        drop_curve: &[],
        media: Media::default(),
        viscosity_pa_s: None,
        beta_curve: &[],
        particle_size_m: 1.0e-5,
        nominal_efficiency: 0.0,
    };
    filter_with(feed, pressure_drop, &settings, None)
}

/// Drop a stream's pressure by the selected clean-drop model, at a constant temperature.
///
/// # Errors
/// [`azoth_core::AzothError::OutOfRange`] or the flash's own error for a state that does not
/// exist; an `ERGUN` geometry that is not a bed is refused by the model before this call.
pub fn filter_with(
    feed: &Stream,
    pressure_drop: Pressure,
    settings: &FilterSettings<'_>,
    inlet_concentration: Option<f64>,
) -> Result<FilterOutcome> {
    let actual_flow_m3h = system_volumetric_flow(feed)? * 3600.0;
    let clean_pa = clean_drop(feed, pressure_drop.value, actual_flow_m3h, settings)?;
    let clean = pascals(clean_pa);

    // `Math.max(0.0, inletPressure - 1.0e-6)` in the class's bara: a millionth of a bar.
    let ceiling = (feed.p.value - 0.1).max(0.0);
    let applied = clean.value.clamp(0.0, ceiling);
    let p_out = pascals(feed.p.value - applied);

    // The class flashes only when the drop is above its own `1e-10` bar threshold; below it
    // the outlet is the cloned inlet untouched. A TP flash at the same state is that state,
    // so the two agree and this needs no second branch.
    let outlet = Stream::from_pt(
        feed.components.clone(),
        feed.z.clone(),
        feed.n,
        p_out,
        feed.t,
    )?;

    let efficiency = nominal_removal_efficiency(
        settings.beta_curve,
        settings.particle_size_m,
        settings.nominal_efficiency,
    );
    let (outlet_concentration, captured) = match inlet_concentration {
        Some(inlet) => {
            let mass_flow = feed.n * feed.molar_mass()?.value;
            (
                inlet * (1.0 - efficiency),
                kilograms_per_hour(mass_flow * inlet * efficiency * 3600.0),
            )
        }
        None => (0.0, kilograms_per_hour(0.0)),
    };

    Ok(FilterOutcome {
        outlet,
        clean_drop: clean,
        applied_drop: pascals(applied),
        nominal_removal_efficiency: efficiency,
        outlet_particle_concentration: outlet_concentration,
        captured_rate: captured,
    })
}

/// `calculateCleanPressureDrop`, in Pa.
fn clean_drop(
    feed: &Stream,
    configured: f64,
    actual_flow_m3h: f64,
    settings: &FilterSettings<'_>,
) -> Result<f64> {
    match settings.model {
        PressureDropModel::Fixed => Ok(configured),
        PressureDropModel::FlowScaled => {
            let reference = settings.reference_flow_m3h.unwrap_or(actual_flow_m3h);
            if reference.is_nan() || reference <= 0.0 {
                return Ok(0.0);
            }
            Ok(configured * (actual_flow_m3h / reference).powf(settings.flow_exponent))
        }
        PressureDropModel::Tabulated => {
            if settings.drop_curve.is_empty() {
                return Ok(configured);
            }
            Ok(curve_drop(settings.drop_curve, actual_flow_m3h))
        }
        PressureDropModel::Ergun => {
            if actual_flow_m3h <= 0.0 {
                return Ok(0.0);
            }
            let density = system_mass_density(feed)?;
            let Some(viscosity) = settings.viscosity_pa_s else {
                return Ok(configured);
            };
            if !density.is_finite() || density <= 0.0 || !viscosity.is_finite() || viscosity <= 0.0
            {
                return Ok(configured);
            }
            let media = settings.media;
            let u = actual_flow_m3h / 3600.0 / media.area_m2;
            let solid = 1.0 - media.void_fraction;
            let void_cube = media.void_fraction.powi(3);
            let viscous = 150.0 * viscosity * solid * solid
                / (void_cube * media.particle_diameter_m * media.particle_diameter_m)
                * u;
            let inertial = 1.75 * density * solid / (void_cube * media.particle_diameter_m) * u * u;
            Ok((media.bed_depth_m * (viscous + inertial)).max(0.0))
        }
    }
}

/// `FilterPressureDropCurve.getPressureDrop`, Pa.
fn curve_drop(curve: &[DropPoint], flow_m3h: f64) -> f64 {
    if flow_m3h <= 0.0 {
        return 0.0;
    }
    let first = curve[0];
    if flow_m3h <= first.flow_m3h {
        return (first.drop_pa * flow_m3h / first.flow_m3h).max(0.0);
    }
    for window in curve.windows(2) {
        if flow_m3h <= window[1].flow_m3h {
            return interpolate(window[0], window[1], flow_m3h);
        }
    }
    let last = curve[curve.len() - 1];
    if curve.len() == 1 {
        return (last.drop_pa * flow_m3h / last.flow_m3h).max(0.0);
    }
    interpolate(curve[curve.len() - 2], last, flow_m3h)
}

fn interpolate(lower: DropPoint, upper: DropPoint, flow_m3h: f64) -> f64 {
    let fraction = (flow_m3h - lower.flow_m3h) / (upper.flow_m3h - lower.flow_m3h);
    (lower.drop_pa + fraction * (upper.drop_pa - lower.drop_pa)).max(0.0)
}

/// `Filter.getNominalRemovalEfficiency`: the beta curve's value, else the scalar.
fn nominal_removal_efficiency(curve: &[BetaPoint], particle_size_m: f64, scalar: f64) -> f64 {
    if curve.is_empty() {
        return scalar.clamp(0.0, 1.0);
    }
    let beta = beta_ratio(curve, particle_size_m);
    (1.0 - 1.0 / beta).clamp(0.0, 1.0)
}

/// `FilterPerformanceCurve.getBetaRatio`: linear in size, logarithmic in the ratio.
fn beta_ratio(curve: &[BetaPoint], particle_size_m: f64) -> f64 {
    let first = curve[0];
    if particle_size_m <= first.particle_size_m {
        return first.beta_ratio;
    }
    let last = curve[curve.len() - 1];
    if particle_size_m >= last.particle_size_m {
        return last.beta_ratio;
    }
    for window in curve.windows(2) {
        if particle_size_m <= window[1].particle_size_m {
            let fraction = (particle_size_m - window[0].particle_size_m)
                / (window[1].particle_size_m - window[0].particle_size_m);
            let lower = window[0].beta_ratio.ln();
            let upper = window[1].beta_ratio.ln();
            return (lower + fraction * (upper - lower)).exp();
        }
    }
    last.beta_ratio
}

#[cfg(test)]
mod tests {
    use super::*;

    fn curve() -> [DropPoint; 2] {
        [
            DropPoint {
                flow_m3h: 0.5,
                drop_pa: 10_000.0,
            },
            DropPoint {
                flow_m3h: 1.5,
                drop_pa: 50_000.0,
            },
        ]
    }

    #[test]
    fn the_tabulated_curve_interpolates_and_extends() {
        assert_eq!(curve_drop(&curve(), 1.0), 30_000.0);
        assert_eq!(curve_drop(&curve(), 0.25), 5_000.0);
        assert_eq!(curve_drop(&curve(), 2.0), 70_000.0);
        assert_eq!(curve_drop(&curve(), 0.0), 0.0);
    }

    #[test]
    fn the_beta_curve_interpolates_logarithmically() {
        let points = [
            BetaPoint {
                particle_size_m: 5.0e-6,
                beta_ratio: 2.0,
            },
            BetaPoint {
                particle_size_m: 2.0e-5,
                beta_ratio: 1000.0,
            },
        ];
        assert!((beta_ratio(&points, 5.0e-6) - 2.0).abs() < 1.0e-12);
        assert!((beta_ratio(&points, 2.0e-5) - 1000.0).abs() < 1.0e-12);
        // The geometric mean at the midpoint, which is what log-linear interpolation is.
        assert!((beta_ratio(&points, 1.25e-5) - (2.0_f64 * 1000.0).sqrt()).abs() < 1.0e-9);
        assert!((beta_ratio(&points, 1.0e-6) - 2.0).abs() < 1.0e-12);
    }

    #[test]
    fn the_scalar_efficiency_is_the_curve_s_absence() {
        assert_eq!(nominal_removal_efficiency(&[], 1.0e-5, 0.99), 0.99);
        assert_eq!(nominal_removal_efficiency(&[], 1.0e-5, 1.5), 1.0);
    }
}
