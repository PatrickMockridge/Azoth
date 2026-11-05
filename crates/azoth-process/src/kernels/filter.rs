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
//!
//! **`run` reads a carried loading state, and the loading is not a dynamic state the record
//! lacks.** A finite capacity turns the solids loading into a fraction, the increase at
//! capacity adds to the clean drop, breakthrough retains the efficiency and a bypass above its
//! cracking drop lowers the applied drop. The class clamps the drop to the inlet ceiling
//! **after** the bypass choice, and this reproduces that order. The accretion over `dt` is not
//! here; it is the separate step the model id raises by `time_step`.

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

/// `Filter`'s carried loading state, as the model id's six optional inputs state it.
///
/// **The class keeps these as fields on the equipment and `run` reads them every time**, so the
/// request is resolved here rather than at the boundary: absent is the class's initial value,
/// which is what makes an existing case's answer unchanged.
#[derive(Debug, Clone, Copy, Default)]
pub struct LoadingRequest {
    /// The accumulated loading, kg; absent is 0.
    pub solids_loading_kg: Option<f64>,
    /// The loading at which breakthrough reaches one, kg; absent is `+Infinity`.
    pub loading_capacity_kg: Option<f64>,
    /// The extra drop at one capacity equivalent, Pa; absent is 0.
    pub pressure_drop_increase_at_capacity_pa: Option<f64>,
    /// The loading fraction where breakthrough starts; absent is 0.8.
    pub breakthrough_start_fraction: Option<f64>,
    /// The bypass's cracking differential pressure, Pa; absent leaves it closed.
    pub bypass_cracking_pressure_pa: Option<f64>,
    /// The differential pressure calling for replacement, Pa; absent is 5 bar.
    pub terminal_pressure_drop_pa: Option<f64>,
}

/// The class's `DEFAULT_TERMINAL_DELTA_P` of 5 bar, in Pa.
pub const DEFAULT_TERMINAL_DROP_PA: f64 = 500_000.0;
/// `setBreakthroughStartFraction`'s ceiling: the divisor `1 - start` stays positive.
pub const MAX_BREAKTHROUGH_START: f64 = 0.999_999;

/// `updateHydraulicState`'s arithmetic, and the capture's current efficiency.
#[derive(Debug, Clone, Copy)]
pub struct LoadingState {
    /// The carried loading, kg, clamped at zero.
    pub solids_loading_kg: f64,
    /// `solids_loading / capacity`, or zero without a finite positive capacity.
    pub loading_fraction: f64,
    /// The fraction past the breakthrough start, capped at one.
    pub breakthrough_fraction: f64,
    /// The nominal efficiency retained by breakthrough and bypass.
    pub efficiency: f64,
    /// The clean drop plus the loading increase, unclamped and pre-bypass, Pa.
    pub unrestricted_drop: f64,
    /// The drop after the bypass choice and before the inlet-pressure clamp, Pa.
    pub restricted_drop: f64,
    /// The bypassed flow fraction.
    pub bypass_fraction: f64,
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
    /// The drop applied, Pa - the restricted one unless it was clamped.
    pub applied_drop: Pressure,
    /// The drop before the inlet-pressure clamp, Pa.
    pub restricted_drop: Pressure,
    /// The clean drop plus the loading increase, before bypass, Pa.
    pub unrestricted_drop: Pressure,
    /// The nominal removal efficiency, a fraction.
    pub nominal_removal_efficiency: f64,
    /// The nominal efficiency retained by breakthrough and bypass.
    pub current_removal_efficiency: f64,
    /// The outlet particle concentration as a mass fraction, zero when the model is off.
    pub outlet_particle_concentration: f64,
    /// The captured mass rate, kg/s.
    pub captured_rate: MassRate,
    /// The carried loading, kg.
    pub solids_loading_kg: f64,
    /// `solids_loading / capacity`, or zero without a finite capacity.
    pub loading_fraction: f64,
    /// The fraction past the breakthrough start, capped at one.
    pub breakthrough_fraction: f64,
    /// The bypassed flow fraction.
    pub bypass_fraction: f64,
    /// `unrestricted / terminal`.
    pub differential_pressure_utilization: f64,
    /// Either the terminal drop or the loading capacity has been reached.
    pub replacement_required: bool,
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
    filter_with(
        feed,
        pressure_drop,
        &settings,
        None,
        &LoadingRequest::default(),
    )
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
    loading: &LoadingRequest,
) -> Result<FilterOutcome> {
    let actual_flow_m3h = system_volumetric_flow(feed)? * 3600.0;
    let clean_pa = clean_drop(feed, pressure_drop.value, actual_flow_m3h, settings)?;
    let nominal = nominal_removal_efficiency(
        settings.beta_curve,
        settings.particle_size_m,
        settings.nominal_efficiency,
    );
    let state = loading_state(clean_pa, nominal, loading);

    // `Math.max(0.0, inletPressure - 1.0e-6)` in the class's bara: a millionth of a bar. The
    // clamp is the class's **last** step: bypass is chosen from the unrestricted drop above.
    let ceiling = (feed.p.value - 0.1).max(0.0);
    let applied = state.restricted_drop.clamp(0.0, ceiling);
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

    let (outlet_concentration, captured) = match inlet_concentration {
        Some(inlet) => {
            let mass_flow = feed.n * feed.molar_mass()?.value;
            (
                inlet * (1.0 - state.efficiency),
                kilograms_per_hour(mass_flow * inlet * state.efficiency * 3600.0),
            )
        }
        None => (0.0, kilograms_per_hour(0.0)),
    };

    let terminal = loading
        .terminal_pressure_drop_pa
        .unwrap_or(DEFAULT_TERMINAL_DROP_PA);
    let (utilization, replacement_required) =
        utilization(state.unrestricted_drop, terminal, state.loading_fraction);

    Ok(FilterOutcome {
        outlet,
        clean_drop: pascals(clean_pa),
        applied_drop: pascals(applied),
        restricted_drop: pascals(state.restricted_drop),
        unrestricted_drop: pascals(state.unrestricted_drop),
        nominal_removal_efficiency: nominal,
        current_removal_efficiency: state.efficiency,
        outlet_particle_concentration: outlet_concentration,
        captured_rate: captured,
        solids_loading_kg: state.solids_loading_kg,
        loading_fraction: state.loading_fraction,
        breakthrough_fraction: state.breakthrough_fraction,
        bypass_fraction: state.bypass_fraction,
        differential_pressure_utilization: utilization,
        replacement_required,
    })
}

/// `getDifferentialPressureUtilization` and `isReplacementRequired`, on the unrestricted drop.
fn utilization(unrestricted_drop: f64, terminal_drop: f64, loading_fraction: f64) -> (f64, bool) {
    let utilization = unrestricted_drop / terminal_drop;
    (utilization, utilization >= 1.0 || loading_fraction >= 1.0)
}

/// `updateHydraulicState` and the capture's current efficiency, from the carried request.
///
/// **The clamp is not here**: this is the class's state before `run` applies the inlet ceiling,
/// and `restricted_drop` is the cracking pressure or the unrestricted drop. The bypass decision
/// reads the unrestricted drop, and the comparison is the class's strict `>`.
fn loading_state(clean_pa: f64, nominal: f64, request: &LoadingRequest) -> LoadingState {
    let solids_loading_kg = request.solids_loading_kg.unwrap_or(0.0).max(0.0);
    let capacity = request.loading_capacity_kg.unwrap_or(f64::INFINITY);
    let loading_fraction = if capacity.is_finite() && capacity > 0.0 {
        solids_loading_kg / capacity
    } else {
        0.0
    };
    // numerics-ok: `setBreakthroughStartFraction` clamps the fraction to [0, 0.999999].
    let start = request
        .breakthrough_start_fraction
        .unwrap_or(0.8)
        .clamp(0.0, MAX_BREAKTHROUGH_START);
    let breakthrough_fraction = if loading_fraction <= start {
        0.0
    } else {
        ((loading_fraction - start) / (1.0 - start)).min(1.0)
    };
    let increase = request
        .pressure_drop_increase_at_capacity_pa
        .unwrap_or(0.0)
        .max(0.0);
    let unrestricted_drop = clean_pa + increase * loading_fraction.max(0.0);
    let (bypass_fraction, restricted_drop) = match request.bypass_cracking_pressure_pa {
        Some(cracking) if unrestricted_drop > cracking => (
            // numerics-ok: `updateHydraulicState` clamps the bypass with Math.max/Math.min.
            (1.0 - (cracking / unrestricted_drop).sqrt()).clamp(0.0, 1.0),
            cracking,
        ),
        _ => (0.0, unrestricted_drop),
    };
    // numerics-ok: `getCurrentRemovalEfficiency` clamps with Math.max(0, Math.min(1, ...)).
    let efficiency =
        (nominal * (1.0 - breakthrough_fraction) * (1.0 - bypass_fraction)).clamp(0.0, 1.0);
    LoadingState {
        solids_loading_kg,
        loading_fraction,
        breakthrough_fraction,
        efficiency,
        unrestricted_drop,
        restricted_drop,
        bypass_fraction,
    }
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

    fn loaded(loading: f64, capacity: f64, increase: f64, start: f64) -> LoadingRequest {
        LoadingRequest {
            solids_loading_kg: Some(loading),
            loading_capacity_kg: Some(capacity),
            pressure_drop_increase_at_capacity_pa: Some(increase),
            breakthrough_start_fraction: Some(start),
            ..LoadingRequest::default()
        }
    }

    #[test]
    fn the_loaded_state_matches_the_capture() {
        let state = loading_state(20_000.0, 0.0, &loaded(9.0, 10.0, 100_000.0, 0.5));
        assert!((state.loading_fraction - 0.9).abs() < 1.0e-12);
        assert!((state.breakthrough_fraction - 0.8).abs() < 1.0e-12);
        assert!((state.unrestricted_drop - 110_000.0).abs() < 1.0e-9);
        assert!((state.restricted_drop - 110_000.0).abs() < 1.0e-9);
        assert_eq!(state.bypass_fraction, 0.0);
    }

    #[test]
    fn the_bypass_state_matches_the_capture() {
        let mut request = loaded(5.0, 10.0, 360_000.0, 0.5);
        request.bypass_cracking_pressure_pa = Some(50_000.0);
        let state = loading_state(20_000.0, 0.99, &request);
        assert!((state.unrestricted_drop - 200_000.0).abs() < 1.0e-9);
        assert!((state.restricted_drop - 50_000.0).abs() < 1.0e-9);
        assert!((state.bypass_fraction - 0.5).abs() < 1.0e-12);
        assert!((state.efficiency - 0.495).abs() < 1.0e-12);
    }

    #[test]
    fn breakthrough_starts_at_the_boundary_and_is_capped() {
        let at = loading_state(20_000.0, 0.0, &loaded(5.0, 10.0, 100_000.0, 0.5));
        assert_eq!(at.breakthrough_fraction, 0.0);
        let above = loading_state(20_000.0, 0.0, &loaded(7.5, 10.0, 100_000.0, 0.5));
        assert!((above.breakthrough_fraction - 0.5).abs() < 1.0e-12);
        let past = loading_state(20_000.0, 0.0, &loaded(20.0, 10.0, 100_000.0, 0.5));
        assert_eq!(past.breakthrough_fraction, 1.0);
        let clamped = loading_state(20_000.0, 0.0, &loaded(20.0, 10.0, 100_000.0, 5.0));
        assert_eq!(clamped.breakthrough_fraction, 1.0);
    }

    #[test]
    fn the_bypass_comparison_is_strict() {
        let mut request = loaded(5.0, 10.0, 360_000.0, 0.5);
        request.bypass_cracking_pressure_pa = Some(200_000.0);
        let state = loading_state(20_000.0, 0.99, &request);
        assert_eq!(state.bypass_fraction, 0.0);
        assert!((state.restricted_drop - 200_000.0).abs() < 1.0e-9);
    }

    #[test]
    fn a_zero_capacity_is_the_sentinel() {
        let state = loading_state(20_000.0, 0.0, &loaded(9.0, 0.0, 100_000.0, 0.5));
        assert_eq!(state.loading_fraction, 0.0);
        assert_eq!(state.breakthrough_fraction, 0.0);
        assert!((state.unrestricted_drop - 20_000.0).abs() < 1.0e-9);
    }

    #[test]
    fn negatives_are_clamped_and_the_bypass_reads_the_unclamped_drop() {
        let state = loading_state(20_000.0, 0.0, &loaded(-3.0, 10.0, -1.0, 0.5));
        assert_eq!(state.solids_loading_kg, 0.0);
        assert_eq!(state.unrestricted_drop, 20_000.0);

        let mut clamp = loaded(0.0, 10.0, 0.0, 0.5);
        clamp.bypass_cracking_pressure_pa = Some(2_500_000.0);
        let state = loading_state(2_600_000.0, 0.0, &clamp);
        assert!((state.unrestricted_drop - 2_600_000.0).abs() < 1.0e-6);
        assert!((state.restricted_drop - 2_500_000.0).abs() < 1.0e-6);
        assert!((state.bypass_fraction - 0.019_419_324_309_079_777).abs() < 1.0e-12);
    }

    #[test]
    fn the_replacement_check_reads_either_threshold() {
        let (utilization_at_dp, required_at_dp) = utilization(200_000.0, 500_000.0, 0.5);
        assert!((utilization_at_dp - 0.4).abs() < 1.0e-12);
        assert!(!required_at_dp);
        let (_, required_at_loading) = utilization(200_000.0, 500_000.0, 1.1);
        assert!(required_at_loading);
        let (utilization_at_terminal, required_at_terminal) =
            utilization(500_000.0, 500_000.0, 0.5);
        assert_eq!(utilization_at_terminal, 1.0);
        assert!(required_at_terminal);
    }
}
