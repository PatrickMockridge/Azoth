//! `thermal.incremental_zone_analysis` - the HTRI-style incremental rating.
//!
//! Transcribed from `IncrementalZoneAnalysis.java` at NeqSim `118d9aff`. The class is a
//! single pass over the caller's zones in insertion order: no subdivision, no iteration and no
//! convergence test. The tube side dispatches on the regime - Shah for `CONDENSING`,
//! Gungor-Winterton for `EVAPORATING`, the single-phase path for everything else, `TWO_PHASE`
//! included - and the shell side is always Kern and ignores its regime.
//!
//! **A zero or non-positive flow is refused.** The class answers one - a zero tube flow gives
//! the laminar coefficient - but a coefficient computed from no flow is not an answer this
//! port publishes.
//!
//! Spec: `specs/calcs/thermal/incremental_zone_analysis.toml`.

use azoth_core::units::{
    Length, ThermalConductivity, kelvin_intervals, pascals, square_meters, watts,
    watts_per_square_meter_kelvin,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::bell_delaware;
use crate::boiling_heat_transfer;
use crate::results::{IncrementalZoneAnalysisResult, IncrementalZoneInput, IncrementalZoneRow};
use crate::shah_condensation;
use crate::shell_and_tube_rating::darcy_friction;
use crate::spec_gen;
use crate::two_phase_pressure_drop;

/// The class's `IncrementalZone` defaults for the optional two-phase fields.
const LIQUID_DENSITY_DEFAULT: f64 = 800.0;
const VAPOR_DENSITY_DEFAULT: f64 = 5.0;
const LIQUID_VISCOSITY_DEFAULT: f64 = 5.0e-4;
const VAPOR_VISCOSITY_DEFAULT: f64 = 1.0e-5;
const LIQUID_CP_DEFAULT: f64 = 2500.0;
const LIQUID_CONDUCTIVITY_DEFAULT: f64 = 0.13;
const QUALITY_IN_DEFAULT: f64 = 0.9;
const QUALITY_OUT_DEFAULT: f64 = 0.7;
const REDUCED_PRESSURE_DEFAULT: f64 = 0.1;
const SURFACE_TENSION_DEFAULT: f64 = 0.01;
const HEAT_OF_VAPORIZATION_DEFAULT: f64 = 200_000.0;

/// The zone's phase regime, as the spec spells it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseRegime {
    /// Single-phase vapour.
    Vapor,
    /// Single-phase liquid.
    Liquid,
    /// Two-phase condensing.
    Condensing,
    /// Two-phase evaporating.
    Evaporating,
    /// Two-phase mixed; the class falls through to the single-phase path.
    TwoPhase,
}

impl std::str::FromStr for PhaseRegime {
    type Err = AzothError;

    fn from_str(value: &str) -> Result<Self> {
        match value {
            "vapor" => Ok(Self::Vapor),
            "liquid" => Ok(Self::Liquid),
            "condensing" => Ok(Self::Condensing),
            "evaporating" => Ok(Self::Evaporating),
            "two_phase" => Ok(Self::TwoPhase),
            other => Err(AzothError::invalid_input(
                "tube_phase_regime",
                format!(
                    "{other} is not one of vapor, liquid, condensing, evaporating or two_phase"
                ),
            )),
        }
    }
}

impl PhaseRegime {
    /// The class's own enum spelling, which the result rows echo.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Vapor => "VAPOR",
            Self::Liquid => "LIQUID",
            Self::Condensing => "CONDENSING",
            Self::Evaporating => "EVAPORATING",
            Self::TwoPhase => "TWO_PHASE",
        }
    }
}

/// The geometry and the zones, in SI magnitudes.
#[derive(Debug, Clone)]
pub struct IncrementalZoneAnalysisRequest {
    pub tube_inside_diameter: f64,
    pub tube_outside_diameter: f64,
    pub tube_length: f64,
    pub tube_count: f64,
    pub tube_passes: f64,
    pub shell_inside_diameter: f64,
    pub baffle_spacing: f64,
    pub tube_pitch: f64,
    pub triangular_pitch: bool,
    pub tube_wall_conductivity: f64,
    pub fouling_tube: f64,
    pub fouling_shell: f64,
    pub zones: Vec<IncrementalZoneInput>,
}

/// Rates every zone in the caller's order.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a geometry the class cannot use, a non-integral count, or
/// a zone whose regime-read properties are not positive finite.
pub fn rate(request: &IncrementalZoneAnalysisRequest) -> Result<IncrementalZoneAnalysisResult> {
    let spec = &spec_gen::INCREMENTAL_ZONE_ANALYSIS_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "tube_inside_diameter" => Some(request.tube_inside_diameter),
            "tube_outside_diameter" => Some(request.tube_outside_diameter),
            "tube_length" => Some(request.tube_length),
            "tube_count" => Some(request.tube_count),
            "tube_passes" => Some(request.tube_passes),
            "shell_inside_diameter" => Some(request.shell_inside_diameter),
            "baffle_spacing" => Some(request.baffle_spacing),
            "tube_pitch" => Some(request.tube_pitch),
            "tube_wall_conductivity" => Some(request.tube_wall_conductivity),
            "fouling_tube" => Some(request.fouling_tube),
            "fouling_shell" => Some(request.fouling_shell),
            _ => None,
        },
        &mut warnings,
    )?;
    validate_geometry(request)?;

    let tube_count = whole(request.tube_count, "tube_count")?;
    let tube_passes = whole(request.tube_passes, "tube_passes")?;
    let n_zones = request.zones.len();

    let mut total_area = 0.0;
    let mut total_duty = 0.0;
    let mut total_tube_dp = 0.0;
    let mut total_shell_dp = 0.0;
    let mut total_u_times_area = 0.0;
    let mut minimum_approach = f64::MAX;
    let mut rows = Vec::with_capacity(n_zones);

    for (index, zone) in request.zones.iter().enumerate() {
        let tube_regime: PhaseRegime = zone.tube_phase_regime.parse().map_err(|_| {
            AzothError::invalid_input(
                "tube_phase_regime",
                format!("{} is not a declared regime", zone.tube_phase_regime),
            )
        })?;
        let shell_regime: PhaseRegime = zone.shell_phase_regime.parse().map_err(|_| {
            AzothError::invalid_input(
                "shell_phase_regime",
                format!("{} is not a declared regime", zone.shell_phase_regime),
            )
        })?;
        let mass_flux = tube_mass_flux(request, tube_count, tube_passes);
        check_zone(zone, tube_regime, mass_flux)?;

        let h_tube = tube_side_htc(request, zone, tube_regime, mass_flux, n_zones);
        let h_shell = shell_side_htc(request, zone);
        let u = overall_u(request, h_tube, h_shell);
        let lmtd = zone_lmtd(zone);
        let area = if u > 0.0 && lmtd > 0.0 {
            zone.duty.value / (u * lmtd)
        } else {
            0.0
        };
        let tube_dp = tube_side_dp(request, zone, tube_regime, mass_flux, n_zones);
        let shell_dp = shell_side_dp(request, zone, n_zones);

        total_area += area;
        total_duty += zone.duty.value;
        total_tube_dp += tube_dp;
        total_shell_dp += shell_dp;
        total_u_times_area += u * area;
        minimum_approach = minimum_approach.min(
            (zone.hot_inlet_temp.value - zone.cold_outlet_temp.value)
                .abs()
                .min((zone.hot_outlet_temp.value - zone.cold_inlet_temp.value).abs()),
        );

        rows.push(IncrementalZoneRow {
            zone_index: (index + 1) as f64,
            zone_name: zone.zone_name.clone(),
            tube_phase_regime: tube_regime.as_str().to_string(),
            shell_phase_regime: shell_regime.as_str().to_string(),
            tube_side_heat_transfer_coefficient: watts_per_square_meter_kelvin(h_tube),
            shell_side_heat_transfer_coefficient: watts_per_square_meter_kelvin(h_shell),
            overall_heat_transfer_coefficient: watts_per_square_meter_kelvin(u),
            lmtd: kelvin_intervals(lmtd),
            duty: zone.duty,
            required_area: square_meters(area),
            tube_side_pressure_drop: pascals(tube_dp),
            shell_side_pressure_drop: pascals(shell_dp),
            hot_inlet_temp: zone.hot_inlet_temp,
            hot_outlet_temp: zone.hot_outlet_temp,
            cold_inlet_temp: zone.cold_inlet_temp,
            cold_outlet_temp: zone.cold_outlet_temp,
        });
    }

    let weighted_u = if total_area > 0.0 {
        total_u_times_area / total_area
    } else {
        0.0
    };

    Ok(IncrementalZoneAnalysisResult {
        number_of_zones: n_zones as f64,
        total_required_area: square_meters(total_area),
        total_duty: watts(total_duty),
        weighted_average_u: watts_per_square_meter_kelvin(weighted_u),
        total_tube_side_pressure_drop: pascals(total_tube_dp),
        total_shell_side_pressure_drop: pascals(total_shell_dp),
        minimum_approach_temperature: kelvin_intervals(minimum_approach),
        zones: rows,
        warnings,
    })
}

fn whole(value: f64, field: &str) -> Result<i64> {
    if !value.is_finite() || value <= 0.0 || (value - value.round()).abs() > 0.0 {
        return Err(AzothError::invalid_input(
            field,
            format!("{value} is not a positive whole count"),
        ));
    }
    Ok(value.round() as i64)
}

fn validate_geometry(request: &IncrementalZoneAnalysisRequest) -> Result<()> {
    let positives = [
        request.tube_inside_diameter,
        request.tube_outside_diameter,
        request.tube_length,
        request.shell_inside_diameter,
        request.baffle_spacing,
        request.tube_pitch,
        request.tube_wall_conductivity,
    ];
    if positives
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
    {
        return Err(AzothError::invalid_input(
            "geometry",
            "the incremental analysis requires positive finite geometry",
        ));
    }
    if !request.fouling_tube.is_finite()
        || request.fouling_tube < 0.0
        || !request.fouling_shell.is_finite()
        || request.fouling_shell < 0.0
    {
        return Err(AzothError::invalid_input(
            "fouling",
            "the incremental analysis requires finite nonnegative fouling resistances",
        ));
    }
    Ok(())
}

fn tube_mass_flux(
    request: &IncrementalZoneAnalysisRequest,
    tube_count: i64,
    tube_passes: i64,
) -> f64 {
    let tubes_per_pass = tube_count / tube_passes;
    let area_per_tube =
        std::f64::consts::PI / 4.0 * request.tube_inside_diameter * request.tube_inside_diameter;
    let total_area = tubes_per_pass as f64 * area_per_tube;
    if total_area > 0.0 {
        // The flux is per zone and the zone's own mass flow divides it.
        1.0 / total_area
    } else {
        0.0
    }
}

fn positive(field: &str, value: f64) -> Result<()> {
    if !value.is_finite() || value <= 0.0 {
        return Err(AzothError::invalid_input(
            field,
            format!("{value} is not positive finite"),
        ));
    }
    Ok(())
}

fn check_zone(
    zone: &IncrementalZoneInput,
    tube_regime: PhaseRegime,
    flux_per_kg: f64,
) -> Result<()> {
    positive("tube_mass_flow", zone.tube_mass_flow.value)?;
    positive("shell_mass_flow", zone.shell_mass_flow.value)?;
    positive("shell_density", zone.shell_density.value)?;
    positive("shell_viscosity", zone.shell_viscosity.value)?;
    positive("shell_cp", zone.shell_cp.value)?;
    positive("shell_conductivity", zone.shell_conductivity.value)?;
    if flux_per_kg <= 0.0 {
        return Err(AzothError::invalid_input(
            "tube_count",
            "the tube flow area is not positive",
        ));
    }
    if !zone.duty.value.is_finite()
        || !zone.hot_inlet_temp.value.is_finite()
        || !zone.hot_outlet_temp.value.is_finite()
        || !zone.cold_inlet_temp.value.is_finite()
        || !zone.cold_outlet_temp.value.is_finite()
    {
        return Err(AzothError::invalid_input(
            "zones",
            format!(
                "the zone {} states a non-finite duty or temperature",
                zone.zone_name
            ),
        ));
    }
    match tube_regime {
        PhaseRegime::Condensing | PhaseRegime::Evaporating => {
            let liquid_density = zone
                .tube_liquid_density
                .map_or(LIQUID_DENSITY_DEFAULT, |v| v.value);
            let liquid_viscosity = zone
                .tube_liquid_viscosity
                .map_or(LIQUID_VISCOSITY_DEFAULT, |v| v.value);
            let liquid_cp = zone.tube_liquid_cp.map_or(LIQUID_CP_DEFAULT, |v| v.value);
            let liquid_k = zone
                .tube_liquid_conductivity
                .map_or(LIQUID_CONDUCTIVITY_DEFAULT, |v| v.value);
            positive("tube_liquid_density", liquid_density)?;
            positive("tube_liquid_viscosity", liquid_viscosity)?;
            positive("tube_liquid_cp", liquid_cp)?;
            positive("tube_liquid_conductivity", liquid_k)?;
            if tube_regime == PhaseRegime::Condensing {
                let reduced = zone
                    .tube_reduced_pressure
                    .unwrap_or(REDUCED_PRESSURE_DEFAULT);
                positive("tube_reduced_pressure", reduced)?;
            } else {
                positive(
                    "tube_vapor_density",
                    zone.tube_vapor_density
                        .map_or(VAPOR_DENSITY_DEFAULT, |v| v.value),
                )?;
                positive(
                    "tube_vapor_viscosity",
                    zone.tube_vapor_viscosity
                        .map_or(VAPOR_VISCOSITY_DEFAULT, |v| v.value),
                )?;
                positive(
                    "tube_heat_of_vaporization",
                    zone.tube_heat_of_vaporization
                        .unwrap_or(HEAT_OF_VAPORIZATION_DEFAULT),
                )?;
            }
        }
        PhaseRegime::Vapor | PhaseRegime::Liquid | PhaseRegime::TwoPhase => {
            positive("tube_density", zone.tube_density.value)?;
            positive("tube_viscosity", zone.tube_viscosity.value)?;
            positive("tube_cp", zone.tube_cp.value)?;
            positive("tube_conductivity", zone.tube_conductivity.value)?;
        }
    }
    Ok(())
}

fn single_phase_htc(mass_flux: f64, viscosity: f64, cp: f64, k: f64, id: f64) -> f64 {
    if viscosity <= 0.0 || k <= 0.0 || id <= 0.0 {
        return 0.0;
    }
    let re = mass_flux * id / viscosity;
    let pr = cp * viscosity / k;
    let nu = if re > 3000.0 {
        let f = darcy_friction(re);
        let value = (f / 8.0) * (re - 1000.0) * pr
            / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
        value.max(0.0)
    } else if re > 2300.0 {
        let f = darcy_friction(3000.0);
        let turbulent =
            (f / 8.0) * 2000.0 * pr / (1.0 + 12.7 * (f / 8.0).sqrt() * (pr.powf(2.0 / 3.0) - 1.0));
        let fraction = (re - 2300.0) / 700.0;
        3.66 + fraction * (turbulent - 3.66)
    } else {
        3.66
    };
    nu * k / id
}

fn tube_side_htc(
    request: &IncrementalZoneAnalysisRequest,
    zone: &IncrementalZoneInput,
    regime: PhaseRegime,
    flux_per_kg: f64,
    n_zones: usize,
) -> f64 {
    let mass_flux = flux_per_kg * zone.tube_mass_flow.value;
    let liquid_density = zone
        .tube_liquid_density
        .map_or(LIQUID_DENSITY_DEFAULT, |v| v.value);
    let liquid_viscosity = zone
        .tube_liquid_viscosity
        .map_or(LIQUID_VISCOSITY_DEFAULT, |v| v.value);
    let liquid_cp = zone.tube_liquid_cp.map_or(LIQUID_CP_DEFAULT, |v| v.value);
    let liquid_k = zone
        .tube_liquid_conductivity
        .map_or(LIQUID_CONDUCTIVITY_DEFAULT, |v| v.value);
    match regime {
        PhaseRegime::Condensing => {
            let h_lo = shah_condensation::liquid_only_htc(
                mass_flux,
                request.tube_inside_diameter,
                liquid_density,
                liquid_viscosity,
                liquid_cp,
                liquid_k,
            );
            shah_condensation::average_htc(
                h_lo,
                zone.tube_reduced_pressure
                    .unwrap_or(REDUCED_PRESSURE_DEFAULT),
                zone.tube_quality_in.unwrap_or(QUALITY_IN_DEFAULT),
                zone.tube_quality_out.unwrap_or(QUALITY_OUT_DEFAULT),
                10,
            )
        }
        PhaseRegime::Evaporating => {
            let heat_flux =
                if request.tube_count > 0.0 && request.tube_inside_diameter > 0.0 && n_zones > 0 {
                    zone.duty.value
                        / (request.tube_count
                            * std::f64::consts::PI
                            * request.tube_inside_diameter
                            * request.tube_length
                            / n_zones as f64)
                } else {
                    0.0
                };
            boiling_heat_transfer::average_htc(
                mass_flux,
                request.tube_inside_diameter,
                liquid_density,
                zone.tube_vapor_density
                    .map_or(VAPOR_DENSITY_DEFAULT, |v| v.value),
                liquid_viscosity,
                liquid_cp,
                liquid_k,
                heat_flux,
                zone.tube_heat_of_vaporization
                    .unwrap_or(HEAT_OF_VAPORIZATION_DEFAULT),
                zone.tube_quality_in.unwrap_or(QUALITY_IN_DEFAULT),
                zone.tube_quality_out.unwrap_or(QUALITY_OUT_DEFAULT),
                10,
            )
        }
        PhaseRegime::Vapor | PhaseRegime::Liquid | PhaseRegime::TwoPhase => single_phase_htc(
            mass_flux,
            zone.tube_viscosity.value,
            zone.tube_cp.value,
            zone.tube_conductivity.value,
            request.tube_inside_diameter,
        ),
    }
}

fn shell_side_htc(request: &IncrementalZoneAnalysisRequest, zone: &IncrementalZoneInput) -> f64 {
    let equivalent = bell_delaware::shell_equivalent_diameter(
        request.tube_outside_diameter,
        request.tube_pitch,
        request.triangular_pitch,
    );
    let area = bell_delaware::crossflow_area(
        request.shell_inside_diameter,
        request.baffle_spacing,
        request.tube_outside_diameter,
        request.tube_pitch,
    );
    if equivalent <= 0.0 || area <= 0.0 {
        return 0.0;
    }
    let mass_flux = zone.shell_mass_flow.value / area;
    let wall_viscosity = zone.shell_viscosity.value;
    bell_delaware::kern_shell_side_htc(
        mass_flux,
        equivalent,
        zone.shell_viscosity.value,
        zone.shell_cp.value,
        zone.shell_conductivity.value,
        wall_viscosity,
    )
}

fn overall_u(request: &IncrementalZoneAnalysisRequest, h_tube: f64, h_shell: f64) -> f64 {
    if h_tube <= 0.0 || h_shell <= 0.0 {
        return 0.0;
    }
    let ratio = request.tube_outside_diameter / request.tube_inside_diameter;
    let total = 1.0 / h_shell
        + request.fouling_shell
        + request.tube_outside_diameter * ratio.ln() / (2.0 * request.tube_wall_conductivity)
        + ratio * request.fouling_tube
        + ratio / h_tube;
    if total > 0.0 { 1.0 / total } else { 0.0 }
}

fn zone_lmtd(zone: &IncrementalZoneInput) -> f64 {
    let dt1 = zone.hot_inlet_temp.value - zone.cold_outlet_temp.value;
    let dt2 = zone.hot_outlet_temp.value - zone.cold_inlet_temp.value;
    if dt1 <= 0.0 || dt2 <= 0.0 {
        return dt1.abs().max(dt2.abs());
    }
    if (dt1 - dt2).abs() < 0.01 {
        return (dt1 + dt2) / 2.0;
    }
    (dt1 - dt2) / (dt1 / dt2).ln()
}

fn tube_side_dp(
    request: &IncrementalZoneAnalysisRequest,
    zone: &IncrementalZoneInput,
    regime: PhaseRegime,
    flux_per_kg: f64,
    n_zones: usize,
) -> f64 {
    let mass_flux = flux_per_kg * zone.tube_mass_flow.value;
    let zone_length = request.tube_length / (n_zones.max(1) as f64);
    if regime == PhaseRegime::Condensing || regime == PhaseRegime::Evaporating {
        let quality = (zone.tube_quality_in.unwrap_or(QUALITY_IN_DEFAULT)
            + zone.tube_quality_out.unwrap_or(QUALITY_OUT_DEFAULT))
            / 2.0;
        return two_phase_pressure_drop::friedel_pressure_drop(
            mass_flux,
            quality,
            request.tube_inside_diameter,
            zone_length,
            zone.tube_liquid_density
                .map_or(LIQUID_DENSITY_DEFAULT, |v| v.value),
            zone.tube_vapor_density
                .map_or(VAPOR_DENSITY_DEFAULT, |v| v.value),
            zone.tube_liquid_viscosity
                .map_or(LIQUID_VISCOSITY_DEFAULT, |v| v.value),
            zone.tube_vapor_viscosity
                .map_or(VAPOR_VISCOSITY_DEFAULT, |v| v.value),
            zone.tube_surface_tension
                .map_or(SURFACE_TENSION_DEFAULT, |v| v.value),
        );
    }
    if zone.tube_viscosity.value <= 0.0 || zone.tube_density.value <= 0.0 || mass_flux <= 0.0 {
        return 0.0;
    }
    let velocity = mass_flux / zone.tube_density.value;
    let re = mass_flux * request.tube_inside_diameter / zone.tube_viscosity.value;
    let friction = darcy_friction(re);
    friction * zone_length / request.tube_inside_diameter
        * zone.tube_density.value
        * velocity
        * velocity
        / 2.0
}

fn shell_side_dp(
    request: &IncrementalZoneAnalysisRequest,
    zone: &IncrementalZoneInput,
    n_zones: usize,
) -> f64 {
    let equivalent = bell_delaware::shell_equivalent_diameter(
        request.tube_outside_diameter,
        request.tube_pitch,
        request.triangular_pitch,
    );
    let area = bell_delaware::crossflow_area(
        request.shell_inside_diameter,
        request.baffle_spacing,
        request.tube_outside_diameter,
        request.tube_pitch,
    );
    if equivalent <= 0.0 || area <= 0.0 {
        return 0.0;
    }
    let mass_flux = zone.shell_mass_flow.value / area;
    let zone_baffles = (request.tube_length / request.baffle_spacing / n_zones.max(1) as f64)
        .round()
        .max(1.0) as usize;
    bell_delaware::kern_shell_side_pressure_drop(
        mass_flux,
        equivalent,
        request.shell_inside_diameter,
        zone_baffles,
        zone.shell_density.value,
        zone.shell_viscosity.value,
        zone.shell_viscosity.value,
    )
}

/// The flat entry point the binding and the spec's case runner call.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn incremental_zone_analysis(
    tube_inside_diameter: Length,
    tube_outside_diameter: Length,
    tube_length: Length,
    tube_count: f64,
    tube_passes: f64,
    shell_inside_diameter: Length,
    baffle_spacing: Length,
    tube_pitch: Length,
    triangular_pitch: bool,
    tube_wall_conductivity: ThermalConductivity,
    fouling_tube: f64,
    fouling_shell: f64,
    zones: &[IncrementalZoneInput],
) -> Result<IncrementalZoneAnalysisResult> {
    rate(&IncrementalZoneAnalysisRequest {
        tube_inside_diameter: tube_inside_diameter.value,
        tube_outside_diameter: tube_outside_diameter.value,
        tube_length: tube_length.value,
        tube_count,
        tube_passes,
        shell_inside_diameter: shell_inside_diameter.value,
        baffle_spacing: baffle_spacing.value,
        tube_pitch: tube_pitch.value,
        triangular_pitch,
        tube_wall_conductivity: tube_wall_conductivity.value,
        fouling_tube,
        fouling_shell,
        zones: zones.to_vec(),
    })
}
