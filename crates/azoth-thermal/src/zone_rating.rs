//! `thermal.zone_rating` - the shell-and-tube rating run zone by zone.
//!
//! Transcribed from `ThermalDesignCalculator.calculateZones(ZoneDefinition[])` at NeqSim
//! `f0c7436c69`. Each zone overrides the eight transport properties and reruns the class's
//! permissive `calculate()`; mass flows, fouling, the deposit, the method and every geometry
//! stay the base request's. The required area is
//! `duty_fraction * total_duty / (U * lmtd)`, zero when either denominator is not positive.
//!
//! **The permissive path's unanswered states are refused.** A zone whose transport properties
//! are not positive finite is an error, not a zero or NaN coefficient.
//!
//! Spec: `specs/calcs/thermal/zone_rating.toml`.

use azoth_core::units::{
    Area, DynamicViscosity, Length, MassRate, ThermalConductivity, square_meters,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::{ZoneDefinitionRow, ZoneRatingResult, ZoneRatingRow};
use crate::shell_and_tube_rating::{RatingRequest, ShellSideMethod, evaluate, validate_strict};
use crate::spec_gen;

/// The base geometry, both mass flows, the fouling pair and the optional group, with the
/// zones. The base's eight transport properties are the class defaults; no zone reads them.
#[derive(Debug, Clone)]
pub struct ZoneRatingRequest {
    pub base: RatingRequest,
    pub zones: Vec<ZoneDefinitionRow>,
}

/// Rates every zone in the caller's order.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a base state the strict rating refuses and for a zone
/// whose stated transport properties are not positive finite, or whose duty or LMTD is not
/// finite.
pub fn rate(request: &ZoneRatingRequest) -> Result<ZoneRatingResult> {
    let spec = &spec_gen::ZONE_RATING_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "tube_outside_diameter" => Some(request.base.tube_outside_diameter),
            "tube_inside_diameter" => Some(request.base.tube_inside_diameter),
            "tube_length" => Some(request.base.tube_length),
            "tube_count" => Some(request.base.tube_count),
            "tube_passes" => Some(request.base.tube_passes),
            "tube_pitch" => Some(request.base.tube_pitch),
            "shell_inside_diameter" => Some(request.base.shell_inside_diameter),
            "baffle_spacing" => Some(request.base.baffle_spacing),
            "baffle_count" => Some(request.base.baffle_count),
            "baffle_cut" => Some(request.base.baffle_cut),
            "tube_wall_conductivity" => Some(request.base.tube_wall_conductivity),
            "tube_mass_flow" => Some(request.base.tube_mass_flow),
            "shell_mass_flow" => Some(request.base.shell_mass_flow),
            "fouling_tube" => Some(request.base.fouling_tube),
            "fouling_shell" => Some(request.base.fouling_shell),
            "tube_fouling_layer_thickness" => request.base.tube_fouling_layer_thickness,
            "tube_fouling_layer_conductivity" => request.base.tube_fouling_layer_conductivity,
            "sealing_pairs" => request.base.sealing_pairs,
            "tube_to_baffle_clearance" => request.base.tube_to_baffle_clearance,
            "shell_to_baffle_clearance" => request.base.shell_to_baffle_clearance,
            "bypass_area" => request.base.bypass_area,
            _ => None,
        },
        &mut warnings,
    )?;
    validate_strict(&request.base)?;

    let mut zones = Vec::with_capacity(request.zones.len());
    for zone in &request.zones {
        check_zone(zone)?;
        let mut zone_request = request.base.clone();
        zone_request.tube_density = zone.tube_density.value;
        zone_request.tube_viscosity = zone.tube_viscosity.value;
        zone_request.tube_cp = zone.tube_cp.value;
        zone_request.tube_conductivity = zone.tube_conductivity.value;
        zone_request.shell_density = zone.shell_density.value;
        zone_request.shell_viscosity = zone.shell_viscosity.value;
        zone_request.shell_cp = zone.shell_cp.value;
        zone_request.shell_conductivity = zone.shell_conductivity.value;

        let out = evaluate(&zone_request)?;
        let u = out.overall_heat_transfer_coefficient.value;
        let lmtd = zone.lmtd.value;
        let area = if u > 0.0 && lmtd > 0.0 {
            zone.duty_fraction * zone.total_duty.value / (u * lmtd)
        } else {
            0.0
        };
        zones.push(ZoneRatingRow {
            zone_name: zone.zone_name.clone(),
            duty_fraction: zone.duty_fraction,
            tube_side_heat_transfer_coefficient: out.tube_side_heat_transfer_coefficient,
            shell_side_heat_transfer_coefficient: out.shell_side_heat_transfer_coefficient,
            overall_heat_transfer_coefficient: out.overall_heat_transfer_coefficient,
            lmtd: zone.lmtd,
            required_area: square_meters(area),
        });
    }

    Ok(ZoneRatingResult {
        zone_count: zones.len() as f64,
        zones,
        warnings,
    })
}

fn check_zone(zone: &ZoneDefinitionRow) -> Result<()> {
    let properties = [
        ("tube_density", zone.tube_density.value),
        ("tube_viscosity", zone.tube_viscosity.value),
        ("tube_cp", zone.tube_cp.value),
        ("tube_conductivity", zone.tube_conductivity.value),
        ("shell_density", zone.shell_density.value),
        ("shell_viscosity", zone.shell_viscosity.value),
        ("shell_cp", zone.shell_cp.value),
        ("shell_conductivity", zone.shell_conductivity.value),
    ];
    for (field, value) in properties {
        if !value.is_finite() || value <= 0.0 {
            return Err(AzothError::invalid_input(
                field,
                format!(
                    "the zone {} states {value}, which is not positive finite",
                    zone.zone_name
                ),
            ));
        }
    }
    if !zone.duty_fraction.is_finite()
        || !zone.total_duty.value.is_finite()
        || !zone.lmtd.value.is_finite()
    {
        return Err(AzothError::invalid_input(
            "zones",
            format!(
                "the zone {} states a duty or LMTD that is not finite",
                zone.zone_name
            ),
        ));
    }
    Ok(())
}

/// The flat entry point the binding and the spec's case runner call.
///
/// One parameter per declared input, in the spec's own order; the optional group precedes
/// the zones.
///
/// # Errors
/// As [`rate`].
#[allow(clippy::too_many_arguments)]
pub fn zone_rating(
    tube_outside_diameter: Length,
    tube_inside_diameter: Length,
    tube_length: Length,
    tube_count: f64,
    tube_passes: f64,
    tube_pitch: Length,
    triangular_pitch: bool,
    shell_inside_diameter: Length,
    baffle_spacing: Length,
    baffle_count: f64,
    baffle_cut: f64,
    tube_wall_conductivity: ThermalConductivity,
    tube_mass_flow: MassRate,
    shell_mass_flow: MassRate,
    fouling_tube: f64,
    fouling_shell: f64,
    shell_side_method: ShellSideMethod,
    tube_fouling_layer_thickness: Option<Length>,
    tube_fouling_layer_conductivity: Option<ThermalConductivity>,
    shell_viscosity_wall: Option<DynamicViscosity>,
    has_sealing: Option<bool>,
    sealing_pairs: Option<f64>,
    tube_to_baffle_clearance: Option<Length>,
    shell_to_baffle_clearance: Option<Length>,
    bypass_area: Option<Area>,
    zones: &[ZoneDefinitionRow],
) -> Result<ZoneRatingResult> {
    let defaults = RatingRequest::default();
    rate(&ZoneRatingRequest {
        base: RatingRequest {
            tube_outside_diameter: tube_outside_diameter.value,
            tube_inside_diameter: tube_inside_diameter.value,
            tube_length: tube_length.value,
            tube_count,
            tube_passes,
            tube_pitch: tube_pitch.value,
            triangular_pitch,
            shell_inside_diameter: shell_inside_diameter.value,
            baffle_spacing: baffle_spacing.value,
            baffle_count,
            baffle_cut,
            tube_wall_conductivity: tube_wall_conductivity.value,
            tube_mass_flow: tube_mass_flow.value,
            shell_mass_flow: shell_mass_flow.value,
            fouling_tube,
            fouling_shell,
            shell_side_method,
            tube_fouling_layer_thickness: tube_fouling_layer_thickness.map(|value| value.value),
            tube_fouling_layer_conductivity: tube_fouling_layer_conductivity
                .map(|value| value.value),
            shell_viscosity_wall: shell_viscosity_wall.map(|value| value.value),
            has_sealing,
            sealing_pairs,
            tube_to_baffle_clearance: tube_to_baffle_clearance.map(|value| value.value),
            shell_to_baffle_clearance: shell_to_baffle_clearance.map(|value| value.value),
            bypass_area: bypass_area.map(|value| value.value),
            // No zone reads the base's transport properties; the class defaults are carried
            // only to keep the request well formed.
            tube_density: defaults.tube_density,
            tube_viscosity: defaults.tube_viscosity,
            tube_cp: defaults.tube_cp,
            tube_conductivity: defaults.tube_conductivity,
            shell_density: defaults.shell_density,
            shell_viscosity: defaults.shell_viscosity,
            shell_cp: defaults.shell_cp,
            shell_conductivity: defaults.shell_conductivity,
        },
        zones: zones.to_vec(),
    })
}
