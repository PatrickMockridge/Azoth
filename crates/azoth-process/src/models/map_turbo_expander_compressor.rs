//! `process.map_turbo_expander_compressor` - the map-based single shaft as a registered id.
//!
//! Spec: `specs/models/process/map_turbo_expander_compressor.toml`. The arithmetic is
//! [`crate::kernels::map_turbo_expander_compressor`]; what is here is the boundary a case and a
//! cross-impl test address.

use azoth_core::units::{
    Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, VolumeRate,
    cubic_meter_per_hour, joules_per_mole, meter, watts,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::compressor_chart::{CompressorChart, CurvePoint};
use crate::kernels::map_turbo_expander_compressor::{
    DEFAULT_BEARING_LOSS, DEFAULT_COMPRESSOR_SPEED, DEFAULT_MECHANICAL_EFFICIENCY,
    MapMachineSettings, ShaftMode, map_turbo_expander,
};
use crate::model_gen;
use crate::models::compressor::{build_chart, chart_kind};
use crate::stream::Stream;

/// One row of `curve`: a measured point, in the units the spec declares.
#[derive(Debug, Clone, PartialEq)]
pub struct MapCurvePoint {
    /// The speed the point was measured at, rpm.
    pub speed: f64,
    /// The measured actual flow.
    pub flow: VolumeRate,
    /// The measured polytropic head, in metres.
    pub head: Length,
    /// The measured polytropic efficiency, as a fraction.
    pub efficiency: f64,
}

/// Result of `process.map_turbo_expander_compressor`.
#[derive(Debug, Clone, PartialEq)]
pub struct MapTurboExpanderResult {
    /// The expander outlet's molar flow.
    pub expander_outlet_n: f64,
    /// The expander outlet's composition.
    pub expander_outlet_z: Vec<f64>,
    /// The expander outlet's pressure.
    pub expander_outlet_p: Pressure,
    /// The expander outlet's temperature.
    pub expander_outlet_t: ThermodynamicTemperature,
    /// The expander outlet's molar enthalpy.
    pub expander_outlet_h: MolarEnergy,
    /// The compressor outlet's molar flow.
    pub compressor_outlet_n: f64,
    /// The compressor outlet's composition.
    pub compressor_outlet_z: Vec<f64>,
    /// The compressor outlet's pressure.
    pub compressor_outlet_p: Pressure,
    /// The compressor outlet's temperature.
    pub compressor_outlet_t: ThermodynamicTemperature,
    /// The compressor outlet's molar enthalpy.
    pub compressor_outlet_h: MolarEnergy,
    /// The shaft speed, rpm.
    pub shaft_speed: f64,
    /// The expander's power that reaches the shaft, W.
    pub available_shaft_power: Power,
    /// The compressor's absorbed power, W.
    pub consumed_compressor_power: Power,
    /// `available - consumed`, W.
    pub power_balance_residual: Power,
    /// Whether the speed hit a bound.
    pub speed_limited: bool,
    /// Which outcome the shaft reached: `balanced`, `under_power_surge` or
    /// `over_power_max_speed`.
    pub operating_status: String,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Balance a map-driven compressor against an expander on one shaft.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a declaration violation or a balanced-speed
/// run with no map; and whatever the expander, the compressor or the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn map_turbo_expander_compressor(
    expander_components: &[String],
    expander_n: f64,
    expander_z: &[f64],
    expander_p: Pressure,
    expander_t: ThermodynamicTemperature,
    compressor_components: &[String],
    compressor_n: f64,
    compressor_z: &[f64],
    compressor_p: Pressure,
    compressor_t: ThermodynamicTemperature,
    expander_out_pressure: Pressure,
    expander_isentropic_efficiency: f64,
    shaft_mode: Option<&str>,
    mechanical_efficiency: Option<f64>,
    bearing_loss: Option<Power>,
    compressor_outlet_pressure: Option<Pressure>,
    compressor_isentropic_efficiency: Option<f64>,
    compressor_speed: Option<f64>,
    min_shaft_speed: Option<f64>,
    max_shaft_speed: Option<f64>,
    curve: Option<&[MapCurvePoint]>,
    chart_type: Option<&str>,
    chart_conditions: Option<&[f64]>,
    impeller_diameter: Option<Length>,
) -> Result<MapTurboExpanderResult> {
    let spec = &model_gen::MAP_TURBO_EXPANDER_COMPRESSOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "expander_isentropic_efficiency" => Some(expander_isentropic_efficiency),
            "mechanical_efficiency" => mechanical_efficiency,
            "expander_t" => Some(expander_t.value),
            "expander_p" => Some(expander_p.value),
            "compressor_t" => Some(compressor_t.value),
            "compressor_p" => Some(compressor_p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let mode = match shaft_mode.unwrap_or("specified_pressures") {
        "specified_pressures" => ShaftMode::SpecifiedPressures,
        "balanced_speed" => ShaftMode::BalancedSpeed,
        other => {
            return Err(AzothError::invalid_input(
                "shaft_mode",
                format!("'{other}' is not one of the two declared shaft modes"),
            ));
        }
    };

    let expander_feed = Stream::from_pt(
        expander_components.to_vec(),
        expander_z.to_vec(),
        expander_n,
        expander_p,
        expander_t,
    )?;
    let compressor_feed = Stream::from_pt(
        compressor_components.to_vec(),
        compressor_z.to_vec(),
        compressor_n,
        compressor_p,
        compressor_t,
    )?;

    let points: Vec<CurvePoint> = curve
        .unwrap_or(&[])
        .iter()
        .map(|point| CurvePoint {
            speed: point.speed,
            flow: point.flow.get::<cubic_meter_per_hour>(),
            head: point.head.get::<meter>(),
            efficiency: point.efficiency,
        })
        .collect();
    let kind = chart_kind(chart_type)?;
    let chart: Option<CompressorChart> = if points.is_empty() {
        None
    } else {
        Some(build_chart(
            compressor_components,
            compressor_t,
            compressor_p,
            compressor_z,
            &points,
            kind,
            chart_conditions,
            impeller_diameter,
        )?)
    };

    let settings = MapMachineSettings {
        expander_out_pressure,
        expander_efficiency: expander_isentropic_efficiency,
        mode,
        mechanical_efficiency: mechanical_efficiency.unwrap_or(DEFAULT_MECHANICAL_EFFICIENCY),
        bearing_loss: bearing_loss.map_or(DEFAULT_BEARING_LOSS, |power| power.value),
        min_speed: min_shaft_speed,
        max_speed: max_shaft_speed,
        chart: chart.as_ref(),
        compressor_speed: compressor_speed.unwrap_or(DEFAULT_COMPRESSOR_SPEED),
        compressor_outlet_pressure,
        compressor_isentropic_efficiency: compressor_isentropic_efficiency.unwrap_or(0.75),
    };
    let outcome = map_turbo_expander(&expander_feed, &compressor_feed, &settings)?;

    Ok(MapTurboExpanderResult {
        expander_outlet_n: outcome.expander_outlet.n,
        expander_outlet_z: outcome.expander_outlet.z.clone(),
        expander_outlet_p: outcome.expander_outlet.p,
        expander_outlet_t: outcome.expander_outlet.t,
        expander_outlet_h: joules_per_mole(outcome.expander_outlet.h.value),
        compressor_outlet_n: outcome.compressor_outlet.n,
        compressor_outlet_z: outcome.compressor_outlet.z.clone(),
        compressor_outlet_p: outcome.compressor_outlet.p,
        compressor_outlet_t: outcome.compressor_outlet.t,
        compressor_outlet_h: joules_per_mole(outcome.compressor_outlet.h.value),
        shaft_speed: outcome.shaft_speed,
        available_shaft_power: watts(outcome.available_shaft_power),
        consumed_compressor_power: watts(outcome.consumed_compressor_power),
        power_balance_residual: watts(outcome.power_balance_residual),
        speed_limited: outcome.speed_limited,
        operating_status: outcome.status.as_str().to_string(),
        warnings,
    })
}
