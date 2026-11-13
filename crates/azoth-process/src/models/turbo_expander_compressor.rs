//! `process.turbo_expander_compressor` - the coupled machine as a registered id.
//!
//! Spec: `specs/models/process/turbo_expander_compressor.toml`. The arithmetic is
//! [`crate::kernels::turbo_expander_compressor`]; what is here is the boundary a case and a cross-impl
//! test address.

use azoth_core::units::{
    Area, Length, Mass, MolarEnergy, Power, Pressure, ThermodynamicTemperature, Velocity,
    VolumeRate, cubic_meters_per_hour, joules_per_mole, meters, pascals, watts,
};
use azoth_core::{Result, Warning, apply_checks};

use crate::expander_chart::{CurvePoint as MapPoint, ExpanderChart};
use crate::kernels::turbo_expander_compressor::{
    self as kernel, DEFAULT_BEARING_LOSS, DEFAULT_COMPRESSOR_DESIGN_EFFICIENCY,
    DEFAULT_COMPRESSOR_DESIGN_HEAD_M, DEFAULT_EXPANDER_DESIGN_EFFICIENCY,
    DEFAULT_IGV_AREA_INCREASE_FACTOR, DEFAULT_IGV_OPENING, DEFAULT_MAXIMUM_IGV_AREA, TurboSettings,
};
use crate::model_gen;
use crate::rotating_design::turbo_expander_compressor::design as design_tec;
use crate::stream::Stream;

/// One row of `uc_curve`.
#[derive(Debug, Clone, PartialEq)]
pub struct TurboUcPoint {
    /// The velocity ratio.
    pub uc: f64,
    /// The efficiency there.
    pub efficiency: f64,
}

/// One row of the two Q/N curves.
#[derive(Debug, Clone, PartialEq)]
pub struct TurboQnPoint {
    /// The flow coefficient ratio.
    pub qn: f64,
    /// The correction factor there.
    pub factor: f64,
}

/// One row of `chart`: an IGV curve point, efficiency alone.
#[derive(Debug, Clone, PartialEq)]
pub struct ExpanderMapPoint {
    /// The IGV opening.
    pub igv: f64,
    /// The velocity ratio.
    pub uc: f64,
    /// The isentropic efficiency.
    pub efficiency: f64,
}

/// One row of `igv_penalty_curve`.
#[derive(Debug, Clone, PartialEq)]
pub struct TurboIgvPenaltyPoint {
    /// The IGV opening.
    pub igv: f64,
    /// The multiplicative factor there.
    pub factor: f64,
}

/// Result of `process.turbo_expander_compressor`.
#[derive(Debug, Clone, PartialEq)]
pub struct TurboExpanderResult {
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
    /// The solved shaft speed, rpm.
    pub speed: f64,
    /// The expander's power, W.
    pub power_expander: Power,
    /// The compressor's absorbed power, W.
    pub power_compressor: Power,
    /// The expander efficiency the balance used.
    pub expander_isentropic_efficiency: f64,
    /// The design efficiency, back-calculated in the inverse.
    pub expander_design_isentropic_efficiency: f64,
    /// The compressor's polytropic efficiency.
    pub compressor_polytropic_efficiency: f64,
    /// The compressor's polytropic head, m.
    pub compressor_polytropic_head: Length,
    /// The expander's `U/C/design_uc`.
    pub uc_ratio_expander: f64,
    /// The compressor's Q/N ratio.
    pub qn_ratio_compressor: f64,
    /// The IGV opening the machine reports.
    pub igv_opening: f64,
    /// Mechanical design output.
    pub design_expander_wheel_diameter: Option<Length>,
    /// Mechanical design output.
    pub design_expander_tip_speed: Option<Velocity>,
    /// Mechanical design output.
    pub design_expander_type: Option<String>,
    /// Mechanical design output.
    pub design_expander_stages: Option<f64>,
    /// Mechanical design output.
    pub design_expander_design_inlet_pressure: Option<Pressure>,
    /// Mechanical design output.
    pub design_expander_design_outlet_pressure: Option<Pressure>,
    /// Mechanical design output.
    pub design_expander_casing_pressure: Option<Pressure>,
    /// Mechanical design output.
    pub design_expander_casing_temperature: Option<ThermodynamicTemperature>,
    /// Mechanical design output.
    pub design_expander_power: Option<Power>,
    /// Mechanical design output.
    pub design_expander_enthalpy_drop: Option<f64>,
    /// Mechanical design output.
    pub design_compressor_impeller_diameter: Option<Length>,
    /// Mechanical design output.
    pub design_compressor_tip_speed: Option<Velocity>,
    /// Mechanical design output.
    pub design_compressor_stages: Option<f64>,
    /// Mechanical design output.
    pub design_compressor_polytropic_head_design: Option<f64>,
    /// Mechanical design output.
    pub design_compressor_head_per_stage: Option<f64>,
    /// Mechanical design output.
    pub design_compressor_casing_type: Option<String>,
    /// Mechanical design output.
    pub design_compressor_design_suction_pressure: Option<Pressure>,
    /// Mechanical design output.
    pub design_compressor_design_discharge_pressure: Option<Pressure>,
    /// Mechanical design output.
    pub design_compressor_casing_pressure: Option<Pressure>,
    /// Mechanical design output.
    pub design_compressor_casing_temperature: Option<ThermodynamicTemperature>,
    /// Mechanical design output.
    pub design_compressor_power: Option<Power>,
    /// Mechanical design output.
    pub design_shaft_diameter: Option<Length>,
    /// Mechanical design output.
    pub design_shaft_weight: Option<Mass>,
    /// Mechanical design output.
    pub design_operating_speed: Option<f64>,
    /// Mechanical design output.
    pub design_max_continuous_speed: Option<f64>,
    /// Mechanical design output.
    pub design_trip_speed: Option<f64>,
    /// Mechanical design output.
    pub design_first_critical_speed: Option<f64>,
    /// Mechanical design output.
    pub design_second_critical_speed: Option<f64>,
    /// Mechanical design output.
    pub design_bearing_span: Option<Length>,
    /// Mechanical design output.
    pub design_bearing_type: Option<String>,
    /// Mechanical design output.
    pub design_expander_seal_type: Option<String>,
    /// Mechanical design output.
    pub design_compressor_seal_type: Option<String>,
    /// Mechanical design output.
    pub design_gear_ratio: Option<f64>,
    /// Mechanical design output.
    pub design_torque: Option<f64>,
    /// Mechanical design output.
    pub design_shaft_stress: Option<Pressure>,
    /// Mechanical design output.
    pub design_surge_flow: Option<VolumeRate>,
    /// Mechanical design output.
    pub design_stonewall_flow: Option<VolumeRate>,
    /// Mechanical design output.
    pub design_expander_nozzle_throat_area: Option<Area>,
    /// Mechanical design output.
    pub design_expander_wheel_weight: Option<Mass>,
    /// Mechanical design output.
    pub design_expander_casing_weight: Option<Mass>,
    /// Mechanical design output.
    pub design_compressor_impeller_weight: Option<Mass>,
    /// Mechanical design output.
    pub design_compressor_casing_weight: Option<Mass>,
    /// Mechanical design output.
    pub design_weight_vessel: Option<Mass>,
    /// Mechanical design output.
    pub design_weight_internals: Option<Mass>,
    /// Mechanical design output.
    pub design_weight_piping: Option<Mass>,
    /// Mechanical design output.
    pub design_weight_electro_instrument: Option<Mass>,
    /// Mechanical design output.
    pub design_weight_structural_steel: Option<Mass>,
    /// Mechanical design output.
    pub design_weight_total: Option<Mass>,
    /// Mechanical design output.
    pub design_module_length: Option<Length>,
    /// Mechanical design output.
    pub design_module_width: Option<Length>,
    /// Mechanical design output.
    pub design_module_height: Option<Length>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

/// Balance an expander and a compressor on one shaft.
///
/// # Errors
/// [`azoth_core::AzothError::InvalidInput`] for a non-positive design value, an efficiency
/// outside `(0, 1]`, a non-positive isentropic drop or a map that is not a map; and whatever
/// the flash refuses.
#[allow(clippy::too_many_arguments)] // one parameter per declared input
pub fn turbo_expander_compressor(
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
    impeller_diameter: Length,
    design_speed: f64,
    design_uc: f64,
    design_qn: f64,
    igv_opening: Option<f64>,
    igv_control_mode: Option<bool>,
    maximum_igv_area: Option<Area>,
    igv_area_increase_factor: Option<f64>,
    bearing_loss: Option<Power>,
    expander_design_isentropic_efficiency: Option<f64>,
    compressor_design_polytropic_efficiency: Option<f64>,
    compressor_design_polytropic_head: Option<Length>,
    expander_out_temperature: Option<ThermodynamicTemperature>,
    uc_curve: Option<&[TurboUcPoint]>,
    qn_efficiency_curve: Option<&[TurboQnPoint]>,
    qn_head_curve: Option<&[TurboQnPoint]>,
    chart: Option<&[ExpanderMapPoint]>,
    igv_penalty_curve: Option<&[TurboIgvPenaltyPoint]>,
    mechanical_design: Option<bool>,
) -> Result<TurboExpanderResult> {
    let spec = &model_gen::TURBO_EXPANDER_COMPRESSOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "design_speed" => Some(design_speed),
            "design_uc" => Some(design_uc),
            "design_qn" => Some(design_qn),
            "impeller_diameter" => Some(impeller_diameter.value),
            "expander_design_isentropic_efficiency" => expander_design_isentropic_efficiency,
            "compressor_design_polytropic_efficiency" => compressor_design_polytropic_efficiency,
            _ => None,
        },
        &mut warnings,
    )?;

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

    let uc: Vec<(f64, f64)> = uc_curve
        .unwrap_or(&[])
        .iter()
        .map(|row| (row.uc, row.efficiency))
        .collect();
    let qn_efficiency: Vec<(f64, f64)> = qn_efficiency_curve
        .unwrap_or(&[])
        .iter()
        .map(|row| (row.qn, row.factor))
        .collect();
    let qn_head: Vec<(f64, f64)> = qn_head_curve
        .unwrap_or(&[])
        .iter()
        .map(|row| (row.qn, row.factor))
        .collect();
    let penalty: Vec<(f64, f64)> = igv_penalty_curve
        .unwrap_or(&[])
        .iter()
        .map(|row| (row.igv, row.factor))
        .collect();
    let map = match chart {
        None => None,
        Some(rows) => {
            let points: Vec<MapPoint> = rows
                .iter()
                .map(|row| MapPoint {
                    igv: row.igv,
                    uc: row.uc,
                    efficiency: row.efficiency,
                    head: 0.0,
                })
                .collect();
            Some(ExpanderChart::new(&points, 1.0)?)
        }
    };

    let settings = TurboSettings {
        expander_out_pressure,
        impeller_diameter: impeller_diameter.value,
        design_speed,
        design_uc,
        design_qn,
        maximum_igv_area: maximum_igv_area.map_or(DEFAULT_MAXIMUM_IGV_AREA, |area| area.value),
        igv_area_increase_factor: igv_area_increase_factor
            .unwrap_or(DEFAULT_IGV_AREA_INCREASE_FACTOR),
        igv_opening: igv_opening.unwrap_or(DEFAULT_IGV_OPENING),
        igv_control_mode: igv_control_mode.unwrap_or(false),
        bearing_loss: bearing_loss.map_or(DEFAULT_BEARING_LOSS, |power| power.value),
        expander_design_efficiency: expander_design_isentropic_efficiency
            .unwrap_or(DEFAULT_EXPANDER_DESIGN_EFFICIENCY),
        compressor_design_efficiency: compressor_design_polytropic_efficiency
            .unwrap_or(DEFAULT_COMPRESSOR_DESIGN_EFFICIENCY),
        compressor_design_head_m: compressor_design_polytropic_head
            .map_or(DEFAULT_COMPRESSOR_DESIGN_HEAD_M, |head| head.value),
        uc_curve: &uc,
        qn_efficiency_curve: &qn_efficiency,
        qn_head_curve: &qn_head,
        chart: map.as_ref(),
        igv_penalty_curve: &penalty,
        expander_out_temperature,
    };

    let outcome = kernel::turbo_expander(&expander_feed, &compressor_feed, &settings)?;

    let design = if mechanical_design.unwrap_or(false) {
        let compressor_head_kj_kg = outcome.compressor_head_m * 9.81 / 1000.0;
        Some(design_tec(
            &expander_feed,
            &outcome.expander_outlet,
            &compressor_feed,
            &outcome.compressor_outlet,
            watts(outcome.power_expander),
            watts(outcome.power_compressor),
            compressor_head_kj_kg,
            outcome.speed,
            1.0,
        )?)
    } else {
        None
    };
    let d = design.as_ref();

    Ok(TurboExpanderResult {
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
        speed: outcome.speed,
        power_expander: watts(outcome.power_expander),
        power_compressor: watts(outcome.power_compressor),
        expander_isentropic_efficiency: outcome.expander_efficiency,
        expander_design_isentropic_efficiency: outcome.expander_design_efficiency,
        compressor_polytropic_efficiency: outcome.compressor_efficiency,
        compressor_polytropic_head: meters(outcome.compressor_head_m),
        uc_ratio_expander: outcome.uc_ratio,
        qn_ratio_compressor: outcome.qn_ratio,
        igv_opening: outcome.igv_opening,
        design_expander_wheel_diameter: d.map(|value| value.expander_wheel_diameter),
        design_expander_tip_speed: d.map(|value| value.expander_tip_speed),
        design_expander_type: d.map(|value| value.expander_type.as_str().to_owned()),
        design_expander_stages: d.map(|value| f64::from(value.expander_stages)),
        design_expander_design_inlet_pressure: d.map(|value| value.expander_design_inlet_pressure),
        design_expander_design_outlet_pressure: d
            .map(|value| value.expander_design_outlet_pressure),
        design_expander_casing_pressure: d.map(|value| value.expander_casing_pressure),
        design_expander_casing_temperature: d.map(|value| value.expander_casing_temperature),
        design_expander_power: d.map(|value| value.expander_power),
        design_expander_enthalpy_drop: d.map(|value| value.expander_enthalpy_drop),
        design_compressor_impeller_diameter: d.map(|value| value.compressor_impeller_diameter),
        design_compressor_tip_speed: d.map(|value| value.compressor_tip_speed),
        design_compressor_stages: d.map(|value| f64::from(value.compressor_stages)),
        design_compressor_polytropic_head_design: d.map(|value| value.compressor_polytropic_head),
        design_compressor_head_per_stage: d.map(|value| value.compressor_head_per_stage),
        design_compressor_casing_type: d
            .map(|value| value.compressor_casing_type.as_str().to_owned()),
        design_compressor_design_suction_pressure: d
            .map(|value| value.compressor_design_suction_pressure),
        design_compressor_design_discharge_pressure: d
            .map(|value| value.compressor_design_discharge_pressure),
        design_compressor_casing_pressure: d.map(|value| value.compressor_casing_pressure),
        design_compressor_casing_temperature: d.map(|value| value.compressor_casing_temperature),
        design_compressor_power: d.map(|value| value.compressor_power),
        design_shaft_diameter: d.map(|value| value.shaft_diameter),
        design_shaft_weight: d.map(|value| value.shaft_weight),
        design_operating_speed: d.map(|value| value.operating_speed),
        design_max_continuous_speed: d.map(|value| value.max_continuous_speed),
        design_trip_speed: d.map(|value| value.trip_speed),
        design_first_critical_speed: d.map(|value| value.first_critical_speed),
        design_second_critical_speed: d.map(|value| value.second_critical_speed),
        design_bearing_span: d.map(|value| value.bearing_span),
        design_bearing_type: d.map(|value| value.bearing_type.to_owned()),
        design_expander_seal_type: d.map(|value| value.expander_seal_type.to_owned()),
        design_compressor_seal_type: d.map(|value| value.compressor_seal_type.to_owned()),
        design_gear_ratio: d.map(|value| value.gear_ratio),
        design_torque: d.map(|value| value.design_torque),
        design_shaft_stress: d.map(|value| pascals(value.design_shaft_stress * 1.0e6)),
        design_surge_flow: d.map(|value| cubic_meters_per_hour(value.design_surge_flow)),
        design_stonewall_flow: d.map(|value| cubic_meters_per_hour(value.design_stonewall_flow)),
        design_expander_nozzle_throat_area: d
            .map(|value| azoth_core::units::square_meters(value.expander_nozzle_throat_area)),
        design_expander_wheel_weight: d.map(|value| value.expander_wheel_weight),
        design_expander_casing_weight: d.map(|value| value.expander_casing_weight),
        design_compressor_impeller_weight: d.map(|value| value.compressor_impeller_weight),
        design_compressor_casing_weight: d.map(|value| value.compressor_casing_weight),
        design_weight_vessel: d.map(|value| value.weight_vessel),
        design_weight_internals: d.map(|value| value.weight_internals),
        design_weight_piping: d.map(|value| value.weight_piping),
        design_weight_electro_instrument: d.map(|value| value.weight_electro_instrument),
        design_weight_structural_steel: d.map(|value| value.weight_structural_steel),
        design_weight_total: d.map(|value| value.weight_total),
        design_module_length: d.map(|value| value.module_length),
        design_module_width: d.map(|value| value.module_width),
        design_module_height: d.map(|value| value.module_height),
        warnings,
    })
}
