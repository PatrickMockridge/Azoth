//! `process.turbo_expander_compressor` - the coupled machine as a registered id.
//!
//! Spec: `specs/models/process/turbo_expander_compressor.toml`. The arithmetic is
//! [`crate::kernels::turbo_expander_compressor`]; what is here is the boundary a case and a cross-impl
//! test address.

use azoth_core::units::{
    Area, Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole, meters,
    watts,
};
use azoth_core::{Result, Warning, apply_checks};

use crate::expander_chart::{CurvePoint as MapPoint, ExpanderChart};
use crate::kernels::turbo_expander_compressor::{
    self as kernel, DEFAULT_BEARING_LOSS, DEFAULT_COMPRESSOR_DESIGN_EFFICIENCY,
    DEFAULT_COMPRESSOR_DESIGN_HEAD_M, DEFAULT_EXPANDER_DESIGN_EFFICIENCY,
    DEFAULT_IGV_AREA_INCREASE_FACTOR, DEFAULT_IGV_OPENING, DEFAULT_MAXIMUM_IGV_AREA, TurboSettings,
};
use crate::model_gen;
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
        warnings,
    })
}
