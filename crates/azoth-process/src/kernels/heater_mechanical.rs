//! `HeatExchangerMechanicalDesign`'s thermal estimate.
//!
//! **Two branches, selected by the equipment rather than by an input.** The heater path takes
//! the utility specification's two temperatures and the process inlet and outlet; the
//! two-stream path takes the exchanger's four stream temperatures. Either way the pairs give
//! an LMTD and an approach, `UA = duty / LMTD`, and each of the five candidate types turns
//! that into a required area with its own allowable approach and typical overall coefficient
//! and its own geometry model. **The sizing coefficient differs**: the heater path divides by
//! the larger of a used coefficient and the type's, the exchanger path by the type's alone.
//! The detailed `ShellAndTubeDesignCalculator`, the fouling and velocity setters and the
//! vessel pressure boundary are outside.

use std::f64::consts::PI;

/// `UtilityStreamSpecification`'s five values and their predicates.
#[derive(Debug, Clone, Copy, Default)]
pub struct UtilitySpecification {
    /// The utility supply temperature, K.
    pub supply_temperature: Option<f64>,
    /// The utility return temperature, K.
    pub return_temperature: Option<f64>,
    /// The minimum approach temperature, K.
    pub approach_temperature: Option<f64>,
    /// The heat-capacity rate, W/K.
    pub heat_capacity_rate: Option<f64>,
    /// The overall heat-transfer coefficient, W/(m2 K).
    pub overall_heat_transfer_coefficient: Option<f64>,
}

impl UtilitySpecification {
    /// `hasSupplyTemperature`.
    #[must_use]
    pub fn has_supply_temperature(&self) -> bool {
        self.supply_temperature.is_some()
    }

    /// `hasReturnTemperature`.
    #[must_use]
    pub fn has_return_temperature(&self) -> bool {
        self.return_temperature.is_some()
    }

    /// `hasApproachTemperature`.
    #[must_use]
    pub fn has_approach_temperature(&self) -> bool {
        self.approach_temperature.is_some()
    }

    /// `hasHeatCapacityRate`: stated **and positive**.
    #[must_use]
    pub fn has_heat_capacity_rate(&self) -> bool {
        self.heat_capacity_rate.is_some_and(|value| value > 0.0)
    }

    /// `hasOverallHeatTransferCoefficient`: stated **and positive**.
    #[must_use]
    pub fn has_overall_heat_transfer_coefficient(&self) -> bool {
        self.overall_heat_transfer_coefficient
            .is_some_and(|value| value > 0.0)
    }
}

/// The five candidate configurations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MechanicalType {
    /// Shell-and-tube with carbon-steel tubes.
    ShellAndTube,
    /// Plate-and-frame.
    PlateAndFrame,
    /// Air-cooled fin fan.
    AirCooler,
    /// Double pipe.
    DoublePipe,
    /// Plate-fin / brazed aluminium.
    PlateFin,
}

impl MechanicalType {
    /// The class's string spelling.
    #[must_use]
    pub fn parse(spelling: &str) -> Option<Self> {
        Some(match spelling {
            "shell_and_tube" => Self::ShellAndTube,
            "plate_and_frame" => Self::PlateAndFrame,
            "air_cooler" => Self::AirCooler,
            "double_pipe" => Self::DoublePipe,
            "plate_fin" => Self::PlateFin,
            _ => return None,
        })
    }

    /// The enum's spelling.
    #[must_use]
    pub fn spelling(self) -> &'static str {
        match self {
            Self::ShellAndTube => "shell_and_tube",
            Self::PlateAndFrame => "plate_and_frame",
            Self::AirCooler => "air_cooler",
            Self::DoublePipe => "double_pipe",
            Self::PlateFin => "plate_fin",
        }
    }

    /// `getTypicalOverallHeatTransferCoefficient`.
    #[must_use]
    pub fn typical_u(self) -> f64 {
        match self {
            Self::ShellAndTube => 580.0,
            Self::PlateAndFrame => 1150.0,
            Self::AirCooler => 90.0,
            Self::DoublePipe => 350.0,
            Self::PlateFin => 600.0,
        }
    }

    /// `getAllowableApproachTemperature`.
    #[must_use]
    pub fn allowable_approach(self) -> f64 {
        match self {
            Self::ShellAndTube => 5.0,
            Self::PlateAndFrame => 2.0,
            Self::AirCooler => 12.0,
            Self::DoublePipe => 8.0,
            Self::PlateFin => 1.5,
        }
    }
}

/// `SelectionCriterion`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SelectionCriterion {
    /// Smallest required area.
    MinArea,
    /// Lowest estimated weight.
    MinWeight,
    /// Lowest estimated pressure drop.
    MinPressureDrop,
}

impl SelectionCriterion {
    /// The class's string spelling.
    #[must_use]
    pub fn parse(spelling: &str) -> Option<Self> {
        Some(match spelling {
            "min_area" => Self::MinArea,
            "min_weight" => Self::MinWeight,
            "min_pressure_drop" => Self::MinPressureDrop,
            _ => return None,
        })
    }
}

/// The design's selection inputs.
#[derive(Debug, Clone, Copy, Default)]
pub struct MechanicalRequest {
    /// The manual type selection; absent evaluates all five and ranks them.
    pub manual_selection: Option<MechanicalType>,
    /// The automatic ranking's metric; absent is the class's `MIN_AREA`.
    pub selection_criterion: Option<SelectionCriterion>,
}

/// One candidate's `HeatExchangerSizingResult`.
#[derive(Debug, Clone, Copy)]
pub struct TypeSizing {
    /// The configuration.
    pub kind: MechanicalType,
    /// The required area, m2.
    pub required_area_m2: f64,
    /// The required UA, W/K.
    pub required_ua_w_per_k: f64,
    /// The type's typical overall coefficient, W/(m2 K).
    pub overall_heat_transfer_coefficient: f64,
    /// The type's approach, K.
    pub approach_temperature_k: f64,
    /// The estimated length, m.
    pub estimated_length_m: f64,
    /// The equivalent inside diameter, m.
    pub inner_diameter_m: f64,
    /// The outside diameter, m.
    pub outer_diameter_m: f64,
    /// The wall thickness, m.
    pub wall_thickness_m: f64,
    /// The tube or plate count.
    pub tube_count: f64,
    /// The tube passes or plate pairs.
    pub tube_passes: f64,
    /// The fin surface area, m2; zero for the unfinned types.
    pub fin_surface_area_m2: f64,
    /// The class's qualitative pressure-drop proxy, bar.
    pub estimated_pressure_drop_bar: f64,
    /// The estimated weight, kg.
    pub estimated_weight_kg: f64,
    /// The module length, m.
    pub module_length_m: f64,
    /// The module width, m.
    pub module_width_m: f64,
    /// The module height, m.
    pub module_height_m: f64,
}

/// What the thermal estimate computed.
#[derive(Debug, Clone)]
pub struct MechanicalDesign {
    /// The log-mean temperature difference, K.
    pub log_mean_temperature_difference: f64,
    /// The approach temperature, K.
    pub approach_temperature: f64,
    /// The overall coefficient the estimate used, W/(m2 K).
    pub used_overall_heat_transfer_coefficient: f64,
    /// The calculated UA, W/K.
    pub calculated_ua: f64,
    /// The five candidates, in the class's enum order.
    pub types: Vec<TypeSizing>,
    /// The selected row.
    pub selected: TypeSizing,
}

fn calculate_lmtd(delta_t1: f64, delta_t2: f64) -> f64 {
    let safe1 = delta_t1.max(1.0e-6);
    let safe2 = delta_t2.max(1.0e-6);
    if (safe1 - safe2).abs() < 1.0e-6 {
        return (safe1 + safe2) / 2.0;
    }
    ((safe1 - safe2) / (safe1 / safe2).ln()).abs()
}

/// `determineSingleStreamTemperatures`, returning `(utility_supply, utility_return,
/// process_in, process_out)` in the class's heating order, or `None` when the utility side is
/// not sufficiently specified.
fn single_stream_temperatures(
    process_in: f64,
    process_out: f64,
    duty: f64,
    utility: &UtilitySpecification,
) -> Option<[f64; 4]> {
    let heating = process_out >= process_in;
    let mut supply = f64::NAN;
    let mut ret = f64::NAN;
    if utility.has_supply_temperature() {
        supply = utility.supply_temperature.unwrap_or(f64::NAN);
    } else if utility.has_approach_temperature() {
        let approach = utility.approach_temperature.unwrap_or(0.0);
        supply = if heating {
            process_out + approach
        } else {
            process_out - approach
        };
    }
    if utility.has_return_temperature() {
        ret = utility.return_temperature.unwrap_or(f64::NAN);
    }
    if ret.is_nan() && utility.has_heat_capacity_rate() && duty > 0.0 && !supply.is_nan() {
        let delta = duty / utility.heat_capacity_rate.unwrap_or(1.0);
        ret = if heating {
            supply - delta
        } else {
            supply + delta
        };
    }
    if ret.is_nan() && !supply.is_nan() {
        ret = supply;
    }
    if supply.is_nan() || ret.is_nan() {
        return None;
    }
    Some(if heating {
        [supply, ret, process_in, process_out]
    } else {
        [process_in, process_out, supply, ret]
    })
}

fn size_type(
    kind: MechanicalType,
    required_area: f64,
    required_ua: f64,
    approach: f64,
) -> TypeSizing {
    let mut row = TypeSizing {
        kind,
        required_area_m2: required_area,
        required_ua_w_per_k: required_ua,
        overall_heat_transfer_coefficient: kind.typical_u(),
        approach_temperature_k: approach,
        estimated_length_m: 0.0,
        inner_diameter_m: 0.0,
        outer_diameter_m: 0.0,
        wall_thickness_m: 0.0,
        tube_count: 0.0,
        tube_passes: 0.0,
        fin_surface_area_m2: 0.0,
        estimated_pressure_drop_bar: 0.0,
        estimated_weight_kg: 0.0,
        module_length_m: 0.0,
        module_width_m: 0.0,
        module_height_m: 0.0,
    };
    match kind {
        MechanicalType::ShellAndTube => {
            let tube_length = 6.0_f64.max((required_area / (PI * 0.019 * 40.0)).min(12.0));
            let count = (required_area / (PI * 0.019 * tube_length)).ceil().max(2.0);
            let passes = if count > 100.0 { 4.0 } else { 2.0 };
            let pitch = 0.019 * 1.25;
            let bundle = count.sqrt() * pitch;
            let inner = bundle.max(0.3);
            let outer = inner + 2.0 * 0.012;
            let wetted = PI * inner * tube_length;
            row.estimated_length_m = tube_length;
            row.inner_diameter_m = inner;
            row.outer_diameter_m = outer;
            row.wall_thickness_m = 0.012;
            row.tube_count = count;
            row.tube_passes = passes;
            row.estimated_pressure_drop_bar = 0.05 * passes * tube_length;
            row.estimated_weight_kg = wetted * 0.012 * 7850.0;
            row.module_length_m = tube_length + 1.5;
            row.module_width_m = outer + 0.8;
            row.module_height_m = outer + 0.8;
        }
        MechanicalType::PlateAndFrame => {
            let pairs = (required_area / 2.2).ceil().max(5.0);
            let plates = pairs * 2.0 + 1.0;
            let stack = 0.6 + plates * 0.003;
            let equivalent = (required_area / PI).sqrt();
            row.estimated_length_m = stack;
            row.inner_diameter_m = equivalent;
            row.outer_diameter_m = equivalent + 0.016;
            row.wall_thickness_m = 0.008;
            row.tube_count = plates;
            row.tube_passes = pairs;
            row.estimated_pressure_drop_bar = 0.02 * pairs;
            row.estimated_weight_kg = plates * 12.0;
            row.module_length_m = stack;
            row.module_width_m = 1.0;
            row.module_height_m = 1.8;
        }
        MechanicalType::AirCooler => {
            let plan = (required_area / 2.5).max(5.0);
            let width = plan.sqrt();
            let length = plan / width;
            let equivalent = (4.0 * plan / PI).sqrt();
            row.estimated_length_m = length;
            row.inner_diameter_m = equivalent;
            row.outer_diameter_m = equivalent + 0.012;
            row.wall_thickness_m = 0.006;
            row.fin_surface_area_m2 = required_area * 2.5;
            row.estimated_pressure_drop_bar = 0.01 * required_area;
            row.estimated_weight_kg = plan * 45.0;
            row.module_length_m = length;
            row.module_width_m = width;
            row.module_height_m = 3.2;
        }
        MechanicalType::DoublePipe => {
            let perimeter = PI * 0.05;
            let length = (required_area / perimeter).max(2.5);
            let outer = 0.11;
            let metal = PI * outer * length * 0.005;
            row.estimated_length_m = length;
            row.inner_diameter_m = 0.05;
            row.outer_diameter_m = outer;
            row.wall_thickness_m = 0.005;
            row.tube_count = 2.0;
            row.tube_passes = 2.0;
            row.estimated_pressure_drop_bar = 0.08 * length;
            row.estimated_weight_kg = metal * 7850.0;
            row.module_length_m = length + 0.5;
            row.module_width_m = outer + 0.6;
            row.module_height_m = outer + 0.6;
        }
        MechanicalType::PlateFin => {
            let volume = required_area / 1000.0;
            let wh = (volume / 6.0).powf(1.0 / 3.0);
            let length = 6.0 * wh;
            row.estimated_length_m = length;
            row.inner_diameter_m = wh;
            row.outer_diameter_m = wh + 2.0 * 0.0016;
            row.wall_thickness_m = 0.0016;
            row.fin_surface_area_m2 = required_area;
            row.estimated_pressure_drop_bar = 0.03 * length;
            row.estimated_weight_kg = volume * 0.30 * 2730.0;
            row.module_length_m = length + 2.0;
            row.module_width_m = wh + 0.8;
            row.module_height_m = wh + 0.5;
        }
    }
    row
}

/// `finalizeThermalEstimates`, shared by both branches.
fn finalize(
    mut approach: f64,
    mut lmtd: f64,
    mut calculated_ua: f64,
    duty: f64,
    used_u: f64,
) -> (f64, f64, f64) {
    if approach.is_nan() || approach < 0.0 {
        approach = 0.0;
    }
    if lmtd.is_nan() || lmtd <= 0.0 {
        if approach > 0.0 {
            lmtd = approach;
        } else if duty > 0.0 && used_u > 0.0 {
            lmtd = duty / used_u;
        } else {
            lmtd = 1.0;
        }
    }
    if (calculated_ua.is_nan() || calculated_ua <= 0.0) && duty > 0.0 && lmtd > 0.0 {
        calculated_ua = duty / lmtd;
    }
    if calculated_ua.is_nan() || calculated_ua <= 0.0 {
        calculated_ua = used_u;
    }
    (lmtd, approach, calculated_ua)
}

/// The five candidates and the selection, with the branch's own sizing coefficient.
fn size_candidates(
    calculated_ua: f64,
    lmtd: f64,
    approach: f64,
    duty: f64,
    used_u: f64,
    use_type_specific_coefficient: bool,
    request: &MechanicalRequest,
) -> MechanicalDesign {
    let kinds = [
        MechanicalType::ShellAndTube,
        MechanicalType::PlateAndFrame,
        MechanicalType::AirCooler,
        MechanicalType::DoublePipe,
        MechanicalType::PlateFin,
    ];
    let mut types = Vec::with_capacity(kinds.len());
    for kind in kinds {
        let type_approach = if approach > 0.0 {
            approach.max(kind.allowable_approach())
        } else {
            kind.allowable_approach()
        };
        let mut ua_for_type = calculated_ua;
        if (ua_for_type.is_nan() || ua_for_type <= 0.0) && duty > 0.0 {
            ua_for_type = duty / type_approach.max(1.0);
        }
        if ua_for_type.is_nan() || ua_for_type <= 0.0 {
            ua_for_type = kind.typical_u();
        }
        let sizing_coefficient = if use_type_specific_coefficient {
            kind.typical_u()
        } else {
            used_u.max(kind.typical_u())
        };
        let area = (ua_for_type / sizing_coefficient.max(1.0e-6)).max(0.1);
        types.push(size_type(kind, area, ua_for_type, type_approach));
    }

    let selected = match request.manual_selection {
        Some(kind) => types
            .iter()
            .find(|row| row.kind == kind)
            .copied()
            .unwrap_or(types[0]),
        None => {
            let criterion = request
                .selection_criterion
                .unwrap_or(SelectionCriterion::MinArea);
            let metric = |row: &TypeSizing| match criterion {
                SelectionCriterion::MinArea => row.required_area_m2,
                SelectionCriterion::MinWeight => row.estimated_weight_kg,
                SelectionCriterion::MinPressureDrop => row.estimated_pressure_drop_bar,
            };
            types
                .iter()
                .copied()
                .min_by(|a, b| metric(a).total_cmp(&metric(b)))
                .unwrap_or(types[0])
        }
    };

    MechanicalDesign {
        log_mean_temperature_difference: lmtd,
        approach_temperature: approach,
        used_overall_heat_transfer_coefficient: used_u,
        calculated_ua,
        types,
        selected,
    }
}

/// `HeatExchangerMechanicalDesign.calcDesign`'s heater branch.
#[must_use]
pub fn design(
    duty_w: f64,
    process_in_k: f64,
    process_out_k: f64,
    utility: &UtilitySpecification,
    request: &MechanicalRequest,
) -> MechanicalDesign {
    let duty = duty_w.abs();
    let mut used_u = 500.0;
    if utility.has_overall_heat_transfer_coefficient() {
        used_u = utility.overall_heat_transfer_coefficient.unwrap_or(500.0);
    }

    let mut lmtd = f64::NAN;
    let mut approach = f64::NAN;
    if let Some(temperatures) =
        single_stream_temperatures(process_in_k, process_out_k, duty, utility)
    {
        let delta_t1 = (temperatures[0] - temperatures[3]).abs();
        let delta_t2 = (temperatures[1] - temperatures[2]).abs();
        lmtd = calculate_lmtd(delta_t1, delta_t2);
        approach = delta_t1.min(delta_t2);
    } else {
        let delta = (process_out_k - process_in_k).abs();
        if delta > 0.0 {
            lmtd = delta;
            approach = delta;
        }
    }
    let mut calculated_ua = f64::NAN;
    if duty > 0.0 && lmtd > 0.0 {
        calculated_ua = duty / lmtd;
    }

    let (lmtd, approach, calculated_ua) = finalize(approach, lmtd, calculated_ua, duty, used_u);
    size_candidates(calculated_ua, lmtd, approach, duty, used_u, false, request)
}

/// `HeatExchangerMechanicalDesign.calcDesign`'s two-stream branch.
///
/// **The branch a `HeatExchanger` takes**, and two things differ from the heater's: the
/// sizing coefficient is the configuration's own typical U rather than the larger of it and a
/// used U, and the used coefficient is the class's 500 whatever the exchanger states.
///
/// The temperatures are the exchanger's own four stream temperatures, ordered by their mean;
/// the UA is `getUAvalue()`, zero meaning absent. `duty_w` is the class's `getDuty()` basis,
/// which the caller resolves: the rating and the delta-T closure write the field, a pinned
/// outlet does not. `effectiveness` is the class's `thermalEffectiveness`, zero where the run
/// did not rate.
#[must_use]
#[allow(clippy::too_many_arguments)] // one parameter per temperature and state it reads
pub fn two_stream_design(
    hot_in_k: f64,
    hot_out_k: f64,
    cold_in_k: f64,
    cold_out_k: f64,
    duty_w: f64,
    ua_w_per_k: f64,
    effectiveness: f64,
    request: &MechanicalRequest,
) -> MechanicalDesign {
    let duty = duty_w.abs();
    let used_u = 500.0;

    // `determineTwoStreamTemperatures`: the two sides ordered by mean, not by which the
    // construction calls hot.
    let hot_mean = (hot_in_k + hot_out_k) / 2.0;
    let cold_mean = (cold_in_k + cold_out_k) / 2.0;
    let mut lmtd = f64::NAN;
    let mut approach = f64::NAN;
    if !hot_mean.is_nan() && !cold_mean.is_nan() {
        let (t0, t1, t2, t3) = if hot_mean >= cold_mean {
            (hot_in_k, hot_out_k, cold_in_k, cold_out_k)
        } else {
            (cold_in_k, cold_out_k, hot_in_k, hot_out_k)
        };
        let delta_t1 = (t0 - t3).abs();
        let delta_t2 = (t1 - t2).abs();
        lmtd = calculate_lmtd(delta_t1, delta_t2);
        approach = delta_t1.min(delta_t2);
    }

    let mut calculated_ua = if ua_w_per_k > 0.0 {
        ua_w_per_k
    } else {
        f64::NAN
    };
    // **Carried although it cannot win through the declared inputs**: the rating that wrote
    // the effectiveness already produced the UA this compares against, so the derived value
    // is never larger; a pinned or delta-T run leaves the field at zero.
    if duty > 0.0 && effectiveness > 0.0 && effectiveness < 0.999 {
        let delta_max = (hot_in_k - cold_in_k).abs();
        if delta_max > 0.0 {
            let c_min = duty / (effectiveness * delta_max);
            let ntu = -(1.0 - effectiveness).ln();
            let ua_from_effectiveness = c_min * ntu;
            if ua_from_effectiveness.is_finite()
                && ua_from_effectiveness > 0.0
                && (calculated_ua.is_nan()
                    || calculated_ua <= 0.0
                    || ua_from_effectiveness > calculated_ua)
            {
                calculated_ua = ua_from_effectiveness;
            }
        }
    }

    let (lmtd, approach, calculated_ua) = finalize(approach, lmtd, calculated_ua, duty, used_u);
    size_candidates(calculated_ua, lmtd, approach, duty, used_u, true, request)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(actual: f64, expected: f64, tolerance: f64, what: &str) {
        let error = (actual - expected).abs() / expected.abs().max(f64::MIN_POSITIVE);
        assert!(
            error <= tolerance,
            "{what}: {actual} against {expected}, relative error {error:.3e} > {tolerance:.3e}"
        );
    }

    fn no_utility() -> UtilitySpecification {
        UtilitySpecification::default()
    }

    #[test]
    fn the_no_utility_estimate_is_the_capture() {
        let design = design(
            100_000.0,
            320.0,
            340.863_381_404_017_6,
            &no_utility(),
            &MechanicalRequest::default(),
        );
        close(
            design.log_mean_temperature_difference,
            20.863_381_404_017_616,
            1.0e-9,
            "lmtd",
        );
        close(design.calculated_ua, 4_793.086_894_221_247, 1.0e-9, "ua");
        close(
            design.used_overall_heat_transfer_coefficient,
            500.0,
            1.0e-12,
            "u",
        );
        close(
            design.types[0].required_area_m2,
            8.263_942_921_071_116,
            1.0e-9,
            "shell area",
        );
        close(
            design.types[1].required_area_m2,
            4.167_901_647_148_911,
            1.0e-9,
            "plate area",
        );
        close(
            design.types[2].required_area_m2,
            9.586_173_788_442_494,
            1.0e-9,
            "air area",
        );
        assert_eq!(design.selected.kind, MechanicalType::PlateAndFrame);
        close(
            design.selected.required_area_m2,
            4.167_901_647_148_911,
            1.0e-9,
            "selected",
        );
    }

    #[test]
    fn one_stated_supply_temperature_collapses_the_return() {
        let utility = UtilitySpecification {
            supply_temperature: Some(450.0),
            ..UtilitySpecification::default()
        };
        let design = design(
            100_000.0,
            320.0,
            340.863_381_404_017_6,
            &utility,
            &MechanicalRequest::default(),
        );
        close(
            design.log_mean_temperature_difference,
            119.264_321_460_286_68,
            1.0e-9,
            "lmtd",
        );
        close(
            design.approach_temperature,
            109.136_618_595_982_38,
            1.0e-9,
            "approach",
        );
        close(design.calculated_ua, 838.473_725_858_027_8, 1.0e-9, "ua");
        close(
            design.selected.required_area_m2,
            0.729_107_587_702_632_8,
            1.0e-9,
            "area",
        );
    }

    #[test]
    fn a_heat_capacity_rate_moves_the_return() {
        let utility = UtilitySpecification {
            supply_temperature: Some(450.0),
            heat_capacity_rate: Some(5000.0),
            ..UtilitySpecification::default()
        };
        let design = design(
            100_000.0,
            320.0,
            340.863_381_404_017_6,
            &utility,
            &MechanicalRequest::default(),
        );
        close(
            design.log_mean_temperature_difference,
            109.567_742_355_256_06,
            1.0e-9,
            "lmtd",
        );
        close(design.calculated_ua, 912.677_379_556_676, 1.0e-9, "ua");
    }

    #[test]
    fn a_stated_u_is_the_used_one_but_not_the_type_s_coefficient() {
        let utility = UtilitySpecification {
            supply_temperature: Some(450.0),
            overall_heat_transfer_coefficient: Some(300.0),
            ..UtilitySpecification::default()
        };
        let design = design(
            100_000.0,
            320.0,
            340.863_381_404_017_6,
            &utility,
            &MechanicalRequest::default(),
        );
        close(
            design.used_overall_heat_transfer_coefficient,
            300.0,
            1.0e-12,
            "used u",
        );
        close(
            design.types[1].required_area_m2,
            0.729_107_587_702_632_8,
            1.0e-9,
            "the plate area still uses its 1150",
        );
    }

    #[test]
    fn an_approach_temperature_supplies_the_utility_temperature() {
        let utility = UtilitySpecification {
            approach_temperature: Some(20.0),
            heat_capacity_rate: Some(6000.0),
            ..UtilitySpecification::default()
        };
        let design = design(
            100_000.0,
            320.0,
            340.863_381_404_017_6,
            &utility,
            &MechanicalRequest::default(),
        );
        close(design.approach_temperature, 20.0, 1.0e-12, "approach");
        close(
            design.log_mean_temperature_difference,
            22.031_780_162_764_73,
            1.0e-9,
            "lmtd",
        );
        close(design.calculated_ua, 4_538.897_866_534_779_5, 1.0e-9, "ua");
    }

    #[test]
    fn the_cooling_direction_reorders_the_pairs() {
        let utility = UtilitySpecification {
            supply_temperature: Some(250.0),
            ..UtilitySpecification::default()
        };
        let design = design(
            -100_000.0,
            320.0,
            298.904_911_855_824_54,
            &utility,
            &MechanicalRequest::default(),
        );
        close(
            design.approach_temperature,
            48.904_911_855_824_54,
            1.0e-9,
            "approach",
        );
        close(
            design.log_mean_temperature_difference,
            58.823_380_930_957_015,
            1.0e-9,
            "lmtd",
        );
        close(design.calculated_ua, 1_700.004_291_072_858_2, 1.0e-9, "ua");
    }

    #[test]
    fn a_manual_selection_wins_and_the_criterion_ranks() {
        let manual = MechanicalRequest {
            manual_selection: Some(MechanicalType::DoublePipe),
            selection_criterion: None,
        };
        let chosen = design(
            100_000.0,
            320.0,
            340.863_381_404_017_6,
            &no_utility(),
            &manual,
        );
        assert_eq!(chosen.selected.kind, MechanicalType::DoublePipe);

        let by_weight = MechanicalRequest {
            manual_selection: None,
            selection_criterion: Some(SelectionCriterion::MinWeight),
        };
        let ranked = design(
            100_000.0,
            320.0,
            340.863_381_404_017_6,
            &no_utility(),
            &by_weight,
        );
        let lightest = ranked
            .types
            .iter()
            .min_by(|a, b| a.estimated_weight_kg.total_cmp(&b.estimated_weight_kg))
            .expect("a row");
        assert_eq!(ranked.selected.kind, lightest.kind);
    }

    /// The capture's `rating_ua_100`: the rating state, read by the two-stream branch.
    fn rating_state() -> MechanicalDesign {
        two_stream_design(
            400.0,
            313.757_079_900_407_86,
            300.0,
            330.597_085_014_208_4,
            4_861.386_048_470_768,
            100.0,
            0.761_089_933_102_921_7,
            &MechanicalRequest::default(),
        )
    }

    #[test]
    fn the_two_stream_rating_is_the_capture() {
        let design = rating_state();
        close(
            design.log_mean_temperature_difference,
            34.383_764_949_499_78,
            1.0e-9,
            "lmtd",
        );
        close(
            design.approach_temperature,
            13.757_079_900_407_86,
            1.0e-9,
            "approach",
        );
        close(
            design.used_overall_heat_transfer_coefficient,
            500.0,
            1.0e-12,
            "the class's 500, not the exchanger's 100",
        );
        close(design.calculated_ua, 100.0, 1.0e-12, "ua");
        assert_eq!(design.selected.kind, MechanicalType::PlateAndFrame);
        close(
            design.selected.required_area_m2,
            0.1,
            1.0e-12,
            "the 0.1 floor",
        );
    }

    #[test]
    fn the_two_stream_area_divides_by_the_type_s_own_coefficient() {
        let design = rating_state();
        let shell = design
            .types
            .iter()
            .find(|row| row.kind == MechanicalType::ShellAndTube)
            .expect("the shell row");
        close(shell.required_area_m2, 100.0 / 580.0, 1.0e-12, "shell area");
        let air = design
            .types
            .iter()
            .find(|row| row.kind == MechanicalType::AirCooler)
            .expect("the air row");
        close(air.required_area_m2, 100.0 / 90.0, 1.0e-12, "air area");
    }

    #[test]
    fn the_effectiveness_arm_is_carried_and_cannot_beat_a_stated_ua() {
        close(
            rating_state().calculated_ua,
            100.0,
            1.0e-12,
            "the stated UA stands",
        );
        // With no usable UA the derived one answers: `cmin * -ln(1 - effectiveness)`.
        let derived = two_stream_design(
            400.0,
            313.757_079_900_407_86,
            300.0,
            330.597_085_014_208_4,
            4_861.386_048_470_768,
            0.0,
            0.761_089_933_102_921_7,
            &MechanicalRequest::default(),
        );
        close(
            derived.calculated_ua,
            91.446_371_334_873_32,
            1.0e-9,
            "derived ua",
        );
    }

    #[test]
    fn an_effectiveness_above_0999_skips_the_derived_arm() {
        // The capture's `rating_ua_1000`: the arm's own `effectiveness < 0.999` guard refuses.
        let design = two_stream_design(
            400.0,
            300.000_085_440_696_64,
            300.0,
            339.654_504_041_168_7,
            6_387.386_513_645_394,
            1000.0,
            0.999_997_845_450_292_1,
            &MechanicalRequest::default(),
        );
        close(
            design.calculated_ua,
            1000.0,
            1.0e-12,
            "the stated UA stands",
        );
    }

    #[test]
    fn the_sides_order_by_mean_not_by_name() {
        let swapped = two_stream_design(
            300.0,
            330.597_085_014_208_4,
            400.0,
            313.757_079_900_407_86,
            4_861.386_048_470_768,
            100.0,
            0.761_089_933_102_921_7,
            &MechanicalRequest::default(),
        );
        close(
            swapped.log_mean_temperature_difference,
            34.383_764_949_499_78,
            1.0e-9,
            "lmtd",
        );
        close(
            swapped.approach_temperature,
            13.757_079_900_407_86,
            1.0e-9,
            "approach",
        );
    }

    #[test]
    fn a_pinned_state_reads_the_classes_500_and_a_zero_duty() {
        // `runSpecifiedStream` never writes the `duty` field, so `determineDuty` sees zero.
        let design = two_stream_design(
            400.0,
            350.0,
            300.0,
            318.438_003_460_016_26,
            0.0,
            500.0,
            0.0,
            &MechanicalRequest::default(),
        );
        close(
            design.log_mean_temperature_difference,
            64.499_058_916_849_1,
            1.0e-9,
            "lmtd",
        );
        close(design.approach_temperature, 50.0, 1.0e-12, "approach");
        close(design.calculated_ua, 500.0, 1.0e-12, "ua");
        close(
            design.selected.required_area_m2,
            0.434_782_608_695_652_16,
            1.0e-9,
            "area",
        );
    }

    #[test]
    fn the_delta_t_state_reads_its_computed_conductance() {
        let design = two_stream_design(
            400.0,
            299.0,
            300.0,
            401.735_496_362_576_6,
            15_901.284_453_854_947,
            157.438_459_939_157_9,
            0.0,
            &MechanicalRequest::default(),
        );
        close(
            design.log_mean_temperature_difference,
            1.334_128_581_919_146_2,
            1.0e-9,
            "lmtd",
        );
        close(design.approach_temperature, 1.0, 1.0e-12, "approach");
        close(design.calculated_ua, 157.438_459_939_157_9, 1.0e-9, "ua");
        close(
            design.selected.required_area_m2,
            0.136_903_008_642_746,
            1.0e-9,
            "area",
        );
    }
}
