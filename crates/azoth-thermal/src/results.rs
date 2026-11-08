//! Result types for the thermal calculations.
//!
//! The same contract as the hydraulics crate's results: one struct per
//! calculation, field names identical to the Python result dataclass and listed in
//! [`CalcResult::FIELDS`], with a test asserting the three agree. See
//! `crates/azoth-core/src/result.rs` for why the duplication is deliberate.

use azoth_core::units::{
    Area, DynamicViscosity, Frequency, HeatTransfer, Length, MassDensity, Pressure,
    SpecificHeatCapacity, TemperatureInterval, ThermodynamicTemperature, Time, Velocity,
};
use azoth_core::{CalcResult, Warning};
use uom::si::f64::Power;

/// Result of `thermal.conduction_plane_wall`.
#[derive(Debug, Clone, PartialEq)]
pub struct ConductionPlaneWallResult {
    /// Heat flow rate through the wall. Signed: it follows the sign of the
    /// temperature difference, and there is no warning for a negative value
    /// because a negative value is a direction rather than a problem.
    pub q: Power,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ConductionPlaneWallResult {
    const CALC_ID: &'static str = "thermal.conduction_plane_wall";
    const FIELDS: &'static [&'static str] = &["q", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.shell_and_tube_rating`.
#[derive(Debug, Clone, PartialEq)]
pub struct ShellAndTubeRatingResult {
    /// The tube-side film coefficient.
    pub tube_side_heat_transfer_coefficient: HeatTransfer,
    /// The shell-side film coefficient.
    pub shell_side_heat_transfer_coefficient: HeatTransfer,
    /// The overall coefficient on the outside area.
    pub overall_heat_transfer_coefficient: HeatTransfer,
    /// The tube-side pressure drop over all passes.
    pub tube_side_pressure_drop: Pressure,
    /// The shell-side pressure drop.
    pub shell_side_pressure_drop: Pressure,
    /// The tube-side velocity in one pass.
    pub tube_side_velocity: Velocity,
    /// The shell-side crossflow velocity.
    pub shell_side_velocity: Velocity,
    /// The tube-side Reynolds number.
    pub tube_side_reynolds: f64,
    /// The shell-side Reynolds number.
    pub shell_side_reynolds: f64,
    /// The effective tube flow diameter.
    pub effective_tube_inside_diameter: Length,
    /// The nominal outside heat-transfer area.
    pub outside_heat_transfer_area: Area,
    /// The deposit's cylindrical resistance on the outside area, m**2*K/W.
    pub tube_fouling_layer_resistance_outside: f64,
    /// The baffle-cut correction factor.
    pub jc: f64,
    /// The baffle-leakage correction factor.
    pub jl: f64,
    /// The bundle-bypass correction factor.
    pub jb: f64,
    /// The unequal-spacing correction factor.
    pub js: f64,
    /// The adverse-temperature-gradient correction factor.
    pub jr: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ShellAndTubeRatingResult {
    const CALC_ID: &'static str = "thermal.shell_and_tube_rating";
    const FIELDS: &'static [&'static str] = &[
        "tube_side_heat_transfer_coefficient",
        "shell_side_heat_transfer_coefficient",
        "overall_heat_transfer_coefficient",
        "tube_side_pressure_drop",
        "shell_side_pressure_drop",
        "tube_side_velocity",
        "shell_side_velocity",
        "tube_side_reynolds",
        "shell_side_reynolds",
        "effective_tube_inside_diameter",
        "outside_heat_transfer_area",
        "tube_fouling_layer_resistance_outside",
        "jc",
        "jl",
        "jb",
        "js",
        "jr",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.boiling_heat_transfer`.
#[derive(Debug, Clone, PartialEq)]
pub struct BoilingHeatTransferResult {
    /// The Chen local coefficient.
    pub chen_htc: HeatTransfer,
    /// The Gungor-Winterton local coefficient.
    pub gungor_winterton_htc: HeatTransfer,
    /// The same coefficient with the horizontal correction.
    pub gungor_winterton_corrected_htc: HeatTransfer,
    /// The Simpson's-rule average over the range.
    pub average_htc: HeatTransfer,
    /// The Lockhart-Martinelli parameter.
    pub martinelli_parameter: f64,
    /// The Chen enhancement factor.
    pub chen_enhancement_factor: f64,
    /// The Chen suppression factor.
    pub chen_suppression_factor: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for BoilingHeatTransferResult {
    const CALC_ID: &'static str = "thermal.boiling_heat_transfer";
    const FIELDS: &'static [&'static str] = &[
        "chen_htc",
        "gungor_winterton_htc",
        "gungor_winterton_corrected_htc",
        "average_htc",
        "martinelli_parameter",
        "chen_enhancement_factor",
        "chen_suppression_factor",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.shah_condensation`.
#[derive(Debug, Clone, PartialEq)]
pub struct ShahCondensationResult {
    /// The all-liquid Dittus-Boelter coefficient.
    pub liquid_only_htc: HeatTransfer,
    /// `G*D/mu_l`.
    pub liquid_only_reynolds: f64,
    /// The local Shah coefficient at the stated quality.
    pub local_htc: HeatTransfer,
    /// The Simpson's-rule average over the range.
    pub average_htc: HeatTransfer,
    /// The vertical-tube coefficient.
    pub vertical_htc: HeatTransfer,
    /// Whether the stated values are inside the validated bands.
    pub in_valid_range: bool,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ShahCondensationResult {
    const CALC_ID: &'static str = "thermal.shah_condensation";
    const FIELDS: &'static [&'static str] = &[
        "liquid_only_htc",
        "liquid_only_reynolds",
        "local_htc",
        "average_htc",
        "vertical_htc",
        "in_valid_range",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.lmtd_correction_factor`.
#[derive(Debug, Clone, PartialEq)]
pub struct LmtdCorrectionFactorResult {
    /// The capacity ratio.
    pub r: f64,
    /// The thermal effectiveness.
    pub p: f64,
    /// The correction factor for the stated shell passes.
    pub ft: f64,
    /// The correction factor for one shell pass.
    pub ft_1_shell_pass: f64,
    /// The correction factor for two shell passes.
    pub ft_2_shell_passes: f64,
    /// The first pass count reaching 0.75, or -1.
    pub required_shell_passes: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for LmtdCorrectionFactorResult {
    const CALC_ID: &'static str = "thermal.lmtd_correction_factor";
    const FIELDS: &'static [&'static str] = &[
        "r",
        "p",
        "ft",
        "ft_1_shell_pass",
        "ft_2_shell_passes",
        "required_shell_passes",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.shell_and_tube_vibration`.
#[derive(Debug, Clone, PartialEq)]
pub struct ShellAndTubeVibrationResult {
    /// The first natural frequency of the pinned span, Hz.
    pub natural_frequency: Frequency,
    /// The Strouhal shedding frequency, Hz.
    pub vortex_shedding_frequency: Frequency,
    /// `f_vs / f_n`.
    pub vortex_shedding_ratio: f64,
    /// Whether the shedding ratio is inside 0.8-1.2.
    pub vortex_shedding_critical: bool,
    /// Connors' critical velocity.
    pub critical_velocity: Velocity,
    /// `V / (0.5 * V_crit)`.
    pub velocity_ratio: f64,
    /// Whether the velocity ratio exceeds one.
    pub fluid_elastic_critical: bool,
    /// The first acoustic resonance frequency, Hz.
    pub acoustic_frequency: Frequency,
    /// `f_vs / f_ac`.
    pub acoustic_ratio: f64,
    /// Whether the acoustic ratio is inside 0.8-1.2.
    pub acoustic_critical: bool,
    /// True only when none of the three mechanisms is critical.
    pub passed: bool,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ShellAndTubeVibrationResult {
    const CALC_ID: &'static str = "thermal.shell_and_tube_vibration";
    const FIELDS: &'static [&'static str] = &[
        "natural_frequency",
        "vortex_shedding_frequency",
        "vortex_shedding_ratio",
        "vortex_shedding_critical",
        "critical_velocity",
        "velocity_ratio",
        "fluid_elastic_critical",
        "acoustic_frequency",
        "acoustic_ratio",
        "acoustic_critical",
        "passed",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.two_phase_pressure_drop`.
#[derive(Debug, Clone, PartialEq)]
pub struct TwoPhasePressureDropResult {
    /// The Friedel frictional gradient at the stated quality, Pa/m.
    pub friedel_gradient: Pressure,
    /// The Friedel drop over the tube length at the stated quality.
    pub friedel_pressure_drop: Pressure,
    /// The Simpson's-rule average over the quality range, times the length.
    pub friedel_average_pressure_drop: Pressure,
    /// The Muller-Steinhagen-Heck frictional gradient, Pa/m.
    pub muller_steinhagen_heck_gradient: Pressure,
    /// The gravitational gradient of the homogeneous mixture, Pa/m.
    pub gravitational_gradient: Pressure,
    /// The momentum-flux change over the quality range.
    pub acceleration_pressure_drop: Pressure,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TwoPhasePressureDropResult {
    const CALC_ID: &'static str = "thermal.two_phase_pressure_drop";
    const FIELDS: &'static [&'static str] = &[
        "friedel_gradient",
        "friedel_pressure_drop",
        "friedel_average_pressure_drop",
        "muller_steinhagen_heck_gradient",
        "gravitational_gradient",
        "acceleration_pressure_drop",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.fouling_model`.
#[derive(Debug, Clone, PartialEq)]
pub struct FoulingModelResult {
    /// The resistance at the stated time, m**2*K/W.
    pub fouling_resistance_at_time: f64,
    /// The resistance after the Euler steps, m**2*K/W.
    pub current_fouling_resistance: f64,
    /// The model's asymptote, m**2*K/W.
    pub asymptotic_fouling_resistance: f64,
    /// The threshold wall temperature; `MAX_VALUE` when no wall exceeds the threshold.
    pub threshold_wall_temperature: ThermodynamicTemperature,
    /// The velocity at which the asymptote crosses the fixed resistance.
    pub threshold_velocity: Velocity,
    /// The predicted time to the target resistance; `-1` past the asymptote.
    pub time_to_target: Time,
    /// The time advanced by the Euler steps.
    pub operating_time: Time,
    /// The Ebert-Panchal rate at the current resistance, m**2*K/(W*h).
    pub fouling_rate: f64,
    /// The Reynolds number of the stated conditions.
    pub reynolds_number: f64,
    /// The wall shear stress of the stated conditions.
    pub wall_shear_stress: Pressure,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for FoulingModelResult {
    const CALC_ID: &'static str = "thermal.fouling_model";
    const FIELDS: &'static [&'static str] = &[
        "fouling_resistance_at_time",
        "current_fouling_resistance",
        "asymptotic_fouling_resistance",
        "threshold_wall_temperature",
        "threshold_velocity",
        "time_to_target",
        "operating_time",
        "fouling_rate",
        "reynolds_number",
        "wall_shear_stress",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.tube_insert_model`.
#[derive(Debug, Clone, PartialEq)]
pub struct TubeInsertModelResult {
    /// `h_enhanced/h_plain`.
    pub heat_transfer_enhancement_ratio: f64,
    /// `f_enhanced/f_plain`.
    pub pressure_drop_penalty_ratio: f64,
    /// The performance evaluation criterion at constant pumping power.
    pub performance_evaluation_criteria: f64,
    /// The enhanced coefficient.
    pub enhanced_heat_transfer_coefficient: HeatTransfer,
    /// The enhanced pressure drop.
    pub enhanced_pressure_drop: Pressure,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TubeInsertModelResult {
    const CALC_ID: &'static str = "thermal.tube_insert_model";
    const FIELDS: &'static [&'static str] = &[
        "heat_transfer_enhancement_ratio",
        "pressure_drop_penalty_ratio",
        "performance_evaluation_criteria",
        "enhanced_heat_transfer_coefficient",
        "enhanced_pressure_drop",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Result of `thermal.tema_standard`.
#[derive(Debug, Clone, PartialEq)]
pub struct TemaStandardResult {
    /// The configuration's relative cost factor.
    pub configuration_cost_factor: f64,
    /// The shell letter's pressure drop factor.
    pub shell_pressure_drop_factor: f64,
    /// Whether the bundle is removable.
    pub bundle_removable: bool,
    /// Whether the configuration handles thermal expansion.
    pub good_thermal_expansion: bool,
    /// The construction class's cost factor.
    pub class_cost_factor: f64,
    /// The class's minimum tube wall.
    pub min_tube_wall_thickness: Length,
    /// The class's minimum corrosion allowance.
    pub min_corrosion_allowance: Length,
    /// The tube size's outside diameter.
    pub tube_outer_diameter: Length,
    /// The same diameter, transported in inches.
    pub tube_outer_diameter_inch: Length,
    /// The layout angle, in degrees.
    pub layout_angle: f64,
    /// The pattern's minimum pitch ratio.
    pub min_pitch_ratio: f64,
    /// The pattern's heat-transfer factor.
    pub pitch_heat_transfer_factor: f64,
    /// The baffle type's heat-transfer factor.
    pub baffle_heat_transfer_factor: f64,
    /// The baffle type's pressure-drop factor.
    pub baffle_pressure_drop_factor: f64,
    /// The minimum tube pitch.
    pub min_tube_pitch: Length,
    /// The truncated tube-count estimate.
    pub estimated_tube_count: f64,
    /// The maximum unsupported tube span.
    pub max_unsupported_span: Length,
    /// The minimum baffle spacing.
    pub min_baffle_spacing: Length,
    /// The maximum baffle spacing.
    pub max_baffle_spacing: Length,
    /// The recommended designation.
    pub recommended_designation: String,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for TemaStandardResult {
    const CALC_ID: &'static str = "thermal.tema_standard";
    const FIELDS: &'static [&'static str] = &[
        "configuration_cost_factor",
        "shell_pressure_drop_factor",
        "bundle_removable",
        "good_thermal_expansion",
        "class_cost_factor",
        "min_tube_wall_thickness",
        "min_corrosion_allowance",
        "tube_outer_diameter",
        "tube_outer_diameter_inch",
        "layout_angle",
        "min_pitch_ratio",
        "pitch_heat_transfer_factor",
        "baffle_heat_transfer_factor",
        "baffle_pressure_drop_factor",
        "min_tube_pitch",
        "estimated_tube_count",
        "max_unsupported_span",
        "min_baffle_spacing",
        "max_baffle_spacing",
        "recommended_designation",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// One zone's stated transport properties, for `thermal.zone_rating`.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneDefinitionRow {
    /// The zone's name.
    pub zone_name: String,
    /// The zone's fraction of the total duty.
    pub duty_fraction: f64,
    /// The exchanger's total duty.
    pub total_duty: Power,
    /// The zone's stated log-mean temperature difference.
    pub lmtd: TemperatureInterval,
    /// The zone's tube-side density.
    pub tube_density: MassDensity,
    /// The zone's tube-side viscosity.
    pub tube_viscosity: DynamicViscosity,
    /// The zone's tube-side specific heat capacity.
    pub tube_cp: SpecificHeatCapacity,
    /// The zone's tube-side thermal conductivity.
    pub tube_conductivity: azoth_core::units::ThermalConductivity,
    /// The zone's shell-side density.
    pub shell_density: MassDensity,
    /// The zone's shell-side viscosity.
    pub shell_viscosity: DynamicViscosity,
    /// The zone's shell-side specific heat capacity.
    pub shell_cp: SpecificHeatCapacity,
    /// The zone's shell-side thermal conductivity.
    pub shell_conductivity: azoth_core::units::ThermalConductivity,
}

/// One zone's rating, for `thermal.zone_rating`.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneRatingRow {
    /// The zone's name, echoed.
    pub zone_name: String,
    /// The zone's duty fraction, echoed.
    pub duty_fraction: f64,
    /// The zone's tube-side film coefficient.
    pub tube_side_heat_transfer_coefficient: HeatTransfer,
    /// The zone's shell-side film coefficient.
    pub shell_side_heat_transfer_coefficient: HeatTransfer,
    /// The zone's overall coefficient on the outside area.
    pub overall_heat_transfer_coefficient: HeatTransfer,
    /// The zone's stated LMTD, echoed.
    pub lmtd: TemperatureInterval,
    /// `duty_fraction * total_duty / (U * lmtd)`, zero when either denominator is not positive.
    pub required_area: Area,
}

/// Result of `thermal.zone_rating`.
#[derive(Debug, Clone, PartialEq)]
pub struct ZoneRatingResult {
    /// The number of zones rated.
    pub zone_count: f64,
    /// One row per zone, in the input order.
    pub zones: Vec<ZoneRatingRow>,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for ZoneRatingResult {
    const CALC_ID: &'static str = "thermal.zone_rating";
    const FIELDS: &'static [&'static str] = &["zone_count", "zones", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}
