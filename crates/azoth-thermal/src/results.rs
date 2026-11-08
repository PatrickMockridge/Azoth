//! Result types for the thermal calculations.
//!
//! The same contract as the hydraulics crate's results: one struct per
//! calculation, field names identical to the Python result dataclass and listed in
//! [`CalcResult::FIELDS`], with a test asserting the three agree. See
//! `crates/azoth-core/src/result.rs` for why the duplication is deliberate.

use azoth_core::units::{Area, Frequency, HeatTransfer, Length, Pressure, Velocity};
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
