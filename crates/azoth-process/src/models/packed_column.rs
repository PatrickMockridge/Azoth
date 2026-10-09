//! `process.packed_column` - the packed column as a registered id.
//!
//! Spec: `specs/models/process/packed_column.toml`. **The arithmetic is
//! [`crate::models::distillation_column`]'s**, because the class's is: `PackedColumn extends
//! `DistillationColumn`, its `run` is `super.run(id)` and the packing is read by a hydraulics
//! report afterwards. What this module adds is the stage count the class's constructor derives
//! from a packed height, the refusal of the one packing parameter the class itself refuses, and
//! the seven quantities `ColumnInternalsDesigner` reports on the far side of the solve.

use azoth_core::units::{
    Length, MolarEnergy, Power, Pressure, ThermodynamicTemperature, joules_per_mole,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::column::capacity::{DEFAULT_MAX_ALLOWABLE_FS_FACTOR, FsLimits, fs_limits};
use crate::kernels::distillation_column as kernel;
use crate::kernels::packed_column::{PackedReport, PackingInputs, packing_report, stage_count};
use crate::models::distillation_column::distillation_column_outcome;

/// Result of `process.packed_column`.
///
/// **The base column's record, plus the seven quantities the packing's report adds.** The packing
/// parameters reach the *separation* only through the stage count, which is why the profile, both
/// products, both duties and the three residuals are the base's unchanged; everything the packing
/// contributes is `ColumnInternalsDesigner`'s report on the far side of the solve, and those are
/// the seven at the end. **The trayed half of that designer is not here**: this class builds it
/// with `internalsType = "packed"`, so `calculateTrayed` never runs for it.
#[derive(Debug, Clone, PartialEq)]
pub struct PackedColumnResult {
    /// Each tray's temperature, K.
    pub tray_temperature: Vec<ThermodynamicTemperature>,
    /// Each tray's pressure, Pa.
    pub tray_pressure: Vec<Pressure>,
    /// Each tray's vapour traffic, mol/s.
    pub tray_gas_n: Vec<f64>,
    /// Each tray's liquid traffic, mol/s.
    pub tray_liquid_n: Vec<f64>,
    /// Distillate molar flow, mol/s.
    pub distillate_n: f64,
    /// Distillate composition.
    pub distillate_z: Vec<f64>,
    /// Distillate pressure.
    pub distillate_p: Pressure,
    /// Distillate temperature.
    pub distillate_t: ThermodynamicTemperature,
    /// Distillate molar enthalpy.
    pub distillate_h: MolarEnergy,
    /// Bottoms molar flow, mol/s.
    pub bottoms_n: f64,
    /// Bottoms composition.
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    pub bottoms_p: Pressure,
    /// Bottoms temperature.
    pub bottoms_t: ThermodynamicTemperature,
    /// Bottoms molar enthalpy.
    pub bottoms_h: MolarEnergy,
    /// **The vapour each tray withdrew**, zero where a tray drew none. The packing does not change
    /// a draw, and the three fractions are this id's own declared inputs - inherited from the base
    /// and overridden in no part - so a stated one draws here as it does there.
    pub gas_side_draw_n: Vec<f64>,
    /// The liquid each tray withdrew as a liquid side draw.
    pub liquid_side_draw_n: Vec<f64>,
    /// The liquid each tray withdrew as a pumparound.
    pub pumparound_n: Vec<f64>,
    /// The condenser's duty, W.
    pub condenser_duty: Power,
    /// The reboiler's duty, W.
    pub reboiler_duty: Power,
    /// Iterations taken.
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    pub temperature_residual: f64,
    /// The products' worst component imbalance against the feed, relative.
    pub mass_residual: f64,
    /// The enthalpy closure.
    pub energy_residual: f64,
    /// `getFsFactor`, at the diameter the report resolved.
    pub fs_factor: f64,
    /// `getFsFactorUtilization`, against this id's own `max_allowable_fs_factor`.
    pub fs_factor_utilization: f64,
    /// `isFsFactorWithinDesignLimit`.
    pub fs_factor_within_design_limit: bool,
    /// `getMinimumDiameterForFsLimit`.
    pub minimum_diameter_for_fs_limit: Length,
    /// The height equivalent to a theoretical plate, m.
    pub hetp: Length,
    /// The packed height over the HETP.
    ///
    /// **It describes a different column from the one solved, and it is published that way.** The
    /// constructor fixes the stage count on a `0.5` m HETP guess and the report computes the real
    /// one afterwards: the capture's 2.3 m case solves five middle trays and reports `1.597`
    /// theoretical stages. Asserted rather than corrected, and filed upstream.
    pub theoretical_stages: f64,
    /// The load, in per cent of flood.
    pub percent_flood: f64,
    /// The velocity at which the bed floods, m/s.
    pub flooding_velocity: f64,
    /// The bed's total pressure drop, Pa.
    pub packing_pressure_drop: Pressure,
    /// Whether the bed is inside the design window.
    ///
    /// **`PackedColumn.isHydraulicsOk`, which the jar shows is the calculator's `isDesignOk`** -
    /// there is one predicate and not two, so this is the hydraulics calculation's own verdict.
    pub hydraulics_ok: bool,
    /// The column's internal diameter, m - the stated one, or the sized one where none was stated.
    pub internal_diameter: Length,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl PackedColumnResult {
    /// The record of one solve: the base column's answer, plus the packing's report.
    ///
    /// **`of` rather than `From`, because the report does not follow from the base record.** It is
    /// read at the middle tray's state on the far side of the solve, so it is a second input and
    /// not a re-projection of the first - the shape `DistillationColumnResult::of` and
    /// `RateBasedPackedColumnResult::of` already have.
    ///
    /// **It takes the *outcome* rather than a base `DistillationColumnResult`, and that is the
    /// class's own split.** `PackedColumn.calcPackingHydraulics` builds one designer with
    /// `internalsType = "packed"`, so `calculateTrayed` **never runs** for this class and there is
    /// no trayed internals report to copy: a packed column that published one would be publishing
    /// a tree NeqSim does not build for it. The base's own quantities are re-read from the outcome
    /// here instead, which is what its 29 fields were.
    #[must_use]
    pub fn of(
        outcome: &kernel::ColumnOutcome,
        limits: &FsLimits,
        report: PackedReport,
        warnings: Vec<Warning>,
    ) -> Self {
        Self {
            tray_temperature: outcome.trays.iter().map(|t| t.temperature).collect(),
            tray_pressure: outcome.trays.iter().map(|t| t.pressure).collect(),
            tray_gas_n: outcome.trays.iter().map(|t| t.gas_n).collect(),
            tray_liquid_n: outcome.trays.iter().map(|t| t.liquid_n).collect(),
            distillate_n: outcome.distillate.n,
            distillate_z: outcome.distillate.z.clone(),
            distillate_p: outcome.distillate.p,
            distillate_t: outcome.distillate.t,
            distillate_h: joules_per_mole(outcome.distillate.h.value),
            bottoms_n: outcome.bottoms.n,
            bottoms_z: outcome.bottoms.z.clone(),
            bottoms_p: outcome.bottoms.p,
            bottoms_t: outcome.bottoms.t,
            bottoms_h: joules_per_mole(outcome.bottoms.h.value),
            gas_side_draw_n: outcome
                .gas_side_draws
                .iter()
                .map(|draw| draw.as_ref().map_or(0.0, |draw| draw.n))
                .collect(),
            liquid_side_draw_n: outcome
                .liquid_side_draws
                .iter()
                .map(|draw| draw.as_ref().map_or(0.0, |draw| draw.n))
                .collect(),
            pumparound_n: outcome
                .pumparounds
                .iter()
                .map(|draw| draw.as_ref().map_or(0.0, |draw| draw.n))
                .collect(),
            condenser_duty: outcome.condenser_duty,
            reboiler_duty: outcome.reboiler_duty,
            iterations: outcome.iterations,
            temperature_residual: outcome.temperature_residual,
            mass_residual: outcome.mass_residual,
            energy_residual: outcome.energy_residual,
            fs_factor: limits.fs_factor,
            fs_factor_utilization: limits.fs_factor_utilization,
            fs_factor_within_design_limit: limits.fs_factor_within_design_limit,
            minimum_diameter_for_fs_limit: limits.minimum_diameter_for_fs_limit,
            hetp: report.hetp,
            theoretical_stages: report.theoretical_stages,
            percent_flood: report.percent_flood,
            flooding_velocity: report.flooding_velocity,
            packing_pressure_drop: report.packing_pressure_drop,
            hydraulics_ok: report.hydraulics_ok,
            internal_diameter: report.internal_diameter,
            warnings,
        }
    }
}

/// Solve a packed column.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a `packing_hydraulic_capacity_factor` that is not positive
/// and finite, which is `setPackingHydraulicCapacityFactor`'s own check and the one parameter of
/// the packing group the class refuses rather than reads; every error
/// [`distillation_column_outcome`] raises, which this id inherits because the solve is the base's;
/// and every error the packing's report raises, which is [`packing_report`]'s - the middle tray
/// failing to rebuild as a two-phase fluid, or the hydraulics refusing the state it is given.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are twenty-six
pub fn packed_column(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    packed_height: f64,
    feed_stage: usize,
    has_reboiler: bool,
    has_condenser: bool,
    top_pressure: Pressure,
    bottom_pressure: Pressure,
    reboiler_temperature: Option<ThermodynamicTemperature>,
    condenser_temperature: Option<ThermodynamicTemperature>,
    temperature_tolerance: f64,
    max_iterations: usize,
    packing_type: Option<&str>,
    structured_packing: Option<bool>,
    design_flood_fraction: Option<f64>,
    packing_hydraulic_capacity_factor: Option<f64>,
    column_diameter: Option<f64>,
    murphree_efficiency: Option<f64>,
    reactive: Option<bool>,
    reactive_start_tray: Option<usize>,
    reactive_end_tray: Option<usize>,
    solver_type: Option<&str>,
    top_specification_type: Option<&str>,
    top_specification_target: Option<f64>,
    top_specification_component: Option<&str>,
    bottom_specification_type: Option<&str>,
    bottom_specification_target: Option<f64>,
    bottom_specification_component: Option<&str>,
    tray_murphree_efficiency: Option<&[f64]>,
    gas_side_draw_fractions: Option<&[f64]>,
    liquid_side_draw_fractions: Option<&[f64]>,
    pumparound_fractions: Option<&[f64]>,
    side_draw_flow_tray: Option<usize>,
    side_draw_flow_phase: Option<&str>,
    side_draw_flow_target: Option<f64>,
    side_draw_flow_tolerance: Option<f64>,
    side_draw_flow_max_iterations: Option<usize>,
    pumparound_return_tray: Option<usize>,
    pumparound_draw_tray: Option<usize>,
    pumparound_draw_fraction: Option<f64>,
    pumparound_temperature_drop: Option<f64>,
    pumparound_tolerance: Option<f64>,
    pumparound_max_iterations: Option<usize>,
    max_allowable_fs_factor: Option<f64>,
) -> Result<PackedColumnResult> {
    // **The class's own check, on the one packing parameter it does not simply read.**
    // `setPackingHydraulicCapacityFactor` throws for a value that is not positive and finite, so
    // the class refuses a factor its own setter would reject rather than carrying it inertly as
    // it carries the other four.
    if let Some(factor) = packing_hydraulic_capacity_factor {
        if !factor.is_finite() || factor <= 0.0 {
            return Err(AzothError::invalid_input(
                "packing_hydraulic_capacity_factor",
                format!(
                    "a packing hydraulic capacity factor of {factor} is not positive and finite, \
                     which `PackedColumn.setPackingHydraulicCapacityFactor` refuses with an \
                     `IllegalArgumentException`"
                ),
            ));
        }
    }
    // **The declared bounds, which nothing was running.** The base column applies its own spec's
    // four, and this id repeats them because it declares them: until this call existed, the packed
    // spec's `top_pressure`, `bottom_pressure`, `temperature_tolerance` and `feed_t` rows were a
    // declaration no code read.
    let spec = &crate::model_gen::PACKED_COLUMN_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "top_pressure" => Some(top_pressure.value),
            "bottom_pressure" => Some(bottom_pressure.value),
            "temperature_tolerance" => Some(temperature_tolerance),
            "feed_t" => Some(feed_t.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // **The other four are declarations the solve is indifferent to, and the report is not.** The
    // class accepts them, every one is read by `ColumnInternalsDesigner` after the column has
    // converged, and this id publishes seven of its quantities - so they are carried to the report
    // instead of being discarded. That is the split `PackedColumn.run` makes: `super.run(id)`
    // first, `calcPackingHydraulics()` afterwards.
    let packing = PackingInputs {
        name: packing_type.map(str::to_string),
        structured: structured_packing,
        design_flood_fraction,
        hydraulic_capacity_factor: packing_hydraulic_capacity_factor,
        column_diameter,
    };

    let stages = stage_count(packed_height);
    let (outcome, warnings) = distillation_column_outcome(
        None,
        components,
        feed_n,
        feed_z,
        feed_p,
        feed_t,
        stages,
        feed_stage,
        has_reboiler,
        has_condenser,
        top_pressure,
        bottom_pressure,
        reboiler_temperature,
        condenser_temperature,
        temperature_tolerance,
        max_iterations,
        murphree_efficiency,
        // **A notional stage is still a stage.** The height derives the tray count, so these are
        // one entry per theoretical stage with the same `NaN` fall-through the column's own
        // vector uses.
        tray_murphree_efficiency,
        solver_type,
        top_specification_type,
        top_specification_target,
        top_specification_component,
        bottom_specification_type,
        bottom_specification_target,
        bottom_specification_component,
        // **The section is the base column's own**, which `PackedColumn` inherits and does not
        // override: its middle trays are the ones the packing's height derives.
        reactive,
        reactive_start_tray,
        reactive_end_tray,
        // **The draws are the base column's and `PackedColumn` inherits all of them**: it
        // overrides `run` and `toJson` and no part of the separation, so a draw on a notional
        // stage is the same split the column makes. The ends are still the ends, so a fraction
        // on one is still refused.
        gas_side_draw_fractions,
        liquid_side_draw_fractions,
        pumparound_fractions,
        side_draw_flow_tray,
        side_draw_flow_phase,
        side_draw_flow_target,
        side_draw_flow_tolerance,
        side_draw_flow_max_iterations,
        pumparound_return_tray,
        pumparound_draw_tray,
        pumparound_draw_fraction,
        pumparound_temperature_drop,
        pumparound_tolerance,
        pumparound_max_iterations,
    )?;

    let report = packing_report(&outcome, components, packed_height, &packing)?;
    // **The diameter the Fs family reads is the report's resolved one, and that is the class's
    // own order**: `calcPackingHydraulics` ends by writing the inherited `internalDiameter` from
    // the stated `columnDiameter` where it is positive and from the designer's sizing otherwise,
    // so a caller that stated none gets its factor at the diameter the bed was sized to - the
    // capture's unstated row reads `0.3` m back.
    let limits = fs_limits(
        &outcome.distillate,
        report.internal_diameter,
        max_allowable_fs_factor.unwrap_or(DEFAULT_MAX_ALLOWABLE_FS_FACTOR),
    )?;
    Ok(PackedColumnResult::of(&outcome, &limits, report, warnings))
}

#[cfg(test)]
mod report {
    use azoth_core::units::{kelvins, pascals};

    use crate::models::packed_column::packed_column;

    /// **The wired report reproduces the capture's own seven, on the row the class solved.**
    ///
    /// The capture's `packed_distillation_binary_2m3` block, end to end: the column is solved from
    /// the row's own inputs and the seven quantities `ColumnInternalsDesigner` reports are held to
    /// what the class printed. `packing_pressure_drop` is the one that discriminates - it reads the
    /// packing factor, the void fraction and the tray's vapour density through Leva's form, so a
    /// port that read the wrong tray or the wrong row's geometry would move it.
    ///
    /// `internal_diameter` goes through the **sizing** branch, because the class's own
    /// `columnDiameter` is `-1.0` on this row: the diameter is not stated, so it is sized and the
    /// sized value is what the report then runs at.
    ///
    /// **`hydraulics_ok` is held here rather than in the case file.** `TestCase` carries a boolean
    /// for an *input* and not for an expectation - `gen_registry`'s `collect_numbers` excludes
    /// flags deliberately, so that a reader cannot mistake one for a quantity - and no spec-driven
    /// runner exists for this model to read one back through `case.flag`. The three cases hold the
    /// other six; this holds the flag, against the same capture.
    #[test]
    #[allow(clippy::too_many_lines)] // one call, one parameter per declared input
    fn the_report_reproduces_the_captures_seven() {
        let out = packed_column(
            &["methane".to_string(), "n-butane".to_string()],
            7.490_704_036_290_964,
            &[0.5, 0.5],
            pascals(2.0e6),
            kelvins(300.0),
            2.3,
            2,
            true,
            true,
            pascals(1.9e6),
            pascals(2.0e6),
            Some(kelvins(373.15)),
            Some(kelvins(253.15)),
            1.0e-6,
            200,
            Some("Pall-Ring-50"),
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            // The Fs limit, unstated, so the class's own `2.5` answers.
            None,
        )
        .expect("the capture's 2.3 m row solves");

        let close = |got: f64, want: f64| (got - want).abs() / want.abs() < 1e-5;
        assert!(
            close(out.hetp.value, 1.44),
            "hetp: {} against the capture's 1.44",
            out.hetp.value
        );
        assert!(
            close(out.theoretical_stages, 1.597_222_222_222_222),
            "theoretical_stages: {}",
            out.theoretical_stages
        );
        assert!(
            close(out.percent_flood, 15.053_642_072_827_294),
            "percent_flood: {}",
            out.percent_flood
        );
        assert!(
            close(out.flooding_velocity, 0.508_073_023_146_439_5),
            "flooding_velocity: {}",
            out.flooding_velocity
        );
        assert!(
            close(out.packing_pressure_drop.value, 0.640_829_527_984_597_7),
            "packing_pressure_drop: {}",
            out.packing_pressure_drop.value
        );
        assert!(
            !out.hydraulics_ok,
            "the class's own row is outside the 40-to-80-per-cent window"
        );
        assert!(
            close(out.internal_diameter.value, 0.3),
            "internal_diameter: {}",
            out.internal_diameter.value
        );
        // **The Fs family, at the diameter the sizing resolved**, from the same capture's
        // `packed_distillation_binary_unstated` row of `process_column_capacity.tsv`. The
        // diameter matters here: this row states none, so the class writes `0.3` m back and the
        // factor is read against that rather than against the caller's absent number.
        assert!(
            close(out.fs_factor, 0.226_628_651_724_175_35),
            "fs_factor: {}",
            out.fs_factor
        );
        assert!(
            close(out.fs_factor_utilization, 0.090_651_460_689_670_14),
            "fs_factor_utilization: {}",
            out.fs_factor_utilization
        );
        assert!(
            close(
                out.minimum_diameter_for_fs_limit.value,
                0.090_325_143_022_695_02
            ),
            "minimum_diameter_for_fs_limit: {}",
            out.minimum_diameter_for_fs_limit.value
        );
        // **The verdict is held here rather than in the case file**, as `hydraulics_ok` is and
        // for the same reason: `TestCase` carries a boolean for an input and not for an
        // expectation, so a flag asserted nowhere would read as verified.
        assert!(
            out.fs_factor_within_design_limit,
            "the row is inside the class's own 2.5 limit"
        );
    }
}

#[cfg(test)]
mod surface_tension {
    use azoth_core::units::{
        kelvins, kilograms_per_cubic_meter, kilograms_per_second, meters, newtons_per_meter,
        pascal_seconds, pascals,
    };
    use azoth_hydraulics::packing_hydraulics::{PackingState, packing_hydraulics};

    use crate::kernels::packed_column::{DESIGNER_SURFACE_TENSION_N_PER_M, stage_count};
    use crate::models::distillation_column::distillation_column_outcome;
    use crate::segment::phase::{Pick, phase_view};
    use crate::stream::Stream;

    /// **Which surface tension the packed column's report takes, settled against the capture.**
    ///
    /// `ColumnInternalsDesigner.getTrayProperties` reads
    /// `fluid.getInterphaseProperties().getSurfaceTension(0, 1)` and answers `0.02` when that is
    /// not finite and positive, while `RateBasedPackedColumn.estimateSurfaceTension` answers
    /// `0.025` on the same class of pair. Two routes, two fallbacks, one quantity, and the
    /// capture decides: `wetted_area = 72.58281335470625` on the 2.3 m row.
    ///
    /// **The control runs first and it is what makes the reading trustworthy.** Onda's
    /// wetted-area expression reads sigma and the liquid side and no diffusivity, so it is the
    /// one captured output that separates the candidates; but flooding velocity, percent flood
    /// and pressure drop read **no sigma at all**, so if they did not reproduce the capture the
    /// tray or its state would be wrong and the sigma reading would be meaningless. They agree to
    /// `1e-6`, which is what says the middle tray is index 3 of 7 and its state is the capture's.
    ///
    /// The answer is `0.02` - the designer's own fallback, and a third constant distinct from the
    /// rate-based path's. `0.025` gives `63.944` and `0.05` gives `39.804`, both asserted *not* to
    /// reproduce, so the capture is shown to decide rather than merely to agree.
    ///
    /// The residual `3e-7` is the port's own tray density and viscosity rather than NeqSim's: the
    /// capture prints this row's per-tray temperature, pressure and two flows and nothing else,
    /// so those two properties are the port's numbers and are the limit on this reading's
    /// precision. Extending the probe to print them is what would tighten it.
    #[test]
    #[allow(clippy::too_many_lines)] // one measurement, kept linear so it reads as one
    fn the_designers_fallback_reproduces_the_capture_and_the_rate_based_constant_does_not() {
        let components = vec!["methane".to_string(), "n-butane".to_string()];
        let (outcome, _) = distillation_column_outcome(
            None,
            &components,
            7.490704036290964,
            &[0.5, 0.5],
            pascals(2.0e6),
            kelvins(300.0),
            stage_count(2.3),
            2,
            true,
            true,
            pascals(1.9e6),
            pascals(2.0e6),
            Some(kelvins(373.15)),
            Some(kelvins(253.15)),
            1.0e-6,
            200,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
            None,
        )
        .expect("the capture's 2.3 m state solves");

        println!("\ntrays = {}", outcome.trays.len());
        for (index, tray) in outcome.trays.iter().enumerate() {
            println!(
                "  tray {index}: T={:.4} P={:.1} gas_n={:.6} liquid_n={:.6}",
                tray.temperature.value, tray.pressure.value, tray.gas_n, tray.liquid_n
            );
        }

        let middle = outcome.trays.len() / 2;
        let tray = &outcome.trays[middle];
        println!("\nMIDDLE = tray {middle}");

        let gas = Stream::from_pt(
            components.clone(),
            tray.gas_z.clone(),
            tray.gas_n,
            tray.pressure,
            tray.temperature,
        )
        .expect("the tray's vapour rebuilds");
        let liquid = Stream::from_pt(
            components.clone(),
            tray.liquid_z.clone(),
            tray.liquid_n,
            tray.pressure,
            tray.temperature,
        )
        .expect("the tray's liquid rebuilds");
        let gas_view = phase_view(&gas, Pick::Gas).expect("a gas phase");
        let liquid_view = phase_view(&liquid, Pick::Liquid).expect("a liquid phase");

        let gas_mass = tray.gas_n * gas_view.molar_mass;
        let liquid_mass = tray.liquid_n * liquid_view.molar_mass;
        println!(
            "  vapour: rho={:.4} kg/m3  M={:.6} kg/mol  mu={:.6e} Pa.s  m={:.6} kg/s",
            gas_view.density, gas_view.molar_mass, gas_view.transport.mu.value, gas_mass
        );
        println!(
            "  liquid: rho={:.4} kg/m3  M={:.6} kg/mol  mu={:.6e} Pa.s  m={:.6} kg/s",
            liquid_view.density,
            liquid_view.molar_mass,
            liquid_view.transport.mu.value,
            liquid_mass
        );

        println!(
            "\n  sigma   wetted_area          u_flood              percent_flood        dp                   hetp"
        );
        for (label, sigma) in [("0.02 ", 0.02), ("0.025", 0.025), ("0.05 ", 0.05)] {
            let out = packing_hydraulics(
                "Pall-Ring-50",
                PackingState {
                    column_diameter: meters(0.3),
                    packed_height: 2.3,
                    vapor_mass_flow: kilograms_per_second(gas_mass),
                    liquid_mass_flow: kilograms_per_second(liquid_mass),
                    vapor_density: kilograms_per_cubic_meter(gas_view.density),
                    liquid_density: kilograms_per_cubic_meter(liquid_view.density),
                    vapor_viscosity: pascal_seconds(gas_view.transport.mu.value),
                    liquid_viscosity: pascal_seconds(liquid_view.transport.mu.value),
                    surface_tension: newtons_per_meter(sigma),
                    vapor_diffusivity: 3.3e-7,
                    liquid_diffusivity: 1.9e-9,
                    hydraulic_capacity_factor: 1.0,
                },
            )
            .expect("the packing hydraulics runs");
            println!(
                "  {label}  {:.12}  {:.12}  {:.12}  {:.12}  {:.6}",
                out.wetted_area,
                out.flooding_velocity,
                out.percent_flood,
                out.total_pressure_drop.value,
                out.hetp.value
            );
        }

        println!(
            "\n  CAPTURE 72.58281335470625  0.5080730231464395  15.053642072827294  0.6408295279845977  1.44"
        );
        println!(
            "  (u_flood, percent_flood and dp read no sigma: they are the control on the tray state)"
        );

        // **The control first: three outputs that read no sigma pin the tray and the state.**
        // They agree to ~3e-7, so the middle tray is index 3 of 7 and the state it is read at is
        // the capture's. Only then is the sigma reading meaningful.
        let hydraulics = |sigma: f64| {
            packing_hydraulics(
                "Pall-Ring-50",
                PackingState {
                    column_diameter: meters(0.3),
                    packed_height: 2.3,
                    vapor_mass_flow: kilograms_per_second(gas_mass),
                    liquid_mass_flow: kilograms_per_second(liquid_mass),
                    vapor_density: kilograms_per_cubic_meter(gas_view.density),
                    liquid_density: kilograms_per_cubic_meter(liquid_view.density),
                    vapor_viscosity: pascal_seconds(gas_view.transport.mu.value),
                    liquid_viscosity: pascal_seconds(liquid_view.transport.mu.value),
                    surface_tension: newtons_per_meter(sigma),
                    vapor_diffusivity: 3.3e-7,
                    liquid_diffusivity: 1.9e-9,
                    hydraulic_capacity_factor: 1.0,
                },
            )
            .expect("the packing hydraulics runs")
        };
        let close = |got: f64, want: f64| (got - want).abs() / want.abs() < 1e-5;
        let control = hydraulics(DESIGNER_SURFACE_TENSION_N_PER_M);
        assert!(
            close(control.flooding_velocity, 0.508_073_023_146_439_5),
            "the flooding velocity pins the tray state: {}",
            control.flooding_velocity
        );
        assert!(
            close(control.percent_flood, 15.053_642_072_827_294),
            "the load pins the tray state: {}",
            control.percent_flood
        );
        assert!(
            close(control.total_pressure_drop.value, 0.640_829_527_984_597_7),
            "the pressure drop pins the tray state: {}",
            control.total_pressure_drop.value
        );

        // **And then the reading: `wetted_area` is the one output that separates the candidates,
        // because Onda's expression reads sigma and the liquid side and no diffusivity.**
        assert!(
            close(control.wetted_area, 72.582_813_354_706_25),
            "the designer's own 0.02 fallback reproduces the capture's wetted area: {}",
            control.wetted_area
        );
        for (rival, label) in [
            (0.025, "the rate-based path's DEFAULT_SURFACE_TENSION"),
            (0.05, "a mid-range value"),
        ] {
            let out = hydraulics(rival);
            assert!(
                !close(out.wetted_area, 72.582_813_354_706_25),
                "{label} must NOT reproduce it, or the capture cannot decide: {}",
                out.wetted_area
            );
        }
    }
}
