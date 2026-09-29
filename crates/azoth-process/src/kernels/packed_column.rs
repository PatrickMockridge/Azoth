//! `PackedColumn`'s own arithmetic.
//!
//! Spec: `specs/models/process/packed_column.toml`. The class extends `DistillationColumn` and
//! its `run` is `super.run(id)` followed by a packing-hydraulics report, so the separation is
//! [`super::distillation_column`]'s. What this module adds is the stage count a constructor
//! derives from a packed height, and [`packing_report`] - the seven quantities
//! `ColumnInternalsDesigner` reads on the far side of the solve.

use azoth_core::Result;
use azoth_core::units::{
    Length, Pressure, ThermodynamicTemperature, kilograms_per_cubic_meter, kilograms_per_second,
    meters, pascal_seconds, pascals,
};
use azoth_hydraulics::packing_hydraulics::{PackingState, packing_hydraulics_with};
use azoth_hydraulics::packing_sizing::{PackingSizingState, packing_sizing};

use super::distillation_column::{
    ColumnOutcome, ColumnSetup, ReactiveSection, SolverType, Specification, distillation_column,
    tray_streams,
};
use crate::column::murphree::Murphree;
use crate::segment::phase::{Pick, phase_view};
use crate::stream::Stream;

/// The HETP guess `PackedColumn`'s constructors divide a packed height by.
///
/// `estimateStages(packedHeight, 0.5)` states it at each call site, and the class's javadoc
/// calls it "the default HETP of ~0.5 m". The hydraulics report computes a real HETP afterwards
/// and never feeds it back, so this guess is what the solved column's stage count is.
pub const HETP_GUESS_M: f64 = 0.5;

/// The surface tension `ColumnInternalsDesigner.getTrayProperties` falls back to, in N/m.
///
/// **A third constant, and it is not the rate-based path's.** The designer asks
/// `fluid.getInterphaseProperties().getSurfaceTension(0, 1)` and answers `0.02` when that is not
/// finite and positive - while `RateBasedPackedColumn.estimateSurfaceTension` answers
/// [`crate::segment::fallbacks::DEFAULT_SURFACE_TENSION`], `0.025`, on the same class of pair.
/// Two routes, two fallbacks, one physical quantity.
///
/// **Settled by measurement, not by reading.** The capture's 2.3 m row pins
/// `wetted_area = 72.58281335470625`, and Onda's wetted-area expression reads the surface tension
/// and the liquid side and no diffusivity - so it is the one captured output that separates the
/// candidates. At the middle tray this constant reproduces it to `3e-7`; `0.025` answers
/// `63.944` and `0.05` answers `39.804`, both outside any tolerance that admits the first. The
/// three outputs that read no surface tension at all - flooding velocity, percent flood and
/// pressure drop - agree with the capture to `1e-6`, which is what says the tray and its state
/// are the right ones before the reading is taken.
///
/// The triangle is asserted in `models::packed_column::surface_tension`.
pub const DESIGNER_SURFACE_TENSION_N_PER_M: f64 = 0.02;

/// The stage count `PackedColumn.estimateStages` derives from a packed height.
///
/// `ceil(packed_height / 0.5)` floored at two. **The floor is reached rather than refused, and
/// that is measured**: the capture's stage-count rows report two middle trays for a height of
/// `0.2` m and for `0.0` m and for `-1.0` m alike, because `estimateStages` has no domain check
/// and neither does `setPackedHeight`.
pub fn stage_count(packed_height_m: f64) -> usize {
    // `Math.ceil(height / 0.5)` then Java's `(int)` cast, which maps a negative or a `NaN` to
    // zero and saturates a positive overflow at the type's maximum - the same three cases
    // Rust's `as usize` saturating cast maps, so the floor decides the first two.
    let stages = (packed_height_m / HETP_GUESS_M).ceil() as usize;
    stages.max(2)
}

/// Everything a packed column's solve takes.
///
/// **The four packing parameters beside the height are not here, and that is the class's own
/// shape rather than an omission.** `PackedColumn.run` is `super.run(id)` followed by
/// `calcPackingHydraulics()`, which goes through `ColumnInternalsDesigner` on the far side of
/// the converged column - so `packing_type`, `structured_packing`, `design_flood_fraction` and
/// `column_diameter` are read only by that report, and the separation is the base's. They travel
/// to it in [`PackingInputs`] rather than here.
#[derive(Debug, Clone)]
pub struct PackedSetup {
    /// The feed, which the column's own `feed_stage` places.
    pub feed: Stream,
    /// The packed bed's height, which the constructor turns into a stage count.
    pub packed_height: f64,
    pub feed_stage: usize,
    pub has_reboiler: bool,
    pub has_condenser: bool,
    pub top_pressure: Pressure,
    pub bottom_pressure: Pressure,
    pub condenser_temperature: Option<ThermodynamicTemperature>,
    pub reboiler_temperature: Option<ThermodynamicTemperature>,
    pub temperature_tolerance: f64,
    pub max_iterations: usize,
    /// The Murphree efficiency, or `None` for no correction - the base column's own
    /// correction, which `PackedColumn` inherits and overrides no part of.
    pub murphree_efficiency: Option<Murphree>,
    /// **Which trays flash reactively**, which `PackedColumn` inherits from the base and does
    /// not override: its middle trays are the ones the packing's height derives.
    pub reactive: ReactiveSection,
    pub top_specification: Option<Specification>,
    pub bottom_specification: Option<Specification>,
    pub solver_type: SolverType,
}

/// `PackedColumn.run`: **`super.run(id)` at the stage count a constructor derives.**
///
/// The whole of this class's arithmetic is the one line below. `estimateStages` is called by
/// the constructor, so the packing changes *which column* is solved and not how it is solved -
/// and the stage count is a function of the height alone, which is why the four other packing
/// parameters are absent from [`PackedSetup`].
///
/// # Errors
/// Whatever the base column's solve raises.
pub fn packed_column(setup: &PackedSetup) -> Result<ColumnOutcome> {
    distillation_column(&ColumnSetup {
        feed: setup.feed.clone(),
        feed_stage: setup.feed_stage,
        number_of_stages: stage_count(setup.packed_height),
        has_reboiler: setup.has_reboiler,
        has_condenser: setup.has_condenser,
        top_pressure: setup.top_pressure,
        bottom_pressure: setup.bottom_pressure,
        condenser_temperature: setup.condenser_temperature,
        reboiler_temperature: setup.reboiler_temperature,
        temperature_tolerance: setup.temperature_tolerance,
        max_iterations: setup.max_iterations,
        murphree_efficiency: setup.murphree_efficiency.clone(),
        // A packed column is a distillation column and takes the base's correction.
        absorber_murphree: None,
        // `PackedColumn` declares no warm start; the base column takes its cold seed.
        initial_state: None,
        top_specification: setup.top_specification.clone(),
        bottom_specification: setup.bottom_specification.clone(),
        top_feed: None,
        tray_temperatures: None,
        solver_type: setup.solver_type,
        reactive: setup.reactive,
        // The draws are the base column's too, and this entry declares none of them.
        gas_side_draw_fractions: None,
        side_draw_flows: Vec::new(),
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
        pumparound_returns: Vec::new(),
        pumparound_inlets: Vec::new(),
        pumparound_tolerance: None,
        pumparound_max_iterations: None,
    })
}

// ---- The report `ColumnInternalsDesigner` produces on the far side of the solve.

/// The packing name `PackedColumn` carries when none is stated.
///
/// The class's own constructor argument, and the name `PackingSpecificationLibrary` falls back to.
pub const DEFAULT_PACKING_NAME: &str = "Pall-Ring-50";
/// The fraction of flood the class sizes a bed to when none is stated.
pub const DEFAULT_DESIGN_FLOOD_FRACTION: f64 = 0.70;
/// The relative hydraulic capacity factor the class carries when none is stated.
pub const DEFAULT_HYDRAULIC_CAPACITY_FACTOR: f64 = 1.0;
/// The internal diameter `PackedColumn` carries when none is stated.
///
/// **At or below zero rather than absent**, which is what makes it the sizing branch:
/// `setInternalDiameter` is not called, the field keeps this value, and `calcPackingHydraulics`
/// sizes the bed instead of reading a diameter.
pub const UNSIZED_COLUMN_DIAMETER_M: f64 = -1.0;

/// The vapour viscosity `PackingHydraulicsCalculator` holds when nobody sets one, in Pa*s.
///
/// **`ColumnInternalsDesigner.calculatePacked` does not set it, and that is read from the jar
/// rather than guessed.** Disassembling the pinned `neqsim-f0c7436.jar` shows the designer calling
/// exactly `setPackingPreset`, `setDesignFloodFraction`, `setVaporDensity`, `setLiquidDensity`,
/// `setLiquidViscosity`, `setSurfaceTension` and `setColumnDiameter` on the calculator - so this
/// field keeps its constructor initialiser. Handing it the tray's own viscosity instead would be a
/// silent improvement on the class.
///
/// **No quantity `process.packed_column` publishes reads it.** The sensitivity would show in
/// `kga`, `kla` and `htu_*`, which this id does not publish; the one published quantity downstream
/// of them, `hetp`, is clamped to the estimate's own band on every captured row.
pub const CALCULATOR_VAPOR_VISCOSITY_PA_S: f64 = 1.0e-5;
/// The vapour diffusivity the calculator holds when nobody sets one, in m**2/s.
///
/// Set by nobody, for the reason [`CALCULATOR_VAPOR_VISCOSITY_PA_S`] records.
pub const CALCULATOR_VAPOR_DIFFUSIVITY_M2_S: f64 = 1.0e-5;
/// The liquid diffusivity the calculator holds when nobody sets one, in m**2/s.
///
/// Set by nobody, for the reason [`CALCULATOR_VAPOR_VISCOSITY_PA_S`] records.
pub const CALCULATOR_LIQUID_DIFFUSIVITY_M2_S: f64 = 1.0e-9;

/// The five packing parameters the report reads, none of which enters the solve.
///
/// **They are absent from [`PackedSetup`] because the separation does not read them, and they are
/// here because the report does.** That split is the class's own: `PackedColumn.run` is
/// `super.run(id)` followed by `calcPackingHydraulics()`, so the column is solved without them and
/// the report is assembled from them afterwards.
#[derive(Debug, Clone, Default)]
pub struct PackingInputs {
    /// The packing's name. `None` is the class's own [`DEFAULT_PACKING_NAME`].
    pub name: Option<String>,
    /// Whether the packing is structured. **`None` leaves the resolved row's own category**, which
    /// is what a column that never calls `setStructuredPacking` does.
    pub structured: Option<bool>,
    /// The fraction of flood the bed is sized to. `None` is the class's own `0.70`.
    pub design_flood_fraction: Option<f64>,
    /// The relative hydraulic capacity factor. `None` is the class's own `1.0`.
    pub hydraulic_capacity_factor: Option<f64>,
    /// The stated internal diameter. **`None`, or a value at or below zero, takes the sizing
    /// branch**, which is the class's own `-1.0`.
    pub column_diameter: Option<f64>,
}

/// `ColumnInternalsDesigner.calculatePacked`, as far as `process.packed_column` publishes it.
///
/// **Seven quantities, and the split is measured rather than chosen.** The capture prints fourteen
/// hydraulics keys in two blocks of seven: the first from `PackedColumn`'s own accessors and the
/// second from `column.getHydraulics()`. Only the first seven are this model's to publish - the
/// second seven belong to `PackingHydraulicsCalculator` and `hydraulics.packing_hydraulics` already
/// publishes them, so the wiring pins them by calling that calculation at the same state.
#[derive(Debug, Clone, PartialEq)]
pub struct PackedReport {
    /// The height equivalent to a theoretical plate, m.
    pub hetp: Length,
    /// The packed height over the HETP.
    ///
    /// **It describes a different column from the one solved, and it is published that way.** The
    /// constructor fixes the stage count on a `0.5` m HETP guess and this is the real one computed
    /// afterwards: the capture's own 2.3 m case solves five middle trays and reports `1.597`
    /// theoretical stages. The class's own inconsistency, asserted rather than corrected, and filed
    /// upstream.
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
    /// the field is assigned from it and no second predicate exists. `PackingProbe` prints the same
    /// key from `PackingHydraulicsCalculator` directly, so the two ids name one predicate.
    pub hydraulics_ok: bool,
    /// The column's internal diameter, m - the stated one, or the sized one where none was stated.
    pub internal_diameter: Length,
}

/// `ColumnInternalsDesigner.calculatePacked`: the bed's hydraulics at the middle tray's state.
///
/// **The state is the middle tray's, at `trays.len() / 2`**, which is the class's own
/// `trays.size() / 2` and not "the middle one" in any other sense. On the capture's seven-tray row
/// it is tray 3, and the three report outputs that read no surface tension - flooding velocity,
/// percent flood and pressure drop - reproduce the capture there to `1e-6`, which is what says the
/// index is right before the surface-tension reading is taken.
///
/// # Errors
/// * whatever rebuilding the middle tray's vapour and liquid refuses, which is the class's own
///   requirement that the tray's fluid be two-phase;
/// * whatever [`packing_sizing`] or [`packing_hydraulics_with`] raises for the state they are given.
pub fn packing_report(
    outcome: &ColumnOutcome,
    components: &[String],
    packed_height: f64,
    packing: &PackingInputs,
) -> Result<PackedReport> {
    let index = outcome.trays.len() / 2;
    let tray = &outcome.trays[index];
    let (vapour, liquid) = tray_streams(tray, components)?;
    let vapor_view = phase_view(&vapour, Pick::Gas)?;
    let liquid_view = phase_view(&liquid, Pick::Liquid)?;
    let vapor_mass_flow = tray.gas_n * vapor_view.molar_mass;
    let liquid_mass_flow = tray.liquid_n * liquid_view.molar_mass;

    let design_flood_fraction = packing
        .design_flood_fraction
        .unwrap_or(DEFAULT_DESIGN_FLOOD_FRACTION);
    let capacity_factor = packing
        .hydraulic_capacity_factor
        .unwrap_or(DEFAULT_HYDRAULIC_CAPACITY_FACTOR);
    let stated_diameter = packing.column_diameter.unwrap_or(UNSIZED_COLUMN_DIAMETER_M);
    let name = packing.name.as_deref().unwrap_or(DEFAULT_PACKING_NAME);

    // **The sizing branch first, then the hydraulics at the diameter it resolved.** This is the
    // composition `calcPackingHydraulics` performs: a diameter at or below zero is sized by
    // `sizeColumnDiameter`, and the sized value is what the full report then runs at.
    let internal_diameter = if stated_diameter > 0.0 {
        meters(stated_diameter)
    } else {
        packing_sizing(
            name,
            PackingSizingState {
                design_flood_fraction,
                vapor_mass_flow: kilograms_per_second(vapor_mass_flow),
                liquid_mass_flow: kilograms_per_second(liquid_mass_flow),
                vapor_density: kilograms_per_cubic_meter(vapor_view.density),
                liquid_density: kilograms_per_cubic_meter(liquid_view.density),
                liquid_viscosity: pascal_seconds(liquid_view.transport.mu.value),
                hydraulic_capacity_factor: capacity_factor,
            },
        )?
        .column_diameter
    };

    let out = packing_hydraulics_with(
        name,
        PackingState {
            column_diameter: internal_diameter,
            packed_height,
            vapor_mass_flow: kilograms_per_second(vapor_mass_flow),
            liquid_mass_flow: kilograms_per_second(liquid_mass_flow),
            // **The tray's own two densities and its liquid viscosity**, which are what the designer
            // passes and what the capture prints per tray.
            vapor_density: kilograms_per_cubic_meter(vapor_view.density),
            liquid_density: kilograms_per_cubic_meter(liquid_view.density),
            vapor_viscosity: pascal_seconds(CALCULATOR_VAPOR_VISCOSITY_PA_S),
            liquid_viscosity: pascal_seconds(liquid_view.transport.mu.value),
            // **The designer's own fallback, because the interface route answers nothing here.**
            // The capture prints both readings as `0.0` on every tray, so the constant is what
            // answers - and it is a third constant, distinct from the rate-based path's `0.025`.
            surface_tension: azoth_core::units::newtons_per_meter(DESIGNER_SURFACE_TENSION_N_PER_M),
            vapor_diffusivity: CALCULATOR_VAPOR_DIFFUSIVITY_M2_S,
            liquid_diffusivity: CALCULATOR_LIQUID_DIFFUSIVITY_M2_S,
            hydraulic_capacity_factor: capacity_factor,
        },
        packing.structured,
    )?;

    Ok(PackedReport {
        hetp: out.hetp,
        theoretical_stages: out.theoretical_stages,
        percent_flood: out.percent_flood,
        flooding_velocity: out.flooding_velocity,
        packing_pressure_drop: pascals(out.total_pressure_drop.value),
        hydraulics_ok: out.design_ok,
        internal_diameter,
    })
}
