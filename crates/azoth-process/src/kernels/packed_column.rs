//! `PackedColumn`'s own arithmetic.
//!
//! Spec: `specs/models/process/packed_column.toml`. The class extends `DistillationColumn` and
//! its `run` is `super.run(id)` followed by a packing-hydraulics report, so the separation is
//! [`super::distillation_column`]'s and the one number this module adds is the stage count a
//! constructor derives from a packed height.

use azoth_core::Result;
use azoth_core::units::{Pressure, ThermodynamicTemperature};

use super::distillation_column::{
    ColumnOutcome, ColumnSetup, ReactiveSection, SolverType, Specification, distillation_column,
};
use crate::stream::Stream;

/// The HETP guess `PackedColumn`'s constructors divide a packed height by.
///
/// `estimateStages(packedHeight, 0.5)` states it at each call site, and the class's javadoc
/// calls it "the default HETP of ~0.5 m". The hydraulics report computes a real HETP afterwards
/// and never feeds it back, so this guess is what the solved column's stage count is.
pub const HETP_GUESS_M: f64 = 0.5;

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
/// `column_diameter` are read only by a report this port does not carry, and the separation is
/// the base's.
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
        top_specification: setup.top_specification.clone(),
        bottom_specification: setup.bottom_specification.clone(),
        top_feed: None,
        tray_temperatures: None,
        solver_type: setup.solver_type,
        // The class inherits `setReactive` and this entry does not declare the section: the
        // packings's own parameters are a report, and the reactive section is the base's.
        reactive: ReactiveSection::None,
        // The draws are the base column's too, and this entry declares none of them.
        gas_side_draw_fractions: None,
        liquid_side_draw_fractions: None,
        pumparound_fractions: None,
    })
}
