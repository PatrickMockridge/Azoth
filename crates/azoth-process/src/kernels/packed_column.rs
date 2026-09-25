//! `PackedColumn`'s own arithmetic.
//!
//! Spec: `specs/models/process/packed_column.toml`. The class extends `DistillationColumn` and
//! its `run` is `super.run(id)` followed by a packing-hydraulics report, so the separation is
//! [`super::distillation_column`]'s and the one number this module adds is the stage count a
//! constructor derives from a packed height.

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
