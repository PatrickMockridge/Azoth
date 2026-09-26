//! The rate-based packed column's segment model: the internals of
//! `process.rate_based_packed_column`.
//!
//! **Not a palette surface, and this paragraph is why.** `unit_ops.rate_based_packed_column`
//! is one entry and its four streams are the whole interface; a segment is not a unit
//! operation a flowsheet connects. So a reader looking for `unit_ops.segment` finds nothing,
//! which is the correct answer, and this is where it is stated rather than left as an absence.
//!
//! **A segment cannot be oracled alone, and that is a difference from `column/`.** A tray is a
//! `SimulationInterface` and runs standalone, which is why `column_tray.tsv` pins its
//! arithmetic before a column exists to contain it. `calculateSegment` is private, takes two
//! flashed systems and writes into them, and there is no probe that drives one - so the
//! segment's arithmetic is pinned through the column's own capture, and
//! `crates/azoth-process/tests/rate_based_packed_column.rs` is where that happens.
//!
//! **The packing is inside the equations here.** Unlike `unit_ops.packed_column`, whose
//! packing is a report read after the solve, this model's wetted area and film coefficients
//! *are* the transfer: [`transport`] calls `hydraulics.packing_hydraulics` every segment and
//! the fluxes are built from what it answers.

pub mod equilibrium;
pub mod fallbacks;
pub mod film;
pub mod heat;
pub mod phase;
pub mod profile;
pub mod step;
pub mod transport;

pub use fallbacks::{Fallbacks, Property};
pub use phase::{PhaseView, Pick, phase_view};
pub use profile::{ColumnOutcome, ProfileSettings, solve_fixed_point_profile};
pub use step::{
    MAX_HEAT_TRANSFER_FRACTION, MAX_TRANSFER_FRACTION, SegmentComputation, SegmentResult,
    calculate_segment, component_moles,
};
pub use transport::{SnapshotSettings, TransportSnapshot, calculate_transport_snapshot};
