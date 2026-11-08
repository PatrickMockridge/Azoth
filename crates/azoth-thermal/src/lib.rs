//! Heat transfer calculations.
//!
//! The first namespace beside `hydraulics`, and the first proof that the spec
//! pipeline is domain-agnostic rather than hydraulics-shaped: this crate's specs,
//! generated range checks, test cases and documentation all flow from the same
//! `specs/calcs/` tree and the same tooling.
//!
//! * [`conduction_plane_wall`] - steady conduction through a slab

pub mod bell_delaware;
pub mod boiling_heat_transfer;
pub mod conduction_plane_wall;
pub mod lmtd_correction_factor;
pub mod results;
pub mod shah_condensation;
pub mod shell_and_tube_rating;
pub mod shell_and_tube_vibration;
pub mod spec_gen;
pub mod two_phase_pressure_drop;

pub use boiling_heat_transfer::{BoilingRequest, boiling_heat_transfer};
pub use conduction_plane_wall::conduction_plane_wall;
pub use lmtd_correction_factor::{LmtdRequest, lmtd_correction_factor, rate as lmtd_rate};
pub use results::{
    BoilingHeatTransferResult, ConductionPlaneWallResult, LmtdCorrectionFactorResult,
    ShahCondensationResult, ShellAndTubeRatingResult, ShellAndTubeVibrationResult,
    TwoPhasePressureDropResult,
};
pub use shah_condensation::{ShahRequest, shah_condensation};
pub use shell_and_tube_rating::{RatingRequest, ShellSideMethod, rate, shell_and_tube_rating};
pub use shell_and_tube_vibration::{
    VibrationRequest, rate as vibration_rate, shell_and_tube_vibration,
};
pub use two_phase_pressure_drop::{
    TwoPhasePressureDropRequest, rate as two_phase_pressure_drop_rate, two_phase_pressure_drop,
};
