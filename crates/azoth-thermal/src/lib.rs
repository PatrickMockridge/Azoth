//! Heat transfer calculations.
//!
//! The first namespace beside `hydraulics`, and the first proof that the spec
//! pipeline is domain-agnostic rather than hydraulics-shaped: this crate's specs,
//! generated range checks, test cases and documentation all flow from the same
//! `specs/calcs/` tree and the same tooling.
//!
//! * [`conduction_plane_wall`] - steady conduction through a slab

pub mod bell_delaware;
pub mod conduction_plane_wall;
pub mod results;
pub mod shell_and_tube_rating;
pub mod spec_gen;

pub use conduction_plane_wall::conduction_plane_wall;
pub use results::{ConductionPlaneWallResult, ShellAndTubeRatingResult};
pub use shell_and_tube_rating::{RatingRequest, ShellSideMethod, rate, shell_and_tube_rating};
