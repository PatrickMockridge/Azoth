//! Heat transfer calculations.
//!
//! The first namespace beside `hydraulics`, and the first proof that the spec
//! pipeline is domain-agnostic rather than hydraulics-shaped: this crate's specs,
//! generated range checks, test cases and documentation all flow from the same
//! `specs/calcs/` tree and the same tooling.
//!
//! * [`conduction_plane_wall`] - steady conduction through a slab

pub mod conduction_plane_wall;
pub mod results;
pub mod spec_gen;

pub use conduction_plane_wall::conduction_plane_wall;
pub use results::ConductionPlaneWallResult;
