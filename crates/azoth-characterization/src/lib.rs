//! Petroleum-fraction characterisation.
//!
//! The front door for a fluid that has no databank row. A crude is an *assay* — a list of
//! cuts, each a molar mass and a density — and everything a cubic needs about a cut follows
//! from those two numbers by correlation. NeqSim's `thermo/characterization/` is that
//! subsystem, and this crate ports it.
//!
//! It is a namespace of its own rather than a family under `eos`, mirroring NeqSim's own
//! package split: what is here is not an equation of state and is not reached from one. A cut
//! is characterised first and *then* handed to a mixture, which is the caller's step.
//!
//! * [`tbp_cut_properties`] — a cut's critical properties, by any of NeqSim's ten models

pub mod model_gen;
pub mod pedersen_plus_split;
pub mod results;
pub mod spec_gen;
pub mod tbp_closure;
pub mod tbp_cut_properties;
pub mod tbp_density;
pub mod unported;
pub mod unported_gen;
pub mod whitson_gamma_split;

pub use pedersen_plus_split::pedersen_plus_split;
pub use results::{
    PedersenPlusSplitResult, TbpClosureResult, TbpCutPropertiesResult, TbpDensityResult,
    WhitsonGammaSplitResult,
};
pub use tbp_closure::{TbpClosureKind, tbp_closure};
pub use tbp_cut_properties::{TbpModel, tbp_cut_properties};
pub use tbp_density::tbp_density;
pub use whitson_gamma_split::{WhitsonDensityModel, gamma, p0_p1, whitson_gamma_split};
