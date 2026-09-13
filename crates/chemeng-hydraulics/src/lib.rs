//! Hydraulics calculations.
//!
//! The slice implemented here runs from the Reynolds number through to the
//! Darcy-Weisbach pressure drop for a straight pipe, with fitting losses
//! available separately:
//!
//! * [`reynolds_number`] - the flow regime, and the input every other calc needs
//! * [`friction_factor_colebrook`] - the implicit, accurate friction factor
//! * [`friction_factor_swamee_jain`] - the explicit approximation to it
//! * [`crane_k_factors`] - fitting losses by the equivalent-length method
//! * [`darcy_weisbach`] - pressure drop over a straight pipe
//!
//! Pipe *with fittings* is a composition of the last two, performed by the
//! `chemeng pipe` CLI rather than by a calc of its own, because the two losses
//! are computed by different methods and adding them is a modelling decision the
//! caller should be able to see.

pub mod results;

pub use results::{
    ColebrookResult, DarcyWeisbachResult, KComponent, KFactorsResult, ReynoldsNumberResult,
    SwameeJainResult,
};
