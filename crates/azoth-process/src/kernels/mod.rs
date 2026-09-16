//! The kernels: each unit operation's arithmetic, as a pure function of its
//! inlet streams and parameters.
//!
//! A kernel reads the streams on its inlets, computes the streams on its outlets,
//! and performs no other interaction - which is the process calculus's definition
//! of what a unit operation is. These are the first few, composing the calcs in
//! `azoth-eos` and `azoth-hydraulics` rather than adding new physics.

pub mod mixer;
pub mod separator;
pub mod splitter;

pub use mixer::mixer;
pub use separator::separator;
pub use splitter::splitter;
