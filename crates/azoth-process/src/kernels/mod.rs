//! The kernels: each unit operation's arithmetic, as a pure function of its
//! inlet streams and parameters.
//!
//! A kernel reads the streams on its inlets, computes the streams on its outlets,
//! and performs no other interaction - which is the process calculus's definition
//! of what a unit operation is. These are the first few, composing the calcs in
//! `azoth-eos` rather than adding new physics.

pub mod cooler;
pub mod heat_exchanger;
pub mod heater;
pub mod mixer;
pub mod pump;
pub mod separator;
pub mod splitter;
pub mod throttling_valve;

pub use cooler::cooler;
pub use heat_exchanger::heat_exchanger;
pub use heater::heater;
pub use mixer::mixer;
pub use pump::pump;
pub use separator::separator;
pub use splitter::splitter;
pub use throttling_valve::throttling_valve;
