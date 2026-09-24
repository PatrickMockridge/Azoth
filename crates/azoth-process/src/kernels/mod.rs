//! The kernels: each unit operation's arithmetic, as a pure function of its
//! inlet streams and parameters.
//!
//! A kernel reads the streams on its inlets, computes the streams on its outlets,
//! and performs no other interaction - which is the process calculus's definition
//! of what a unit operation is. These are the first few, composing the calcs in
//! `azoth-eos` rather than adding new physics.

pub mod compressor;
pub mod cooler;
pub mod expander;
pub mod filter;
pub mod gas_scrubber;
pub mod heat_exchanger;
pub mod heater;
pub mod manifold;
pub mod mixer;
pub mod pipe;
pub mod pump;
pub mod separator;
pub mod shortcut_distillation_column;
pub mod splitter;
pub mod throttling_valve;

pub use compressor::compressor;
pub use cooler::cooler;
pub use expander::expander;
pub use filter::filter;
pub use gas_scrubber::gas_scrubber;
pub use heat_exchanger::heat_exchanger;
pub use heater::heater;
pub use mixer::mixer;
pub use pipe::{FlowRegime, PipeOutcome, pipe};
pub use pump::pump;
pub use separator::separator;
pub use shortcut_distillation_column::{ShortcutColumn, shortcut_distillation_column};
pub use splitter::splitter;
pub use throttling_valve::throttling_valve;
