//! The kernels: each unit operation's arithmetic, as a pure function of its
//! inlet streams and parameters.
//!
//! A kernel reads the streams on its inlets, computes the streams on its outlets,
//! and performs no other interaction - which is the process calculus's definition
//! of what a unit operation is. These are the first few, composing the calcs in
//! `azoth-eos` rather than adding new physics.

pub mod absorption_column;
pub mod component_splitter;
pub mod compressor;
pub mod cooler;
pub mod distillation_column;
pub mod ejector;
pub mod expander;
pub mod filter;
pub mod flare;
pub mod gas_scrubber;
pub mod gibbs_reactor;
pub mod heat_exchanger;
pub mod heater;
pub mod manifold;
pub mod mixer;
pub mod packed_column;
pub mod pipe;
pub mod plug_flow_reactor;
pub mod pump;
pub mod separator;
pub mod shortcut_distillation_column;
pub mod splitter;
pub mod stirred_tank_reactor;
pub mod stripping_column;
pub mod tank;
pub mod three_phase_separator;
pub mod throttling_valve;

pub use component_splitter::component_splitter;
pub use compressor::compressor;
pub use cooler::cooler;
pub use distillation_column::{ColumnOutcome, ColumnSetup, TrayProfile, distillation_column};
pub use ejector::{EjectorSetup, ejector};
pub use expander::expander;
pub use filter::filter;
pub use flare::{FlareNumbers, flare};
pub use gas_scrubber::gas_scrubber;
pub use gibbs_reactor::{
    EnergyMode as GibbsEnergyMode, ReactorNumbers as GibbsNumbers, ReactorSetup as GibbsSetup,
    gibbs_reactor,
};
pub use heat_exchanger::heat_exchanger;
pub use heater::heater;
pub use mixer::mixer;
pub use pipe::{FlowRegime, PipeOutcome, pipe};
pub use pump::pump;
pub use separator::separator;
pub use shortcut_distillation_column::{ShortcutColumn, shortcut_distillation_column};
pub use splitter::splitter;
pub use stirred_tank_reactor::{ReactorSetup, stirred_tank_reactor};
pub use tank::tank;
pub use three_phase_separator::{Entrainment, three_phase_separator};
pub use throttling_valve::throttling_valve;
