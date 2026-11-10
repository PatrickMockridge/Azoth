//! Calculation sheets: the standalone calculations, as a catalogue a front end can read.
//!
//! **A second document, not a second rule set.** A sheet names one registered calculation and
//! the inputs it was given; the calculation's own spec is the authority on what those inputs
//! mean and what the answer is. This crate projects that registry for a form and runs the
//! library's own kernels by id.
//!
//! The catalogue is generated from `specs/calcs/`. The UI projection deliberately drops prose
//! that names the upstream implementation: the library, CLI and docs keep the full record.

pub mod book;
pub mod catalogue;
pub mod catalogue_gen;
pub mod dispatch;
pub mod dispatch_gen;
pub mod dispatch_support;
pub mod form;
pub mod workspace;
