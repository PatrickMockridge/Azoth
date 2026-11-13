//! The rotating machines' mechanical designs.
//!
//! `ExpanderMechanicalDesign`, `PumpMechanicalDesign` and `CompressorMechanicalDesign`
//! are the classes whose `calcDesign` this crate carries. They are separate from the
//! steady kernels because none of them is part of `run`: a caller asks for the design.

pub mod compressor;
pub mod expander;
pub mod pump;
