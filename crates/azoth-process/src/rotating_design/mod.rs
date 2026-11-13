//! The rotating machines' mechanical designs.
//!
//! `ExpanderMechanicalDesign` and `PumpMechanicalDesign` are the two classes whose
//! `calcDesign` this crate carries. They are separate from the steady kernels because
//! neither is part of `run`: a caller asks for the design.

pub mod expander;
