//! The `azoth` command-line tool, as a library.
//!
//! Split into a lib and a thin binary so everything except the wiring is reachable from
//! a test: the arguments in [`cli`], the composition in [`pipe`] - a flow rate and a bore
//! into a velocity, a friction factor, the fitting loss, the warnings that survive to the
//! report - and the rendering in [`report`]. What is left in `main.rs` is resolving the
//! arguments into a calculation and printing what comes back.

pub mod cli;
pub mod pipe;
pub mod report;
