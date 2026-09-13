//! The `chemeng` command-line tool, as a library.
//!
//! Split into a lib and a thin binary so the composition in [`pipe`] can be
//! tested. The interesting logic is not the argument parsing - it is turning a
//! flow rate and a bore into a velocity, choosing a friction factor, adding the
//! fitting loss, and deciding which warnings survive to the report. None of that
//! is reachable from an integration test against a binary.

pub mod pipe;
pub mod report;
