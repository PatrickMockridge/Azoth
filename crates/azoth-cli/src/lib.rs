//! The `azoth` command-line tool, as a library.
//!
//! Split into a lib and a thin binary so everything except the wiring is reachable from
//! a test: the arguments in [`cli`], the composition in [`pipe`] - a flow rate and a bore
//! into a velocity, a friction factor, the fitting loss, the warnings that survive to the
//! report - and the rendering in [`report`]. What is left in `main.rs` is resolving the
//! arguments into a calculation and printing what comes back.
//!
//! [`edit`] and [`forms`] are the middleware's two doors here: one command against a document,
//! and the palette as the forms a front-end renders.

pub mod check;
pub mod cli;
pub mod edit;
pub mod forms;
pub mod mcp;
pub mod mcp_http;
pub mod pipe;
pub mod report;
pub mod run;
pub mod serve;
pub mod session;
