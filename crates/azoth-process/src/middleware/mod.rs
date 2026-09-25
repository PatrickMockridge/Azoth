//! The interoperation surface: what a front-end reads and writes.
//!
//! `docs/src/architecture/middleware.md` states what sits between the schema and the things
//! that drive it - a flowsheet editor, a notebook, an agent. This module is that surface, and
//! every one of those is a **binding** of it rather than a second implementation: the session
//! is `session::Workspace`, the graph is `graph`, a unit op's declaration is `form`, an edit is
//! a `command::Command`, and one call's answer is the `envelope::Envelope`.
//!
//! **Nothing here is a second rule set.** The checker decides what a document may say, the
//! executor decides what it computes, and `executor::json` writes the values; this layer
//! projects those and adds no verdict of its own.

pub mod diagnostic;
pub mod form;
pub mod graph;
