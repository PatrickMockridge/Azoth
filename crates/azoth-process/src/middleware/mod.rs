//! The interoperation surface: what a front-end reads and writes.
//!
//! `docs/src/architecture/middleware.md` states what sits between the schema and the things
//! that drive it - a flowsheet editor, a notebook, an agent. This module is that surface, and
//! every one of those is a **binding** of it rather than a second implementation: the session
//! is [`session::Workspace`], the graph is [`graph`], a unit op's declaration is [`form`], an
//! edit is a [`command::Command`], and one call's answer is the [`envelope::Envelope`].
//!
//! **Nothing here is a second rule set.** The checker decides what a document may say, the
//! executor decides what it computes, and `executor::json` writes the values; this layer
//! projects those and adds no verdict of its own.

pub mod command;
pub mod diagnostic;
pub mod envelope;
pub mod form;
pub mod graph;
pub mod session;
pub mod tools;

use azoth_core::Result;

use crate::unit_op::UnitOpSpec;

/// The palette as a front-end reads it: one form per entry, and optionally the agent's tools.
///
/// **One assembly, so a binding cannot spell the document differently.** The CLI and the wasm
/// module differ in where the palette comes from and in nothing else: both call this with the
/// specs they loaded, and `tools` is `false` for a widget generator and `true` for an agent's
/// tool list - the same declaration read two ways rather than two documents.
///
/// # Errors
/// A write that fails, which for this shape it cannot.
pub fn catalogue(palette: &[UnitOpSpec], with_tools: bool) -> Result<String> {
    let mut document = serde_json::json!({ "unit_ops": form::forms(palette) });
    if with_tools {
        let object = document
            .as_object_mut()
            .expect("the object this function just built");
        object.insert(
            "tools".to_string(),
            serde_json::json!(tools::tools(palette)),
        );
    }
    serde_json::to_string_pretty(&document).map_err(|error| {
        azoth_core::AzothError::invalid_input(
            "catalogue",
            format!("it could not be written: {error}"),
        )
    })
}
