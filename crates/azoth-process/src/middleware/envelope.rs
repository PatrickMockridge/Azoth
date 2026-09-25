//! One call's answer: everything a front-end needs, in one document.
//!
//! **A canvas asks once per edit and gets one object back**, rather than four calls whose answers
//! can interleave with a fifth edit. What is in it is what the layers below produce and nothing
//! re-derived here: the document is `Flowsheet::to_toml`, the graph is `graph::graph`, the
//! diagnostics are `DiagnosticRecord`s, the paths are `Session::paths`, and `session` is
//! `executor::json`'s `SessionReport` embedded as a value. There is no second stream encoder and
//! no second spelling of a path.

use azoth_core::Result;
use serde::Serialize;

use crate::executor::json::SessionReport;
use crate::middleware::diagnostic::{DiagnosticRecord, records};
use crate::middleware::graph::Graph;
use crate::middleware::session::Workspace;

/// Everything one call answers with.
#[derive(Debug, Clone, Serialize)]
pub struct Envelope {
    /// Whether the document can run: no diagnostic of error severity.
    pub ok: bool,
    /// Whether the values are older than the document.
    pub dirty: bool,
    /// Which of the class's two orders the next run takes: `insertion` or `topological`.
    ///
    /// **Here because it is a property of the session and not of a widget.** A top bar's control
    /// reads it from the same object as everything else rather than holding a copy, so a control
    /// and the session cannot disagree after an edit that was refused.
    pub execution_order: &'static str,
    pub flowsheet: FlowsheetView,
    pub diagnostics: Vec<DiagnosticRecord>,
    /// Every value's path, empty until the document has run.
    pub paths: Vec<String>,
    /// The run's report, or `null` where the document has not run or the last run failed.
    pub session: Option<SessionReport>,
    /// Why the last run failed, where one did.
    ///
    /// **A failure the checker cannot see.** The checker holds a document to what is decidable
    /// from the document; whether a named substance exists is the databank's answer, and it
    /// arrives here rather than as a diagnostic, because it is a fact about the world and not
    /// about the document.
    pub run_error: Option<String>,
}

/// The document itself, in the three forms a front-end reads it in.
#[derive(Debug, Clone, Serialize)]
pub struct FlowsheetView {
    pub id: String,
    pub name: String,
    /// The TOML text, which is what a save writes and what a load reads.
    pub document: String,
    pub graph: Graph,
}

/// Build the envelope.
///
/// # Errors
/// Whatever writing the document or projecting the graph refuses - a parameter value JSON cannot
/// carry, which the checker refuses as a `ParameterKind` too.
pub fn envelope(workspace: &Workspace) -> Result<Envelope> {
    Ok(Envelope {
        ok: workspace.ok(),
        dirty: workspace.dirty(),
        execution_order: workspace.order().name(),
        flowsheet: FlowsheetView {
            id: workspace.flowsheet().id.clone(),
            name: workspace.flowsheet().name.clone(),
            document: workspace.document()?,
            graph: workspace.graph()?,
        },
        diagnostics: records(workspace.diagnostics()),
        paths: workspace.paths(),
        session: workspace.report().transpose()?,
        run_error: workspace.run_error().map(str::to_string),
    })
}

/// The envelope as JSON.
///
/// # Errors
/// Whatever [`envelope`] refuses, or a write that fails - which for this shape it cannot.
pub fn to_json(workspace: &Workspace) -> Result<String> {
    serde_json::to_string(&envelope(workspace)?).map_err(|error| {
        azoth_core::AzothError::invalid_input(
            "json",
            format!("the envelope could not be written: {error}"),
        )
    })
}
