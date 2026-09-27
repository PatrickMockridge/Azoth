//! The resources a client may read: the document, its diagnostics, its report and the palette.
//!
//! `docs/src/calculus/session.md` states the claim this module is the operator for - **a read is
//! not an edit** - and `docs/src/architecture/middleware.md` is where the surface is described.
//! Every function here takes a `&Workspace` and returns a string, so there is no call in this
//! file that can change a document, and the theorem's model in `lean/Azoth/Session.lean` is the
//! shape of the signatures rather than a property somebody has to keep.
//!
//! **A resource is a second surface for a document every call already returns**, which is why
//! this module is a handful of readings rather than a store. The four are the four things a
//! client asks for by name rather than by the envelope it just received:
//!
//! | uri | what it is | volatile |
//! |---|---|---|
//! | `azoth://document` | the flowsheet, as the TOML a save writes | yes |
//! | `azoth://diagnostics` | what the checker last said | yes |
//! | `azoth://report` | what the last run computed, or why there is nothing | yes |
//! | `azoth://catalogue` | the palette, the units and the substances | no |
//!
//! **The fourth is the one that is not volatile**, and the difference is not tidiness: the
//! catalogue is a fact about the *library*, so a client may cache it for the session, and the
//! other three are facts about a *document*, which an edit changes. That is what `volatile`
//! tells a transport, and `Session::read` beside this module answers with it so the transport
//! does not have to know which is which.

use azoth_core::Result;
use serde::Serialize;

use crate::middleware::diagnostic::records;
use crate::middleware::session::Workspace;
use crate::middleware::{catalogue, envelope};
use crate::unit_op::UnitOpSpec;

/// One resource this server publishes.
#[derive(Debug, Clone, Serialize)]
pub struct Resource {
    /// The URI a client reads it by. The scheme is this server's and the rest is a path.
    pub uri: &'static str,
    /// A human name, and the one a client displays.
    pub name: &'static str,
    /// One sentence on what reading it answers.
    pub description: &'static str,
    /// The media type of the text, which a client may use to choose a reader.
    pub mime_type: &'static str,
    /// **Whether an edit changes what this answers.** A transport turns it into the cache hint
    /// the specification wants, and `false` is a property of the library rather than of a
    /// document: the catalogue is the same for every session this palette serves.
    pub volatile: bool,
}

/// The resources, in the order a client lists them.
pub const RESOURCES: [Resource; 4] = [
    Resource {
        uri: "azoth://document",
        name: "The flowsheet",
        description: "the document as TOML, which is what `azoth save` writes",
        mime_type: "application/toml",
        volatile: true,
    },
    Resource {
        uri: "azoth://diagnostics",
        name: "The diagnostics",
        description: "what the checker last said about the document, as structured records",
        mime_type: "application/json",
        volatile: true,
    },
    Resource {
        uri: "azoth://report",
        name: "The report",
        description: "what the last run computed, or a sentence saying there is nothing current",
        mime_type: "application/json",
        volatile: true,
    },
    Resource {
        uri: "azoth://catalogue",
        name: "The catalogue",
        description: "the palette, the units and the substances: a fact about the library",
        mime_type: "application/json",
        volatile: false,
    },
];

/// One resource's answer.
#[derive(Debug, Clone, Serialize)]
pub struct Reading {
    pub uri: &'static str,
    pub mime_type: &'static str,
    /// The text of the resource.
    pub text: String,
    /// **Whether the document has changed since this text was computed.** The stamp a client's
    /// cache is answerable by: a report read while `dirty` is true describes a document the
    /// session no longer holds, and a client that caches it is caching a stale number. It is the
    /// same flag the envelope carries and not a second mechanism.
    pub dirty: bool,
}

/// Whether this server publishes a URI.
#[must_use]
pub fn known(uri: &str) -> bool {
    RESOURCES.iter().any(|resource| resource.uri == uri)
}

/// **Read one resource, or `None` for a URI this server does not publish.**
///
/// The `None` is deliberate rather than an error type: an unknown URI is a client's mistake
/// about a name, which the transport answers as a bad request; every *known* URI answers with
/// text, including `azoth://report` before the document has run, where the text is the sentence
/// saying so rather than an absence. A resource that failed to read would be one a client has to
/// treat as an error, and "this document has not run" is not one.
#[must_use]
pub fn read(workspace: &Workspace, palette: &[UnitOpSpec], uri: &str) -> Option<Reading> {
    let resource = RESOURCES.iter().find(|resource| resource.uri == uri)?;
    let text = match uri {
        "azoth://document" => document(workspace),
        "azoth://diagnostics" => diagnostics(workspace),
        "azoth://report" => report(workspace),
        "azoth://catalogue" => catalogue(palette, false).ok()?,
        _ => return None,
    };
    Some(Reading {
        uri: resource.uri,
        mime_type: resource.mime_type,
        text,
        dirty: workspace.dirty(),
    })
}

/// The document as TOML. A write that fails is a sentence, because a resource has no error lane:
/// the only way `to_toml` refuses is a value JSON cannot carry, which the edit that made it would
/// have refused first.
fn document(workspace: &Workspace) -> String {
    workspace
        .document()
        .unwrap_or_else(|error| format!("the document cannot be written: {error}"))
}

/// The checker's records, as JSON.
fn diagnostics(workspace: &Workspace) -> String {
    let records = records(workspace.diagnostics());
    serde_json::to_string_pretty(&records)
        .unwrap_or_else(|error| format!("the diagnostics cannot be written: {error}"))
}

/// **What the last run computed, or a sentence saying there is nothing current.**
///
/// **The staleness sentence is the point of this resource.** The values a session holds after an
/// edit describe the document as it was *before* that edit, so answering with them and no word
/// about it would be a client reading a number the session has already disowned - which is the
/// mistake `dirty` exists to prevent, one layer out.
fn report(workspace: &Workspace) -> String {
    match workspace.report() {
        // **Never run, which is not the same as gone stale.** A session opens `dirty` - there are
        // no values, so they are older than any document - and answering "the document has
        // changed since the last run" for a document that has never run would send a client
        // looking for a run it thinks it missed.
        None => "the document has not run".to_string(),
        Some(_) if workspace.dirty() => {
            "the document has changed since the last run, so there is nothing current: \
             read this resource again after a run"
                .to_string()
        }
        Some(Ok(report)) => serde_json::to_string_pretty(&report)
            .unwrap_or_else(|error| format!("the report cannot be written: {error}")),
        Some(Err(error)) => format!("the report cannot be written: {error}"),
    }
}

/// The envelope's own reader, kept here so a caller that wants *everything* has one call.
///
/// **Not a resource**, deliberately: it is what every tool call already answers with, and a
/// resource for it would be a second surface for the thing the first one is made of.
///
/// # Errors
/// Whatever the envelope refuses.
pub fn envelope_json(workspace: &Workspace) -> Result<String> {
    envelope::to_json(workspace)
}
