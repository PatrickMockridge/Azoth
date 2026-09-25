//! The azoth middleware, compiled for a browser.
//!
//! **This crate holds no logic.** Every function below is a call into
//! `azoth_process::middleware` and a string out, which is what keeps a browser a consumer of the
//! wire layer rather than a second implementation of it - and what keeps the whole surface
//! covered by the tests that run on the host.
//!
//! **The boundary is strings**, deliberately: a document in, a JSON document out. Nothing here
//! needs `serde-wasm-bindgen` or a JS object graph, and a caller that would rather speak HTTP
//! would call the same four functions.
//!
//! **Two measurements the crate rests on**, both taken rather than assumed:
//!
//! * `rustc --print cfg --target wasm32-unknown-unknown` reports `panic="abort"`, so a panic
//!   traps the module instance rather than unwinding. Every fallible path is a `Result` in the
//!   layer below, and this crate maps each one to a JS exception instead of letting it reach a
//!   panic. There is no separate wasm profile: the workspace's `panic = "unwind"` is inert on
//!   this target, and a profile that changed nothing would be a setting with no measurement
//!   behind it.
//! * The databank is already embedded (`include_str!` throughout `azoth-eos`,
//!   `azoth-reactions`, `azoth-standards`), so the module carries its own data with no fetch.
//!   Measured: the release module is 2.1 MB with the full component, UNIFAC, Pitzer, MBWR and
//!   reaction tables in it.
//!
//! There is deliberately no `[lints] workspace = true`: `#[wasm_bindgen]` expands to `unsafe
//! extern` blocks, and the workspace forbids `unsafe_code`, which an item-level `allow` cannot
//! relax. Keeping this crate free of logic is what makes that acceptable - `cargo clippy
//! --workspace` still lints it, and the wire layer keeps the `forbid`.

use wasm_bindgen::prelude::*;

use azoth_process::UnitOpSpec;
use azoth_process::load_palette_text;
use azoth_process::middleware::command::Command;
use azoth_process::middleware::session::Workspace;
use azoth_process::middleware::{catalogue, envelope};
use azoth_process::palette_gen::PALETTE;

/// The palette the module carries.
fn palette() -> Vec<UnitOpSpec> {
    // The bundle is generated from the shipped specs and held to them byte for byte by
    // `the_embedded_palette_is_the_directory`, so a failure here is a broken build rather than a
    // caller's mistake.
    load_palette_text(PALETTE).expect("the embedded palette parses")
}

/// The palette as one form per unit operation, and - with `with_tools` - the agent's tools.
///
/// # Errors
/// A JS exception if the document cannot be written, which the embedded palette cannot cause.
#[wasm_bindgen]
pub fn palette_json(with_tools: bool) -> Result<String, JsValue> {
    catalogue(&palette(), with_tools).map_err(failed)
}

/// A live document: the session, as the browser holds it.
///
/// **The one stateful thing on this side of the boundary.** A canvas opens one, applies commands
/// to it, and asks it to run; each call answers with the whole envelope, so what a front-end
/// draws is always the answer to the call it just made rather than a document it reassembled.
#[wasm_bindgen]
pub struct Editor {
    workspace: Workspace,
}

#[wasm_bindgen]
impl Editor {
    /// Open a document from its TOML.
    ///
    /// # Errors
    /// A JS exception on a document this schema cannot read. A document that reads but does not
    /// check opens with its diagnostics: a canvas has to be able to show a broken flowsheet.
    #[wasm_bindgen(constructor)]
    pub fn new(document: &str) -> Result<Editor, JsValue> {
        Workspace::open(document, palette())
            .map(|workspace| Editor { workspace })
            .map_err(failed)
    }

    /// Apply one command, and answer with the envelope it left.
    ///
    /// The command is the JSON object `middleware::command::Command` reads, e.g.
    /// `{"command":"remove_instance","id":"hx1"}`.
    ///
    /// # Errors
    /// A JS exception for a command that is not one, or one that cannot take effect. A command
    /// that names something absent is *not* an error: the envelope that comes back carries the
    /// checker's own diagnostic for it.
    pub fn apply(&mut self, command: &str) -> Result<String, JsValue> {
        let command: Command = serde_json::from_str(command).map_err(failed)?;
        self.workspace.apply(&command).map_err(failed)?;
        self.envelope()
    }

    /// Run the document, and answer with the envelope it left.
    ///
    /// # Errors
    /// A JS exception if the run refuses the document - an unresolved substance, a port that
    /// cannot be wired - which the envelope also reports, because a run that failed is a state a
    /// canvas draws rather than an exception it cannot read.
    pub fn run(&mut self) -> Result<String, JsValue> {
        let _ = self.workspace.run();
        self.envelope()
    }

    /// The envelope as it stands, without running anything.
    ///
    /// # Errors
    /// A JS exception if the document or the graph cannot be written.
    pub fn envelope(&self) -> Result<String, JsValue> {
        envelope::to_json(&self.workspace).map_err(failed)
    }

    /// The document as TOML, which is what a save writes.
    ///
    /// # Errors
    /// A JS exception if the document cannot be written.
    pub fn document(&self) -> Result<String, JsValue> {
        self.workspace.document().map_err(failed)
    }

    /// The connection graph alone, for a canvas that re-draws without re-reading the envelope.
    ///
    /// # Errors
    /// A JS exception if a parameter value cannot be written as JSON.
    pub fn graph(&self) -> Result<String, JsValue> {
        serde_json::to_string(&self.workspace.graph().map_err(failed)?).map_err(failed)
    }

    /// One value by path, e.g. `p1.outlet.P`.
    ///
    /// # Errors
    /// A JS exception for a path that names no value, or before the document has run.
    pub fn value(&self, path: &str) -> Result<f64, JsValue> {
        self.workspace.value(path).map_err(failed)
    }

    /// Whether the document can run.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn ok(&self) -> bool {
        self.workspace.ok()
    }

    /// Whether the values are older than the document.
    #[wasm_bindgen(getter)]
    #[must_use]
    pub fn dirty(&self) -> bool {
        self.workspace.dirty()
    }
}

/// Every refusal, as a JS exception carrying the library's own sentence.
///
/// One place, so a binding cannot lose the reason a call failed: "calculation failed" is not
/// actionable and the layer below always says what was wrong instead.
fn failed(error: impl std::fmt::Display) -> JsValue {
    JsValue::from_str(&error.to_string())
}
