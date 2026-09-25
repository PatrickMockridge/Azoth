//! One session, and the calls a client may make on it.
//!
//! **This is the shared half of every transport.** `azoth mcp` speaks JSON-RPC over stdio and
//! `azoth serve` speaks HTTP over a socket, and neither of them knows what an edit is: both hand a
//! tool name and its arguments to [`Session::call`] and pass the answer on. A second dispatch, a
//! second place to read a command or a second notion of what a refusal is, is what this module
//! exists to make impossible.
//!
//! **The session is one document, held for the life of the process.** A call sees the effect of
//! the call before it, which is what makes a client's turns cumulative rather than independent -
//! and what makes the document on disk a *starting point* rather than an input: nothing here
//! writes it, because a session is not a file editor. `flowsheet.document` is in every answer, so
//! a client that wants the result already has it.

use std::path::Path;

use azoth_core::{AzothError, Result};
use azoth_process::middleware::command::Command;
use azoth_process::middleware::session::Workspace;
use azoth_process::middleware::{envelope, tools};
use azoth_process::{ExecutionOrder, UnitOpSpec, load_palette};
use serde_json::{Value, json};

/// A document, the palette it is checked against, and whether a call runs it.
pub struct Session {
    workspace: Workspace,
    palette: Vec<UnitOpSpec>,
    run: bool,
}

impl Session {
    /// Open a document, once, for every client this process will serve.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a palette that does not load or a document that does not
    /// read. A document that reads and does not *check* opens: a client has to be able to see a
    /// broken flowsheet, because showing it is how it gets fixed.
    pub fn open(flowsheet: &Path, palette_dir: &Path, run_after_each_call: bool) -> Result<Self> {
        let palette = load_palette(palette_dir)
            .map_err(|error| AzothError::invalid_input("palette", error))?;
        let text = std::fs::read_to_string(flowsheet).map_err(|error| {
            AzothError::invalid_input("flowsheet", format!("{}: {error}", flowsheet.display()))
        })?;
        let workspace = Workspace::open(&text, palette.clone())
            .map_err(|error| AzothError::invalid_input("flowsheet", error.to_string()))?;
        Ok(Self {
            workspace,
            palette,
            run: run_after_each_call,
        })
    }

    /// The tools, which are `middleware::tools`'s own list and nothing else.
    #[must_use]
    pub fn tools(&self) -> Vec<tools::Tool> {
        tools::tools(&self.palette)
    }

    /// One tool call: the envelope, or the sentence that refused the call.
    ///
    /// **`Err` is about the *call*, and `ok: false` is about the *document***. A name that is not a
    /// tool, arguments a command cannot be read from, or a command that cannot take effect are
    /// refusals the caller can correct. A call that landed and left the document refused by the
    /// checker is an `Ok` whose envelope says `ok: false` with the diagnostics, because a broken
    /// document is a state a client shows and refusing it would say the call had not happened.
    ///
    /// **`run` overrides what the process was started with.** `None` is the invocation's own
    /// policy - `azoth mcp` runs unless it was told not to - and `Some(false)` is a client saying
    /// that *this* edit should not pay the physics. It is an `Option` rather than a `bool` because
    /// the two are different statements: "as usual" and "no".
    pub fn call(
        &mut self,
        name: &str,
        arguments: &Value,
        run: Option<bool>,
    ) -> std::result::Result<Value, String> {
        if !self.tools().iter().any(|tool| tool.name == name) {
            return Err(format!("Unknown tool: {name}"));
        }

        // **The tool's name is the command's tag.** `tools.rs` writes it as the schema's `const`,
        // a transport carries the tool in its own field, and the command enum reads it from the
        // object - so this is the single point where the spellings meet, and it is a copy of one
        // string rather than a fourth place it is written down.
        let mut arguments = arguments.clone();
        if let Some(object) = arguments.as_object_mut() {
            object.insert("command".to_string(), json!(name));
        }
        let command: Command =
            serde_json::from_value(arguments).map_err(|error| error.to_string())?;
        self.workspace
            .apply(&command)
            .map_err(|error| error.to_string())?;

        if run.unwrap_or(self.run) {
            // A run the edit made impossible is reported in the envelope rather than here, so the
            // error is deliberately dropped: it is the same fact, said in its own field.
            let _ = self.workspace.run();
        }
        self.envelope()
    }

    /// The envelope as it stands, running the document only if a caller asks for one.
    ///
    /// **A look does not run.** A client polling for somebody else's edit wants the document and
    /// the diagnostics, and a run on every poll would be the physics paid for a fact that had not
    /// changed - so `Some(true)` is the only thing here that runs, and `None` and `Some(false)` are
    /// the same request.
    ///
    /// # Errors
    /// A sentence if the document or the graph cannot be written.
    pub fn settle(&mut self, run: Option<bool>) -> std::result::Result<Value, String> {
        if run == Some(true) {
            let _ = self.workspace.run();
        }
        self.envelope()
    }

    /// Run the document, and answer with the envelope.
    ///
    /// # Errors
    /// A sentence if the document cannot be written, which for a document that opened cannot
    /// happen.
    pub fn run(&mut self) -> std::result::Result<Value, String> {
        let _ = self.workspace.run();
        self.envelope()
    }

    /// Set the order the next run takes: `insertion` or `topological`.
    ///
    /// # Errors
    /// A sentence for a name that is neither.
    pub fn set_order(&mut self, order: &str) -> std::result::Result<Value, String> {
        self.workspace
            .set_order(ExecutionOrder::parse(order).map_err(|error| error.to_string())?);
        self.envelope()
    }

    /// The envelope as it stands.
    ///
    /// # Errors
    /// A sentence if the document or the graph cannot be written.
    pub fn envelope(&self) -> std::result::Result<Value, String> {
        let document = self.envelope_json()?;
        serde_json::from_str(&document).map_err(|error| error.to_string())
    }

    /// The envelope as the codec wrote it.
    ///
    /// **The bytes and not a re-serialisation of them.** A transport that wants to hand a client
    /// the document should hand it *this*: `serde_json::Value`'s object holds its keys sorted, so
    /// a value that went through [`Session::envelope`] would come back in a different order from
    /// the one `executor::json` writes - and the order is the schema's.
    ///
    /// # Errors
    /// A sentence if the document or the graph cannot be written.
    pub fn envelope_json(&self) -> std::result::Result<String, String> {
        envelope::to_json(&self.workspace).map_err(|error| error.to_string())
    }
}
