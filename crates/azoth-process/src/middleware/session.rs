//! The live document: a flowsheet, its diagnostics, and its last run.
//!
//! `middleware.md`'s session/document layer — "the live flowsheet plus its named stream results
//! and a dirty/clean flag". It is the one stateful type in this layer, and every front door is a
//! binding of it rather than a second implementation: a notebook holds one, a browser holds one,
//! and the CLI builds one per invocation because a command line has nowhere to keep it.
//!
//! **Two states, and the difference between them is what a widget reads.** A document is
//! *checked* whenever it changes, which is cheap and local, and *run* only when a caller asks,
//! because a run is the physics. So the diagnostics are always current and the values are not:
//! `dirty` is true from the moment an edit lands until a run succeeds, and a run that failed
//! leaves the previous values gone rather than stale — they described a document that no longer
//! exists, and showing them beside the current one is the mistake this flag exists to prevent.

use std::collections::BTreeMap;

use azoth_core::{AzothError, Result};

use crate::check::{Diagnostic, Severity, validate};
use crate::executor::json::SessionReport;
use crate::executor::{Session, session_report};
use crate::flowsheet::Flowsheet;
use crate::middleware::command::{self, Command};
use crate::middleware::graph::{self, Graph};
use crate::order::ExecutionOrder;
use crate::stream::Stream;
use crate::unit_op::UnitOpSpec;

/// A flowsheet, the palette it is checked against, and what the last check and run said.
pub struct Workspace {
    flowsheet: Flowsheet,
    palette: Vec<UnitOpSpec>,
    diagnostics: Vec<Diagnostic>,
    session: Option<Session>,
    run_error: Option<String>,
    order: ExecutionOrder,
    dirty: bool,
}

impl Workspace {
    /// Open a document from its TOML, and check it.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] on a document this schema cannot read. A document that reads
    /// but does not check opens with its diagnostics rather than failing: a canvas has to be able
    /// to show a broken flowsheet, because showing it is how a user fixes it.
    pub fn open(text: &str, palette: Vec<UnitOpSpec>) -> Result<Self> {
        Ok(Self::of(Flowsheet::from_toml(text)?, palette))
    }

    /// The same, from a flowsheet already in hand.
    #[must_use]
    pub fn of(flowsheet: Flowsheet, palette: Vec<UnitOpSpec>) -> Self {
        let diagnostics = validate(&flowsheet, &palette);
        Self {
            flowsheet,
            palette,
            diagnostics,
            session: None,
            run_error: None,
            order: ExecutionOrder::default(),
            dirty: true,
        }
    }

    /// The document as TOML, which is what a front-end holds and saves.
    ///
    /// # Errors
    /// Whatever `Flowsheet::to_toml` refuses.
    pub fn document(&self) -> Result<String> {
        self.flowsheet.to_toml()
    }

    /// The connection graph, as the canvas draws it.
    ///
    /// # Errors
    /// Whatever the projection refuses - a parameter value JSON cannot carry.
    pub fn graph(&self) -> Result<Graph> {
        graph::graph(&self.flowsheet, &self.palette)
    }

    /// Apply one edit, and re-check.
    ///
    /// **The check happens on every edit, which is why it is the cheap half.** A canvas that had
    /// to ask for diagnostics could draw a red arrow a frame late; a canvas that gets them back
    /// from the edit cannot.
    ///
    /// # Errors
    /// Whatever [`command::apply`] refuses about the command itself. A command that names
    /// something absent is not refused here - it leaves the document as it was and the check that
    /// follows is what says so.
    pub fn apply(&mut self, command: &Command) -> Result<()> {
        command::apply(&mut self.flowsheet, &self.palette, command)?;
        self.recheck();
        self.dirty = true;
        Ok(())
    }

    /// Re-run the checker over the document as it stands.
    pub fn recheck(&mut self) {
        self.diagnostics = validate(&self.flowsheet, &self.palette);
    }

    /// Run the document at its own declared inputs.
    ///
    /// # Errors
    /// Whatever the run refuses: a fluid the databank does not resolve, a port it cannot wire, a
    /// unit that refuses its arguments. The failure is also kept, so a front-end that only reads
    /// the envelope can show it.
    pub fn run(&mut self) -> Result<()> {
        self.run_with(&BTreeMap::new())
    }

    /// Run with boundary values overriding the document's own, keyed by name.
    ///
    /// # Errors
    /// Whatever [`Session::run`] refuses, including a name the document does not declare.
    pub fn run_with(&mut self, feeds: &BTreeMap<String, Stream>) -> Result<()> {
        match Session::run(self.flowsheet.clone(), &self.palette, feeds, self.order) {
            Ok(session) => {
                self.session = Some(session);
                self.run_error = None;
                self.dirty = false;
                Ok(())
            }
            Err(error) => {
                // The values that were showing described the document as it was before the edit
                // that made this run necessary, so they go rather than staying to be misread.
                self.session = None;
                self.run_error = Some(error.to_string());
                self.dirty = true;
                Err(error)
            }
        }
    }

    /// What the checker last said, over the document as it stands.
    #[must_use]
    pub fn diagnostics(&self) -> &[Diagnostic] {
        &self.diagnostics
    }

    /// Whether the document can run: no diagnostic of error severity.
    ///
    /// **A warning does not make a document un-runnable**, which is the line the two severities
    /// are drawn on: an unused feed is declared and consumed by nothing, and the document runs.
    #[must_use]
    pub fn ok(&self) -> bool {
        !self
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.severity() == Severity::Error)
    }

    /// Whether the values are older than the document.
    #[must_use]
    pub fn dirty(&self) -> bool {
        self.dirty
    }

    /// The run's report, where the values are current.
    ///
    /// # Errors
    /// Whatever the codec refuses - a magnitude that is not finite.
    pub fn report(&self) -> Option<Result<SessionReport>> {
        self.session.as_ref().map(session_report)
    }

    /// Every value's path, where the document has run.
    #[must_use]
    pub fn paths(&self) -> Vec<String> {
        self.session
            .as_ref()
            .map(Session::paths)
            .unwrap_or_default()
    }

    /// One value by path, where the document has run.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a path that names no value, naming what it knows.
    pub fn value(&self, path: &str) -> Result<f64> {
        match &self.session {
            Some(session) => session.value(path),
            None => Err(AzothError::invalid_input(
                path,
                "the document has not run, so it has no values",
            )),
        }
    }

    /// Why the last run failed, where one did.
    #[must_use]
    pub fn run_error(&self) -> Option<&str> {
        self.run_error.as_deref()
    }

    #[must_use]
    pub fn flowsheet(&self) -> &Flowsheet {
        &self.flowsheet
    }

    #[must_use]
    pub fn palette(&self) -> &[UnitOpSpec] {
        &self.palette
    }

    /// Which of the class's two orders the next run takes.
    #[must_use]
    pub fn order(&self) -> ExecutionOrder {
        self.order
    }

    /// Set the execution order. **Setting it does not re-run**: an order is a property of the next
    /// run, and a document whose values are current stays current until something changes it.
    pub fn set_order(&mut self, order: ExecutionOrder) {
        self.order = order;
    }
}
