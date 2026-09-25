//! A diagnostic as a front-end receives it.
//!
//! [`crate::check::Diagnostic`] is a Rust enum and reaches a widget as nothing: it does not
//! serialise, and its `Debug` form is a line for a terminal. [`DiagnosticRecord`] is the wire
//! form - a stable `code` to switch on, a `target` to put a mark on, a `message` to show and a
//! `detail` for whatever else the variant carries - and `code`'s and `target`'s derivations live
//! in `check.rs` beside `severity` and `location`, so one file knows what a diagnostic is about
//! however it is encoded.

use serde::Serialize;
use serde_json::json;

use crate::check::{Diagnostic, Target};

/// One diagnostic, as a widget or an agent reads it.
#[derive(Debug, Clone, Serialize)]
pub struct DiagnosticRecord {
    /// The variant's own name, snake-cased and stable.
    pub code: &'static str,
    /// `"error"` or `"warning"`.
    pub severity: &'static str,
    /// The table it is about: `palette`, `instances`, `connections`, `inputs`, `products`,
    /// `flowsheet`.
    pub section: &'static str,
    /// The id within that table, dotted as the document spells it; empty for the table itself.
    pub path: String,
    /// What to put a mark on.
    pub target: Target,
    /// The human line, which is the same sentence the CLI prints.
    pub message: String,
    /// The variant's other fields, for a tooltip or an agent.
    ///
    /// **An object always, and often empty.** The subject a diagnostic names is in `target`, so
    /// what is left here is the variant's own additions - a count, or the unit that is not in the
    /// vocabulary - and a variant with none carries `{}` rather than a missing key, because a
    /// reader switching on `code` should not also have to test for the field's presence.
    pub detail: serde_json::Value,
}

impl From<&Diagnostic> for DiagnosticRecord {
    fn from(diagnostic: &Diagnostic) -> Self {
        let detail = match diagnostic {
            Diagnostic::UnknownDimension { dimension, .. } => json!({ "dimension": dimension }),
            Diagnostic::UnknownParameterUnit { unit, .. } => json!({ "unit": unit }),
            Diagnostic::UnknownUnitOp { unit, .. } => json!({ "unit": unit }),
            Diagnostic::OverfedPort { count, .. } | Diagnostic::UnderfedPort { count, .. } => {
                json!({ "count": count })
            }
            _ => json!({}),
        };
        let location = diagnostic.location();
        Self {
            code: diagnostic.code(),
            severity: diagnostic.severity().name(),
            section: location.section,
            path: location.path,
            target: diagnostic.target(),
            message: diagnostic.message(),
            detail,
        }
    }
}

/// Every diagnostic of a check, in order.
#[must_use]
pub fn records(diagnostics: &[Diagnostic]) -> Vec<DiagnosticRecord> {
    diagnostics.iter().map(DiagnosticRecord::from).collect()
}
