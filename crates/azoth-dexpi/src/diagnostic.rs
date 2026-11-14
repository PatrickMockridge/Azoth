//! What an exchange found, as data a front-end can switch on.
//!
//! **These are not the checker's diagnostics.** [`azoth_process::check`] answers whether a
//! flowsheet is a valid process; these answer whether a DEXPI document and its projection are
//! sound. The two are reported separately and neither is translated into the other: a
//! `missing_parameter` here means the DEXPI step carried no value for a palette parameter, and
//! the checker's own `MissingParameter` is what finally refuses the document.

use serde::Serialize;

/// How much a diagnostic matters.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Severity {
    /// An exchange that cannot be made, or a document that cannot be read as one.
    Error,
    /// An exchange that can be made, with something named that it lost.
    Warning,
    /// Evidence, not a loss.
    Info,
}

/// One finding about a document or a projection.
#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct DexpiDiagnostic {
    pub severity: Severity,
    /// A stable machine name, such as `unmapped_unit`; the field a client matches on.
    pub code: &'static str,
    /// The DEXPI object id, the azoth endpoint, or the file, whichever the finding names.
    pub subject: String,
    pub message: String,
}

impl DexpiDiagnostic {
    #[must_use]
    pub fn error(
        code: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            severity: Severity::Error,
            code,
            subject: subject.into(),
            message: message.into(),
        }
    }

    #[must_use]
    pub fn warning(
        code: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            severity: Severity::Warning,
            code,
            subject: subject.into(),
            message: message.into(),
        }
    }

    #[must_use]
    pub fn info(
        code: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) -> Self {
        Self {
            severity: Severity::Info,
            code,
            subject: subject.into(),
            message: message.into(),
        }
    }

    /// Whether this finding says the exchange cannot be made.
    #[must_use]
    pub fn is_error(&self) -> bool {
        self.severity == Severity::Error
    }
}
