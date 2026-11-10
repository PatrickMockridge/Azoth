//! `azoth sheets` — the calculation book, headless.
//!
//! The same document the editor's Calculations workspace opens. It answers with the sheet
//! envelope, so a script and a browser read the same object.

use std::path::Path;

use azoth_core::{AzothError, Result};
use azoth_sheets::workspace::{Command, SheetWorkspace};

/// Open a book, apply one command and/or calculate, and answer with the envelope.
///
/// # Errors
/// A book that does not read, a command that cannot take effect, or a sheet id the book does not
/// carry.
pub fn run(book: &Path, command: Option<&str>, calculate: Option<&str>) -> Result<String> {
    let document = std::fs::read_to_string(book).map_err(|error| {
        AzothError::invalid_input("book", format!("{}: {error}", book.display()))
    })?;
    let mut workspace = SheetWorkspace::open(&document)?;
    if let Some(text) = command {
        let command: Command = serde_json::from_str(text).map_err(|error| {
            AzothError::invalid_input("command", format!("it is not a sheet command: {error}"))
        })?;
        workspace.apply(&command)?;
    }
    match calculate {
        Some("all") => workspace.calculate(None)?,
        Some(id) => workspace.calculate(Some(id))?,
        None => {}
    }
    workspace.to_json()
}
