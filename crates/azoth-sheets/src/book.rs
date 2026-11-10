//! A saved calculation book.
//!
//! **A document of its own, and not a flowsheet.** A sheet names one registered calculation and
//! the inputs it was given; results are derived and never saved. The inputs are SI magnitudes for
//! quantities, which is the same convention the flowsheet's stream records cross with.

use std::collections::BTreeMap;

use azoth_core::{AzothError, Result};
use serde::{Deserialize, Serialize};

/// The book's id when a document does not give one.
const DEFAULT_ID: &str = "calcbook_1";

/// The book's name when a document does not give one.
const DEFAULT_NAME: &str = "Untitled calculation book";

/// One calculation sheet.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Sheet {
    pub id: String,
    pub name: String,
    pub calculation: String,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inputs: BTreeMap<String, toml::Value>,
}

/// A collection of sheets, as it is saved.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct CalcBook {
    #[serde(default = "default_id")]
    pub id: String,
    #[serde(default = "default_name")]
    pub name: String,
    #[serde(default)]
    pub sheets: Vec<Sheet>,
}

fn default_id() -> String {
    DEFAULT_ID.to_string()
}

fn default_name() -> String {
    DEFAULT_NAME.to_string()
}

impl Default for CalcBook {
    fn default() -> Self {
        Self {
            id: default_id(),
            name: default_name(),
            sheets: Vec::new(),
        }
    }
}

impl CalcBook {
    /// Read a book from its TOML.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a document the schema does not read. An empty document
    /// is a blank book, which is what New starts from.
    pub fn parse(document: &str) -> Result<Self> {
        if document.trim().is_empty() {
            return Ok(Self::default());
        }
        toml::from_str(document).map_err(|error| {
            AzothError::invalid_input("calculation book", format!("it could not be read: {error}"))
        })
    }

    /// The book as TOML, which is what a save writes.
    ///
    /// # Errors
    /// A value the TOML writer cannot express, which this shape does not produce.
    pub fn to_document(&self) -> Result<String> {
        toml::to_string(self).map_err(|error| {
            AzothError::invalid_input(
                "calculation book",
                format!("it could not be written: {error}"),
            )
        })
    }

    /// One sheet by id.
    #[must_use]
    pub fn sheet(&self, id: &str) -> Option<&Sheet> {
        self.sheets.iter().find(|sheet| sheet.id == id)
    }

    /// Whether the book already names a sheet.
    #[must_use]
    pub fn takes(&self, id: &str) -> bool {
        self.sheet(id).is_some()
    }
}
