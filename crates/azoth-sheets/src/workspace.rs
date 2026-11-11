//! The calculation book as a live document.
//!
//! **One envelope, as the flowsheet has.** A command edits the book, a calculate runs the sheets
//! that are ready, and every answer is the whole book: what a front end draws is the answer to
//! the call it just made, not a document it reassembled. Results are derived and never saved.

use std::collections::{BTreeMap, BTreeSet};

use azoth_core::{AzothError, Result};
use serde::{Deserialize, Serialize};
use serde_json::{Value, json};

use crate::book::{CalcBook, Sheet};
use crate::catalogue;
use crate::dispatch;
use crate::form::{CalcForm, CalcInput};

/// One edit to a calculation book.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "command", rename_all = "snake_case", deny_unknown_fields)]
pub enum Command {
    AddSheet {
        id: String,
        calculation: String,
        #[serde(default)]
        name: Option<String>,
    },
    RemoveSheet {
        id: String,
    },
    RenameSheet {
        id: String,
        name: String,
    },
    SetInput {
        sheet: String,
        input: String,
        value: Value,
    },
    UnsetInput {
        sheet: String,
        input: String,
    },
    SetBookName {
        name: String,
    },
}

/// A calculation's answer.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CalcOutcome {
    pub outputs: BTreeMap<String, Value>,
    pub warnings: Vec<Value>,
}

/// A calculation the kernel refused.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct CalcError {
    pub message: String,
    pub field: Option<String>,
}

/// What a sheet's last calculate left.
#[derive(Debug, Clone, Default)]
struct Slot {
    result: Option<CalcOutcome>,
    error: Option<CalcError>,
}

/// One structural fact about a sheet, before any calculation runs.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct SheetDiagnostic {
    pub code: &'static str,
    pub message: String,
    pub field: Option<String>,
}

/// The book, its results, and which results are older than their inputs.
#[derive(Debug, Clone)]
pub struct SheetWorkspace {
    book: CalcBook,
    slots: BTreeMap<String, Slot>,
    stale: BTreeSet<String>,
}

impl SheetWorkspace {
    /// Open a book from its TOML.
    ///
    /// # Errors
    /// A document the schema does not read, or one that names a sheet twice.
    pub fn open(document: &str) -> Result<Self> {
        let book = CalcBook::parse(document)?;
        let mut seen = BTreeSet::new();
        for sheet in &book.sheets {
            if !seen.insert(sheet.id.as_str()) {
                return Err(AzothError::invalid_input(
                    "calculation book",
                    format!("sheet id `{}` appears twice", sheet.id),
                ));
            }
        }
        let stale = book.sheets.iter().map(|sheet| sheet.id.clone()).collect();
        Ok(Self {
            book,
            slots: BTreeMap::new(),
            stale,
        })
    }

    /// The book as TOML.
    ///
    /// # Errors
    /// A value the TOML writer cannot express.
    pub fn document(&self) -> Result<String> {
        self.book.to_document()
    }

    /// Apply one edit.
    ///
    /// # Errors
    /// A command that cannot take effect: a duplicate sheet id, an absent sheet, an input the
    /// calculation does not declare, or a value of the wrong kind.
    pub fn apply(&mut self, command: &Command) -> Result<()> {
        match command {
            Command::AddSheet {
                id,
                calculation,
                name,
            } => {
                let form = form_of(calculation).ok_or_else(|| {
                    AzothError::invalid_input(
                        "calculation",
                        format!("`{calculation}` is not a calculation in this build"),
                    )
                })?;
                if !form.available {
                    return Err(AzothError::invalid_input(
                        "calculation",
                        form.refusal
                            .unwrap_or("this calculation is not available in a sheet"),
                    ));
                }
                if self.book.takes(id) {
                    return Err(AzothError::invalid_input(
                        "id",
                        format!("the book already has a sheet `{id}`"),
                    ));
                }
                let name = name
                    .as_deref()
                    .map(str::trim)
                    .filter(|name| !name.is_empty())
                    .unwrap_or(form.name)
                    .to_string();
                self.book.sheets.push(Sheet {
                    id: id.clone(),
                    name,
                    calculation: calculation.clone(),
                    inputs: BTreeMap::new(),
                });
                self.stale.insert(id.clone());
            }
            Command::RemoveSheet { id } => {
                self.book.sheets.retain(|sheet| sheet.id != *id);
                self.slots.remove(id);
                self.stale.remove(id);
            }
            Command::RenameSheet { id, name } => {
                let name = name.trim();
                if name.is_empty() {
                    return Err(AzothError::invalid_input("name", "a sheet needs a name"));
                }
                let sheet = self
                    .book
                    .sheets
                    .iter_mut()
                    .find(|sheet| sheet.id == *id)
                    .ok_or_else(|| AzothError::invalid_input("id", format!("no sheet `{id}`")))?;
                sheet.name = name.to_string();
            }
            Command::SetInput {
                sheet,
                input,
                value,
            } => {
                let declaration = self.input_of(sheet, input)?;
                let stored = stored_value(declaration, input, value)?;
                let sheet = self
                    .book
                    .sheets
                    .iter_mut()
                    .find(|candidate| candidate.id == *sheet)
                    .ok_or_else(|| {
                        AzothError::invalid_input("sheet", format!("no sheet `{sheet}`"))
                    })?;
                sheet.inputs.insert(input.clone(), stored);
                self.stale.insert(sheet.id.clone());
            }
            Command::UnsetInput { sheet, input } => {
                self.input_of(sheet, input)?;
                let sheet = self
                    .book
                    .sheets
                    .iter_mut()
                    .find(|candidate| candidate.id == *sheet)
                    .ok_or_else(|| {
                        AzothError::invalid_input("sheet", format!("no sheet `{sheet}`"))
                    })?;
                sheet.inputs.remove(input);
                self.stale.insert(sheet.id.clone());
            }
            Command::SetBookName { name } => {
                let name = name.trim();
                if name.is_empty() {
                    return Err(AzothError::invalid_input("name", "a book needs a name"));
                }
                self.book.name = name.to_string();
            }
        }
        Ok(())
    }

    /// Calculate one sheet, or every ready sheet when `sheet` is `None`.
    ///
    /// # Errors
    /// A sheet id the book does not carry.
    pub fn calculate(&mut self, sheet: Option<&str>) -> Result<()> {
        if let Some(id) = sheet {
            if self.book.sheet(id).is_none() {
                return Err(AzothError::invalid_input(
                    "sheet",
                    format!("no sheet `{id}`"),
                ));
            }
        }
        let ids: Vec<String> = match sheet {
            Some(id) => vec![id.to_string()],
            None => self
                .book
                .sheets
                .iter()
                .filter(|candidate| self.stale.contains(&candidate.id))
                .map(|candidate| candidate.id.clone())
                .collect(),
        };
        for id in ids {
            let Some(sheet) = self.book.sheet(&id) else {
                continue;
            };
            if !self.diagnostics(sheet).is_empty() {
                continue;
            }
            let inputs = serde_json::to_value(&sheet.inputs).map_err(|error| {
                AzothError::invalid_input(
                    "inputs",
                    format!("they could not be read for the calculation: {error}"),
                )
            })?;
            let slot = match dispatch::run(&sheet.calculation, &inputs) {
                dispatch::Outcome::Ok(outcome) => Slot {
                    result: Some(outcome),
                    error: None,
                },
                dispatch::Outcome::Error { message, field } => Slot {
                    result: None,
                    error: Some(CalcError { message, field }),
                },
            };
            self.slots.insert(id.clone(), slot);
            self.stale.remove(&id);
        }
        Ok(())
    }

    /// The answer a front end draws: the whole book, with per-sheet state.
    ///
    /// # Errors
    /// A book whose document cannot be written.
    pub fn to_json(&self) -> Result<String> {
        let document = self.document()?;
        let mut sheets = Vec::with_capacity(self.book.sheets.len());
        let mut ok = true;
        for sheet in &self.book.sheets {
            let diagnostics = self.diagnostics(sheet);
            let ready = diagnostics.is_empty();
            ok &= ready;
            let missing = self.missing(sheet);
            let slot = self.slots.get(&sheet.id);
            sheets.push(json!({
                "id": sheet.id,
                "name": sheet.name,
                "calculation": sheet.calculation,
                "inputs": serde_json::to_value(&sheet.inputs).unwrap_or_else(|_| json!({})),
                "ready": ready,
                "missing": missing,
                "stale": self.stale.contains(&sheet.id),
                "result": slot.and_then(|slot| slot.result.as_ref()),
                "error": slot.and_then(|slot| slot.error.as_ref()),
                "diagnostics": diagnostics,
            }));
        }
        let envelope = json!({
            "ok": ok,
            "dirty": !self.stale.is_empty(),
            "book": {
                "id": self.book.id,
                "name": self.book.name,
                "document": document,
                "sheets": sheets,
            },
        });
        serde_json::to_string(&envelope).map_err(|error| {
            AzothError::invalid_input(
                "calculation book",
                format!("its envelope could not be written: {error}"),
            )
        })
    }

    /// Why a sheet cannot run yet, if it cannot.
    fn diagnostics(&self, sheet: &Sheet) -> Vec<SheetDiagnostic> {
        let Some(form) = form_of(&sheet.calculation) else {
            return vec![SheetDiagnostic {
                code: "unknown_calculation",
                message: format!("`{}` is not a calculation in this build", sheet.calculation),
                field: None,
            }];
        };
        if !form.available {
            return vec![SheetDiagnostic {
                code: "unavailable_calculation",
                message: form
                    .refusal
                    .unwrap_or("this calculation is not available in a sheet")
                    .to_string(),
                field: None,
            }];
        }
        let mut out = Vec::new();
        for name in sheet.inputs.keys() {
            if !form.inputs.iter().any(|input| input.name == name) {
                out.push(SheetDiagnostic {
                    code: "unknown_input",
                    message: format!("`{name}` is not an input of {}", form.name),
                    field: Some(name.clone()),
                });
            }
        }
        for input in form.inputs {
            let Some(value) = sheet.inputs.get(input.name) else {
                if input.required {
                    out.push(SheetDiagnostic {
                        code: "missing_input",
                        message: format!("`{}` is required", input.name),
                        field: Some(input.name.to_string()),
                    });
                }
                continue;
            };
            if let Err(reason) = kind_of(input, value) {
                out.push(SheetDiagnostic {
                    code: "input_kind",
                    message: format!("`{}` {reason}", input.name),
                    field: Some(input.name.to_string()),
                });
            }
        }
        out
    }

    /// The required inputs a sheet has not set.
    fn missing(&self, sheet: &Sheet) -> Vec<&'static str> {
        let Some(form) = form_of(&sheet.calculation) else {
            return Vec::new();
        };
        form.inputs
            .iter()
            .filter(|input| input.required && !sheet.inputs.contains_key(input.name))
            .map(|input| input.name)
            .collect()
    }

    /// The declaration of one input of one sheet.
    fn input_of(&self, sheet: &str, input: &str) -> Result<&'static CalcInput> {
        let Some(book_sheet) = self.book.sheet(sheet) else {
            return Err(AzothError::invalid_input(
                "sheet",
                format!("no sheet `{sheet}`"),
            ));
        };
        let Some(form) = form_of(&book_sheet.calculation) else {
            return Err(AzothError::invalid_input(
                "calculation",
                format!(
                    "`{}` is not a calculation in this build",
                    book_sheet.calculation
                ),
            ));
        };
        form.inputs
            .iter()
            .find(|candidate| candidate.name == input)
            .ok_or_else(|| {
                AzothError::invalid_input(
                    "input",
                    format!("`{input}` is not an input of {}", form.name),
                )
            })
    }
}

/// The form a sheet names, if this build carries one.
fn form_of(calculation: &str) -> Option<&'static CalcForm> {
    catalogue::forms()
        .iter()
        .find(|form| form.id == calculation)
}

/// A value shaped by the declaration, or a refusal.
fn stored_value(declaration: &CalcInput, name: &str, value: &Value) -> Result<toml::Value> {
    match declaration.kind {
        "quantity" => {
            let Some(number) = value.as_f64() else {
                return Err(AzothError::invalid_input(name, "takes a number"));
            };
            if !number.is_finite() {
                return Err(AzothError::invalid_input(name, "takes a finite number"));
            }
            Ok(toml::Value::Float(number))
        }
        "boolean" => value
            .as_bool()
            .map(toml::Value::Boolean)
            .ok_or_else(|| AzothError::invalid_input(name, "takes true or false")),
        "enum" => {
            let Some(text) = value.as_str() else {
                return Err(AzothError::invalid_input(
                    name,
                    "takes one of its declared names",
                ));
            };
            if !declaration.values.contains(&text) {
                return Err(AzothError::invalid_input(
                    name,
                    format!("takes one of {}", declaration.values.join(", ")),
                ));
            }
            Ok(toml::Value::String(text.to_string()))
        }
        _ => value
            .as_str()
            .map(|text| toml::Value::String(text.to_string()))
            .ok_or_else(|| AzothError::invalid_input(name, "takes text")),
    }
}

/// Why a stored value is the wrong kind for its input, if it is.
fn kind_of(declaration: &CalcInput, value: &toml::Value) -> std::result::Result<(), &'static str> {
    match (declaration.kind, value) {
        ("quantity", toml::Value::Float(number)) if number.is_finite() => Ok(()),
        // A TOML author writes `rho = 998`, and TOML keeps that an integer. Refusing it would
        // make a hand-written book incomplete for a number that is the quantity it names.
        ("quantity", toml::Value::Integer(_)) => Ok(()),
        ("quantity", _) => Err("takes a finite number"),
        ("boolean", toml::Value::Boolean(_)) => Ok(()),
        ("boolean", _) => Err("takes true or false"),
        ("enum", toml::Value::String(text)) if declaration.values.contains(&text.as_str()) => {
            Ok(())
        }
        ("enum", _) => Err("takes one of its declared names"),
        ("string", toml::Value::String(_)) => Ok(()),
        ("string", _) => Err("takes text"),
        _ => Err("has a kind the sheet editor cannot store"),
    }
}
