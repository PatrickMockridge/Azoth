//! The catalogue, as one JSON array.
//!
//! **Provenance is joined here, not generated a second time.** The spec projection carries what a
//! form draws; the verification status and the hashes are the core's generated provenance table,
//! looked up by id. A calculation with no entry is shown as unverified rather than omitted.

use azoth_core::provenance_gen::ALL_PROVENANCE;
use azoth_core::{AzothError, Result};
use serde_json::{Map, Value, json};

use crate::catalogue_gen::CALCS;
use crate::form::{CalcForm, ExampleValue};

/// Every registered calculation, in id order.
#[must_use]
pub fn forms() -> &'static [CalcForm] {
    CALCS
}

/// The catalogue as a front end reads it.
///
/// # Errors
/// A projection that cannot be written, which for this shape it cannot.
pub fn to_value() -> Result<Value> {
    let calcs = CALCS
        .iter()
        .map(form_json)
        .collect::<Result<Vec<Value>>>()?;
    Ok(Value::Array(calcs))
}

/// The catalogue as the string a wire carries.
///
/// # Errors
/// A projection that cannot be written, which for this shape it cannot.
pub fn to_json() -> Result<String> {
    serde_json::to_string(&to_value()?).map_err(|error| {
        AzothError::invalid_input(
            "calcs",
            format!("the catalogue could not be written: {error}"),
        )
    })
}

/// One calculation, with its worked example and provenance attached.
fn form_json(form: &CalcForm) -> Result<Value> {
    let mut object = serde_json::to_value(form)
        .map_err(|error| {
            AzothError::invalid_input(
                form.id,
                format!("the calculation's form could not be written: {error}"),
            )
        })?
        .as_object()
        .cloned()
        .unwrap_or_else(Map::new);
    object.insert("worked_example".to_string(), example_json(form));
    object.insert("provenance".to_string(), provenance_json(form.id));
    Ok(Value::Object(object))
}

/// The worked example as an object, because it is an input map and not a list of pairs.
fn example_json(form: &CalcForm) -> Value {
    let mut object = Map::new();
    for entry in form.worked_example {
        object.insert(
            entry.name.to_string(),
            match entry.value {
                ExampleValue::Number(number) => json!(number),
                ExampleValue::Flag(flag) => json!(flag),
                ExampleValue::Name(name) => json!(name),
            },
        );
    }
    Value::Object(object)
}

/// The provenance block a report shows; no upstream class path is in it.
fn provenance_json(id: &str) -> Value {
    let Some(provenance) = ALL_PROVENANCE.iter().find(|entry| entry.calc_id == id) else {
        return json!({
            "verification": "unverified",
            "spec_sha256": null,
            "rust_sha256": null,
            "python_sha256": null,
            "tests_active": 0,
            "tests_skipped": 0,
            "validation_cases": 0,
        });
    };
    json!({
        "verification": provenance.verification.as_str(),
        "spec_sha256": provenance.spec_sha256,
        "rust_sha256": provenance.rust_sha256,
        "python_sha256": provenance.python_sha256,
        "tests_active": provenance.tests_active,
        "tests_skipped": provenance.tests_skipped,
        "validation_cases": provenance.validation_cases,
    })
}
