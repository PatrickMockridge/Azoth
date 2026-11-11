//! Reading a sheet's inputs for the generated dispatcher.
//!
//! Every helper answers an [`Outcome`] error rather than panicking: a missing or mistyped input
//! is a sentence beside the field, the same shape a kernel's refusal takes.

use azoth_core::Warning;
use serde_json::{Value, json};

use crate::dispatch::Outcome;

fn required(name: &str, wanted: &str) -> Outcome {
    Outcome::Error {
        message: format!("`{name}` is required ({wanted})"),
        field: Some(name.to_string()),
    }
}

fn wrong(name: &str, wanted: &str) -> Outcome {
    Outcome::Error {
        message: format!("`{name}` takes {wanted}"),
        field: Some(name.to_string()),
    }
}

/// A required number.
pub fn number(inputs: &Value, name: &str) -> Result<f64, Outcome> {
    let Some(value) = inputs.get(name) else {
        return Err(required(name, "a number"));
    };
    let Some(number) = value.as_f64() else {
        return Err(wrong(name, "a number"));
    };
    if number.is_finite() {
        Ok(number)
    } else {
        Err(wrong(name, "a finite number"))
    }
}

/// A number that may be absent.
pub fn optional_number(inputs: &Value, name: &str) -> Result<Option<f64>, Outcome> {
    match inputs.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(value) => {
            let Some(number) = value.as_f64() else {
                return Err(wrong(name, "a number"));
            };
            if number.is_finite() {
                Ok(Some(number))
            } else {
                Err(wrong(name, "a finite number"))
            }
        }
    }
}

/// A required list of numbers.
pub fn numbers(inputs: &Value, name: &str) -> Result<Vec<f64>, Outcome> {
    let Some(value) = inputs.get(name) else {
        return Err(required(name, "an array of numbers"));
    };
    let Some(entries) = value.as_array() else {
        return Err(wrong(name, "an array of numbers"));
    };
    entries
        .iter()
        .map(|entry| {
            entry
                .as_f64()
                .ok_or_else(|| wrong(name, "an array of numbers"))
        })
        .collect()
}

/// A list of numbers that may be absent.
pub fn optional_numbers(inputs: &Value, name: &str) -> Result<Option<Vec<f64>>, Outcome> {
    match inputs.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(_) => numbers(inputs, name).map(Some),
    }
}

/// A required flag.
pub fn flag(inputs: &Value, name: &str) -> Result<bool, Outcome> {
    let Some(value) = inputs.get(name) else {
        return Err(required(name, "true or false"));
    };
    value.as_bool().ok_or_else(|| wrong(name, "true or false"))
}

/// A flag that may be absent.
pub fn optional_flag(inputs: &Value, name: &str) -> Result<Option<bool>, Outcome> {
    match inputs.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(value) => value
            .as_bool()
            .map(Some)
            .ok_or_else(|| wrong(name, "true or false")),
    }
}

/// A required name.
pub fn text<'a>(inputs: &'a Value, name: &str) -> Result<&'a str, Outcome> {
    let Some(value) = inputs.get(name) else {
        return Err(required(name, "a name"));
    };
    value.as_str().ok_or_else(|| wrong(name, "a name"))
}

/// A name that may be absent.
pub fn optional_text<'a>(inputs: &'a Value, name: &str) -> Result<Option<&'a str>, Outcome> {
    match inputs.get(name) {
        None | Some(Value::Null) => Ok(None),
        Some(value) => value
            .as_str()
            .map(Some)
            .ok_or_else(|| wrong(name, "a name")),
    }
}

/// A required whole number.
pub fn unsigned(inputs: &Value, name: &str) -> Result<u64, Outcome> {
    let Some(value) = inputs.get(name) else {
        return Err(required(name, "a whole number"));
    };
    if let Some(number) = value.as_u64() {
        return Ok(number);
    }
    // A sheet stores every number as a float, so an integer arrives as `1.0` in JSON. It is still
    // a whole number, and refusing it would make a spec's `model = 1` unsettable.
    match value.as_f64() {
        Some(number) if number >= 0.0 && number.fract() == 0.0 => Ok(number as u64),
        _ => Err(wrong(name, "a whole number")),
    }
}

/// A kernel's refusal, as an outcome.
#[must_use]
pub fn refused(error: azoth_core::AzothError) -> Outcome {
    Outcome::Error {
        message: error.to_string(),
        field: error.field().map(str::to_string),
    }
}

/// The kernel's warnings, as the wire carries them.
#[must_use]
pub fn warnings(list: &[Warning]) -> Vec<Value> {
    list.iter()
        .map(|warning| {
            json!({
                "code": warning.code.as_str(),
                "message": warning.message,
                "field": warning.field,
            })
        })
        .collect()
}
