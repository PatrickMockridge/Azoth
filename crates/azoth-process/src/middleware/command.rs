//! One edit to a document, as a typed value.
//!
//! `middleware.md`'s command model: every edit is a command that re-runs the checker, which is
//! what makes a canvas's red arrow a structured diagnostic rather than a string. A front-end does
//! not edit TOML text and does not hold a second copy of the document — it sends one command and
//! reads back the document the library wrote.
//!
//! **The checker is the only rule set.** A command is a *structure* edit and nothing more: it
//! adds, removes, or sets a value in the typed [`Flowsheet`], and whatever is wrong with the
//! result is what `validate` says. That is why `remove_instance` on a name nobody declared is a
//! no-op rather than a refusal — the checker already has `UnknownInstance` for it, and a second
//! opinion about the same fact would be a second rule set. Where a command *is* refused, the
//! refusal is about the command and not the document: a value JSON cannot carry, or a field name
//! the schema has no word for.

use std::collections::BTreeMap;

use azoth_core::{AzothError, Result};
use serde::{Deserialize, Serialize};

use crate::flowsheet::{Connection, Flowsheet, Input, Instance, Recycle, split_endpoint};
use crate::unit_op::UnitOpSpec;

/// One edit.
///
/// **Internally tagged on `command`**, so the wire form is one object rather than a tag beside a
/// payload: `{"command": "remove_instance", "id": "hx1"}`. The tag is consumed by the enum's own
/// content deserializer, and the variant structs deny an unknown field, which
/// `a_command_carrying_an_unknown_field_is_refused` measures rather than assumes.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "command", rename_all = "snake_case", deny_unknown_fields)]
pub enum Command {
    /// Add a unit operation, with whatever parameters are known at that moment.
    AddInstance {
        id: String,
        unit: String,
        #[serde(default)]
        parameters: BTreeMap<String, toml::Value>,
    },
    /// Remove a unit operation, and every connection and recycle that names it.
    ///
    /// **Idempotent, and it takes its edges with it.** Removing an instance and leaving its
    /// connections behind would turn one edit into a document full of `UnknownInstance`, so the
    /// command removes what it made unaddressable - which is the one place a command does more
    /// than the name it was given.
    RemoveInstance {
        id: String,
    },
    Connect {
        from: String,
        to: String,
    },
    /// Remove the first connection or recycle whose `from` and `to` are these.
    Disconnect {
        from: String,
        to: String,
    },
    SetParameter {
        instance: String,
        name: String,
        value: serde_json::Value,
    },
    UnsetParameter {
        instance: String,
        name: String,
    },
    AddInput {
        name: String,
        components: Vec<String>,
        n: f64,
        z: Vec<f64>,
        #[serde(rename = "P")]
        p: f64,
        #[serde(rename = "T")]
        t: f64,
    },
    /// Remove a feed, and every connection that draws from it.
    RemoveInput {
        name: String,
    },
    AddProduct {
        name: String,
    },
    /// Remove a product, and every connection that writes to it.
    RemoveProduct {
        name: String,
    },
    AddRecycle {
        stream: String,
        from: String,
        to: String,
    },
    RemoveRecycle {
        stream: String,
    },
    /// Set one of the tear's convergence parameters.
    ///
    /// `field` is one of the seven `[[recycles]]` keys, and `value` is a number for all but
    /// `acceleration_method`, which is a name. The name is not checked here beyond the seven
    /// fields existing: `broyden` is refused by `Recycle::settings` with the measurement that
    /// closes it, and that is the checker's answer rather than a second one.
    SetRecycle {
        stream: String,
        field: String,
        value: serde_json::Value,
    },
    /// Return one of the tear's settings to silence, which is the class's own default.
    ///
    /// **The inverse `set_recycle` needed**, and it needed it because the difference between an
    /// unstated setting and a stated one is real: `[[recycles]]` writes only what it states, and a
    /// form that showed the class's defaults could not tell a declaration that *meant* this number
    /// from one that never mentioned it. A field cleared to empty had nowhere to go, so the
    /// control that looked like it could clear one could not.
    UnsetRecycleField {
        stream: String,
        field: String,
    },
    /// Place a node on the canvas.
    ///
    /// **Addressed by node id**, `instance:sep1`, which is what the projection gives a canvas and
    /// [`crate::split_node_id`] reads back. A `set` needs an addressee where a `remove` does not -
    /// removing something absent is already what was asked for, and placing something absent
    /// cannot take effect at all - so this is the one command that refuses a name the document
    /// does not declare.
    SetPosition {
        node: String,
        x: f64,
        y: f64,
    },
}

/// The seven fields a `[[recycles]]` entry carries: **the one list**, which `set_recycle` and
/// `unset_recycle_field` both check against and both tool schemas publish.
pub const RECYCLE_FIELDS: [&str; 7] = [
    "flow_tolerance",
    "composition_tolerance",
    "temperature_tolerance",
    "pressure_tolerance",
    "max_iterations",
    "minimum_flow",
    "acceleration_method",
];

/// Apply one command to a document.
///
/// **The palette is needed only for the one thing the document does not know**: whether a
/// parameter's value should be a number, a word or a list, which `AddInstance` and
/// `SetParameter` write as the kind the model declares. Every other field is the document's own.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a command that cannot take effect: a JSON value TOML has no
/// shape for, a recycle field the schema does not carry, or a `set_position` naming a node the
/// document does not declare.
pub fn apply(flowsheet: &mut Flowsheet, palette: &[UnitOpSpec], command: &Command) -> Result<()> {
    match command {
        Command::AddInstance {
            id,
            unit,
            parameters,
        } => {
            // **A repeated id replaces rather than duplicating.** The checker has
            // `DuplicateInstance` for a document that holds two, and a canvas that re-sends an
            // add must not be able to create one - so the command is idempotent on its id, which
            // is also what makes it retryable.
            flowsheet.instances.retain(|instance| instance.id != *id);
            flowsheet.instances.push(Instance {
                id: id.clone(),
                unit: unit.clone(),
                parameters: parameters.clone(),
            });
        }
        Command::RemoveInstance { id } => {
            flowsheet.instances.retain(|instance| instance.id != *id);
            let gone = |endpoint: &str| {
                split_endpoint(endpoint).is_some_and(|(instance, _, _)| instance == id)
            };
            flowsheet
                .connections
                .retain(|connection| !gone(&connection.from) && !gone(&connection.to));
            flowsheet.recycles.retain(|recycle| {
                !gone(&recycle.from) && !gone(&recycle.to) && recycle.stream != *id
            });
        }
        Command::Connect { from, to } => {
            if !flowsheet
                .connections
                .iter()
                .any(|connection| connection.from == *from && connection.to == *to)
            {
                flowsheet.connections.push(Connection {
                    from: from.clone(),
                    to: to.clone(),
                });
            }
        }
        Command::Disconnect { from, to } => {
            let mut removed = false;
            flowsheet.connections.retain(|connection| {
                if !removed && connection.from == *from && connection.to == *to {
                    removed = true;
                    return false;
                }
                true
            });
            flowsheet.recycles.retain(|recycle| {
                if !removed && recycle.from == *from && recycle.to == *to {
                    removed = true;
                    return false;
                }
                true
            });
        }
        Command::SetParameter {
            instance,
            name,
            value,
        } => {
            let kind = parameter_kind(flowsheet, palette, instance, name);
            let value = toml_value(value, kind)?;
            if let Some(instance) = flowsheet
                .instances
                .iter_mut()
                .find(|candidate| candidate.id == *instance)
            {
                instance.parameters.insert(name.clone(), value);
            }
        }
        Command::UnsetParameter { instance, name } => {
            if let Some(instance) = flowsheet
                .instances
                .iter_mut()
                .find(|candidate| candidate.id == *instance)
            {
                instance.parameters.remove(name);
            }
        }
        Command::AddInput {
            name,
            components,
            n,
            z,
            p,
            t,
        } => {
            flowsheet.inputs.retain(|input| input.name != *name);
            flowsheet.inputs.push(Input {
                name: name.clone(),
                components: components.clone(),
                n: *n,
                z: z.clone(),
                p: *p,
                t: *t,
            });
        }
        Command::RemoveInput { name } => {
            flowsheet.inputs.retain(|input| input.name != *name);
            flowsheet
                .connections
                .retain(|connection| connection.from != *name && connection.to != *name);
            flowsheet
                .recycles
                .retain(|recycle| recycle.from != *name && recycle.to != *name);
        }
        Command::AddProduct { name } => {
            if !flowsheet.products.contains(name) {
                flowsheet.products.push(name.clone());
            }
        }
        Command::RemoveProduct { name } => {
            flowsheet.products.retain(|product| product != name);
            flowsheet
                .connections
                .retain(|connection| connection.from != *name && connection.to != *name);
            flowsheet
                .recycles
                .retain(|recycle| recycle.from != *name && recycle.to != *name);
        }
        Command::AddRecycle { stream, from, to } => {
            flowsheet
                .recycles
                .retain(|recycle| recycle.stream != *stream);
            flowsheet.recycles.push(Recycle::new(stream, from, to));
        }
        Command::RemoveRecycle { stream } => {
            flowsheet
                .recycles
                .retain(|recycle| recycle.stream != *stream);
        }
        Command::SetRecycle {
            stream,
            field,
            value,
        } => {
            check_recycle_field(field)?;
            let Some(recycle) = flowsheet
                .recycles
                .iter_mut()
                .find(|recycle| recycle.stream == *stream)
            else {
                return Ok(());
            };
            set_recycle(recycle, field, value)?;
        }
        Command::UnsetRecycleField { stream, field } => {
            check_recycle_field(field)?;
            let Some(recycle) = flowsheet
                .recycles
                .iter_mut()
                .find(|recycle| recycle.stream == *stream)
            else {
                return Ok(());
            };
            match field.as_str() {
                "flow_tolerance" => recycle.flow_tolerance = None,
                "composition_tolerance" => recycle.composition_tolerance = None,
                "temperature_tolerance" => recycle.temperature_tolerance = None,
                "pressure_tolerance" => recycle.pressure_tolerance = None,
                "max_iterations" => recycle.max_iterations = None,
                "minimum_flow" => recycle.minimum_flow = None,
                _ => recycle.acceleration_method = None,
            }
        }
        Command::SetPosition { node, x, y } => {
            let Some((role, name)) = crate::split_node_id(node) else {
                return Err(AzothError::invalid_input(
                    "node",
                    format!(
                        "`{node}` is not a node id; they are `instance:sep1`, `input:feed_1` and `product:out`"
                    ),
                ));
            };
            if !declares(flowsheet, role, name) {
                return Err(AzothError::invalid_input(
                    "node",
                    format!("`{node}` names no node this document declares"),
                ));
            }
            let layout = flowsheet.layout.get_or_insert_with(Default::default);
            let placed = match role {
                crate::NodeRole::Instance => &mut layout.instances,
                crate::NodeRole::Input => &mut layout.inputs,
                crate::NodeRole::Product => &mut layout.products,
            };
            placed.insert(name.to_string(), [*x, *y]);
        }
    }
    Ok(())
}

/// The field name a recycle command names, refused by name when it is not one of the seven.
///
/// **One check for `set` and `unset`**, so the two cannot come to different conclusions about what
/// a recycle has - which is the failure a second copy of a seven-name list produces.
fn check_recycle_field(field: &str) -> Result<()> {
    if RECYCLE_FIELDS.contains(&field) {
        return Ok(());
    }
    Err(AzothError::invalid_input(
        "field",
        format!(
            "`{field}` is not a recycle parameter; the seven are {}",
            RECYCLE_FIELDS.join(", ")
        ),
    ))
}

/// Whether a document declares the node a role and a name name.
fn declares(flowsheet: &Flowsheet, role: crate::NodeRole, name: &str) -> bool {
    match role {
        crate::NodeRole::Instance => flowsheet.instances.iter().any(|i| i.id == name),
        crate::NodeRole::Input => flowsheet.inputs.iter().any(|i| i.name == name),
        crate::NodeRole::Product => flowsheet.products.iter().any(|p| p == name),
    }
}

/// The kind a parameter's value has, where a form knows one.
fn parameter_kind(
    flowsheet: &Flowsheet,
    palette: &[UnitOpSpec],
    instance: &str,
    parameter: &str,
) -> Option<&'static str> {
    let unit = flowsheet
        .instances
        .iter()
        .find(|candidate| candidate.id == instance)?
        .unit
        .as_str();
    palette
        .iter()
        .find(|spec| spec.id == unit)
        .and_then(|spec| crate::model_inputs_gen::inputs_for(&spec.id))
        .and_then(|entry| entry.inputs.iter().find(|input| input.name == parameter))
        .map(|input| input.kind)
}

/// A JSON value as the TOML value a document holds.
///
/// **Shaped by the parameter's kind**, because a form's field is not the same fact as a JSON
/// number: a vector parameter's value is a list and a boolean's is a switch, and a front-end that
/// sent `1` where a list belongs would write a document the checker refuses. Where no kind is
/// known the JSON kind is written as itself, which is the same thing.
fn toml_value(value: &serde_json::Value, kind: Option<&str>) -> Result<toml::Value> {
    match kind {
        Some("vector") => {
            let array = value.as_array().ok_or_else(|| {
                AzothError::invalid_input(
                    "value",
                    format!("a list was expected and {value} is not one"),
                )
            })?;
            let numbers = array
                .iter()
                .map(|entry| {
                    entry.as_f64().map(toml::Value::Float).ok_or_else(|| {
                        AzothError::invalid_input(
                            "value",
                            format!("every entry of the list must be a number and {entry} is not"),
                        )
                    })
                })
                .collect::<Result<Vec<_>>>()?;
            Ok(toml::Value::Array(numbers))
        }
        Some("boolean") => value.as_bool().map(toml::Value::Boolean).ok_or_else(|| {
            AzothError::invalid_input(
                "value",
                format!("a switch was expected and {value} is not one"),
            )
        }),
        Some("enum" | "string") => value
            .as_str()
            .map(|text| toml::Value::String(text.to_string()))
            .ok_or_else(|| {
                AzothError::invalid_input(
                    "value",
                    format!("a name was expected and {value} is not one"),
                )
            }),
        // A quantity, or a parameter no model declares: a number where the JSON is a number, and
        // the JSON's own kind otherwise - which the checker's kind rule then rules on.
        _ => match value {
            serde_json::Value::Null => Err(AzothError::invalid_input(
                "value",
                "TOML has no null, so a parameter cannot be set to one",
            )),
            serde_json::Value::Bool(value) => Ok(toml::Value::Boolean(*value)),
            serde_json::Value::Number(number) => Ok(match number.as_i64() {
                Some(integer) => toml::Value::Integer(integer),
                None => toml::Value::Float(number.as_f64().unwrap_or(f64::NAN)),
            }),
            serde_json::Value::String(text) => Ok(toml::Value::String(text.clone())),
            serde_json::Value::Array(items) => items
                .iter()
                .map(|item| toml_value(item, None))
                .collect::<Result<Vec<_>>>()
                .map(toml::Value::Array),
            serde_json::Value::Object(map) => map
                .iter()
                .map(|(key, item)| Ok((key.clone(), toml_value(item, None)?)))
                .collect::<Result<toml::map::Map<_, _>>>()
                .map(toml::Value::Table),
        },
    }
}

/// One of a tear's convergence parameters.
fn set_recycle(recycle: &mut Recycle, field: &str, value: &serde_json::Value) -> Result<()> {
    if field == "acceleration_method" {
        let name = value.as_str().ok_or_else(|| {
            AzothError::invalid_input(
                "value",
                format!("`acceleration_method` takes a name and {value} is not one"),
            )
        })?;
        recycle.acceleration_method = Some(name.to_string());
        return Ok(());
    }
    let number = value.as_f64().ok_or_else(|| {
        AzothError::invalid_input(
            "value",
            format!("`{field}` takes a number and {value} is not one"),
        )
    })?;
    match field {
        "flow_tolerance" => recycle.flow_tolerance = Some(number),
        "composition_tolerance" => recycle.composition_tolerance = Some(number),
        "temperature_tolerance" => recycle.temperature_tolerance = Some(number),
        "pressure_tolerance" => recycle.pressure_tolerance = Some(number),
        "minimum_flow" => recycle.minimum_flow = Some(number),
        "max_iterations" => {
            let count = u32::try_from(number as i64).map_err(|_| {
                AzothError::invalid_input(
                    "value",
                    format!("`max_iterations` takes a whole count and {number} is not one"),
                )
            })?;
            recycle.max_iterations = Some(count);
        }
        _ => unreachable!("the seven fields are checked before this is called"),
    }
    Ok(())
}
