//! A semantic inventory of a Process document, for comparing two exchanges.
//!
//! **Not a byte comparison.** NeqSim's DOM transformer indents differently and its minted ids
//! are its own; what the two implementations must agree on is what the file *says* - the steps,
//! their types, the streams between them, and the quantities. This function reduces a document
//! to that, keyed by the names the file carries rather than by its structural ids, so the
//! comparison does not depend on serialisation order.

use std::collections::BTreeMap;

use crate::model::{DataValue, Document, Object};

/// The quantities the inventory records, with the property that holds them.
const QUANTITIES: &[&str] = &["MassFlow", "Pressure", "Temperature"];

/// Reduce a document to `key=value` pairs.
#[must_use]
pub fn inventory(document: &Document) -> BTreeMap<String, String> {
    let mut pairs = BTreeMap::new();
    pairs.insert("root.name".to_string(), document.model.name.clone());
    pairs.insert("root.uri".to_string(), document.model.uri.clone());
    for import in &document.model.imports {
        pairs.insert(format!("import.{}", import.prefix), import.source.clone());
    }

    let steps: Vec<&Object> = document
        .model
        .children
        .iter()
        .filter_map(|child| match child {
            crate::model::Child::Object(object) => Some(object),
            _ => None,
        })
        .flat_map(|root| step_objects(root))
        .collect();

    let mut port_steps: BTreeMap<String, String> = BTreeMap::new();
    for step in &steps {
        let identifier = step.data_string("Identifier").unwrap_or_default();
        pairs.insert(format!("step.{identifier}.type"), step.ty.clone());
        if let Some(ports) = step.components("Ports") {
            for port in ports.children.iter().filter_map(|child| match child {
                crate::model::ComponentChild::Object(object) => Some(object),
                _ => None,
            }) {
                if let Some(id) = &port.id {
                    port_steps.insert(id.clone(), identifier.to_string());
                }
            }
        }
    }

    let streams = stream_objects(document);
    for stream in &streams {
        let label = stream
            .data_string("Label")
            .or_else(|| stream.data_string("Identifier"))
            .unwrap_or_default();
        let source = stream
            .references("Source")
            .and_then(|references| references.objects.first())
            .map(|target| {
                port_steps
                    .get(target.trim_start_matches('#'))
                    .cloned()
                    .unwrap_or_default()
            })
            .unwrap_or_default();
        let target = stream
            .references("Target")
            .and_then(|references| references.objects.first())
            .map(|target| {
                port_steps
                    .get(target.trim_start_matches('#'))
                    .cloned()
                    .unwrap_or_default()
            })
            .unwrap_or_default();
        pairs.insert(format!("stream.{label}.source"), source);
        pairs.insert(format!("stream.{label}.target"), target);
        for property in QUANTITIES {
            let Some(components) = stream.components(property) else {
                continue;
            };
            for qualified in components.children.iter().filter_map(|child| match child {
                crate::model::ComponentChild::Object(object) => Some(object),
                _ => None,
            }) {
                let Some(value) = qualified.data_value("Value") else {
                    continue;
                };
                if let Some(magnitude) = value.quantity_magnitude() {
                    pairs.insert(
                        format!("stream.{label}.{property}"),
                        crate::writer::format_double(magnitude),
                    );
                }
                if let Some(unit) = value.quantity_unit() {
                    pairs.insert(format!("stream.{label}.{property}.unit"), unit.to_string());
                }
            }
        }
    }

    pairs.insert("step.count".to_string(), steps.len().to_string());
    pairs.insert("stream.count".to_string(), streams.len().to_string());
    pairs
}

fn step_objects(root: &Object) -> Vec<&Object> {
    let mut steps = Vec::new();
    let Some(conceptual) = root.components("ConceptualModel") else {
        return steps;
    };
    for model in conceptual.children.iter().filter_map(|child| match child {
        crate::model::ComponentChild::Object(object) => Some(object),
        _ => None,
    }) {
        let Some(process_steps) = model.components("ProcessSteps") else {
            continue;
        };
        for step in process_steps
            .children
            .iter()
            .filter_map(|child| match child {
                crate::model::ComponentChild::Object(object) => Some(object),
                _ => None,
            })
        {
            steps.push(step);
        }
    }
    steps
}

fn stream_objects(document: &Document) -> Vec<&Object> {
    let mut streams = Vec::new();
    for child in &document.model.children {
        let crate::model::Child::Object(root) = child else {
            continue;
        };
        let Some(conceptual) = root.components("ConceptualModel") else {
            continue;
        };
        for model in conceptual.children.iter().filter_map(|child| match child {
            crate::model::ComponentChild::Object(object) => Some(object),
            _ => None,
        }) {
            let Some(connections) = model.components("ProcessConnections") else {
                continue;
            };
            for stream in connections.children.iter().filter_map(|child| match child {
                crate::model::ComponentChild::Object(object) => Some(object),
                _ => None,
            }) {
                streams.push(stream);
            }
        }
    }
    streams
}

/// A stream's identifier, or its label when it has none.
#[must_use]
pub fn stream_label(object: &Object) -> Option<&str> {
    object
        .data_string("Label")
        .or_else(|| object.data_string("Identifier"))
}

/// The magnitude and unit of one quantity, if the stream carries it.
#[must_use]
pub fn quantity<'a>(object: &'a Object, property: &str) -> Option<(f64, &'a str)> {
    let components = object.components(property)?;
    for child in &components.children {
        let crate::model::ComponentChild::Object(qualified) = child else {
            continue;
        };
        let value = qualified.data_value("Value")?;
        let magnitude = value.quantity_magnitude()?;
        let unit = value.quantity_unit()?;
        return Some((magnitude, unit));
    }
    None
}

/// The string of a `DataValue`, where it is one.
#[must_use]
pub fn as_string(value: &DataValue) -> Option<&str> {
    match value {
        DataValue::String(text) => Some(text),
        _ => None,
    }
}
