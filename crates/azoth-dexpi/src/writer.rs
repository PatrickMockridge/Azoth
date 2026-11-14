//! DEXPI XML out of the typed metamodel.
//!
//! **This is the only place XML syntax is produced.** A caller builds a [`Document`] and this
//! serialises it, so escaping, attribute order and indentation are decided once. The output is
//! deterministic and idempotent: `write(read(write(document)))` is byte-identical to
//! `write(document)`, which is held by a test rather than asserted here.
//!
//! A non-finite `Double` is refused rather than written: XML Schema's `xsd:double` has lexical
//! forms for `NaN` and `INF`, and a DEXPI file carrying either is a quantity nothing should
//! consume.

use quick_xml::Writer;
use quick_xml::events::{BytesDecl, BytesEnd, BytesStart, BytesText, Event};

use crate::model::{
    Child, ComponentChild, Components, Data, DataValue, Document, Import, Model, Object, References,
};

/// Why a document could not be written.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum WriteError {
    /// A quantity's magnitude is not finite.
    NonFinite { property: String, value: String },
    /// The underlying writer refused.
    Io(String),
}

impl std::fmt::Display for WriteError {
    fn fmt(&self, formatter: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            Self::NonFinite { property, value } => {
                write!(
                    formatter,
                    "`{property}` carries `{value}`, which is not finite"
                )
            }
            Self::Io(message) => write!(formatter, "{message}"),
        }
    }
}

impl std::error::Error for WriteError {}

impl From<quick_xml::Error> for WriteError {
    fn from(error: quick_xml::Error) -> Self {
        Self::Io(error.to_string())
    }
}

impl From<std::io::Error> for WriteError {
    fn from(error: std::io::Error) -> Self {
        Self::Io(error.to_string())
    }
}

/// Serialise a document as DEXPI XML.
///
/// # Errors
/// [`WriteError`] for a non-finite quantity or a writer failure.
pub fn write(document: &Document) -> Result<String, WriteError> {
    let mut writer = Writer::new_with_indent(Vec::new(), b' ', 4);
    writer.write_event(Event::Decl(BytesDecl::new("1.0", Some("UTF-8"), None)))?;
    write_model(&mut writer, &document.model)?;
    String::from_utf8(writer.into_inner())
        .map_err(|error| WriteError::Io(format!("the writer produced invalid UTF-8: {error}")))
}

fn write_model(writer: &mut Writer<Vec<u8>>, model: &Model) -> Result<(), WriteError> {
    let mut start = BytesStart::new("Model");
    start.push_attribute(("name", model.name.as_str()));
    start.push_attribute(("uri", model.uri.as_str()));
    writer.write_event(Event::Start(start))?;
    for import in &model.imports {
        write_import(writer, import)?;
    }
    for child in &model.children {
        write_child(writer, child)?;
    }
    writer.write_event(Event::End(BytesEnd::new("Model")))?;
    Ok(())
}

fn write_import(writer: &mut Writer<Vec<u8>>, import: &Import) -> Result<(), WriteError> {
    let mut start = BytesStart::new("Import");
    start.push_attribute(("prefix", import.prefix.as_str()));
    start.push_attribute(("source", import.source.as_str()));
    writer.write_event(Event::Empty(start))?;
    Ok(())
}

fn write_child(writer: &mut Writer<Vec<u8>>, child: &Child) -> Result<(), WriteError> {
    match child {
        Child::Object(object) => write_object(writer, object),
        Child::Data(data) => write_data(writer, data),
        Child::Components(components) => write_components(writer, components),
        Child::References(references) => write_references(writer, references),
    }
}

fn write_object(writer: &mut Writer<Vec<u8>>, object: &Object) -> Result<(), WriteError> {
    let mut start = BytesStart::new("Object");
    if let Some(id) = &object.id {
        start.push_attribute(("id", id.as_str()));
    }
    if let Some(name) = &object.name {
        start.push_attribute(("name", name.as_str()));
    }
    start.push_attribute(("type", object.ty.as_str()));
    if object.children.is_empty() {
        writer.write_event(Event::Empty(start))?;
        return Ok(());
    }
    writer.write_event(Event::Start(start))?;
    for child in &object.children {
        write_child(writer, child)?;
    }
    writer.write_event(Event::End(BytesEnd::new("Object")))?;
    Ok(())
}

fn write_components(
    writer: &mut Writer<Vec<u8>>,
    components: &Components,
) -> Result<(), WriteError> {
    let mut start = BytesStart::new("Components");
    start.push_attribute(("property", components.property.as_str()));
    writer.write_event(Event::Start(start))?;
    for child in &components.children {
        match child {
            ComponentChild::Object(object) => write_object(writer, object)?,
            ComponentChild::Reference { object } => {
                let mut reference = BytesStart::new("ObjectReference");
                reference.push_attribute(("object", object.as_str()));
                writer.write_event(Event::Empty(reference))?;
            }
        }
    }
    writer.write_event(Event::End(BytesEnd::new("Components")))?;
    Ok(())
}

fn write_references(
    writer: &mut Writer<Vec<u8>>,
    references: &References,
) -> Result<(), WriteError> {
    let joined = references.objects.join(" ");
    let mut start = BytesStart::new("References");
    start.push_attribute(("property", references.property.as_str()));
    start.push_attribute(("objects", joined.as_str()));
    writer.write_event(Event::Empty(start))?;
    Ok(())
}

fn write_data(writer: &mut Writer<Vec<u8>>, data: &Data) -> Result<(), WriteError> {
    let mut start = BytesStart::new("Data");
    start.push_attribute(("property", data.property.as_str()));
    writer.write_event(Event::Start(start))?;
    for value in &data.values {
        write_value(writer, value)?;
    }
    writer.write_event(Event::End(BytesEnd::new("Data")))?;
    Ok(())
}

fn write_value(writer: &mut Writer<Vec<u8>>, value: &DataValue) -> Result<(), WriteError> {
    match value {
        DataValue::String(text) => write_text_element(writer, "String", text)?,
        DataValue::Double(number) => {
            if !number.is_finite() {
                return Err(WriteError::NonFinite {
                    property: "Double".to_string(),
                    value: format!("{number}"),
                });
            }
            write_text_element(writer, "Double", &format_double(*number))?;
        }
        DataValue::Integer(number) => write_text_element(writer, "Integer", &number.to_string())?,
        DataValue::Boolean(flag) => {
            write_text_element(writer, "Boolean", if *flag { "true" } else { "false" })?
        }
        DataValue::DateTime(text) => write_text_element(writer, "DateTime", text)?,
        DataValue::Reference(target) => {
            let mut start = BytesStart::new("DataReference");
            start.push_attribute(("data", target.as_str()));
            writer.write_event(Event::Empty(start))?;
        }
        DataValue::Aggregated { ty, data } => {
            let mut start = BytesStart::new("AggregatedDataValue");
            start.push_attribute(("type", ty.as_str()));
            writer.write_event(Event::Start(start))?;
            for entry in data {
                write_data(writer, entry)?;
            }
            writer.write_event(Event::End(BytesEnd::new("AggregatedDataValue")))?;
        }
        DataValue::Undefined => {
            writer.write_event(Event::Empty(BytesStart::new("Undefined")))?;
        }
    }
    Ok(())
}

fn write_text_element(
    writer: &mut Writer<Vec<u8>>,
    name: &str,
    text: &str,
) -> Result<(), WriteError> {
    let start = BytesStart::new(name);
    writer.write_event(Event::Start(start))?;
    writer.write_event(Event::Text(BytesText::new(text)))?;
    writer.write_event(Event::End(BytesEnd::new(name)))?;
    Ok(())
}

/// An `f64` in the shortest form that reads back as the same value.
///
/// `xsd:double` accepts `1000.0`, `0.1` and `1e-7`, which is exactly the set this emits.
#[must_use]
pub fn format_double(value: f64) -> String {
    format!("{value:?}")
}
