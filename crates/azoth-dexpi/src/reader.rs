//! DEXPI XML into the typed metamodel.
//!
//! **Fail-closed, and no DTD.** A document type declaration is refused, so there is no entity
//! expansion to exploit and no external subset to fetch; document size and nesting depth are
//! capped; an element the metamodel does not know is an error rather than a subtree walked past.
//! A reader that silently skipped what it did not understand would report a complete exchange
//! over a partial one, which is the failure this crate exists to avoid.
//!
//! The document is first parsed into a small untyped tree and then converted. Two passes rather
//! than one recursive descent because the metamodel's own structure depends on element names,
//! and a tree makes the limits and the end-tag checking one place instead of one per rule.

use quick_xml::Reader;
use quick_xml::escape::unescape;
use quick_xml::events::Event;

use crate::diagnostic::DexpiDiagnostic;
use crate::model::{
    Child, ComponentChild, Components, Data, DataValue, Document, Import, Model, Object, References,
};

/// The largest document this reader accepts, in bytes.
pub const MAX_BYTES: usize = 64 * 1024 * 1024;

/// The deepest element nesting this reader accepts.
pub const MAX_DEPTH: usize = 128;

/// A read document and anything worth saying about it.
#[derive(Debug, Clone, PartialEq)]
pub struct ReadReport {
    pub document: Document,
    pub diagnostics: Vec<DexpiDiagnostic>,
}

/// Read a DEXPI XML instance document.
///
/// # Errors
/// [`DexpiDiagnostic`] naming the first thing that made the document unreadable: a DTD, a
/// document over [`MAX_BYTES`], nesting over [`MAX_DEPTH`], a malformed document, or an element
/// the metamodel does not know.
pub fn read(xml: &str) -> Result<ReadReport, DexpiDiagnostic> {
    if xml.len() > MAX_BYTES {
        return Err(DexpiDiagnostic::error(
            "document_too_large",
            "",
            format!(
                "the document is {} bytes; this reader accepts {MAX_BYTES}",
                xml.len()
            ),
        ));
    }
    let root = parse_tree(xml)?;
    if root.name != "Model" {
        return Err(DexpiDiagnostic::error(
            "unsupported_root",
            root.name,
            "a DEXPI 2.0 instance exchange has a `Model` root; a Proteus `PlantModel` is a \
             different serialisation and is not read by this profile",
        ));
    }

    let mut imports = Vec::new();
    let mut children = Vec::new();
    for child in root.elements() {
        match child.name.as_str() {
            "Import" => imports.push(Import {
                prefix: child.required("prefix")?,
                source: child.required("source")?,
            }),
            "Object" => children.push(Child::Object(object(child)?)),
            "Data" => children.push(Child::Data(data(child)?)),
            other => {
                return Err(DexpiDiagnostic::error(
                    "unsupported_element",
                    other,
                    "a DEXPI 2.0 instance exchange holds `Import`, `Object` and `Data` at its \
                     root; a model-definition element belongs in the information model, not in \
                     an exchange",
                ));
            }
        }
    }

    Ok(ReadReport {
        document: Document {
            model: Model {
                name: root.attribute("name").unwrap_or_default(),
                uri: root.attribute("uri").unwrap_or_default(),
                imports,
                children,
            },
        },
        diagnostics: Vec::new(),
    })
}

fn object(element: &RawElement) -> Result<Object, DexpiDiagnostic> {
    let mut children = Vec::new();
    for child in element.elements() {
        match child.name.as_str() {
            "Components" => children.push(Child::Components(components(child)?)),
            "Data" => children.push(Child::Data(data(child)?)),
            "References" => children.push(Child::References(References {
                property: child.required("property")?,
                objects: child
                    .required("objects")?
                    .split_whitespace()
                    .map(ToString::to_string)
                    .collect(),
            })),
            other => {
                return Err(DexpiDiagnostic::error(
                    "unsupported_element",
                    other,
                    "an `Object` holds `Components`, `Data` and `References`",
                ));
            }
        }
    }
    Ok(Object {
        id: element.attribute("id"),
        name: element.attribute("name"),
        ty: element.required("type")?,
        children,
    })
}

fn components(element: &RawElement) -> Result<Components, DexpiDiagnostic> {
    let mut children = Vec::new();
    for child in element.elements() {
        match child.name.as_str() {
            "Object" => children.push(ComponentChild::Object(object(child)?)),
            "ObjectReference" => children.push(ComponentChild::Reference {
                object: child.required("object")?,
            }),
            other => {
                return Err(DexpiDiagnostic::error(
                    "unsupported_element",
                    other,
                    "`Components` holds `Object` and `ObjectReference`",
                ));
            }
        }
    }
    Ok(Components {
        property: element.required("property")?,
        children,
    })
}

fn data(element: &RawElement) -> Result<Data, DexpiDiagnostic> {
    let mut values = Vec::new();
    for child in element.elements() {
        values.push(value(child)?);
    }
    Ok(Data {
        property: element.required("property")?,
        values,
    })
}

fn value(element: &RawElement) -> Result<DataValue, DexpiDiagnostic> {
    let text = element.text();
    Ok(match element.name.as_str() {
        "String" => DataValue::String(text),
        "Double" => DataValue::Double(parse_f64(element, &text)?),
        "Integer" => DataValue::Integer(text.trim().parse().map_err(|_| {
            DexpiDiagnostic::error(
                "bad_number",
                "Integer",
                format!("`{text}` is not an integer"),
            )
        })?),
        "Boolean" => match text.trim() {
            "true" | "1" => DataValue::Boolean(true),
            "false" | "0" => DataValue::Boolean(false),
            other => {
                return Err(DexpiDiagnostic::error(
                    "bad_boolean",
                    "Boolean",
                    format!("`{other}` is not a boolean"),
                ));
            }
        },
        "DateTime" => DataValue::DateTime(text),
        "DataReference" => DataValue::Reference(element.required("data")?),
        "Undefined" => DataValue::Undefined,
        "AggregatedDataValue" => {
            let mut entries = Vec::new();
            for child in element.elements() {
                if child.name != "Data" {
                    return Err(DexpiDiagnostic::error(
                        "unsupported_element",
                        child.name.clone(),
                        "`AggregatedDataValue` holds `Data`",
                    ));
                }
                entries.push(data(child)?);
            }
            DataValue::Aggregated {
                ty: element.required("type")?,
                data: entries,
            }
        }
        other => {
            return Err(DexpiDiagnostic::error(
                "unsupported_value",
                other,
                "a `Data` value is one of `String`, `Double`, `Integer`, `Boolean`, `DateTime`, \
                 `DataReference`, `AggregatedDataValue` and `Undefined`",
            ));
        }
    })
}

fn parse_f64(element: &RawElement, text: &str) -> Result<f64, DexpiDiagnostic> {
    text.trim().parse::<f64>().map_err(|_| {
        DexpiDiagnostic::error(
            "bad_number",
            element.name.clone(),
            format!("`{text}` is not a double"),
        )
    })
}

// --- the untyped tree -------------------------------------------------------------------

pub(crate) struct RawElement {
    pub(crate) name: String,
    attributes: Vec<(String, String)>,
    pub(crate) children: Vec<RawNode>,
}

pub(crate) enum RawNode {
    Element(RawElement),
    Text(String),
}

impl RawElement {
    pub(crate) fn attribute(&self, name: &str) -> Option<String> {
        self.attributes
            .iter()
            .find(|(key, _)| key == name)
            .map(|(_, value)| value.clone())
    }

    pub(crate) fn required(&self, name: &str) -> Result<String, DexpiDiagnostic> {
        self.attribute(name).ok_or_else(|| {
            DexpiDiagnostic::error(
                "missing_attribute",
                self.name.clone(),
                format!("`{}` requires its `{name}` attribute", self.name),
            )
        })
    }

    pub(crate) fn elements(&self) -> impl Iterator<Item = &RawElement> {
        self.children.iter().filter_map(|child| match child {
            RawNode::Element(element) => Some(element),
            RawNode::Text(_) => None,
        })
    }

    pub(crate) fn text(&self) -> String {
        let mut text = String::new();
        for child in &self.children {
            if let RawNode::Text(chunk) = child {
                text.push_str(chunk);
            }
        }
        text
    }
}

/// The five predefined XML entities and a numeric character reference.
fn resolve_reference(name: &str) -> Option<String> {
    let predefined = match name {
        "amp" => Some('&'),
        "lt" => Some('<'),
        "gt" => Some('>'),
        "apos" => Some('\''),
        "quot" => Some('"'),
        _ => None,
    };
    if let Some(character) = predefined {
        return Some(character.to_string());
    }
    let digits = name.strip_prefix('#')?;
    let value = match digits
        .strip_prefix('x')
        .or_else(|| digits.strip_prefix('X'))
    {
        Some(hexadecimal) => u32::from_str_radix(hexadecimal, 16).ok()?,
        None => digits.parse().ok()?,
    };
    char::from_u32(value).map(String::from)
}

pub(crate) fn parse_tree(xml: &str) -> Result<RawElement, DexpiDiagnostic> {
    let mut reader = Reader::from_str(xml);
    reader.config_mut().check_end_names = true;
    let mut stack: Vec<RawElement> = Vec::new();
    let mut root: Option<RawElement> = None;

    loop {
        let event = reader.read_event().map_err(|error| {
            DexpiDiagnostic::error(
                "malformed_document",
                "",
                format!("at byte {}: {error}", reader.buffer_position()),
            )
        })?;
        match event {
            Event::DocType(_) => {
                return Err(DexpiDiagnostic::error(
                    "doctype",
                    "",
                    "a document type declaration is refused: this exchange has no entities to \
                     expand and no external subset to fetch",
                ));
            }
            Event::Start(start) => {
                if stack.len() >= MAX_DEPTH {
                    return Err(DexpiDiagnostic::error(
                        "document_too_deep",
                        "",
                        format!("element nesting exceeds the {MAX_DEPTH} this reader accepts"),
                    ));
                }
                stack.push(element(&start)?);
            }
            Event::Empty(start) => {
                let empty = element(&start)?;
                attach(&mut stack, &mut root, empty)?;
            }
            Event::End(_) => {
                let closed = stack.pop().ok_or_else(|| {
                    DexpiDiagnostic::error("malformed_document", "", "an end tag with nothing open")
                })?;
                attach(&mut stack, &mut root, closed)?;
            }
            Event::Text(text) => {
                let decoded = text
                    .decode()
                    .map_err(|error| DexpiDiagnostic::error("bad_text", "", error.to_string()))?;
                let unescaped = unescape(&decoded)
                    .map_err(|error| DexpiDiagnostic::error("bad_text", "", error.to_string()))?;
                if unescaped.trim().is_empty() {
                    continue;
                }
                let Some(parent) = stack.last_mut() else {
                    return Err(DexpiDiagnostic::error(
                        "malformed_document",
                        "",
                        "text outside the root element",
                    ));
                };
                parent.children.push(RawNode::Text(unescaped.into_owned()));
            }
            Event::GeneralRef(reference) => {
                let name = reference
                    .decode()
                    .map_err(|error| DexpiDiagnostic::error("bad_text", "", error.to_string()))?;
                let resolved = resolve_reference(&name).ok_or_else(|| {
                    DexpiDiagnostic::error(
                        "unknown_entity",
                        name.to_string(),
                        "only the five predefined entities and character references resolve; \
                         this exchange expands no others because it declares no DTD",
                    )
                })?;
                let Some(parent) = stack.last_mut() else {
                    return Err(DexpiDiagnostic::error(
                        "malformed_document",
                        "",
                        "a reference outside the root element",
                    ));
                };
                parent.children.push(RawNode::Text(resolved));
            }
            Event::Comment(_) | Event::Decl(_) | Event::PI(_) => {}
            Event::Eof => break,
            other => {
                return Err(DexpiDiagnostic::error(
                    "unsupported_event",
                    "",
                    format!("{other:?}"),
                ));
            }
        }
    }

    if let Some(open) = stack.last() {
        return Err(DexpiDiagnostic::error(
            "unterminated_element",
            open.name.clone(),
            "the document ended before this element closed",
        ));
    }
    root.ok_or_else(|| {
        DexpiDiagnostic::error("empty_document", "", "the document has no root element")
    })
}

fn attach(
    stack: &mut [RawElement],
    root: &mut Option<RawElement>,
    completed: RawElement,
) -> Result<(), DexpiDiagnostic> {
    match stack.last_mut() {
        Some(parent) => parent.children.push(RawNode::Element(completed)),
        None => {
            if root.is_some() {
                return Err(DexpiDiagnostic::error(
                    "malformed_document",
                    "",
                    "the document has more than one root element",
                ));
            }
            *root = Some(completed);
        }
    }
    Ok(())
}

fn element(start: &quick_xml::events::BytesStart<'_>) -> Result<RawElement, DexpiDiagnostic> {
    let name = String::from_utf8_lossy(start.name().as_ref()).into_owned();
    let mut attributes = Vec::new();
    for attribute in start.attributes() {
        let attribute = attribute.map_err(|error| {
            DexpiDiagnostic::error("bad_attribute", name.clone(), error.to_string())
        })?;
        let key = String::from_utf8_lossy(attribute.key.as_ref()).into_owned();
        let value = attribute
            .unescape_value()
            .map_err(|error| {
                DexpiDiagnostic::error("bad_attribute", key.clone(), error.to_string())
            })?
            .into_owned();
        attributes.push((key, value));
    }
    Ok(RawElement {
        name,
        attributes,
        children: Vec::new(),
    })
}
