//! The DEXPI 2.0 XML metamodel, typed.
//!
//! An instance exchange is a `Model` with `Import` elements and one tree of `Object`s. Every
//! object carries a `type` - a name reference such as `Process/Process.Pumping` - and three
//! kinds of child: `Components` (composition, holding objects), `Data` (a property with a typed
//! value) and `References` (a property pointing at object ids or names).
//!
//! **Child order is preserved because the schema is a sequence**, and because the shipped
//! fixture interleaves `Data` and `Components`: a model that regrouped them would serialise a
//! different document than it read.

/// One exchange document: the `Model` root.
#[derive(Debug, Clone, PartialEq)]
pub struct Document {
    pub model: Model,
}

/// The `Model` root: its identity, the information models it imports, and its content.
#[derive(Debug, Clone, PartialEq)]
pub struct Model {
    pub name: String,
    pub uri: String,
    pub imports: Vec<Import>,
    pub children: Vec<Child>,
}

/// One information-model import. The `source` is one of the official model URIs.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Import {
    pub prefix: String,
    pub source: String,
}

/// One child of a `Model` or of an `Object`.
#[derive(Debug, Clone, PartialEq)]
pub enum Child {
    Object(Object),
    Data(Data),
    Components(Components),
    References(References),
}

/// An object: its identity, its type, and its properties.
#[derive(Debug, Clone, PartialEq)]
pub struct Object {
    pub id: Option<String>,
    pub name: Option<String>,
    /// The DEXPI type, a name reference such as `Process/Process.Pumping`.
    pub ty: String,
    pub children: Vec<Child>,
}

/// A composition property: a named group of objects.
#[derive(Debug, Clone, PartialEq)]
pub struct Components {
    pub property: String,
    pub children: Vec<ComponentChild>,
}

/// What a `Components` group can hold.
#[derive(Debug, Clone, PartialEq)]
pub enum ComponentChild {
    Object(Object),
    /// An `ObjectReference`, by id or name.
    Reference {
        object: String,
    },
}

/// A data property and its value or values.
///
/// **A vector, not one value.** The schema allows several values under one `Data`; a model
/// that kept only the first would silently drop the rest.
#[derive(Debug, Clone, PartialEq)]
pub struct Data {
    pub property: String,
    pub values: Vec<DataValue>,
}

/// One value under a `Data` property.
#[derive(Debug, Clone, PartialEq)]
pub enum DataValue {
    String(String),
    Double(f64),
    Integer(i64),
    Boolean(bool),
    DateTime(String),
    /// A `DataReference`, by id or name - the unit references use this.
    Reference(String),
    /// An `AggregatedDataValue` with a type and nested properties.
    Aggregated {
        ty: String,
        data: Vec<Data>,
    },
    Undefined,
}

/// A reference property pointing at object ids or names.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct References {
    pub property: String,
    pub objects: Vec<String>,
}

impl Object {
    /// One data value, by property name.
    #[must_use]
    pub fn data_value(&self, property: &str) -> Option<&DataValue> {
        self.children.iter().find_map(|child| match child {
            Child::Data(data) if data.property == property => data.values.first(),
            _ => None,
        })
    }

    /// Every value of one data property, in document order.
    ///
    /// **A vector, not one value**: a `PolyLine`'s `Points` is one `Data` carrying
    /// repeated aggregated values, and a reader that took the first would draw one point.
    pub fn data_values<'a>(
        &'a self,
        property: &'a str,
    ) -> impl Iterator<Item = &'a DataValue> + 'a {
        self.children
            .iter()
            .filter_map(move |child| match child {
                Child::Data(data) if data.property == property => Some(data.values.iter()),
                _ => None,
            })
            .flatten()
    }

    /// One string property, by name.
    #[must_use]
    pub fn data_string(&self, property: &str) -> Option<&str> {
        match self.data_value(property)? {
            DataValue::String(text) => Some(text),
            _ => None,
        }
    }

    /// One composition group, by property name.
    #[must_use]
    pub fn components(&self, property: &str) -> Option<&Components> {
        self.children.iter().find_map(|child| match child {
            Child::Components(components) if components.property == property => Some(components),
            _ => None,
        })
    }

    /// One reference property, by name.
    #[must_use]
    pub fn references(&self, property: &str) -> Option<&References> {
        self.children.iter().find_map(|child| match child {
            Child::References(references) if references.property == property => Some(references),
            _ => None,
        })
    }

    /// Every object held directly under one composition property.
    pub fn objects(&self, property: &str) -> impl Iterator<Item = &Object> {
        self.components(property)
            .into_iter()
            .flat_map(|components| components.children.iter())
            .filter_map(|child| match child {
                ComponentChild::Object(object) => Some(object),
                ComponentChild::Reference { .. } => None,
            })
    }

    /// One object held under a composition property, by its id.
    #[must_use]
    pub fn object(&self, property: &str, id: &str) -> Option<&Object> {
        self.objects(property)
            .find(|object| object.id.as_deref() == Some(id))
    }
}

impl DataValue {
    /// The value as an `f64`: a double, or an integer widened.
    #[must_use]
    pub fn as_f64(&self) -> Option<f64> {
        match self {
            Self::Double(value) => Some(*value),
            #[allow(clippy::cast_precision_loss)]
            Self::Integer(value) => Some(*value as f64),
            _ => None,
        }
    }

    /// The unit reference of a physical quantity, if this value is one.
    #[must_use]
    pub fn quantity_unit(&self) -> Option<&str> {
        let Self::Aggregated { ty, data } = self else {
            return None;
        };
        if ty != "Core/PhysicalQuantities.PhysicalQuantity" {
            return None;
        }
        let unit = data.iter().find(|entry| entry.property == "Unit")?;
        match unit.values.first()? {
            DataValue::Reference(unit) => Some(unit.as_str()),
            _ => None,
        }
    }

    /// The magnitude of a physical quantity, if this value is one.
    #[must_use]
    pub fn quantity_magnitude(&self) -> Option<f64> {
        let Self::Aggregated { ty, data } = self else {
            return None;
        };
        if ty != "Core/PhysicalQuantities.PhysicalQuantity" {
            return None;
        }
        data.iter()
            .find(|entry| entry.property == "Value")
            .and_then(|entry| entry.values.first())
            .and_then(DataValue::as_f64)
    }
}
