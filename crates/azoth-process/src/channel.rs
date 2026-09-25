//! Channel types: the field record a stream carries, and the port that declares it.
//!
//! `docs/src/calculus/process.md` fixes the shape: a channel type is a finite
//! record of field names each at a dimension, plus a polarity. This module is the
//! concrete form of that — the schema a unit operation declares its ports in, and
//! the thing a flowsheet's connections are checked against.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// A channel type: the finite record of named fields a stream carries.
///
/// Field names are the symbols in the balance — `n` for molar flow, `z` for
/// composition, `P`, `T`, `h`. Each carries a dimension (from the vocabulary) and
/// a shape. Nothing else crosses a channel.
pub type ChannelType = BTreeMap<String, FieldType>;

/// Which way a channel flows through a unit operation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Direction {
    In,
    Out,
}

/// How many streams a port carries. `One` is the default; `Many` is a mixer's
/// inlets or a splitter's outlets, declared per port rather than as a number.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Multiplicity {
    #[default]
    One,
    Many,
}

/// Whether a field is a single quantity or one entry per component.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize, Default)]
#[serde(rename_all = "snake_case")]
pub enum Shape {
    #[default]
    Scalar,
    Vector,
}

/// One field of a channel: a named quantity at a dimension.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct FieldType {
    /// A dimension id from `specs/vocabulary/vocabulary.toml`, e.g. `pressure`.
    pub dimension: String,
    #[serde(default)]
    pub shape: Shape,
}

/// A port: the declared interface of one side of a unit operation.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Port {
    pub name: String,
    pub direction: Direction,
    #[serde(default)]
    pub multiplicity: Multiplicity,
    /// The field record a stream on this port carries.
    pub fields: ChannelType,
}
