//! A flowsheet: instances of unit operations wired together by named streams.
//!
//! The model follows the calculus. A stream is a name; it is produced exactly
//! once and consumed exactly once. A connection joins an outlet (or a feed, the
//! environment's outlet) to an inlet (or a product, the environment's inlet). A
//! recycle is a connection that closes a loop, with its stream named as the tear.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

/// One instance of a unit operation, with the parameter values a user typed on its
/// form.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Instance {
    pub id: String,
    /// The id of a palette spec, e.g. `unit_ops.pump`.
    pub unit: String,
    #[serde(default)]
    pub parameters: BTreeMap<String, toml::Value>,
}

/// A connection: a stream produced at `from` and consumed at `to`.
///
/// `from` is a feed name or `instance.port` (an outlet); `to` is a product name
/// or `instance.port` (an inlet). A bare name in `from` position is a feed, in
/// `to` position a product; the dot separates an instance from its port.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Connection {
    pub from: String,
    pub to: String,
}

/// A torn loop: a connection whose stream closes a cycle and is the tear.
///
/// `stream` is the recycle's name — the fresh name the calculus binds with
/// restriction, used once where it is written and once where it is read.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Recycle {
    pub stream: String,
    pub from: String,
    pub to: String,
}

/// A flowsheet.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Flowsheet {
    pub id: String,
    pub name: String,
    /// Boundary inlets — stream names the environment produces.
    #[serde(default)]
    pub feeds: Vec<String>,
    /// Boundary outlets — stream names the environment consumes.
    #[serde(default)]
    pub products: Vec<String>,
    #[serde(default)]
    pub instances: Vec<Instance>,
    #[serde(default)]
    pub connections: Vec<Connection>,
    #[serde(default)]
    pub recycles: Vec<Recycle>,
}
