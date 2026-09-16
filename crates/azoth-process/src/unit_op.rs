//! A unit-operation spec: the palette entry.
//!
//! One file under `specs/unit_ops/` declares one unit operation: its id, name,
//! where its port shape comes from, the parameters a user fills in, and its ports.
//! The declaration fixes the *interface* — what crosses each channel and what the
//! unit op takes as form input — and is what `check.rs` holds a flowsheet to. The
//! arithmetic is a later tranche; a spec here is a declaration, not a kernel.

use std::collections::BTreeMap;

use serde::{Deserialize, Serialize};

use crate::channel::Port;

/// Where a unit operation's port shape comes from — its NeqSim class, named so a
/// reader can find the port the declaration was read from.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Source {
    pub standard: String,
}

/// One parameter a user types on the unit operation's form.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Param {
    /// The parameter's canonical unit, for a quantity parameter; omitted for an
    /// enum or boolean. Drawn from the same vocabulary as a calc input's unit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    pub description: String,
}

/// A unit-operation spec.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct UnitOpSpec {
    pub id: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<Source>,
    #[serde(default)]
    pub parameters: BTreeMap<String, Param>,
    pub ports: Vec<Port>,
}
