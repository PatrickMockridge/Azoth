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
///
/// **Unknown keys are refused rather than ignored.** TOML scopes a bare key to the most
/// recent table header, so a `notes =` written one line too low lands in `[source]` — and
/// serde's default is to drop it without a word, which is an assertion nothing reads. This
/// was found by writing exactly that line into `unit_ops.ejector`.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Source {
    pub standard: String,
}

/// One parameter a user types on the unit operation's form.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Param {
    /// The parameter's canonical unit, for a quantity parameter; omitted for an
    /// enum or boolean. Drawn from the same vocabulary as a calc input's unit.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    pub description: String,
}

/// A unit-operation spec.
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct UnitOpSpec {
    pub id: String,
    pub name: String,
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub source: Option<Source>,
    #[serde(default)]
    pub parameters: BTreeMap<String, Param>,
    pub ports: Vec<Port>,
    /// What the declaration cannot say by itself, where there is such a thing.
    ///
    /// **The palette entry has no `assumptions` block the way a calc spec does**, and a
    /// unit operation that is *not ported* needs somewhere to say so: the calculus's rule
    /// is that a thing is not "out of scope", it is not ported with the class that would
    /// close it named. `unit_ops.simple_absorber` is the one today — its NeqSim class is a
    /// fixed-point loop over MDEA/CO₂ loading rather than the stage-wise absorber its name
    /// and its ports describe, and the chemistry that would close it is unported.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub notes: Option<String>,
}
