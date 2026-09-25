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
///
/// **The convergence parameters are `Recycle`'s own, and every one is optional.** A
/// declaration that states none takes the class's defaults, which is what keeps
/// `specs/flowsheets/demo.toml` validating unchanged — and a declaration that states one
/// carries it rather than inventing a number. See [`crate::recycle`] for what each means.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Recycle {
    pub stream: String,
    pub from: String,
    pub to: String,
    /// `Recycle.flowTolerance`, defaulting to `1e-2`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub flow_tolerance: Option<f64>,
    /// `Recycle.compositionTolerance`, defaulting to `1e-2`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub composition_tolerance: Option<f64>,
    /// `Recycle.temperatureTolerance`, defaulting to `1e-2`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub temperature_tolerance: Option<f64>,
    /// `Recycle.pressureTolerance`, defaulting to `1e-2`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub pressure_tolerance: Option<f64>,
    /// `Recycle.maxIterations`, defaulting to `10`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub max_iterations: Option<u32>,
    /// `Recycle.minimumFlow`, kg/hr, defaulting to `1e-20`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub minimum_flow: Option<f64>,
    /// `Recycle`'s `AccelerationMethod`: `direct_substitution`, `wegstein` or `broyden`.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub acceleration_method: Option<String>,
}

impl Recycle {
    /// A tear with the class's defaults, for a caller that states none.
    ///
    /// Every convergence parameter starts unstated, which is what a declaration that carries
    /// only `stream`, `from` and `to` means - and what keeps a runtime `Recycle` and a declared
    /// one the same object.
    #[must_use]
    pub fn new(stream: &str, from: &str, to: &str) -> Self {
        Self {
            stream: stream.to_string(),
            from: from.to_string(),
            to: to.to_string(),
            flow_tolerance: None,
            composition_tolerance: None,
            temperature_tolerance: None,
            pressure_tolerance: None,
            max_iterations: None,
            minimum_flow: None,
            acceleration_method: None,
        }
    }

    /// The tear's convergence settings, with the class's defaults where the declaration is
    /// silent.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] on an `acceleration_method` that is none of the
    /// three names, which a declaration can carry and the class's own enum cannot.
    pub fn settings(&self) -> azoth_core::Result<crate::recycle::RecycleSettings> {
        let defaults = crate::recycle::RecycleSettings::default();
        Ok(crate::recycle::RecycleSettings {
            flow_tolerance: self.flow_tolerance.unwrap_or(defaults.flow_tolerance),
            composition_tolerance: self
                .composition_tolerance
                .unwrap_or(defaults.composition_tolerance),
            temperature_tolerance: self
                .temperature_tolerance
                .unwrap_or(defaults.temperature_tolerance),
            pressure_tolerance: self
                .pressure_tolerance
                .unwrap_or(defaults.pressure_tolerance),
            max_iterations: self.max_iterations.unwrap_or(defaults.max_iterations),
            minimum_flow_kg_per_hr: self.minimum_flow.unwrap_or(defaults.minimum_flow_kg_per_hr),
            acceleration: match self.acceleration_method.as_deref() {
                None => defaults.acceleration,
                Some(name) => crate::recycle::Acceleration::named(name)?,
            },
        })
    }
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
