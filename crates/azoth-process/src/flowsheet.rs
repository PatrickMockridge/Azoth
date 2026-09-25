//! A flowsheet: instances of unit operations wired together by named streams.
//!
//! The model follows the calculus. A stream is a name; it is produced exactly
//! once and consumed exactly once. A connection joins an outlet (or a feed, the
//! environment's outlet) to an inlet (or a product, the environment's inlet). A
//! recycle is a connection that closes a loop, with its stream named as the tear.

use std::collections::BTreeMap;

use azoth_core::units::{kelvins, pascals};
use serde::{Deserialize, Serialize};

use crate::stream::Stream;

/// One boundary inlet: the stream the environment produces, and the record its user specifies.
///
/// **This is the input half of a self-contained simulation.** The fields are the port record's
/// own names — `n`, `z`, `P`, `T` out of the five every palette port declares — so the unit is
/// the field's rather than written here, exactly as an instance parameter's bare number takes
/// its unit from the palette entry that declares it.
///
/// **`h` is the one field of the record a user does not write, and that is the whole point of
/// the split.** It is a state function of `(T, P, z)`, so [`Stream::from_pt`] calculates it -
/// which is why `process.pump`'s own spec has an inlet take four fields and an outlet five.
/// `deny_unknown_fields` is what makes a written `h` an error rather than a key serde reads
/// past, because a silently ignored field is the one failure a reader cannot see.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Input {
    /// The stream's name — what a connection's `from`, and a recycle's `from`, address.
    pub name: String,
    /// The fluid's substances, by name, resolved against the component databank.
    pub components: Vec<String>,
    /// Molar flow, `mol/s`.
    pub n: f64,
    /// Composition, one entry per component.
    pub z: Vec<f64>,
    /// Pressure, `Pa`, under the declaration's own name.
    #[serde(rename = "P")]
    pub p: f64,
    /// Temperature, `K`, under the declaration's own name.
    #[serde(rename = "T")]
    pub t: f64,
}

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

impl Flowsheet {
    /// Read a flowsheet from TOML.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] on a document this schema cannot read, naming what
    /// the parser said - the same refusal [`crate::check::validate`] gives a *valid* document it
    /// cannot accept.
    pub fn from_toml(text: &str) -> azoth_core::Result<Self> {
        toml::from_str(text)
            .map_err(|error| azoth_core::AzothError::invalid_input("flowsheet", error.to_string()))
    }

    /// Write a flowsheet back to TOML.
    ///
    /// **This is the shadow of `Azoth.Rho`'s round trip.** `*@P ≅ P` is the claim that printing a
    /// parsed process and reading it again is the identity, and the code's version of it is that
    /// `from_toml(to_toml(f))` is `f`. The test beside this holds both directions: the value comes
    /// back unchanged, and writing it a second time is byte-identical to the first - which is the
    /// half that catches a writer depending on something it does not preserve.
    ///
    /// **The output is not the input's bytes.** TOML has no comments to round-trip and this schema
    /// writes its fields in the struct's order, so a hand-written file comes back tidied. What is
    /// claimed is the *value*, and the claim is checkable because the value is what the checker,
    /// the executor and the registry all read.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] if the value cannot be written, which TOML's own
    /// rule about ordering a table before a scalar makes possible in principle - a `Vec` field
    /// placed before a plain one in the struct would trip it - and which the test over
    /// `specs/flowsheets/` would catch rather than a caller.
    pub fn to_toml(&self) -> azoth_core::Result<String> {
        toml::to_string(self)
            .map_err(|error| azoth_core::AzothError::invalid_input("flowsheet", error.to_string()))
    }
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

/// A flowsheet: a self-contained simulation.
///
/// **The document declares both halves of the run.** Its `inputs` are what the user specifies and
/// its `products` are what the run calculates, and the field order below is not cosmetic: TOML has
/// no way to reopen a key after an array of tables, so `products` — a plain array — is written
/// *before* `[[inputs]]`, and `toml::to_string` refuses the other order outright.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Flowsheet {
    pub id: String,
    pub name: String,
    /// Boundary outlets — streams the environment consumes. **Names only**: an output is
    /// calculated, so it has nothing to state.
    #[serde(default)]
    pub products: Vec<String>,
    /// Boundary inlets — the streams the environment produces, with the record that makes them.
    #[serde(default)]
    pub inputs: Vec<Input>,
    #[serde(default)]
    pub instances: Vec<Instance>,
    #[serde(default)]
    pub connections: Vec<Connection>,
    #[serde(default)]
    pub recycles: Vec<Recycle>,
}

impl Flowsheet {
    /// The input names, in declaration order — what a bare `from` addresses.
    pub fn input_names(&self) -> impl Iterator<Item = &str> {
        self.inputs.iter().map(|input| input.name.as_str())
    }

    /// The boundary inlets as streams, each built from the record the document declares.
    ///
    /// **The fluid is resolved here**, so an input naming a substance the databank does not carry
    /// is refused by name at this call and not silently dropped. A composition whose length does
    /// not match its component list is refused too - the checker states the same rule as a
    /// diagnostic, and this is what a caller who skipped `validate` meets.
    ///
    /// # Errors
    /// [`azoth_core::AzothError::InvalidInput`] for either case above, or wherever the flash that
    /// fixes the inlet enthalpy refuses its state.
    pub fn input_streams(&self) -> azoth_core::Result<BTreeMap<String, Stream>> {
        let mut streams = BTreeMap::new();
        for input in &self.inputs {
            if input.z.len() != input.components.len() {
                return Err(azoth_core::AzothError::invalid_input(
                    format!("inputs.{}", input.name),
                    format!(
                        "`z` has {} entries and `components` names {} substances; a composition \
                         is one mole fraction per substance",
                        input.z.len(),
                        input.components.len()
                    ),
                ));
            }
            let stream = Stream::from_pt(
                input.components.clone(),
                input.z.clone(),
                input.n,
                pascals(input.p),
                kelvins(input.t),
            )?;
            streams.insert(input.name.clone(), stream);
        }
        Ok(streams)
    }
}

/// An input of one substance, for the tests that build a `Flowsheet` in code.
///
/// **Not a public constructor, deliberately.** The checker's and the order's own tests build a
/// graph where the input's record is not the thing under test, and a caller-visible shorthand
/// would be a second way to state a boundary. The record it does carry is a *valid* one, so a
/// test that ignores it does not fail the checker's own rule about it.
#[cfg(test)]
pub fn named_input(name: &str) -> Input {
    Input {
        name: name.to_string(),
        components: vec!["methane".to_string()],
        n: 1.0,
        z: vec![1.0],
        p: 1.0e5,
        t: 300.0,
    }
}
