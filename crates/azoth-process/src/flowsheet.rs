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
/// A connection endpoint, split into the three things it can say.
///
/// **The grammar is `name` or `instance.port` or `instance.port[i]`**, and it lives here rather
/// than in the three modules that read it - the checker, the execution order and the session each
/// had their own `split_once('.')`, which is three places for a grammar to disagree with itself.
///
/// **The index is what a `many` outlet needs and nothing else does.** `unit_ops.splitter` declares
/// its outlets as one port named `products`, so a two-way split returns two streams under one
/// name; `split1.products[0]` is the first of them. An index on an inlet is not a thing a
/// connection has - a `many` *inlet* takes several connections to the same port name, which is
/// how a mixer is wired - so a consumer endpoint never carries one, and the checker refuses one
/// that does.
///
/// A bare name - no dot - is a boundary stream, and returns `None` for the instance.
#[must_use]
pub fn split_endpoint(endpoint: &str) -> Option<(&str, &str, Option<usize>)> {
    let (instance, rest) = endpoint.split_once('.')?;
    let Some((port, index)) = rest.split_once('[') else {
        return Some((instance, rest, None));
    };
    let index = index.strip_suffix(']')?.parse().ok()?;
    Some((instance, port, Some(index)))
}

/// The path a stream is bound to, and the one a session's `value` answers to.
///
/// `instance.port` for a `one` outlet, `instance.port[i]` for the i-th stream of a `many` one -
/// the same grammar [`split_endpoint`] reads, so a path a run produced is a path a connection
/// could have written.
#[must_use]
pub fn stream_path(instance: &str, port: &str, index: Option<usize>) -> String {
    match index {
        Some(index) => format!("{instance}.{port}[{index}]"),
        None => format!("{instance}.{port}"),
    }
}

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
#[serde(deny_unknown_fields)]
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
#[serde(deny_unknown_fields)]
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
#[serde(deny_unknown_fields)]
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

    /// Whether this declaration's acceleration is one the port carries, and why not if it is not.
    ///
    /// **`broyden` is declared and refused, and the reason is measured.** Its arithmetic *is*
    /// transcribed - `BroydenAccelerator` is held to `captures/process_acceleration.tsv` call by
    /// call - but what a `Recycle` does with the answer is not reproducible: `applyStreamValues`
    /// writes the accelerated fractions through `Component.setx` on **both** phases, and a
    /// two-phase system does not read back what it was written (the capture records `[1.0, 0.0]`
    /// written and `[1.0, 0.4]` read). A port whose stream carries one composition has nothing to
    /// model that with, and the difference is not small: on the condensing graph at
    /// `flowTolerance = 1e-8`, NeqSim runs the full hundred passes and does **not** converge, with
    /// `sep1.liquid` at `8.6243` mol/s (`captures/process_flowsheet_accelerated.tsv`), where
    /// applying the transcribed step to the composition gives a loop that reports converged after
    /// eighteen passes at `0.0833`. A wrong answer that reads as converged is the one outcome this
    /// port may not produce, so the declaration is refused rather than run.
    ///
    /// **Wegstein is carried**, and its step is bounded (`q` is clamped to `[-5, 0]`) where
    /// Broyden's is not - so the same write-back gap moves it by `1.9e-11` relative on that same
    /// loop, which the test over the capture holds.
    #[must_use]
    pub fn unsupported_acceleration(&self) -> Option<azoth_core::AzothError> {
        if self.acceleration_method.as_deref() == Some("broyden") {
            return Some(azoth_core::AzothError::invalid_input(
                "acceleration_method",
                "`broyden` is declared and refused: `BroydenAccelerator` is transcribed and held \
                 to the acceleration capture, but `Recycle.applyStreamValues` writes the \
                 accelerated composition through `Component.setx` on both phases, and a two-phase \
                 system does not read back what it was written - so the class's *effect* is not \
                 reproducible from a single-composition stream. Measured on the condensing graph, \
                 NeqSim does not converge in a hundred passes at 8.6243 mol/s where the \
                 transcribed step converges in eighteen at 0.0833",
            ));
        }
        None
    }
}

/// A flowsheet: a self-contained simulation.
///
/// **The document declares both halves of the run.** Its `inputs` are what the user specifies and
/// its `products` are what the run calculates, and the field order below is not cosmetic: TOML has
/// no way to reopen a key after an array of tables, so `products` — a plain array — is written
/// *before* `[[inputs]]`, and `toml::to_string` refuses the other order outright.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
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
    /// Where an editor draws each node, where one has been placed.
    ///
    /// **Last, and a table rather than a field on `Instance`.** `products` is names only — an
    /// output is calculated and has nothing to state — so a per-instance position could not cover
    /// the boundary, and one table covers all three kinds uniformly. Written last because TOML has
    /// no way to reopen a key after an array of tables, so a table declared after them is the only
    /// order the writer can emit; `specs/flowsheets/` holds it to that.
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub layout: Option<Layout>,
}

/// Where the editor draws each node, by name.
///
/// **A view and not physics.** Nothing in the run reads it: the checker ignores it, the executor
/// never sees it, and a flowsheet written by hand carries none — the projection derives a
/// deterministic layout from the connection graph and uses these only where they exist. It lives
/// in the document so that one artifact is the whole editor state, which is what makes the round
/// trip lossless for a figure drawn rather than only for a document typed.
#[derive(Debug, Clone, Default, PartialEq, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Layout {
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub instances: BTreeMap<String, [f64; 2]>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub inputs: BTreeMap<String, [f64; 2]>,
    #[serde(default, skip_serializing_if = "BTreeMap::is_empty")]
    pub products: BTreeMap<String, [f64; 2]>,
}

impl Layout {
    /// The position stored for one node, by role.
    #[must_use]
    pub fn position(&self, role: crate::check::NodeRole, name: &str) -> Option<[f64; 2]> {
        let placed = match role {
            crate::check::NodeRole::Instance => &self.instances,
            crate::check::NodeRole::Input => &self.inputs,
            crate::check::NodeRole::Product => &self.products,
        };
        placed.get(name).copied()
    }
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
