//! The checker: a flowsheet is held to the calculus's rules, made mechanical.
//!
//! `docs/src/calculus/process.md` states a unit operation as a process on typed,
//! directional channels, conservation as linearity, and feedback as restriction.
//! This module turns those into checks a machine can run:
//!
//! - every instance names a palette unit op, and only its declared parameters;
//! - a connection joins an outlet (or feed) to an inlet (or product);
//! - a Port-to-Port connection joins dimension-compatible field records;
//! - an input's declared record could be a stream — one mole fraction per substance;
//! - an endpoint's position, where it wrote one, is a position the port has;
//! - **linearity** — a `one` port is consumed/produced exactly once, a `many` port
//!   at least once, and every feed/product is used exactly once;
//! - every loop in the instance graph passes through a declared recycle.

use std::collections::{HashMap, HashSet};

use azoth_core::unit_vocab_gen::dimension_exponents;
use azoth_core::units::UNIT_NAMES;

use crate::channel::{ChannelType, Direction, FieldType, Multiplicity};
use crate::flowsheet::{Connection, Flowsheet, Instance, Recycle, split_endpoint};
use crate::unit_op::UnitOpSpec;

/// What a `validate` run found wrong.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Diagnostic {
    UnknownDimension {
        unit_op: String,
        port: String,
        field: String,
        dimension: String,
    },
    /// A parameter's declared unit is not one this crate can convert.
    ///
    /// **The sibling of `UnknownDimension`, and it was missing until P11.** A port
    /// field's dimension was held to the vocabulary from the start and a parameter's unit
    /// was not, which is how `unit_ops.tank` carried a `pressure_drop` in `Pa` for a class
    /// whose `run` never reads a pressure drop — a declared quantity with nothing behind it
    /// and nothing to catch it.
    UnknownParameterUnit {
        unit_op: String,
        parameter: String,
        unit: String,
    },
    DuplicateUnitOpId {
        id: String,
    },
    DuplicatePort {
        unit_op: String,
        port: String,
    },
    DuplicateInstance {
        id: String,
    },
    DuplicateFeed {
        name: String,
    },
    DuplicateProduct {
        name: String,
    },
    NameCollision {
        name: String,
    },
    UnknownUnitOp {
        instance: String,
        unit: String,
    },
    UnknownParameter {
        instance: String,
        parameter: String,
    },
    UnknownFeed {
        name: String,
    },
    UnknownProduct {
        name: String,
    },
    UnknownInstance {
        name: String,
    },
    UnknownPort {
        instance: String,
        port: String,
    },
    ProducerNotOutlet {
        instance: String,
        port: String,
    },
    ConsumerNotInlet {
        instance: String,
        port: String,
    },
    TypeMismatch {
        from: String,
        to: String,
        detail: String,
    },
    OverfedPort {
        instance: String,
        port: String,
        count: usize,
    },
    UnderfedPort {
        instance: String,
        port: String,
        count: usize,
    },
    UnusedFeed {
        name: String,
    },
    UnusedProduct {
        name: String,
    },
    /// An input's declared record does not describe a stream.
    ///
    /// **A shape rule, which is this checker's character.** It compares dimensions by exponent
    /// tuple and never resolves the databank, so whether a named substance *exists* is not
    /// decidable here - that refusal is the kernel's, at `Flowsheet::input_streams`. What is
    /// decidable is whether the record could be a stream at all: a composition is one mole
    /// fraction per substance, and a fluid names at least one.
    InputRecord {
        name: String,
        detail: String,
    },
    /// An endpoint names a stream of a port by a position the port does not have.
    ///
    /// **`split1.products[0]` is the only thing this can be about**, and the three ways to get it
    /// wrong are all one mistake: an index on a `one` outlet, an index on an *inlet*, and an index
    /// past the end of a `many` one. The range is not decidable here - how many streams a splitter
    /// returns is its `split_factors`' length, which is a value and not a declaration - so an
    /// index one past the end is refused by the run rather than by this.
    EndpointIndex {
        endpoint: String,
        detail: String,
    },
    UnrecycledLoop {
        detail: String,
    },
}

/// How bad a diagnostic is.
///
/// **Two levels, because a front-end draws two.** The middleware's rule is that a red arrow is a
/// structured `OverfedPort` and not a string, and what makes it red rather than amber is this.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Severity {
    /// The document cannot run as written.
    Error,
    /// The document runs, and something about it is probably not what the author meant.
    Warning,
}

impl Severity {
    /// The name a front-end switches on, which is the variant's own.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Self::Error => "error",
            Self::Warning => "warning",
        }
    }
}

/// Where in a document a diagnostic is about.
///
/// **A section and a path**, which is what a widget needs to put a mark on the right row: the
/// section is the table (`unit_ops`, `instances`, `connections`, `inputs`, `products`) and the path
/// is the id within it, dotted the way the document spells it. **An empty path is the section
/// whole** - a cycle nobody declared is a fact about `connections` and not about one row of it.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Location {
    /// The table the diagnostic is about.
    pub section: &'static str,
    /// The id within it, dotted; empty for the section itself.
    pub path: String,
}

impl std::fmt::Display for Location {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        if self.path.is_empty() {
            f.write_str(self.section)
        } else {
            write!(f, "{}/{}", self.section, self.path)
        }
    }
}

impl Diagnostic {
    /// How bad this is.
    ///
    /// **Only two variants are warnings**, and both are the same statement: the document declares
    /// something nothing uses. It runs; it is probably not what was meant. Everything else makes
    /// the document unrunnable, which is the line the two levels are drawn on.
    #[must_use]
    pub fn severity(&self) -> Severity {
        match self {
            Self::UnusedFeed { .. } | Self::UnusedProduct { .. } => Severity::Warning,
            _ => Severity::Error,
        }
    }

    /// Where in the document this is about.
    #[must_use]
    pub fn location(&self) -> Location {
        let at = |section: &'static str, path: String| Location { section, path };
        match self {
            Self::UnknownDimension {
                unit_op,
                port,
                field,
                ..
            } => at("palette", format!("{unit_op}.ports.{port}.{field}")),
            Self::UnknownParameterUnit {
                unit_op, parameter, ..
            } => at("palette", format!("{unit_op}.parameters.{parameter}")),
            Self::DuplicateUnitOpId { id } => at("palette", id.clone()),
            Self::DuplicatePort { unit_op, port } => {
                at("palette", format!("{unit_op}.ports.{port}"))
            }
            Self::DuplicateInstance { id } => at("instances", id.clone()),
            Self::DuplicateFeed { name } => at("inputs", name.clone()),
            Self::DuplicateProduct { name } => at("products", name.clone()),
            Self::NameCollision { name } => at("flowsheet", name.clone()),
            Self::UnknownUnitOp { instance, .. } => at("instances", format!("{instance}.unit")),
            Self::UnknownParameter {
                instance,
                parameter,
            } => at("instances", format!("{instance}.parameters.{parameter}")),
            Self::UnknownFeed { name } => at("inputs", name.clone()),
            Self::UnknownProduct { name } => at("products", name.clone()),
            Self::UnknownInstance { name } => at("connections", name.clone()),
            Self::UnknownPort { instance, port }
            | Self::ProducerNotOutlet { instance, port }
            | Self::ConsumerNotInlet { instance, port }
            | Self::OverfedPort { instance, port, .. }
            | Self::UnderfedPort { instance, port, .. } => {
                at("instances", format!("{instance}.ports.{port}"))
            }
            Self::TypeMismatch { from, to, .. } => at("connections", format!("{from} -> {to}")),
            Self::UnusedFeed { name } => at("inputs", name.clone()),
            Self::UnusedProduct { name } => at("products", name.clone()),
            Self::InputRecord { name, .. } => at("inputs", name.clone()),
            Self::EndpointIndex { endpoint, .. } => at("connections", endpoint.clone()),
            Self::UnrecycledLoop { .. } => at("connections", String::new()),
        }
    }
}

/// The human line, which is what a log or a terminal gets.
///
/// **`Debug` is left alone.** The Python bridge renders a diagnostic with `{d:?}` today and
/// `python/tests/test_process.py` asserts that a broken flowsheet's lines contain `OverfedPort`,
/// so a `Display` that replaced it would move a published surface for no gain. This is additive:
/// the variant's name is still in the debug form, and a caller that wants prose has one.
impl std::fmt::Display for Diagnostic {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{}: {}", self.location(), self.message())
    }
}

impl Diagnostic {
    /// What went wrong, in a sentence.
    ///
    /// Naming the subject here rather than in `Display` keeps one place that knows how to say a
    /// variant and one place that knows where it is.
    #[must_use]
    pub fn message(&self) -> String {
        match self {
            Self::UnknownDimension { dimension, .. } => {
                format!("`{dimension}` is not a dimension the vocabulary carries")
            }
            Self::UnknownParameterUnit { unit, .. } => {
                format!("`{unit}` is not a unit the vocabulary carries")
            }
            Self::DuplicateUnitOpId { id } => format!("`{id}` is declared twice"),
            Self::DuplicatePort { port, .. } => format!("port `{port}` is declared twice"),
            Self::DuplicateInstance { id } => format!("instance `{id}` is declared twice"),
            Self::DuplicateFeed { name } => format!("feed `{name}` is declared twice"),
            Self::DuplicateProduct { name } => format!("product `{name}` is declared twice"),
            Self::NameCollision { name } => {
                format!("`{name}` is used both as a feed and as a product")
            }
            Self::UnknownUnitOp { unit, .. } => {
                format!("`{unit}` is not a unit operation the palette declares")
            }
            Self::UnknownParameter { parameter, .. } => {
                format!("`{parameter}` is not a parameter the entry declares")
            }
            Self::UnknownFeed { name } => format!("feed `{name}` is not declared"),
            Self::UnknownProduct { name } => format!("product `{name}` is not declared"),
            Self::UnknownInstance { name } => format!("instance `{name}` is not declared"),
            Self::UnknownPort { port, .. } => format!("`{port}` is not a port of this entry"),
            Self::ProducerNotOutlet { port, .. } => {
                format!("`{port}` produces, and is not an outlet")
            }
            Self::ConsumerNotInlet { port, .. } => {
                format!("`{port}` consumes, and is not an inlet")
            }
            Self::TypeMismatch { detail, .. } => detail.clone(),
            Self::OverfedPort { port, count, .. } => {
                format!("`{port}` takes one stream and {count} are connected")
            }
            Self::UnderfedPort { port, count, .. } => {
                format!("`{port}` takes at least one stream and {count} are connected")
            }
            Self::UnusedFeed { name } => format!("feed `{name}` is declared and never consumed"),
            Self::UnusedProduct { name } => {
                format!("product `{name}` is declared and never produced")
            }
            Self::InputRecord { detail, .. } => detail.clone(),
            Self::EndpointIndex { detail, .. } => detail.clone(),
            Self::UnrecycledLoop { detail } => detail.clone(),
        }
    }
}

/// The palette's own well-formedness: unique ids and ports, and every field names
/// a dimension the vocabulary knows.
#[must_use]
pub fn validate_palette(specs: &[UnitOpSpec]) -> Vec<Diagnostic> {
    let mut diags = Vec::new();
    let mut ids: HashSet<&str> = HashSet::new();

    for spec in specs {
        if !ids.insert(&spec.id) {
            diags.push(Diagnostic::DuplicateUnitOpId {
                id: spec.id.clone(),
            });
        }

        let mut ports: HashSet<&str> = HashSet::new();
        for port in &spec.ports {
            if !ports.insert(&port.name) {
                diags.push(Diagnostic::DuplicatePort {
                    unit_op: spec.id.clone(),
                    port: port.name.clone(),
                });
            }
            for (field_name, field) in &port.fields {
                if dimension_exponents(&field.dimension).is_none() {
                    diags.push(Diagnostic::UnknownDimension {
                        unit_op: spec.id.clone(),
                        port: port.name.clone(),
                        field: field_name.clone(),
                        dimension: field.dimension.clone(),
                    });
                }
            }
        }

        for (name, parameter) in &spec.parameters {
            if let Some(unit) = &parameter.unit {
                if !UNIT_NAMES.contains(&unit.as_str()) {
                    diags.push(Diagnostic::UnknownParameterUnit {
                        unit_op: spec.id.clone(),
                        parameter: name.clone(),
                        unit: unit.clone(),
                    });
                }
            }
        }
    }
    diags
}

/// Validate a flowsheet against a palette. Returns every defect found, empty when
/// the flowsheet honours the rules.
#[must_use]
pub fn validate(flowsheet: &Flowsheet, palette: &[UnitOpSpec]) -> Vec<Diagnostic> {
    let specs: HashMap<&str, &UnitOpSpec> = palette.iter().map(|s| (s.id.as_str(), s)).collect();
    let instances: HashMap<&str, &Instance> = flowsheet
        .instances
        .iter()
        .map(|i| (i.id.as_str(), i))
        .collect();

    let mut diags = Vec::new();

    // Names: instance ids and boundary names must be unique and disjoint.
    let mut instance_ids: HashSet<&str> = HashSet::new();
    for instance in &flowsheet.instances {
        if !instance_ids.insert(&instance.id) {
            diags.push(Diagnostic::DuplicateInstance {
                id: instance.id.clone(),
            });
        }
    }
    let mut feed_names: HashSet<&str> = HashSet::new();
    for feed in flowsheet.input_names() {
        if !feed_names.insert(feed) {
            diags.push(Diagnostic::DuplicateFeed {
                name: feed.to_string(),
            });
        }
        if instance_ids.contains(feed) {
            diags.push(Diagnostic::NameCollision {
                name: feed.to_string(),
            });
        }
    }
    // Inputs: the record has to be able to be a stream at all. Whether the substances it names
    // exist is the databank's answer and not this checker's - see `Diagnostic::InputRecord`.
    for input in &flowsheet.inputs {
        let detail = if input.components.is_empty() {
            Some(format!(
                "the input `{}` names no substance, so it is not a fluid",
                input.name
            ))
        } else if input.z.len() != input.components.len() {
            Some(format!(
                "the input `{}` has {} mole fractions for {} substances",
                input.name,
                input.z.len(),
                input.components.len()
            ))
        } else {
            None
        };
        if let Some(detail) = detail {
            diags.push(Diagnostic::InputRecord {
                name: input.name.clone(),
                detail,
            });
        }
    }
    let mut product_names: HashSet<&str> = HashSet::new();
    for product in &flowsheet.products {
        if !product_names.insert(product) {
            diags.push(Diagnostic::DuplicateProduct {
                name: product.clone(),
            });
        }
        if instance_ids.contains(product.as_str()) {
            diags.push(Diagnostic::NameCollision {
                name: product.clone(),
            });
        }
    }

    // Instances: the unit op exists, and the parameters are the ones it declares.
    for instance in &flowsheet.instances {
        match specs.get(instance.unit.as_str()) {
            None => diags.push(Diagnostic::UnknownUnitOp {
                instance: instance.id.clone(),
                unit: instance.unit.clone(),
            }),
            Some(spec) => {
                for parameter in instance.parameters.keys() {
                    if !spec.parameters.contains_key(parameter) {
                        diags.push(Diagnostic::UnknownParameter {
                            instance: instance.id.clone(),
                            parameter: parameter.clone(),
                        });
                    }
                }
            }
        }
    }

    // Connections: resolve endpoints, check direction, check types, and count
    // each port's incoming/outgoing streams for the linearity pass. A recycle is
    // a connection that also closes a loop; both join a producer to a consumer.
    let edges: Vec<(&str, &str)> = flowsheet
        .connections
        .iter()
        .map(|c| (c.from.as_str(), c.to.as_str()))
        .chain(
            flowsheet
                .recycles
                .iter()
                .map(|r| (r.from.as_str(), r.to.as_str())),
        )
        .collect();

    let mut incoming: HashMap<String, usize> = HashMap::new();
    let mut outgoing: HashMap<String, usize> = HashMap::new();
    let mut feed_uses: HashMap<&str, usize> = HashMap::new();
    let mut product_uses: HashMap<&str, usize> = HashMap::new();

    for (from, to) in edges {
        let producer = resolve_producer(from);
        let consumer = resolve_consumer(to);

        match producer {
            Endpoint::Feed(name) => {
                if !feed_names.contains(name) {
                    diags.push(Diagnostic::UnknownFeed {
                        name: name.to_string(),
                    });
                }
                *feed_uses.entry(name).or_insert(0) += 1;
            }
            Endpoint::Port {
                instance,
                port,
                index,
            } => {
                // **The index is checked where the producer is**, because that is the only side a
                // `many` outlet has: a connection addresses one of `split1.products`' streams by
                // position, while a `many` inlet takes several connections to the same port name.
                check_index(&mut diags, &instances, &specs, instance, port, index);
                check_port(
                    &mut diags,
                    &instances,
                    &specs,
                    instance,
                    port,
                    Direction::Out,
                    &mut outgoing,
                );
            }
            // `resolve_producer` never yields a Product.
            Endpoint::Product(_) => {}
        }

        match consumer {
            Endpoint::Product(name) => {
                if !product_names.contains(name) {
                    diags.push(Diagnostic::UnknownProduct {
                        name: name.to_string(),
                    });
                }
                *product_uses.entry(name).or_insert(0) += 1;
            }
            Endpoint::Port {
                instance,
                port,
                index,
            } => {
                // A consumer endpoint is one port name however many connections reach it, so an
                // index here is not a position - the resolver drops it, and this reports it.
                check_index(&mut diags, &instances, &specs, instance, port, index);
                check_port(
                    &mut diags,
                    &instances,
                    &specs,
                    instance,
                    port,
                    Direction::In,
                    &mut incoming,
                );
            }
            // `resolve_consumer` never yields a Feed.
            Endpoint::Feed(_) => {}
        }

        // Type compatibility, only where both ends declare a channel type — that
        // is, a Port-to-Port edge. A feed or product adopts the port it meets.
        if let (
            Endpoint::Port {
                instance: from_i,
                port: from_p,
                ..
            },
            Endpoint::Port {
                instance: to_i,
                port: to_p,
                ..
            },
        ) = (producer, consumer)
        {
            if let (Some(from_t), Some(to_t)) = (
                channel_type(&instances, &specs, from_i, from_p),
                channel_type(&instances, &specs, to_i, to_p),
            ) {
                if let Err(detail) = compatible(from_t, to_t) {
                    diags.push(Diagnostic::TypeMismatch {
                        from: format!("{from_i}.{from_p}"),
                        to: format!("{to_i}.{to_p}"),
                        detail,
                    });
                }
            }
        }
    }

    // Linearity: each port's stream count must match its multiplicity.
    for instance in &flowsheet.instances {
        let Some(spec) = specs.get(instance.unit.as_str()) else {
            continue;
        };
        for port in &spec.ports {
            let key = format!("{}.{}", instance.id, port.name);
            let count = match port.direction {
                Direction::In => *incoming.get(&key).unwrap_or(&0),
                Direction::Out => *outgoing.get(&key).unwrap_or(&0),
            };
            let ok = match port.multiplicity {
                Multiplicity::One => count == 1,
                Multiplicity::Many => count >= 1,
            };
            if !ok {
                diags.push(if count == 0 {
                    Diagnostic::UnderfedPort {
                        instance: instance.id.clone(),
                        port: port.name.clone(),
                        count,
                    }
                } else {
                    Diagnostic::OverfedPort {
                        instance: instance.id.clone(),
                        port: port.name.clone(),
                        count,
                    }
                });
            }
        }
    }

    // Boundary streams: each feed and product is used exactly once.
    for feed in flowsheet.input_names() {
        if feed_uses.get(feed).copied().unwrap_or(0) != 1 {
            diags.push(Diagnostic::UnusedFeed {
                name: feed.to_string(),
            });
        }
    }
    for product in &flowsheet.products {
        if product_uses.get(product.as_str()).copied().unwrap_or(0) != 1 {
            diags.push(Diagnostic::UnusedProduct {
                name: product.clone(),
            });
        }
    }

    // Feedback is restriction: the instance graph with the recycles removed must
    // be acyclic, or a loop goes undeclared.
    if let Some(detail) = undeclared_loop(&flowsheet.connections, &flowsheet.recycles) {
        diags.push(Diagnostic::UnrecycledLoop { detail });
    }

    diags
}

/// A resolved endpoint of an edge.
#[derive(Debug, Clone, Copy)]
enum Endpoint<'a> {
    Feed(&'a str),
    Product(&'a str),
    Port {
        instance: &'a str,
        port: &'a str,
        /// The stream of a `many` outlet, where the endpoint named one.
        index: Option<usize>,
    },
}

fn resolve_producer(s: &str) -> Endpoint<'_> {
    match split_endpoint(s) {
        Some((instance, port, index)) => Endpoint::Port {
            instance,
            port,
            index,
        },
        None => Endpoint::Feed(s),
    }
}

fn resolve_consumer(s: &str) -> Endpoint<'_> {
    match split_endpoint(s) {
        // **The index is carried rather than dropped, because it is not a position on this side.**
        // A `many` *inlet* takes several connections to the same port name, so an index written
        // here is a mistake - and one this resolver cannot silently discard, or `check_index`
        // would never see it. The rest of the pass reads the port name either way.
        Some((instance, port, index)) => Endpoint::Port {
            instance,
            port,
            index,
        },
        None => Endpoint::Product(s),
    }
}

/// Check that an endpoint's index, where it wrote one, is a position the port has.
///
/// **Three ways to be wrong, one mistake**: an index on a `one` outlet says "the second of one
/// stream", an index on an inlet says it about a port that takes several connections to one name,
/// and an index on a port the entry does not declare duplicates `UnknownPort`. All three are
/// refused rather than dropped, because a dropped index binds the *wrong* stream silently - which
/// is the failure the `many` outlet existed to avoid.
fn check_index(
    diags: &mut Vec<Diagnostic>,
    instances: &HashMap<&str, &Instance>,
    specs: &HashMap<&str, &UnitOpSpec>,
    instance: &str,
    port: &str,
    index: Option<usize>,
) {
    let Some(index) = index else {
        return;
    };
    let endpoint = crate::flowsheet::stream_path(instance, port, Some(index));
    let detail = match instances
        .get(instance)
        .and_then(|inst| specs.get(inst.unit.as_str()))
        .and_then(|spec| spec.ports.iter().find(|p| p.name == port))
    {
        // `UnknownPort` and `UnknownInstance` are already reported by `check_port`.
        None => return,
        Some(p) if p.direction == Direction::In => {
            "an inlet is one port name however many streams reach it, so it has no positions - \
             connect it once per stream instead"
        }
        Some(p) if p.multiplicity != Multiplicity::Many => {
            "the port declares one stream, so it has no second"
        }
        Some(_) => return,
    };
    diags.push(Diagnostic::EndpointIndex {
        endpoint: endpoint.clone(),
        detail: format!("`{endpoint}` names a position: {detail}"),
    });
}

/// Check that `instance.port` exists with the expected direction, and record the
/// stream it carries in `counts`.
fn check_port(
    diags: &mut Vec<Diagnostic>,
    instances: &HashMap<&str, &Instance>,
    specs: &HashMap<&str, &UnitOpSpec>,
    instance: &str,
    port: &str,
    expected: Direction,
    counts: &mut HashMap<String, usize>,
) {
    let Some(inst) = instances.get(instance) else {
        diags.push(Diagnostic::UnknownInstance {
            name: instance.to_string(),
        });
        return;
    };
    let Some(spec) = specs.get(inst.unit.as_str()) else {
        return; // UnknownUnitOp is already reported in the instance pass.
    };
    match spec.ports.iter().find(|p| p.name == port) {
        None => diags.push(Diagnostic::UnknownPort {
            instance: instance.to_string(),
            port: port.to_string(),
        }),
        Some(p) if p.direction != expected => {
            let d = match expected {
                Direction::Out => Diagnostic::ProducerNotOutlet {
                    instance: instance.to_string(),
                    port: port.to_string(),
                },
                Direction::In => Diagnostic::ConsumerNotInlet {
                    instance: instance.to_string(),
                    port: port.to_string(),
                },
            };
            diags.push(d);
        }
        Some(_) => {
            *counts.entry(format!("{instance}.{port}")).or_insert(0) += 1;
        }
    }
}

/// The channel type a port declares, if the instance and port resolve.
fn channel_type<'s>(
    instances: &HashMap<&str, &Instance>,
    specs: &'s HashMap<&str, &UnitOpSpec>,
    instance: &str,
    port: &str,
) -> Option<&'s ChannelType> {
    let spec = specs.get(instances.get(instance)?.unit.as_str())?;
    Some(&spec.ports.iter().find(|p| p.name == port)?.fields)
}

/// Two channel types are compatible when they carry the same field names, each at
/// the same dimension (compared by exponent tuple, not by id string) and shape.
fn compatible(a: &ChannelType, b: &ChannelType) -> Result<(), String> {
    if a.len() != b.len() {
        return Err(format!(
            "field sets differ: {:?} vs {:?}",
            a.keys().collect::<Vec<_>>(),
            b.keys().collect::<Vec<_>>()
        ));
    }
    for (name, field) in a {
        match b.get(name) {
            None => return Err(format!("field {name:?} present on only one side")),
            Some(other) if !fields_compatible(field, other) => {
                return Err(format!("field {name:?} differs: {field:?} vs {other:?}"));
            }
            Some(_) => {}
        }
    }
    Ok(())
}

/// Two fields are compatible at the same dimension (exponents) and shape.
fn fields_compatible(a: &FieldType, b: &FieldType) -> bool {
    a.shape == b.shape
        && dimension_exponents(&a.dimension).is_some()
        && dimension_exponents(&a.dimension) == dimension_exponents(&b.dimension)
}

/// The instance graph's non-recycle edges must be acyclic. Returns a description
/// of the undeclared loop when one exists.
fn undeclared_loop(connections: &[Connection], recycles: &[Recycle]) -> Option<String> {
    let mut edges: Vec<(String, String)> = Vec::new();
    let mut nodes: HashSet<String> = HashSet::new();

    for c in connections {
        if let (Some((a, _)), Some((b, _))) = (split_port(&c.from), split_port(&c.to)) {
            nodes.insert(a.clone());
            nodes.insert(b.clone());
            edges.push((a, b));
        }
    }
    // Recycles are torn by declaration, so they are not edges of the graph that
    // must be acyclic.
    let _ = recycles;

    // Kahn's algorithm: repeatedly drop nodes with no incoming edge.
    let mut indegree: HashMap<&str, usize> = HashMap::new();
    for (_, to) in &edges {
        *indegree.entry(to.as_str()).or_insert(0) += 1;
    }
    let mut ready: Vec<&str> = nodes
        .iter()
        .map(String::as_str)
        .filter(|n| indegree.get(n).copied().unwrap_or(0) == 0)
        .collect();
    let mut processed = 0usize;
    while let Some(node) = ready.pop() {
        processed += 1;
        for (from, to) in &edges {
            if from == node {
                let entry = indegree.entry(to.as_str()).or_insert(0);
                *entry = entry.saturating_sub(1);
                if *entry == 0 {
                    ready.push(to.as_str());
                }
            }
        }
    }

    if processed < nodes.len() {
        Some("an instance cycle is not declared as a recycle".to_string())
    } else {
        None
    }
}

fn split_port(s: &str) -> Option<(String, String)> {
    s.split_once('.')
        .map(|(a, b)| (a.to_string(), b.to_string()))
}

#[cfg(test)]
mod tests {
    use std::collections::BTreeMap;

    use super::*;
    use crate::channel::{Port, Shape};
    use crate::flowsheet::named_input;

    fn field(dimension: &str) -> FieldType {
        FieldType {
            dimension: dimension.to_string(),
            shape: Shape::Scalar,
        }
    }

    fn vec_field(dimension: &str) -> FieldType {
        FieldType {
            dimension: dimension.to_string(),
            shape: Shape::Vector,
        }
    }

    /// The shared material-stream record the palette adopts.
    fn stream() -> ChannelType {
        ChannelType::from([
            ("n".to_string(), field("molar_flow")),
            ("z".to_string(), vec_field("dimensionless")),
            ("P".to_string(), field("pressure")),
            ("T".to_string(), field("thermodynamic_temperature")),
            ("h".to_string(), field("molar_energy")),
        ])
    }

    fn port(name: &str, direction: Direction, multiplicity: Multiplicity) -> Port {
        Port {
            name: name.to_string(),
            direction,
            multiplicity,
            fields: stream(),
        }
    }

    /// A one-inlet / one-outlet unit op, e.g. a pump or pipe.
    fn two_port(id: &str) -> UnitOpSpec {
        UnitOpSpec {
            id: id.to_string(),
            name: id.to_string(),
            source: None,
            parameters: BTreeMap::new(),
            ports: vec![
                port("feed", Direction::In, Multiplicity::One),
                port("discharge", Direction::Out, Multiplicity::One),
            ],
            notes: None,
        }
    }

    /// A many-inlet / one-outlet unit op, e.g. a mixer.
    fn mixer(id: &str) -> UnitOpSpec {
        UnitOpSpec {
            id: id.to_string(),
            name: id.to_string(),
            source: None,
            parameters: BTreeMap::new(),
            ports: vec![
                port("feed", Direction::In, Multiplicity::Many),
                port("product", Direction::Out, Multiplicity::One),
            ],
            notes: None,
        }
    }

    fn splitter(id: &str) -> UnitOpSpec {
        UnitOpSpec {
            id: id.to_string(),
            name: id.to_string(),
            source: None,
            parameters: BTreeMap::new(),
            ports: vec![
                port("feed", Direction::In, Multiplicity::One),
                port("products", Direction::Out, Multiplicity::Many),
            ],
            notes: None,
        }
    }

    fn assert_clean(flowsheet: &Flowsheet, palette: &[UnitOpSpec]) {
        let diags = validate(flowsheet, palette);
        assert!(diags.is_empty(), "expected clean, got {diags:?}");
    }

    fn has(diags: &[Diagnostic], want: Diagnostic) -> bool {
        diags.contains(&want)
    }

    fn has_variant(diags: &[Diagnostic], want: impl Fn(&Diagnostic) -> bool) -> bool {
        diags.iter().any(want)
    }

    #[test]
    fn a_feed_through_one_unit_op_is_clean() {
        let palette = vec![two_port("unit_ops.pump")];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("feed_1")],
            products: vec!["purge".into()],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            }],
            connections: vec![
                Connection {
                    from: "feed_1".into(),
                    to: "p1.feed".into(),
                },
                Connection {
                    from: "p1.discharge".into(),
                    to: "purge".into(),
                },
            ],
            recycles: vec![],
        };
        assert_clean(&flowsheet, &palette);
    }

    /// **An input has to be able to be a stream.** This is the rule the schema gained with
    /// `[[inputs]]`: the record the user writes is checked for shape here, and the substances it
    /// names are checked for existence by `Flowsheet::input_streams`, which is where the databank
    /// is. Splitting it that way is deliberate - this checker compares dimensions by exponent
    /// tuple and resolves no names at all.
    #[test]
    fn an_input_whose_record_is_not_a_stream_is_reported() {
        let palette = vec![two_port("unit_ops.pump")];
        let mut flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("feed_1")],
            products: vec!["purge".into()],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            }],
            connections: vec![
                Connection {
                    from: "feed_1".into(),
                    to: "p1.feed".into(),
                },
                Connection {
                    from: "p1.discharge".into(),
                    to: "purge".into(),
                },
            ],
            recycles: vec![],
        };
        assert_clean(&flowsheet, &palette);

        // Two substances and one mole fraction.
        flowsheet.inputs[0].components = vec!["methane".into(), "n-butane".into()];
        let diags = validate(&flowsheet, &palette);
        assert!(
            has_variant(&diags, |d| matches!(
                d,
                Diagnostic::InputRecord { name, detail }
                    if name == "feed_1" && detail.contains("1 mole fractions for 2 substances")
            )),
            "{diags:?}"
        );
        let record = diags
            .iter()
            .find(|d| matches!(d, Diagnostic::InputRecord { .. }))
            .expect("reported");
        assert_eq!(record.severity(), Severity::Error, "it cannot run");
        assert_eq!(record.location().section, "inputs");
        assert_eq!(record.location().path, "feed_1");

        // And a fluid that names nothing is not a fluid.
        flowsheet.inputs[0].components = vec![];
        flowsheet.inputs[0].z = vec![];
        let diags = validate(&flowsheet, &palette);
        assert!(
            has_variant(&diags, |d| matches!(
                d,
                Diagnostic::InputRecord { detail, .. } if detail.contains("names no substance")
            )),
            "{diags:?}"
        );
    }

    #[test]
    fn an_unknown_unit_op_is_reported() {
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![],
            products: vec![],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.not_real".into(),
                parameters: BTreeMap::new(),
            }],
            connections: vec![],
            recycles: vec![],
        };
        let diags = validate(&flowsheet, &[]);
        assert!(has(
            &diags,
            Diagnostic::UnknownUnitOp {
                instance: "p1".into(),
                unit: "unit_ops.not_real".into()
            }
        ));
    }

    #[test]
    fn an_undeclared_parameter_is_reported() {
        let mut pump = two_port("unit_ops.pump");
        pump.parameters.insert(
            "dp".to_string(),
            crate::unit_op::Param {
                unit: Some("Pa".to_string()),
                description: "rise".to_string(),
            },
        );
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![],
            products: vec![],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::from([("efficiency".to_string(), toml::Value::from(0.8))]),
            }],
            connections: vec![],
            recycles: vec![],
        };
        let diags = validate(&flowsheet, &[pump]);
        assert!(has(
            &diags,
            Diagnostic::UnknownParameter {
                instance: "p1".into(),
                parameter: "efficiency".into()
            }
        ));
    }

    #[test]
    fn an_inlet_used_as_a_producer_is_reported() {
        let palette = vec![two_port("unit_ops.pump")];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("feed_1")],
            products: vec![],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            }],
            // `from` is an inlet: an inlet cannot be a producer.
            connections: vec![Connection {
                from: "p1.feed".into(),
                to: "p1.discharge".into(),
            }],
            recycles: vec![],
        };
        let diags = validate(&flowsheet, &palette);
        assert!(has(
            &diags,
            Diagnostic::ProducerNotOutlet {
                instance: "p1".into(),
                port: "feed".into()
            }
        ));
    }

    #[test]
    fn an_outlet_used_as_a_consumer_is_reported() {
        let palette = vec![two_port("unit_ops.pump")];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("feed_1")],
            products: vec![],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            }],
            // `to` is an outlet: an outlet cannot be a consumer.
            connections: vec![Connection {
                from: "feed_1".into(),
                to: "p1.discharge".into(),
            }],
            recycles: vec![],
        };
        let diags = validate(&flowsheet, &palette);
        assert!(has(
            &diags,
            Diagnostic::ConsumerNotInlet {
                instance: "p1".into(),
                port: "discharge".into()
            }
        ));
    }

    #[test]
    fn a_type_mismatch_between_two_ports_is_reported() {
        let mut weird = two_port("unit_ops.weird");
        // A port whose record is not the material stream: only n and P.
        weird.ports[0].fields = ChannelType::from([
            ("n".to_string(), field("molar_flow")),
            ("P".to_string(), field("pressure")),
        ]);
        let palette = vec![two_port("unit_ops.pump"), weird];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![],
            products: vec![],
            instances: vec![
                Instance {
                    id: "p1".into(),
                    unit: "unit_ops.pump".into(),
                    parameters: BTreeMap::new(),
                },
                Instance {
                    id: "w1".into(),
                    unit: "unit_ops.weird".into(),
                    parameters: BTreeMap::new(),
                },
            ],
            connections: vec![Connection {
                from: "p1.discharge".into(),
                to: "w1.feed".into(),
            }],
            recycles: vec![],
        };
        let diags = validate(&flowsheet, &palette);
        assert!(has_variant(&diags, |d| matches!(
            d,
            Diagnostic::TypeMismatch { .. }
        )));
    }

    #[test]
    fn a_double_fed_inlet_is_reported() {
        let palette = vec![two_port("unit_ops.pump")];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("a"), named_input("b")],
            products: vec![],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            }],
            connections: vec![
                Connection {
                    from: "a".into(),
                    to: "p1.feed".into(),
                },
                Connection {
                    from: "b".into(),
                    to: "p1.feed".into(),
                },
            ],
            recycles: vec![],
        };
        let diags = validate(&flowsheet, &palette);
        assert!(has(
            &diags,
            Diagnostic::OverfedPort {
                instance: "p1".into(),
                port: "feed".into(),
                count: 2
            }
        ));
    }

    #[test]
    fn an_unrecycled_loop_is_reported_and_a_declared_recycle_is_not() {
        let palette = vec![two_port("unit_ops.pump")];
        let instances = vec![
            Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            },
            Instance {
                id: "p2".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            },
        ];

        // p1 -> p2 and p2 -> p1, with no declared recycle.
        let looped = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![],
            products: vec![],
            instances: instances.clone(),
            connections: vec![
                Connection {
                    from: "p1.discharge".into(),
                    to: "p2.feed".into(),
                },
                Connection {
                    from: "p2.discharge".into(),
                    to: "p1.feed".into(),
                },
            ],
            recycles: vec![],
        };
        assert!(has_variant(&validate(&looped, &palette), |d| matches!(
            d,
            Diagnostic::UnrecycledLoop { .. }
        )));

        // The same loop, but the back edge is declared as a recycle.
        let torn = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![],
            products: vec![],
            instances,
            connections: vec![Connection {
                from: "p1.discharge".into(),
                to: "p2.feed".into(),
            }],
            recycles: vec![Recycle::new("r1", "p2.discharge", "p1.feed")],
        };
        assert!(!has_variant(&validate(&torn, &palette), |d| matches!(
            d,
            Diagnostic::UnrecycledLoop { .. }
        )));
    }

    #[test]
    fn an_unknown_parameter_unit_in_the_palette_is_reported() {
        let mut spec = two_port("unit_ops.tank");
        spec.parameters.insert(
            "volume".to_string(),
            crate::unit_op::Param {
                unit: Some("m**3".to_string()),
                description: "the tank's volume".to_string(),
            },
        );
        let diags = validate_palette(&[spec]);
        // `m**3` is not `m**3/s`: the vocabulary has no plain volume unit, so a volume
        // cannot be declared until one is added. That is the state this records rather
        // than a spelling to be worked around.
        assert!(has(
            &diags,
            Diagnostic::UnknownParameterUnit {
                unit_op: "unit_ops.tank".into(),
                parameter: "volume".into(),
                unit: "m**3".into(),
            }
        ));
    }

    #[test]
    fn an_unknown_dimension_in_the_palette_is_reported() {
        let mut spec = two_port("unit_ops.pump");
        spec.ports[0].fields.insert(
            "bad".to_string(),
            FieldType {
                dimension: "no_such_dimension".to_string(),
                shape: Shape::Scalar,
            },
        );
        let diags = validate_palette(&[spec]);
        // **Matched on the field rather than on the whole variant**, because the variant now
        // carries where it is as well as what it is - which is the point of C7.
        let found = diags
            .iter()
            .find(|diag| matches!(diag, Diagnostic::UnknownDimension { dimension, .. } if dimension == "no_such_dimension"))
            .expect("the dimension is reported");
        assert_eq!(found.severity(), Severity::Error);
        assert_eq!(found.location().section, "palette");
        assert_eq!(found.location().path, "unit_ops.pump.ports.feed.bad");
        assert!(
            found.message().contains("no_such_dimension"),
            "the dimension is named in the message: {found}"
        );
        assert_eq!(
            found.to_string(),
            format!("{}: {}", found.location(), found.message())
        );
    }

    /// **Every diagnostic says how bad it is and where, and the two levels are drawn on whether
    /// the document can run at all.**
    ///
    /// The middleware's rule is that a red arrow is a structured `OverfedPort` and not a string,
    /// so the structured half has to be reachable: this walks a flowsheet with four kinds of
    /// defect at once and checks that each carries a section, a path and a level - and that the
    /// one defect that leaves the document *runnable* is the one that comes back a warning.
    #[test]
    fn every_diagnostic_carries_a_severity_and_a_location() {
        let palette = vec![mixer("unit_ops.mixer")];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            // Declared and never consumed, which is the warning.
            inputs: vec![named_input("spare")],
            products: vec!["out".into()],
            instances: vec![
                Instance {
                    id: "m1".into(),
                    unit: "unit_ops.mixer".into(),
                    parameters: BTreeMap::new(),
                },
                // A unit the palette does not carry.
                Instance {
                    id: "m2".into(),
                    unit: "unit_ops.nosuch".into(),
                    parameters: BTreeMap::new(),
                },
            ],
            connections: vec![
                // A port that is not an inlet, fed by a stream nothing produces.
                Connection {
                    from: "nowhere".into(),
                    to: "m2.product".into(),
                },
            ],
            recycles: Vec::new(),
        };

        let diags = validate(&flowsheet, &palette);
        assert!(
            !diags.is_empty(),
            "the flowsheet is broken and nothing was reported"
        );
        for diag in &diags {
            let location = diag.location();
            assert!(!location.section.is_empty(), "`{diag:?}` names no section");
            // The two levels, and no third.
            assert!(matches!(
                diag.severity(),
                Severity::Error | Severity::Warning
            ));
        }

        let warnings: Vec<&Diagnostic> = diags
            .iter()
            .filter(|diag| diag.severity() == Severity::Warning)
            .collect();
        // Two, and both are the same statement: `spare` is declared and nothing consumes it,
        // `out` is declared and nothing produces it. Every other defect here makes the document
        // unrunnable, which is where the line between the two levels is.
        assert_eq!(warnings.len(), 2, "the unused declarations: {warnings:?}");
        assert!(warnings.iter().all(|diag| matches!(
            diag,
            Diagnostic::UnusedFeed { .. } | Diagnostic::UnusedProduct { .. }
        )));
        assert_eq!(warnings[0].location().section, "inputs");
        assert_eq!(warnings[0].location().path, "spare");
        assert_eq!(warnings[1].location().section, "products");
        assert_eq!(warnings[1].location().path, "out");

        // **A red arrow, which is the one the middleware's rule names.** `m1`'s two ports are
        // both unfed, and the mis-wired `m2.product` connection surfaces as the *producer* being
        // unknown (`inputs/nowhere`) rather than as a consumer that is not an inlet - the checker
        // reports the connection's `from` first and the port it lands on is then not reached. Worth
        // knowing if a front-end wants to draw the arrow on the connection rather than on the port.
        let red = diags
            .iter()
            .find(|diag| matches!(diag, Diagnostic::UnderfedPort { .. }))
            .expect("the unfed port is reported");
        assert_eq!(red.severity(), Severity::Error);
        assert_eq!(red.location().section, "instances");
        assert_eq!(red.location().path, "m1.ports.feed");
        assert!(
            diags
                .iter()
                .any(|diag| matches!(diag, Diagnostic::UnknownFeed { name } if name == "nowhere")),
            "the unknown producer is reported: {diags:?}"
        );
        // The display line is the location and the message together.
        assert_eq!(
            red.to_string(),
            format!("{}: {}", red.location(), red.message())
        );
        assert_eq!(
            red.to_string(),
            "instances/m1.ports.feed: `feed` takes at least one stream and 0 are connected"
        );
    }

    #[test]
    fn a_mixer_accepts_many_feeds() {
        let palette = vec![mixer("unit_ops.mixer")];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("a"), named_input("b")],
            products: vec!["out".into()],
            instances: vec![Instance {
                id: "m1".into(),
                unit: "unit_ops.mixer".into(),
                parameters: BTreeMap::new(),
            }],
            connections: vec![
                Connection {
                    from: "a".into(),
                    to: "m1.feed".into(),
                },
                Connection {
                    from: "b".into(),
                    to: "m1.feed".into(),
                },
                Connection {
                    from: "m1.product".into(),
                    to: "out".into(),
                },
            ],
            recycles: vec![],
        };
        assert_clean(&flowsheet, &palette);
    }

    /// A feed into a splitter and its first product out, so an index has somewhere to be written.
    fn split_flowsheet() -> Flowsheet {
        Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("in")],
            products: vec!["out".into()],
            instances: vec![Instance {
                id: "s1".into(),
                unit: "unit_ops.splitter".into(),
                parameters: BTreeMap::new(),
            }],
            connections: vec![
                Connection {
                    from: "in".into(),
                    to: "s1.feed".into(),
                },
                Connection {
                    from: "s1.products[0]".into(),
                    to: "out".into(),
                },
            ],
            recycles: vec![],
        }
    }

    /// **The three ways to write a position a port does not have**, and they are one mistake made
    /// against three different declarations. A dropped index would bind the *wrong* stream
    /// silently, which is the failure the `many` outlet existed to avoid - so all three are
    /// refused rather than ignored.
    #[test]
    fn an_index_is_checked_against_the_port_it_names() {
        let palette = vec![splitter("unit_ops.splitter")];

        // A position on a `many` outlet is what the grammar is *for*, so this is clean.
        assert_clean(&split_flowsheet(), &palette);

        // A `one` outlet has no second stream.
        let two_port_palette = vec![two_port("unit_ops.pump")];
        let mut flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            inputs: vec![named_input("in")],
            products: vec!["out".into()],
            instances: vec![Instance {
                id: "p1".into(),
                unit: "unit_ops.pump".into(),
                parameters: BTreeMap::new(),
            }],
            connections: vec![
                Connection {
                    from: "in".into(),
                    to: "p1.feed".into(),
                },
                Connection {
                    from: "p1.discharge[0]".into(),
                    to: "out".into(),
                },
            ],
            recycles: vec![],
        };
        let diags = validate(&flowsheet, &two_port_palette);
        assert!(
            diags.iter().any(|d| matches!(
                d,
                Diagnostic::EndpointIndex { endpoint, detail }
                    if endpoint == "p1.discharge[0]" && detail.contains("declares one stream")
            )),
            "{diags:?}"
        );

        // An inlet is one name however many connections reach it, so it has no positions. The
        // resolver drops the index, so the port check still passes and this is the only report.
        flowsheet.instances[0].unit = "unit_ops.splitter".into();
        flowsheet.connections[0].to = "p1.feed[0]".into();
        flowsheet.connections[1].from = "p1.products".into();
        let diags = validate(&flowsheet, &palette);
        assert!(
            diags.iter().any(|d| matches!(
                d,
                Diagnostic::EndpointIndex { endpoint, detail }
                    if endpoint == "p1.feed[0]" && detail.contains("no positions")
            )),
            "{diags:?}"
        );
    }

    /// The grammar itself, which three modules used to spell for themselves.
    #[test]
    fn the_endpoint_grammar_is_one_place() {
        assert_eq!(split_endpoint("feed_1"), None);
        assert_eq!(split_endpoint("p1.outlet"), Some(("p1", "outlet", None)));
        assert_eq!(
            split_endpoint("s1.products[3]"),
            Some(("s1", "products", Some(3)))
        );
        // A malformed index is not an endpoint with an index - it is a name the resolver reads as
        // a boundary, which is how `UnknownFeed` reaches it rather than a parse error nobody sees.
        assert_eq!(split_endpoint("s1.products[three]"), None);
        assert_eq!(split_endpoint("s1.products[1"), None);

        assert_eq!(
            crate::flowsheet::stream_path("p1", "outlet", None),
            "p1.outlet"
        );
        assert_eq!(
            crate::flowsheet::stream_path("s1", "products", Some(2)),
            "s1.products[2]"
        );
    }
}
