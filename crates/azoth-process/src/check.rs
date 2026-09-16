//! The checker: a flowsheet is held to the calculus's rules, made mechanical.
//!
//! `docs/src/calculus/process.md` states a unit operation as a process on typed,
//! directional channels, conservation as linearity, and feedback as restriction.
//! This module turns those into checks a machine can run:
//!
//! - every instance names a palette unit op, and only its declared parameters;
//! - a connection joins an outlet (or feed) to an inlet (or product);
//! - a Port-to-Port connection joins dimension-compatible field records;
//! - **linearity** — a `one` port is consumed/produced exactly once, a `many` port
//!   at least once, and every feed/product is used exactly once;
//! - every loop in the instance graph passes through a declared recycle.

use std::collections::{HashMap, HashSet};

use azoth_core::unit_vocab_gen::dimension_exponents;

use crate::channel::{ChannelType, Direction, FieldType, Multiplicity};
use crate::flowsheet::{Connection, Flowsheet, Instance, Recycle};
use crate::unit_op::UnitOpSpec;

/// What a `validate` run found wrong.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Diagnostic {
    UnknownDimension {
        dimension: String,
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
    UnrecycledLoop {
        detail: String,
    },
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
            for field in port.fields.values() {
                if dimension_exponents(&field.dimension).is_none() {
                    diags.push(Diagnostic::UnknownDimension {
                        dimension: field.dimension.clone(),
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
    for feed in &flowsheet.feeds {
        if !feed_names.insert(feed) {
            diags.push(Diagnostic::DuplicateFeed { name: feed.clone() });
        }
        if instance_ids.contains(feed.as_str()) {
            diags.push(Diagnostic::NameCollision { name: feed.clone() });
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
            Endpoint::Port { instance, port } => {
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
            Endpoint::Port { instance, port } => {
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
            },
            Endpoint::Port {
                instance: to_i,
                port: to_p,
            },
        ) = (producer, consumer)
        {
            if let (Some(from_t), Some(to_t)) = (
                channel_type(&instances, &specs, from_i, from_p),
                channel_type(&instances, &specs, to_i, to_p),
            ) {
                if let Err(detail) = compatible(from_t, to_t) {
                    diags.push(Diagnostic::TypeMismatch { detail });
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
    for feed in &flowsheet.feeds {
        if feed_uses.get(feed.as_str()).copied().unwrap_or(0) != 1 {
            diags.push(Diagnostic::UnusedFeed { name: feed.clone() });
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
    Port { instance: &'a str, port: &'a str },
}

fn resolve_producer(s: &str) -> Endpoint<'_> {
    match s.split_once('.') {
        Some((instance, port)) => Endpoint::Port { instance, port },
        None => Endpoint::Feed(s),
    }
}

fn resolve_consumer(s: &str) -> Endpoint<'_> {
    match s.split_once('.') {
        Some((instance, port)) => Endpoint::Port { instance, port },
        None => Endpoint::Product(s),
    }
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
            feeds: vec!["feed_1".into()],
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

    #[test]
    fn an_unknown_unit_op_is_reported() {
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            feeds: vec![],
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
            feeds: vec![],
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
            feeds: vec!["feed_1".into()],
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
            feeds: vec!["feed_1".into()],
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
            feeds: vec![],
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
            feeds: vec!["a".into(), "b".into()],
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
            feeds: vec![],
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
            feeds: vec![],
            products: vec![],
            instances,
            connections: vec![Connection {
                from: "p1.discharge".into(),
                to: "p2.feed".into(),
            }],
            recycles: vec![Recycle {
                stream: "r1".into(),
                from: "p2.discharge".into(),
                to: "p1.feed".into(),
            }],
        };
        assert!(!has_variant(&validate(&torn, &palette), |d| matches!(
            d,
            Diagnostic::UnrecycledLoop { .. }
        )));
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
        assert!(has(
            &diags,
            Diagnostic::UnknownDimension {
                dimension: "no_such_dimension".into()
            }
        ));
    }

    #[test]
    fn a_mixer_accepts_many_feeds() {
        let palette = vec![mixer("unit_ops.mixer")];
        let flowsheet = Flowsheet {
            id: "f".into(),
            name: "f".into(),
            feeds: vec!["a".into(), "b".into()],
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
}
