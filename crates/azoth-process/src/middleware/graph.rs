//! The connection graph, as the editor's own document.
//!
//! A flowsheet's connection graph is already a node/edge set, so the projection uses **xyflow's
//! schema** rather than one invented beside it: `id`, `type`, `position`, `data` on a node,
//! `source`/`target` with `sourceHandle`/`targetHandle` on an edge. The editor's model is then
//! the schema's rather than a translation of it.
//!
//! **Two rules make the join mechanical.** A node's id is `{role}:{name}`, which is
//! [`NodeRole::node_id`], and a **handle's id is the stream's own path** - `feed_1`,
//! `sep1.liquid`, `s1.products[0]` - so an edge's handle, the session's key for the same stream
//! and the edge's `data.path` are one string rather than three that have to be reconciled.
//!
//! **The projection never refuses a document.** A canvas has to draw a broken flowsheet, because
//! drawing it is how a user fixes it: an instance naming a unit op the palette does not carry
//! gets a node with no ports, and an edge naming an instance nobody declared gets the node id it
//! implies and no node to attach to - which xyflow reports, and which the checker's
//! `UnknownInstance` already explains. Nothing is repaired and nothing is dropped.

use std::collections::HashMap;

use azoth_core::{AzothError, Result};
use serde::Serialize;

use crate::NodeRole;
use crate::channel::{Direction, Multiplicity, Port};
use crate::flowsheet::{Flowsheet, Instance, Recycle, split_endpoint, stream_path};
use crate::model_inputs_gen;
use crate::unit_op::UnitOpSpec;

/// The whole graph.
#[derive(Debug, Clone, Serialize)]
pub struct Graph {
    pub id: String,
    pub name: String,
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

/// One node: an instance, a feed or a product.
#[derive(Debug, Clone, Serialize)]
pub struct Node {
    /// `{role}:{name}`, from [`NodeRole::node_id`], and what an edge's `source`/`target` names.
    pub id: String,
    /// The component xyflow renders it with: `unit_op` or `stream`.
    #[serde(rename = "type")]
    pub node_type: &'static str,
    /// Which of the three kinds it is: `instance`, `input` or `product`.
    pub role: &'static str,
    pub position: Position,
    pub data: NodeData,
}

/// Where the canvas puts a node.
#[derive(Debug, Clone, Copy, Serialize)]
pub struct Position {
    pub x: f64,
    pub y: f64,
}

/// What a node carries: its own name, and whatever its kind has to say.
#[derive(Debug, Clone, Serialize)]
pub struct NodeData {
    /// The instance id, the feed's name or the product's name.
    pub name: String,
    /// An instance's unit-op id, e.g. `unit_ops.separator`.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit: Option<String>,
    /// The palette entry's own name, e.g. `Separator`, which is what a title bar shows.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub unit_name: Option<String>,
    /// The instance's parameter values, as the document holds them.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub parameters: Option<serde_json::Value>,
    /// The instance's ports, with the handle ids an edge attaches to.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub ports: Option<Ports>,
    /// A feed's declared record.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub input: Option<InputRecord>,
}

/// The two sides of an instance's declaration.
#[derive(Debug, Clone, Serialize)]
pub struct Ports {
    pub inlets: Vec<Handle>,
    pub outlets: Vec<Handle>,
}

/// One port, and the handle ids a stream on it is addressed by.
///
/// A `one` port has exactly one handle, named by the port; a `many` **outlet** has one per
/// stream it returns, named by the port and the position; a `many` **inlet** has one, however
/// many connections reach it - which is how a mixer is wired.
#[derive(Debug, Clone, Serialize)]
pub struct Handle {
    pub name: String,
    pub multiplicity: &'static str,
    pub handles: Vec<String>,
}

/// A feed's record, in the units the schema declares rather than as `Quantity` pairs.
///
/// **Bare numbers, because the record's own units are the schema's.** `n` is mol/s, `P` is Pa and
/// `T` is K in the document, so what a form writes back is what it read with no conversion in
/// between - which is what keeps the round trip through the editor exact.
#[derive(Debug, Clone, Serialize)]
pub struct InputRecord {
    pub components: Vec<String>,
    pub n: f64,
    pub z: Vec<f64>,
    #[serde(rename = "P")]
    pub p: f64,
    #[serde(rename = "T")]
    pub t: f64,
}

/// One `[[recycles]]` entry's seven settings, as the document states them.
///
/// **`null` is a setting the declaration leaves to the class's own default**, which is not the
/// same fact as the default: `demo.toml` states none of the seven, and a widget that showed the
/// defaults and wrote them back would turn a silence into a pinned number.
#[derive(Debug, Clone, Serialize)]
pub struct Settings {
    pub flow_tolerance: Option<f64>,
    pub composition_tolerance: Option<f64>,
    pub temperature_tolerance: Option<f64>,
    pub pressure_tolerance: Option<f64>,
    pub max_iterations: Option<u32>,
    pub minimum_flow: Option<f64>,
    pub acceleration_method: Option<String>,
}

impl From<&Recycle> for Settings {
    fn from(recycle: &Recycle) -> Self {
        Self {
            flow_tolerance: recycle.flow_tolerance,
            composition_tolerance: recycle.composition_tolerance,
            temperature_tolerance: recycle.temperature_tolerance,
            pressure_tolerance: recycle.pressure_tolerance,
            max_iterations: recycle.max_iterations,
            minimum_flow: recycle.minimum_flow,
            acceleration_method: recycle.acceleration_method.clone(),
        }
    }
}

/// One edge: a connection, or a recycle that tears a loop.
#[derive(Debug, Clone, Serialize)]
pub struct Edge {
    /// `e{index}`, over the connections then the recycles in document order.
    ///
    /// **A `[[connections]]` entry has no id in the schema**, so its position is its identity.
    /// The checker refuses two entries that join the same pair (`DuplicateConnection`), so a
    /// document that passes has one edge per pair and the pair is also a name a command can use.
    pub id: String,
    pub source: String,
    /// The stream's path at the producing end, which is also a session key.
    #[serde(rename = "sourceHandle")]
    pub source_handle: String,
    pub target: String,
    /// The port's path at the consuming end.
    #[serde(rename = "targetHandle")]
    pub target_handle: String,
    pub data: EdgeData,
}

#[derive(Debug, Clone, Serialize)]
pub struct EdgeData {
    /// `connection` or `recycle`.
    pub kind: &'static str,
    /// The endpoints as the document wrote them, so an edge whose node is missing still says
    /// what it referred to.
    pub from: String,
    pub to: String,
    /// The session path of the stream this edge carries: the producer's path for a connection,
    /// and the tear's name for a recycle.
    pub path: String,
    /// A tear's seven convergence settings; absent on a connection, which has none.
    ///
    /// **On the edge because a tear *is* a connection with a name.** A panel for the selected edge
    /// needs the document's own values, and reading them a second time in a front-end would be a
    /// second reader of the format.
    #[serde(skip_serializing_if = "Option::is_none")]
    pub settings: Option<Settings>,
}

/// Project a flowsheet onto the graph a canvas draws.
///
/// # Errors
/// [`AzothError::InvalidInput`] if a parameter's value cannot be written as JSON, which TOML's
/// datetime is and JSON has no shape for. A document carrying one is refused by the checker as a
/// `ParameterKind` too, so this is a refusal of the same value rather than a second opinion.
pub fn graph(flowsheet: &Flowsheet, palette: &[UnitOpSpec]) -> Result<Graph> {
    let specs: HashMap<&str, &UnitOpSpec> = palette
        .iter()
        .map(|spec| (spec.id.as_str(), spec))
        .collect();
    let wired = highest_wired_position(flowsheet);
    let positions = positions(flowsheet);

    let mut nodes = Vec::with_capacity(
        flowsheet.instances.len() + flowsheet.inputs.len() + flowsheet.products.len(),
    );
    for instance in &flowsheet.instances {
        let spec = specs.get(instance.unit.as_str()).copied();
        let role = NodeRole::Instance;
        nodes.push(Node {
            id: role.node_id(&instance.id),
            node_type: "unit_op",
            role: role.name(),
            position: positions
                .get(&instance.id)
                .copied()
                .unwrap_or(Position { x: 0.0, y: 0.0 }),
            data: instance_data(instance, spec, &wired)?,
        });
    }
    for input in &flowsheet.inputs {
        let role = NodeRole::Input;
        nodes.push(Node {
            id: role.node_id(&input.name),
            node_type: "stream",
            role: role.name(),
            position: positions
                .get(&input.name)
                .copied()
                .unwrap_or(Position { x: 0.0, y: 0.0 }),
            data: NodeData {
                name: input.name.clone(),
                unit: None,
                unit_name: None,
                parameters: None,
                ports: None,
                input: Some(InputRecord {
                    components: input.components.clone(),
                    n: input.n,
                    z: input.z.clone(),
                    p: input.p,
                    t: input.t,
                }),
            },
        });
    }
    for product in &flowsheet.products {
        let role = NodeRole::Product;
        nodes.push(Node {
            id: role.node_id(product),
            node_type: "stream",
            role: role.name(),
            position: positions
                .get(product)
                .copied()
                .unwrap_or(Position { x: 0.0, y: 0.0 }),
            data: NodeData {
                name: product.clone(),
                unit: None,
                unit_name: None,
                parameters: None,
                ports: None,
                input: None,
            },
        });
    }

    let mut edges = Vec::with_capacity(flowsheet.connections.len() + flowsheet.recycles.len());
    for (index, connection) in flowsheet.connections.iter().enumerate() {
        edges.push(edge(
            format!("e{index}"),
            "connection",
            None,
            &connection.from,
            &connection.to,
        ));
    }
    for (index, recycle) in flowsheet.recycles.iter().enumerate() {
        edges.push(edge(
            format!("e{}", flowsheet.connections.len() + index),
            "recycle",
            Some(recycle),
            &recycle.from,
            &recycle.to,
        ));
    }

    Ok(Graph {
        id: flowsheet.id.clone(),
        name: flowsheet.name.clone(),
        nodes,
        edges,
    })
}

/// One edge, from the two endpoints as the document wrote them.
///
/// The handle ids are the stream paths the endpoints name: the producer's path at the source, and
/// the port's path at the target - an index belongs to the producing side only, because a `many`
/// inlet is one name however many streams reach it.
fn edge(id: String, kind: &'static str, tear: Option<&Recycle>, from: &str, to: &str) -> Edge {
    // A bare name is a boundary: a feed where the stream is produced, a product where it is
    // consumed. An endpoint with a dot is an instance's port, whatever the role of the side.
    let source_handle = match split_endpoint(from) {
        Some((instance, port, index)) => stream_path(instance, port, index),
        None => from.to_string(),
    };
    let target_handle = match split_endpoint(to) {
        Some((instance, port, _)) => stream_path(instance, port, None),
        None => to.to_string(),
    };
    let (source_role, source_name) = match split_endpoint(from) {
        Some((instance, _, _)) => (NodeRole::Instance, instance),
        None => (NodeRole::Input, from),
    };
    let (target_role, target_name) = match split_endpoint(to) {
        Some((instance, _, _)) => (NodeRole::Instance, instance),
        None => (NodeRole::Product, to),
    };
    Edge {
        id,
        source: source_role.node_id(source_name),
        source_handle: source_handle.clone(),
        target: target_role.node_id(target_name),
        target_handle,
        data: EdgeData {
            kind,
            from: from.to_string(),
            to: to.to_string(),
            path: tear.map_or(source_handle, |recycle| recycle.stream.clone()),
            settings: tear.map(Settings::from),
        },
    }
}

fn instance_data(
    instance: &Instance,
    spec: Option<&UnitOpSpec>,
    wired: &HashMap<(String, String), usize>,
) -> Result<NodeData> {
    let parameters = serde_json::to_value(&instance.parameters).map_err(|error| {
        AzothError::invalid_input(
            &instance.id,
            format!(
                "`{}`'s parameters cannot be written as JSON: {error}",
                instance.id
            ),
        )
    })?;
    let Some(spec) = spec else {
        // A unit op the palette does not carry: the node is drawn with no ports, and the
        // checker's `UnknownUnitOp` is what says so.
        return Ok(NodeData {
            name: instance.id.clone(),
            unit: Some(instance.unit.clone()),
            unit_name: None,
            parameters: Some(parameters),
            ports: None,
            input: None,
        });
    };

    let declared = declared_positions(spec, instance);
    let port_handle = |port: &Port| {
        let multiplicity = match port.multiplicity {
            Multiplicity::One => "one",
            Multiplicity::Many => "many",
        };
        let handles = match (port.direction, port.multiplicity) {
            (Direction::Out, Multiplicity::Many) => {
                let wired = wired
                    .get(&(instance.id.clone(), port.name.clone()))
                    .copied()
                    .unwrap_or(0);
                // **How many streams a `many` outlet returns is a value and not a declaration.**
                // The checker says so itself. So the handles are the union of what the document
                // addresses and what the instance's own vector says: a document mid-edit names a
                // position the vector does not have yet, and a vector that grew has positions no
                // connection names yet. Drawing both is what lets either be fixed.
                let count = wired.max(declared.unwrap_or(0)).max(1);
                (0..count)
                    .map(|index| stream_path(&instance.id, &port.name, Some(index)))
                    .collect()
            }
            (Direction::Out, Multiplicity::One) => {
                vec![stream_path(&instance.id, &port.name, None)]
            }
            // A `many` inlet is one name, however many connections reach it.
            (Direction::In, _) => vec![stream_path(&instance.id, &port.name, None)],
        };
        Handle {
            name: port.name.clone(),
            multiplicity,
            handles,
        }
    };

    Ok(NodeData {
        name: instance.id.clone(),
        unit: Some(instance.unit.clone()),
        unit_name: Some(spec.name.clone()),
        parameters: Some(parameters),
        ports: Some(Ports {
            inlets: spec
                .ports
                .iter()
                .filter(|port| port.direction == Direction::In)
                .map(port_handle)
                .collect(),
            outlets: spec
                .ports
                .iter()
                .filter(|port| port.direction == Direction::Out)
                .map(port_handle)
                .collect(),
        }),
        input: None,
    })
}

/// How many streams a `many` outlet returns, where the instance's own parameters say.
///
/// **The entry's one vector parameter, and only when there is exactly one.** Measured over the
/// shipped palette: the three entries declaring a `many` outlet - `splitter`, `component_splitter`
/// and `manifold` - each declare exactly one vector parameter, `split_factors`, and the count is
/// its length. An entry with two vectors would be ambiguous, so the count is then unknown and the
/// wiring alone decides.
fn declared_positions(spec: &UnitOpSpec, instance: &Instance) -> Option<usize> {
    let vectors: Vec<&str> = model_inputs_gen::inputs_for(&spec.id)?
        .inputs
        .iter()
        .filter(|input| input.kind == "vector")
        .map(|input| input.name)
        .filter(|name| spec.parameters.contains_key(*name))
        .collect();
    let [vector] = vectors.as_slice() else {
        return None;
    };
    instance
        .parameters
        .get(*vector)
        .and_then(toml::Value::as_array)
        .map(Vec::len)
}

/// The highest position a connection names out of each `many` outlet, plus one.
fn highest_wired_position(flowsheet: &Flowsheet) -> HashMap<(String, String), usize> {
    let mut wired: HashMap<(String, String), usize> = HashMap::new();
    let endpoints = flowsheet
        .connections
        .iter()
        .map(|connection| (connection.from.as_str(), connection.to.as_str()))
        .chain(
            flowsheet
                .recycles
                .iter()
                .map(|recycle| (recycle.from.as_str(), recycle.to.as_str())),
        );
    for (from, to) in endpoints {
        for text in [from, to] {
            if let Some((instance, port, Some(index))) = split_endpoint(text) {
                let entry = wired
                    .entry((instance.to_string(), port.to_string()))
                    .or_insert(0);
                *entry = (*entry).max(index + 1);
            }
        }
    }
    wired
}

/// Where every node goes: the document's own layout, else a deterministic layering.
///
/// **The fallback is what a hand-written flowsheet gets**, and it is stable across a reload - a
/// column per step of the longest path, a row per node within it, feeds left and products right.
/// A node the document places keeps its own position, so a document that has been arranged draws
/// as arranged and one that has never been opened draws as something.
fn positions(flowsheet: &Flowsheet) -> HashMap<String, Position> {
    const COLUMN: f64 = 240.0;
    const ROW: f64 = 120.0;

    let index: HashMap<&str, usize> = flowsheet
        .instances
        .iter()
        .enumerate()
        .map(|(i, instance)| (instance.id.as_str(), i))
        .collect();

    // Longest path from a boundary, over the instance-to-instance edges the order walk uses -
    // the connections less the declared recycles. Bounded by the instance count, so a cycle that
    // no recycle tears stops growing rather than diverging.
    let mut depth = vec![0_usize; flowsheet.instances.len()];
    let edges: Vec<(usize, usize)> = flowsheet
        .connections
        .iter()
        .filter_map(|connection| {
            let from = split_endpoint(&connection.from);
            let to = split_endpoint(&connection.to);
            match (from, to) {
                (Some((from, _, _)), Some((to, _, _))) => Some((from, to)),
                _ => None,
            }
        })
        .filter_map(|(from, to)| Some((*index.get(from)?, *index.get(to)?)))
        .collect();
    let torn: Vec<(String, String)> = flowsheet
        .recycles
        .iter()
        .filter_map(|recycle| {
            let from = split_endpoint(&recycle.from);
            let to = split_endpoint(&recycle.to);
            match (from, to) {
                (Some((from, _, _)), Some((to, _, _))) => Some((from.to_string(), to.to_string())),
                _ => None,
            }
        })
        .collect();
    for _ in 0..flowsheet.instances.len() {
        let mut moved = false;
        for (source, target) in &edges {
            let (from, to) = (
                &flowsheet.instances[*source].id,
                &flowsheet.instances[*target].id,
            );
            if torn.iter().any(|(a, b)| a == from && b == to) {
                continue;
            }
            if depth[*target] < depth[*source] + 1 {
                depth[*target] = depth[*source] + 1;
                moved = true;
            }
        }
        if !moved {
            break;
        }
    }

    let mut placed = HashMap::new();
    let mut rows: HashMap<usize, usize> = HashMap::new();

    for (index, input) in flowsheet.inputs.iter().enumerate() {
        placed.insert(
            input.name.clone(),
            Position {
                x: 0.0,
                y: ROW * index as f64,
            },
        );
    }
    for (index, instance) in flowsheet.instances.iter().enumerate() {
        let column = depth[index] + 1;
        let row = {
            let next = rows.entry(column).or_insert(0);
            let this = *next;
            *next += 1;
            this
        };
        placed.insert(
            instance.id.clone(),
            Position {
                x: COLUMN * column as f64,
                y: ROW * row as f64,
            },
        );
    }
    let last = depth.iter().max().copied().unwrap_or(0) + 2;
    for (index, product) in flowsheet.products.iter().enumerate() {
        placed.insert(
            product.clone(),
            Position {
                x: COLUMN * last as f64,
                y: ROW * index as f64,
            },
        );
    }

    // The document's own placements win, node by node.
    if let Some(layout) = &flowsheet.layout {
        for (name, [x, y]) in layout
            .instances
            .iter()
            .chain(layout.inputs.iter())
            .chain(layout.products.iter())
        {
            placed.insert(name.clone(), Position { x: *x, y: *y });
        }
    }
    placed
}
