//! A P&ID as a flowsheet scaffold.
//!
//! **A P&ID cannot become a runnable flowsheet, and this does not pretend otherwise.** It has
//! no composition and no process state. What an import proves is the drawing's topology: which
//! tagged item each stream leaves, by which nozzle, on which line, and where it goes. Feeds are
//! left undeclared and no parameter is invented, so `azoth check` reports the scaffold.
//!
//! **A loop with no marker is inferred, and said so.** DEXPI has no recycle object; a drawing's
//! loop is a ring of segments. The importer finds a cycle, chooses the lexicographically first
//! edge in it as the tear - deterministic, not meaningful - and reports `recycle_inferred`. An
//! azoth export records the real tear in `Azoth/Recycle`, which is read instead.

use std::collections::{BTreeMap, BTreeSet};

use azoth_process::flowsheet::{
    Connection, Flowsheet, Instance, Instrument, InstrumentRange, Recycle, SafetyFunction, Signal,
};
use azoth_process::{Direction, UnitOpSpec};

use crate::diagnostic::{DexpiDiagnostic, Severity};
use crate::extensions;
use crate::ids::{model_name, sanitize_instance_id, uniquify};
use crate::model::{Child, ComponentChild, DataValue, Document, Object};
use crate::plant::{self, DEFAULT_LINE, PASS_THROUGH_CLASSES};
use crate::reader;
use crate::validate::validate;

/// How strict an import is.
#[derive(Debug, Clone, Copy, Default)]
pub struct PlantImportOptions {
    /// Turn every warning into a refusal.
    pub strict: bool,
}

/// The flowsheet an exchange became, and what it cost.
#[derive(Debug, Clone)]
pub struct PlantImportReport {
    pub flowsheet: Flowsheet,
    pub diagnostics: Vec<DexpiDiagnostic>,
    pub document: Document,
}

/// Read a DEXPI Plant exchange as a flowsheet scaffold.
///
/// # Errors
/// The diagnostics that refused the import: a document that will not read, a document that
/// fails the semantic rules, or - under [`PlantImportOptions::strict`] - any warning.
pub fn import_plant(
    xml: &str,
    palette: &[UnitOpSpec],
    options: &PlantImportOptions,
) -> Result<PlantImportReport, Vec<DexpiDiagnostic>> {
    let read = reader::read(xml).map_err(|diagnostic| vec![diagnostic])?;
    import_plant_document(read.document, palette, options)
}

/// The same import over a document already parsed - the door a projection uses.
///
/// # Errors
/// The same diagnostics [`import_plant`] reports.
pub fn import_plant_document(
    document: Document,
    palette: &[UnitOpSpec],
    options: &PlantImportOptions,
) -> Result<PlantImportReport, Vec<DexpiDiagnostic>> {
    let validation = validate(&document);
    if !validation.is_valid() {
        return Err(validation.diagnostics);
    }
    let plant_model = find_plant_model(&document).ok_or_else(|| {
        vec![DexpiDiagnostic::error(
            "missing_plant_model",
            "",
            "the document holds no `Plant/PlantModel`",
        )]
    })?;

    let mut diagnostics = Vec::new();
    let mut taken = BTreeSet::new();
    let steps = read_steps(plant_model, &mut diagnostics, &mut taken);

    let mut owners: BTreeMap<String, NodeOwner> = BTreeMap::new();
    for (index, step) in steps.iter().enumerate() {
        for nozzle in &step.nozzles {
            for node in &nozzle.nodes {
                owners.insert(
                    node.clone(),
                    NodeOwner {
                        step: Some(index),
                        port: nozzle.port.clone(),
                        boundary: None,
                    },
                );
            }
        }
    }
    let (feeds, products) = collect_boundaries(&document, &mut owners);
    let segments = read_segments(plant_model);

    let mut connections = Vec::new();
    let mut marked_recycles: Vec<(String, String, String)> = Vec::new();
    for segment in &segments {
        let Some(source) = owners.get(&segment.source_node) else {
            diagnostics.push(DexpiDiagnostic::warning(
                "unresolved_segment",
                segment.id.clone(),
                "the segment's source node belongs to no item, nozzle or boundary connector - \
                 an inline pass-through fitting is the usual cause, and this port draws none",
            ));
            continue;
        };
        let Some(target) = owners.get(&segment.target_node) else {
            diagnostics.push(DexpiDiagnostic::warning(
                "unresolved_segment",
                segment.id.clone(),
                "the segment's target node belongs to no item, nozzle or boundary connector",
            ));
            continue;
        };
        let Some(from) = endpoint(source, true, &steps, palette, &mut diagnostics) else {
            diagnostics.push(DexpiDiagnostic::warning(
                "connection_dropped",
                segment.id.clone(),
                "the segment's source item is not one this port models",
            ));
            continue;
        };
        let Some(to) = endpoint(target, false, &steps, palette, &mut diagnostics) else {
            diagnostics.push(DexpiDiagnostic::warning(
                "connection_dropped",
                segment.id.clone(),
                "the segment's target item is not one this port models",
            ));
            continue;
        };
        if let Some(recycle) = segment.recycle.clone() {
            marked_recycles.push((recycle, from, to));
            continue;
        }
        connections.push(Connection {
            from,
            to,
            line: segment.line.clone(),
        });
    }

    // **The loops no marker names are found, not assumed.** One tear per cycle - the
    // lexicographically first edge in it - until the instance graph is acyclic.
    let mut recycles: Vec<Recycle> = marked_recycles
        .into_iter()
        .map(|(name, from, to)| Recycle::new(&name, &from, &to))
        .collect();
    infer_recycles(&mut connections, &mut recycles, &mut diagnostics);

    for feed in &feeds {
        diagnostics.push(DexpiDiagnostic::warning(
            "input_composition_missing",
            feed.clone(),
            "DEXPI carries no composition, so this feed is left undeclared: add an [[inputs]] \
             record with components, z, n, P and T, and `azoth check` stops reporting \
             `unknown_feed`",
        ));
    }

    let instances = steps
        .iter()
        .filter_map(|step| match &step.unit {
            Some(unit) => Some(Instance {
                id: step.instance_id.clone()?,
                unit: unit.clone(),
                parameters: BTreeMap::new(),
            }),
            None => None,
        })
        .collect();

    let root = find_engineering_model(&document);
    let item_tags = item_tags(plant_model);
    let instrumentation = read_instruments(plant_model, &item_tags, &mut diagnostics);

    let flowsheet = Flowsheet {
        id: root
            .and_then(|root| root.data_string(extensions::FLOWSHEET_ID))
            .map_or_else(|| model_name(&document.model.name), |id| id.to_string()),
        name: root
            .and_then(|root| root.data_string(extensions::FLOWSHEET_NAME))
            .map_or_else(|| document.model.name.clone(), |name| name.to_string()),
        products: products.into_iter().collect(),
        inputs: Vec::new(),
        instances,
        connections,
        recycles,
        instruments: instrumentation.instruments,
        signals: instrumentation.signals,
        layout: None,
    };

    if options.strict {
        for diagnostic in &mut diagnostics {
            if diagnostic.severity == Severity::Warning {
                diagnostic.severity = Severity::Error;
            }
        }
    }
    if diagnostics.iter().any(DexpiDiagnostic::is_error) {
        return Err(diagnostics);
    }

    Ok(PlantImportReport {
        flowsheet,
        diagnostics,
        document,
    })
}

struct Step {
    identifier: String,
    unit: Option<String>,
    instance_id: Option<String>,
    nozzles: Vec<Nozzle>,
}

struct Nozzle {
    nodes: Vec<String>,
    port: Option<String>,
}

struct NodeOwner {
    step: Option<usize>,
    port: Option<String>,
    boundary: Option<String>,
}

struct Segment {
    id: String,
    source_node: String,
    target_node: String,
    line: Option<String>,
    recycle: Option<String>,
}

fn read_steps(
    plant_model: &Object,
    diagnostics: &mut Vec<DexpiDiagnostic>,
    taken: &mut BTreeSet<String>,
) -> Vec<Step> {
    let mut steps = Vec::new();
    let Some(items) = plant_model.components("TaggedPlantItems") else {
        return steps;
    };
    for item in items.children.iter().filter_map(|child| match child {
        ComponentChild::Object(object) => Some(object),
        _ => None,
    }) {
        let class = plant::class_of(&item.ty);
        if PASS_THROUGH_CLASSES.contains(&class) {
            continue;
        }
        let identifier = item
            .data_string("TagName")
            .or(item.id.as_deref())
            .unwrap_or_default()
            .to_string();
        let unit = match plant::for_class(class) {
            Some(unit) => Some(unit.to_string()),
            None => {
                diagnostics.push(DexpiDiagnostic::warning(
                    "unsupported_plant_class",
                    identifier.clone(),
                    format!(
                        "`{class}` has no reviewed azoth unit; the item is kept as \
                         `dexpi.unmapped`"
                    ),
                ));
                Some("dexpi.unmapped".to_string())
            }
        };
        let sanitized = sanitize_instance_id(&identifier);
        let instance_id = uniquify(&sanitized, taken);
        if instance_id != sanitized {
            diagnostics.push(DexpiDiagnostic::warning(
                "renamed_on_import",
                identifier.clone(),
                format!(
                    "`{identifier}` collides with another name; this instance is `{instance_id}`"
                ),
            ));
        }
        let mut nozzles = Vec::new();
        if let Some(component) = item.components("Nozzles") {
            for nozzle in component.children.iter().filter_map(|child| match child {
                ComponentChild::Object(object) => Some(object),
                _ => None,
            }) {
                nozzles.push(Nozzle {
                    nodes: nozzle_node_ids(nozzle),
                    port: nozzle.data_string("SubTagName").map(ToString::to_string),
                });
            }
        }
        steps.push(Step {
            identifier,
            unit,
            instance_id: Some(instance_id),
            nozzles,
        });
    }
    steps
}

fn nozzle_node_ids(nozzle: &Object) -> Vec<String> {
    nozzle
        .components("Nodes")
        .map(|component| {
            component
                .children
                .iter()
                .filter_map(|child| match child {
                    ComponentChild::Object(object) => object.id.clone(),
                    ComponentChild::Reference { .. } => None,
                })
                .collect()
        })
        .unwrap_or_default()
}

/// The off-page connectors, wherever a segment holds one, split into feeds and products.
fn collect_boundaries(
    document: &Document,
    owners: &mut BTreeMap<String, NodeOwner>,
) -> (BTreeSet<String>, BTreeSet<String>) {
    let mut feeds = BTreeSet::new();
    let mut products = BTreeSet::new();
    walk(&document.model.children, &mut |object| {
        let inbound = object.ty == "Plant/Piping.FlowInPipeOffPageConnector";
        let outbound = object.ty == "Plant/Piping.FlowOutPipeOffPageConnector";
        if !inbound && !outbound {
            return;
        }
        let name = object
            .data_string("PipeConnectorNumber")
            .or(object.id.as_deref())
            .unwrap_or_default()
            .to_string();
        let nodes = object
            .components("Nodes")
            .map(|component| {
                component
                    .children
                    .iter()
                    .filter_map(|child| match child {
                        ComponentChild::Object(node) => node.id.clone(),
                        ComponentChild::Reference { .. } => None,
                    })
                    .collect::<Vec<_>>()
            })
            .unwrap_or_default();
        for node in nodes {
            owners.insert(
                node,
                NodeOwner {
                    step: None,
                    port: None,
                    boundary: Some(name.clone()),
                },
            );
        }
        if inbound {
            feeds.insert(name);
        } else {
            products.insert(name);
        }
    });
    (feeds, products)
}

fn read_segments(plant_model: &Object) -> Vec<Segment> {
    let mut segments = Vec::new();
    let Some(component) = plant_model.components("PipingNetworkSystems") else {
        return segments;
    };
    for system in component.children.iter().filter_map(|child| match child {
        ComponentChild::Object(object) => Some(object),
        _ => None,
    }) {
        // The placeholder line is azoth's own; it reads back as no line at all.
        let line = system
            .data_string("LineNumber")
            .filter(|line| *line != DEFAULT_LINE)
            .map(ToString::to_string);
        let Some(segments_component) = system.components("Segments") else {
            continue;
        };
        for segment in segments_component
            .children
            .iter()
            .filter_map(|child| match child {
                ComponentChild::Object(object) => Some(object),
                _ => None,
            })
        {
            let Some(source_node) = first_reference(segment, "SourceNode") else {
                continue;
            };
            let Some(target_node) = first_reference(segment, "TargetNode") else {
                continue;
            };
            segments.push(Segment {
                id: segment.id.clone().unwrap_or_default(),
                source_node: source_node.trim_start_matches('#').to_string(),
                target_node: target_node.trim_start_matches('#').to_string(),
                line: line.clone(),
                recycle: segment
                    .data_string(extensions::RECYCLE)
                    .map(ToString::to_string),
            });
        }
    }
    segments
}

fn first_reference(object: &Object, property: &str) -> Option<String> {
    object
        .references(property)
        .and_then(|references| references.objects.first())
        .cloned()
}

/// One endpoint of a segment, named the way a connection addresses it.
fn endpoint(
    owner: &NodeOwner,
    source: bool,
    steps: &[Step],
    palette: &[UnitOpSpec],
    diagnostics: &mut Vec<DexpiDiagnostic>,
) -> Option<String> {
    if let Some(name) = &owner.boundary {
        return Some(name.clone());
    }
    let step = &steps[owner.step?];
    let unit = step.unit.as_ref()?;
    let instance = step
        .instance_id
        .clone()
        .unwrap_or_else(|| sanitize_instance_id(&step.identifier));
    let direction = if source {
        Direction::Out
    } else {
        Direction::In
    };
    let ports: Vec<_> = palette
        .iter()
        .find(|spec| spec.id == *unit)
        .map(|spec| {
            spec.ports
                .iter()
                .filter(|port| port.direction == direction)
                .collect::<Vec<_>>()
        })
        .unwrap_or_default();

    let chosen = if let Some(port) = &owner.port {
        // The nozzle names the port. A drawing written by NeqSim names it `INLET`/`OUTLET_n`
        // instead of by the palette's port, so a name the unit does not declare falls back.
        if ports.iter().any(|candidate| candidate.name == *port) {
            port.clone()
        } else {
            // NeqSim's own writer names nozzles `INLET` and `OUTLET_n`; the suffix is the
            // outlet's position, and the first inlet is the only inlet it draws.
            let positional = match port.as_str() {
                "INLET" => Some(0),
                other => other
                    .strip_prefix("OUTLET_")
                    .and_then(|suffix| suffix.parse::<usize>().ok())
                    .and_then(|number| number.checked_sub(1)),
            };
            positional
                .and_then(|index| ports.get(index))
                .map(|port| port.name.clone())
                .unwrap_or_else(|| {
                    diagnostics.push(DexpiDiagnostic::warning(
                        "port_assignment_assumed",
                        step.identifier.clone(),
                        format!(
                            "`{unit}` declares no `{port}` port; the first {} port is used",
                            if source { "outlet" } else { "inlet" }
                        ),
                    ));
                    ports
                        .first()
                        .map_or_else(|| "port".to_string(), |port| port.name.clone())
                })
        }
    } else {
        match ports.as_slice() {
            [] => "port".to_string(),
            [single] => single.name.clone(),
            many => {
                diagnostics.push(DexpiDiagnostic::warning(
                    "port_assignment_assumed",
                    step.identifier.clone(),
                    format!(
                        "`{unit}` declares {} {} ports and the drawing does not name which is \
                         which; assignment is by document order and unverified",
                        many.len(),
                        if source { "outlet" } else { "inlet" }
                    ),
                ));
                many[0].name.clone()
            }
        }
    };
    Some(format!("{instance}.{chosen}"))
}

fn infer_recycles(
    connections: &mut Vec<Connection>,
    recycles: &mut Vec<Recycle>,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) {
    let mut tear = recycles.len();
    while let Some(cycle) = find_cycle(connections) {
        let chosen = cycle.iter().min_by(|left, right| {
            (left.from.as_str(), left.to.as_str()).cmp(&(right.from.as_str(), right.to.as_str()))
        });
        let Some(edge) = chosen else {
            break;
        };
        tear += 1;
        let name = format!("recycle_{tear}");
        diagnostics.push(DexpiDiagnostic::info(
            "recycle_inferred",
            name.clone(),
            format!(
                "the drawing closes a loop and states no tear; `{}` -> `{}` is taken as one, \
                 which is deterministic rather than meaningful",
                edge.from, edge.to
            ),
        ));
        connections.remove(edge.index);
        recycles.push(Recycle::new(&name, &edge.from, &edge.to));
    }
}

struct CycleEdge {
    from: String,
    to: String,
    index: usize,
}

/// A cycle over the instance graph, or `None` when the connections are acyclic.
fn find_cycle(connections: &[Connection]) -> Option<Vec<CycleEdge>> {
    let mut adjacency: BTreeMap<&str, Vec<usize>> = BTreeMap::new();
    for (index, connection) in connections.iter().enumerate() {
        if connection.from.contains('.') && connection.to.contains('.') {
            adjacency
                .entry(instance_of(&connection.from))
                .or_default()
                .push(index);
        }
    }
    let mut state: BTreeMap<&str, u8> = BTreeMap::new();
    let mut path: Vec<usize> = Vec::new();
    for node in adjacency.keys() {
        if state.get(node).copied().unwrap_or(0) == 0 {
            if let Some(cycle) = visit(node, &adjacency, connections, &mut state, &mut path) {
                return Some(cycle);
            }
        }
    }
    None
}

fn visit<'a>(
    node: &'a str,
    adjacency: &BTreeMap<&'a str, Vec<usize>>,
    connections: &'a [Connection],
    state: &mut BTreeMap<&'a str, u8>,
    path: &mut Vec<usize>,
) -> Option<Vec<CycleEdge>> {
    state.insert(node, 1);
    if let Some(edges) = adjacency.get(node) {
        for &index in edges {
            let target: &'a str = instance_of(&connections[index].to);
            // numerics-ok: 0 is the DFS colour for "unvisited", not a magnitude.
            match state.get(target).copied().unwrap_or(0) {
                1 => {
                    // The back edge closes a cycle; the ring is the suffix of the path that
                    // starts at the target.
                    let start = path
                        .iter()
                        .position(|&entry| instance_of(&connections[entry].from) == target)?;
                    let mut cycle = path[start..]
                        .iter()
                        .map(|&entry| CycleEdge {
                            from: connections[entry].from.clone(),
                            to: connections[entry].to.clone(),
                            index: entry,
                        })
                        .collect::<Vec<_>>();
                    cycle.push(CycleEdge {
                        from: connections[index].from.clone(),
                        to: connections[index].to.clone(),
                        index,
                    });
                    return Some(cycle);
                }
                0 => {
                    path.push(index);
                    if let Some(cycle) = visit(target, adjacency, connections, state, path) {
                        return Some(cycle);
                    }
                    path.pop();
                }
                _ => {}
            }
        }
    }
    state.insert(node, 2);
    None
}

/// The instance an endpoint names: `V-101.outlet` is `V-101`.
fn instance_of(endpoint: &str) -> &str {
    endpoint
        .split_once('.')
        .map_or(endpoint, |(instance, _)| instance)
}

fn find_plant_model(document: &Document) -> Option<&Object> {
    let mut found = None;
    walk(&document.model.children, &mut |object| {
        if object.ty == "Plant/PlantModel" && found.is_none() {
            found = Some(object);
        }
    });
    found
}

fn find_engineering_model(document: &Document) -> Option<&Object> {
    document
        .model
        .children
        .iter()
        .find_map(|child| match child {
            Child::Object(root) if root.ty == "Core/EngineeringModel" => Some(root),
            _ => None,
        })
}

fn walk<'a>(children: &'a [Child], visit_object: &mut impl FnMut(&'a Object)) {
    for child in children {
        if let Child::Object(object) = child {
            walk_object(object, visit_object);
        }
    }
}

fn walk_object<'a>(object: &'a Object, visit_object: &mut impl FnMut(&'a Object)) {
    visit_object(object);
    for nested in &object.children {
        if let Child::Components(components) = nested {
            for child in &components.children {
                if let ComponentChild::Object(object) = child {
                    walk_object(object, visit_object);
                }
            }
        }
    }
}

struct Instrumentation {
    instruments: Vec<Instrument>,
    signals: Vec<Signal>,
}

/// The tagged items' object ids to their tags, so an instrument's `SensingLocation` reference
/// becomes the instance id the flowsheet uses.
fn item_tags(plant_model: &Object) -> BTreeMap<String, String> {
    let mut map = BTreeMap::new();
    let Some(items) = plant_model.components("TaggedPlantItems") else {
        return map;
    };
    for item in items.children.iter().filter_map(|child| match child {
        ComponentChild::Object(object) => Some(object),
        _ => None,
    }) {
        if let (Some(id), Some(tag)) = (item.id.clone(), item.data_string("TagName")) {
            map.insert(id, tag.to_string());
        }
    }
    map
}

/// The instruments and the signals between them, as the Plant model states them.
///
/// **The model's fields are the exchange's, and the extensions fill the rest.** Category,
/// modifier and number, the location reference and the safety class are standard; the tag,
/// service, measured and actuated endpoints, range, safe state and intervals are `Azoth/…`
/// properties, because the DEXPI model has no field for them here.
fn read_instruments(
    plant_model: &Object,
    items: &BTreeMap<String, String>,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) -> Instrumentation {
    let mut instruments = Vec::new();
    let mut tags_by_object: BTreeMap<String, String> = BTreeMap::new();
    let mut queued: Vec<(String, String, Option<String>)> = Vec::new();
    let Some(component) = plant_model.components("ProcessInstrumentationFunctions") else {
        return Instrumentation {
            instruments,
            signals: Vec::new(),
        };
    };
    for object in component.children.iter().filter_map(|child| match child {
        ComponentChild::Object(object) => Some(object),
        _ => None,
    }) {
        let tag = object
            .data_string(extensions::TAG)
            .map(ToString::to_string)
            .unwrap_or_else(|| reconstruct_tag(object));
        if let Some(id) = &object.id {
            tags_by_object.insert(id.clone(), tag.clone());
        }
        let location = read_location(object);
        if object.data_value("Location").is_some() && location.is_none() {
            diagnostics.push(DexpiDiagnostic::warning(
                "unsupported_instrument_location",
                tag.clone(),
                "the location reference is not one of Field, ControlPanel or CentralLocation",
            ));
        }
        let safety = read_safety(object, &tag, diagnostics);
        instruments.push(Instrument {
            id: tag.clone(),
            function: function_letters(object),
            service: object
                .data_string(extensions::SERVICE)
                .map(ToString::to_string),
            loop_number: function_number(object),
            measures: object
                .data_string(extensions::MEASURES)
                .map(ToString::to_string)
                .or_else(|| reference_tag(object, "SensingLocation", items)),
            actuates: object
                .data_string(extensions::ACTUATES)
                .map(ToString::to_string)
                .or_else(|| reference_tag(object, "ActuatingLocation", items)),
            location,
            range: read_range(object),
            safety,
        });
        if let Some(signals) = object.components("SignalConveyingFunctions") {
            for signal in signals.children.iter().filter_map(|child| match child {
                ComponentChild::Object(object) => Some(object),
                _ => None,
            }) {
                let Some(target) = first_reference(signal, "Target") else {
                    continue;
                };
                queued.push((
                    tag.clone(),
                    target.trim_start_matches('#').to_string(),
                    read_signal_kind(signal),
                ));
            }
        }
    }

    let mut signals = Vec::new();
    for (from, target, kind) in queued {
        match tags_by_object.get(&target) {
            Some(to) => signals.push(Signal {
                from,
                to: to.clone(),
                kind,
            }),
            None => diagnostics.push(DexpiDiagnostic::warning(
                "unknown_signal_target",
                from,
                format!("the signal's target `{target}` is not a declared instrument"),
            )),
        }
    }
    Instrumentation {
        instruments,
        signals,
    }
}

fn reconstruct_tag(object: &Object) -> String {
    let letters = function_letters(object);
    match function_number(object) {
        Some(number) if !letters.is_empty() => format!("{letters}-{number}"),
        Some(number) => number,
        None if !letters.is_empty() => letters,
        None => object.id.clone().unwrap_or_default(),
    }
}

fn function_letters(object: &Object) -> String {
    let category = object
        .data_string("ProcessInstrumentationFunctionCategory")
        .unwrap_or_default();
    let modifier = object
        .data_string("ProcessInstrumentationFunctionModifier")
        .unwrap_or_default();
    if !category.is_empty() {
        return format!("{category}{modifier}");
    }
    // A file that states no category: the tag's leading letters are the function.
    let tag = object
        .data_string(extensions::TAG)
        .or(object.id.as_deref())
        .unwrap_or_default();
    tag.chars()
        .take_while(char::is_ascii_alphabetic)
        .collect::<String>()
        .to_ascii_uppercase()
}

fn function_number(object: &Object) -> Option<String> {
    if let Some(number) = object.data_string("ProcessInstrumentationFunctionNumber") {
        if !number.is_empty() {
            return Some(number.to_string());
        }
    }
    let tag = object
        .data_string(extensions::TAG)
        .or(object.id.as_deref())
        .unwrap_or_default();
    let digits: String = tag
        .chars()
        .rev()
        .take_while(char::is_ascii_digit)
        .collect::<String>()
        .chars()
        .rev()
        .collect();
    (!digits.is_empty()).then_some(digits)
}

fn read_location(object: &Object) -> Option<String> {
    let DataValue::Reference(reference) = object.data_value("Location")? else {
        return None;
    };
    match reference.rsplit('.').next()? {
        "Field" => Some("field".to_string()),
        "ControlPanel" => Some("panel".to_string()),
        "CentralLocation" => Some("central".to_string()),
        _ => None,
    }
}

fn read_safety(
    object: &Object,
    tag: &str,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) -> Option<SafetyFunction> {
    let class = object.data_string("SafetyRelevanceClass")?;
    let sil = class
        .trim()
        .strip_prefix("SIL")
        .map(str::trim)
        .and_then(|value| value.parse::<u8>().ok());
    let Some(sil) = sil else {
        diagnostics.push(DexpiDiagnostic::warning(
            "unsupported_safety_class",
            tag,
            format!("`{class}` is not a `SIL n` class, so the safety function is not read"),
        ));
        return None;
    };
    Some(SafetyFunction {
        function: object
            .data_string(extensions::SAFETY_FUNCTION)
            .map_or_else(|| tag.to_string(), ToString::to_string),
        sil,
        voting: object
            .data_string("VotingSystemRepresentation")
            .map(ToString::to_string),
        safe_state: object
            .data_string(extensions::SAFE_STATE)
            .map(ToString::to_string),
        proof_test_hours: number(object, extensions::PROOF_TEST_HOURS),
        closure_seconds: number(object, extensions::CLOSURE_SECONDS),
    })
}

fn read_range(object: &Object) -> Option<InstrumentRange> {
    Some(InstrumentRange {
        min: number(object, extensions::RANGE_MIN)?,
        max: number(object, extensions::RANGE_MAX)?,
        unit: object
            .data_string(extensions::RANGE_UNIT)
            .map(ToString::to_string)?,
    })
}

fn read_signal_kind(object: &Object) -> Option<String> {
    let DataValue::Reference(reference) = object.data_value("SignalConveyingType")? else {
        return None;
    };
    match reference.rsplit('.').next()? {
        "ElectricalSignalConveying" => Some("electric".to_string()),
        "PneumaticSignalConveying" => Some("pneumatic".to_string()),
        "HydraulicSignalConveying" => Some("hydraulic".to_string()),
        "CapillarySignalConveying" => Some("capillary".to_string()),
        "ConductedRadiationSignalConveying" => Some("radiation".to_string()),
        _ => None,
    }
}

fn number(object: &Object, property: &str) -> Option<f64> {
    object.data_value(property).and_then(DataValue::as_f64)
}

fn reference_tag(
    object: &Object,
    property: &str,
    map: &BTreeMap<String, String>,
) -> Option<String> {
    let target = first_reference(object, property)?;
    map.get(target.trim_start_matches('#')).cloned()
}
