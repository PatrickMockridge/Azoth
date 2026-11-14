//! The Proteus serialisation - DEXPI 1.x's XML - read into the DEXPI 2.0 Plant model.
//!
//! **A projection, not a second importer.** Proteus carries the same drawing in a different
//! shape: an `Equipment` element with a `ComponentClass` instead of an `Object` with a `type`,
//! `Connection` elements with `FromID`/`ToID` instead of segment node references, and
//! `GenericAttributes` instead of `Data`. This module projects those into the metamodel, and
//! [`crate::plant::import::import_plant_document`] reads the result - so the topology rules,
//! the port assignment, the line numbers and the recycle inference are shared, not written
//! twice.
//!
//! **What the projection drops is named.** Positions, labels, centroids and presentation are
//! graphics; design attributes are not carried into the flowsheet. Inline instrumentation and
//! fittings this port has no unit for are collapsed through, with an
//! `inline_component_skipped` diagnostic per class. Nothing is guessed silently.

use std::collections::{BTreeMap, BTreeSet};

use azoth_process::UnitOpSpec;

use crate::diagnostic::DexpiDiagnostic;
use crate::model::{
    Child, ComponentChild, Components, Data, DataValue, Document, Import, Model, Object, References,
};
use crate::plant::import::{PlantImportOptions, PlantImportReport, import_plant_document};
use crate::plant::{self, PASS_THROUGH_CLASSES};
use crate::profile::{CORE_MODEL, PLANT_MODEL};
use crate::reader::RawElement;

/// Read a Proteus `PlantModel` as a flowsheet scaffold.
///
/// # Errors
/// The diagnostics that refused the import: a file that is not Proteus, one with no plant
/// content, or - through the shared importer - any document the Plant rules refuse.
pub fn import_proteus(
    xml: &str,
    palette: &[UnitOpSpec],
    options: &PlantImportOptions,
) -> Result<PlantImportReport, Vec<DexpiDiagnostic>> {
    let mut diagnostics = Vec::new();
    let projection = project(xml, &mut diagnostics).map_err(|diagnostic| vec![diagnostic])?;
    let mut report = import_plant_document(projection.document, palette, options)?;
    if let Some(drawing) = projection.drawing {
        let positions =
            crate::drawing::layout::apply(&mut report.flowsheet, &drawing, &report.bindings);
        report.layout_positions = positions;
        report.drawing = Some(drawing);
    }
    diagnostics.extend(report.diagnostics);
    report.diagnostics = diagnostics;
    Ok(report)
}

/// The projected document and the drawing it was projected from.
pub(crate) struct Projection {
    pub(crate) document: Document,
    pub(crate) drawing: Option<crate::drawing::Drawing>,
}

/// Whether a file is a Proteus `PlantModel` rather than a DEXPI 2.0 `Model`.
#[must_use]
pub fn is_proteus(xml: &str) -> bool {
    xml.contains("<PlantModel")
}

struct Item {
    object_id: String,
    class: String,
    tag: String,
    nozzles: Vec<ProjectedNozzle>,
}

struct ProjectedNozzle {
    node: String,
    port: Option<String>,
}

/// Project a Proteus document into the DEXPI 2.0 Plant metamodel.
fn project(
    xml: &str,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) -> Result<Projection, DexpiDiagnostic> {
    let root = crate::reader::parse_tree(xml)?;
    if root.name != "PlantModel" {
        return Err(DexpiDiagnostic::error(
            "unsupported_root",
            root.name,
            "a Proteus file has a `PlantModel` root; this reader projects exactly that",
        ));
    }

    let mut elements: Vec<&RawElement> = Vec::new();
    collect(&root, &mut elements);
    let by_id: BTreeMap<String, &RawElement> = elements
        .iter()
        .filter_map(|element| element.attribute("ID").map(|id| (id, *element)))
        .collect();

    // Nodes per element: a nozzle's process nodes, or a piping component's connection points.
    let nodes_of = |element: &RawElement| -> Vec<String> {
        let mut nodes = Vec::new();
        for points in element
            .elements()
            .filter(|child| child.name == "ConnectionPoints")
        {
            for node in points.elements().filter(|child| child.name == "Node") {
                let process = node.attribute("Type").as_deref() == Some("process");
                if let Some(id) = node.attribute("ID") {
                    if process || points.elements().count() == 1 {
                        nodes.push(id);
                    }
                }
            }
        }
        nodes
    };

    let mut items: Vec<Item> = Vec::new();
    let mut ids: BTreeMap<String, String> = BTreeMap::new();
    let mut item_by_element: BTreeMap<String, usize> = BTreeMap::new();
    let mut nozzle_owner: BTreeMap<String, usize> = BTreeMap::new();

    // **Equipment first.** A mapped class is an item; a nested chamber or bundle is not, and
    // its nozzles belong to the enclosing item through their `Association`.
    let mut stack: Vec<(&RawElement, Option<usize>)> = vec![(&root, None)];
    while let Some((element, owner)) = stack.pop() {
        let class = element.attribute("ComponentClass");
        let mut current = owner;
        if let Some(class) = &class {
            if let Some(unit) = plant::for_class(class) {
                if !element.name.starts_with("Nozzle") && element.name != "PipingComponent" {
                    let object_id = format!("Item{}", items.len() + 1);
                    if let Some(source) = element.attribute("ID") {
                        ids.insert(source, object_id.clone());
                    }
                    let tag = tag_of(element)
                        .unwrap_or_else(|| element.attribute("ID").unwrap_or_default());
                    items.push(Item {
                        object_id: object_id.clone(),
                        class: class.clone(),
                        tag,
                        nozzles: Vec::new(),
                    });
                    current = Some(items.len() - 1);
                    item_by_element.insert(
                        element.attribute("ID").unwrap_or_else(|| object_id.clone()),
                        items.len() - 1,
                    );
                }
                let _ = unit;
            }
        }
        if element.name == "Nozzle" {
            if let Some(owner) = current {
                nozzle_owner.insert(element.attribute("ID").unwrap_or_default(), owner);
                let nodes = nodes_of(element);
                if let Some(node) = nodes.first() {
                    items[owner].nozzles.push(ProjectedNozzle {
                        node: node.clone(),
                        port: sub_tag(element),
                    });
                }
            }
        }
        // An `Association` names a nozzle located on this element or its parent.
        for association in element
            .elements()
            .filter(|child| child.name == "Association")
        {
            if association.attribute("Type").as_deref() == Some("is the location of") {
                if let (Some(target), Some(owner)) = (association.attribute("ItemID"), current) {
                    nozzle_owner.insert(target, owner);
                }
            }
        }
        for child in element.elements() {
            stack.push((child, current));
        }
    }
    // Nozzles recorded by association after the item was built.
    for (nozzle_id, owner) in &nozzle_owner {
        if let Some(nozzle) = by_id.get(nozzle_id.as_str()) {
            let nodes = nodes_of(nozzle);
            if let Some(node) = nodes.first() {
                if !items[*owner]
                    .nozzles
                    .iter()
                    .any(|existing| existing.node == *node)
                {
                    items[*owner].nozzles.push(ProjectedNozzle {
                        node: node.clone(),
                        port: sub_tag(nozzle),
                    });
                }
            }
        }
    }

    // **Inline components.** A valve in a segment is an item too; a fitting is collapsed.
    for element in &elements {
        let Some(class) = element.attribute("ComponentClass") else {
            continue;
        };
        if !element.name.starts_with("PipingComponent")
            && element.name != "Nozzle"
            && element.name != "Equipment"
        {
            continue;
        }
        if plant::for_class(&class).is_none() {
            continue;
        }
        if item_by_element.contains_key(&element.attribute("ID").unwrap_or_default()) {
            continue;
        }
        let nodes = nodes_of(element);
        if nodes.len() < 2 {
            continue;
        }
        let object_id = format!("Item{}", items.len() + 1);
        if let Some(source) = element.attribute("ID") {
            ids.insert(source, object_id.clone());
        }
        let tag = tag_of(element).unwrap_or_else(|| element.attribute("ID").unwrap_or_default());
        // **The nozzle names no palette port.** A Proteus connection point does not say which
        // side of a valve a stream enters; the importer's direction rule does, so the name is
        // left unstated rather than guessed from node order.
        items.push(Item {
            object_id,
            class: class.clone(),
            tag,
            nozzles: vec![
                ProjectedNozzle {
                    node: nodes[0].clone(),
                    port: None,
                },
                ProjectedNozzle {
                    node: nodes[1].clone(),
                    port: None,
                },
            ],
        });
        item_by_element.insert(element.attribute("ID").unwrap_or_default(), items.len() - 1);
    }

    // **Boundaries.** Off-page connectors are the drawing's environment inlets and outlets.
    let mut boundaries: Vec<(String, String, String, bool)> = Vec::new();
    for element in &elements {
        let Some(class) = element.attribute("ComponentClass") else {
            continue;
        };
        let inbound = class == "FlowInPipeOffPageConnector";
        let outbound = class == "FlowOutPipeOffPageConnector";
        if !inbound && !outbound {
            continue;
        }
        let nodes = nodes_of(element);
        let Some(node) = nodes.first() else {
            continue;
        };
        let name = connector_name(element);
        boundaries.push((
            element.attribute("ID").unwrap_or_default(),
            node.clone(),
            name,
            inbound,
        ));
    }

    // **The segment graph.** Each `Connection` is an edge between two component nodes; a
    // pass-through's two nodes are joined internally so the contraction below walks through it.
    let mut terminal: BTreeSet<String> = BTreeSet::new();
    for item in &items {
        for nozzle in &item.nozzles {
            terminal.insert(nozzle.node.clone());
        }
    }
    for (_, node, _, _) in &boundaries {
        terminal.insert(node.clone());
    }

    // **A component is terminal if the projection already models it.** An item, an item's
    // nozzle and a boundary connector are all drawn; everything else a `Connection` reaches is
    // an inline fitting or instrument, and is joined through.
    let mut component_terminal: BTreeSet<&str> = BTreeSet::new();
    for id in item_by_element.keys() {
        component_terminal.insert(id.as_str());
    }
    for id in nozzle_owner.keys() {
        component_terminal.insert(id.as_str());
    }
    for (_, _, name, _) in &boundaries {
        component_terminal.insert(name.as_str());
    }
    for element in &elements {
        let class = element.attribute("ComponentClass").unwrap_or_default();
        if class == "FlowInPipeOffPageConnector" || class == "FlowOutPipeOffPageConnector" {
            if let Some(id) = element.attribute("ID") {
                component_terminal
                    .insert(by_id.get_key_value(&id).map_or("", |(key, _)| key.as_str()));
            }
        }
    }

    let mut edges: BTreeMap<(String, String), String> = BTreeMap::new();
    let mut internal: Vec<(String, String)> = Vec::new();
    let mut skipped: BTreeSet<String> = BTreeSet::new();
    for system in elements
        .iter()
        .filter(|element| element.name == "PipingNetworkSystem")
    {
        let line = attribute_value(system, "LineNumber").unwrap_or_default();
        for segment in system
            .elements()
            .filter(|element| element.name == "PipingNetworkSegment")
        {
            for connection in segment
                .elements()
                .filter(|element| element.name == "Connection")
            {
                let Some(from) = connection.attribute("FromID") else {
                    continue;
                };
                let Some(to) = connection.attribute("ToID") else {
                    continue;
                };
                let Some(source) = node_of(
                    &by_id,
                    &from,
                    connection.attribute("FromNode").as_deref(),
                    &nodes_of,
                ) else {
                    diagnostics.push(DexpiDiagnostic::warning(
                        "unresolved_component",
                        from.clone(),
                        "the connection's source component is not in the file",
                    ));
                    continue;
                };
                let Some(target) = node_of(
                    &by_id,
                    &to,
                    connection.attribute("ToNode").as_deref(),
                    &nodes_of,
                ) else {
                    diagnostics.push(DexpiDiagnostic::warning(
                        "unresolved_component",
                        to.clone(),
                        "the connection's target component is not in the file",
                    ));
                    continue;
                };
                edges
                    .entry((source.clone(), target.clone()))
                    .or_insert_with(|| line.clone());
                for id in [&from, &to] {
                    let class = by_id
                        .get(id.as_str())
                        .and_then(|element| element.attribute("ComponentClass"))
                        .unwrap_or_default();
                    if !component_terminal.contains(id.as_str())
                        && plant::for_class(&class).is_none()
                        && !skipped.contains(id.as_str())
                    {
                        skipped.insert(id.clone());
                        diagnostics.push(DexpiDiagnostic::info(
                            "inline_component_skipped",
                            id.clone(),
                            format!(
                                "`{class}` is a fitting or instrument this port has no unit for; \
                                 the stream crosses it and the topology is joined through it"
                            ),
                        ));
                        if PASS_THROUGH_CLASSES.contains(&class.as_str()) {
                            if let Some(element) = by_id.get(id.as_str()) {
                                let nodes = nodes_of(element);
                                if nodes.len() >= 2 {
                                    internal.push((nodes[0].clone(), nodes[1].clone()));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    // Every inline component's nodes are internal; contract them even where a Connection
    // reached only one side.
    for element in &elements {
        let Some(class) = element.attribute("ComponentClass") else {
            continue;
        };
        if plant::for_class(&class).is_some() {
            continue;
        }
        let nodes = nodes_of(element);
        if nodes.len() >= 2 {
            internal.push((nodes[0].clone(), nodes[1].clone()));
        }
    }

    let collapsed = contract(edges, &internal);
    let dropped = collapsed
        .keys()
        .filter(|(source, target)| !terminal.contains(source) || !terminal.contains(target))
        .count();
    if dropped > 0 {
        diagnostics.push(DexpiDiagnostic::info(
            "topology_gap",
            "",
            format!(
                "{dropped} connection(s) end at a component this port does not model as an \
                 item; they are not drawn in the flowsheet"
            ),
        ));
    }

    // **The projected document.**
    let mut tagged_items = Vec::new();
    for item in &items {
        let mut nozzles = Vec::new();
        for (index, nozzle) in item.nozzles.iter().enumerate() {
            let mut children = vec![Child::Components(Components {
                property: "Nodes".to_string(),
                children: vec![ComponentChild::Object(Object {
                    id: Some(nozzle.node.clone()),
                    name: None,
                    ty: "Plant/Piping.PipingNode".to_string(),
                    children: Vec::new(),
                })],
            })];
            if let Some(port) = &nozzle.port {
                children.push(Child::Data(string_data("SubTagName", port)));
            }
            nozzles.push(ComponentChild::Object(Object {
                id: Some(format!("{}Nozzle{}", item.object_id, index + 1)),
                name: None,
                ty: "Plant/ProcessEquipment.Nozzle".to_string(),
                children,
            }));
        }
        tagged_items.push(ComponentChild::Object(Object {
            id: Some(item.object_id.clone()),
            name: None,
            ty: native_type(&item.class),
            children: vec![
                Child::Data(string_data("TagName", &item.tag)),
                Child::Components(Components {
                    property: "Nozzles".to_string(),
                    children: nozzles,
                }),
            ],
        }));
    }

    let mut segments_by_line: BTreeMap<String, Vec<ComponentChild>> = BTreeMap::new();
    for (index, ((source, target), line)) in collapsed.iter().enumerate() {
        if !terminal.contains(source) || !terminal.contains(target) {
            continue;
        }
        let line = line.clone();
        segments_by_line
            .entry(line)
            .or_default()
            .push(ComponentChild::Object(Object {
                id: Some(format!("Segment{}", index + 1)),
                name: None,
                ty: "Plant/Piping.PipingNetworkSegment".to_string(),
                children: vec![
                    Child::References(References {
                        property: "SourceNode".to_string(),
                        objects: vec![format!("#{source}")],
                    }),
                    Child::References(References {
                        property: "TargetNode".to_string(),
                        objects: vec![format!("#{target}")],
                    }),
                ],
            }));
    }
    let mut systems = Vec::new();
    for (index, (line, segments)) in segments_by_line.into_iter().enumerate() {
        let mut children = Vec::new();
        if !line.is_empty() {
            children.push(Child::Data(string_data("LineNumber", &line)));
        }
        children.push(Child::Components(Components {
            property: "Segments".to_string(),
            children: segments,
        }));
        systems.push(ComponentChild::Object(Object {
            id: Some(format!("PipingNetworkSystem{}", index + 1)),
            name: None,
            ty: "Plant/Piping.PipingNetworkSystem".to_string(),
            children,
        }));
    }

    let mut plant_children = vec![Child::Components(Components {
        property: "TaggedPlantItems".to_string(),
        children: tagged_items,
    })];
    if !systems.is_empty() {
        plant_children.push(Child::Components(Components {
            property: "PipingNetworkSystems".to_string(),
            children: systems,
        }));
    }
    // Boundary connectors are placed where the Plant import expects them: inside segments.
    // Here they are held on the plant model and re-emitted into the first system.
    if !boundaries.is_empty() {
        let connectors: Vec<ComponentChild> = boundaries
            .iter()
            .enumerate()
            .map(|(index, (source, node, name, inbound))| {
                ids.insert(source.clone(), format!("BoundaryConnector{}", index + 1));
                ComponentChild::Object(Object {
                    id: Some(format!("BoundaryConnector{}", index + 1)),
                    name: None,
                    ty: if *inbound {
                        "Plant/Piping.FlowInPipeOffPageConnector".to_string()
                    } else {
                        "Plant/Piping.FlowOutPipeOffPageConnector".to_string()
                    },
                    children: vec![
                        Child::Components(Components {
                            property: "Nodes".to_string(),
                            children: vec![ComponentChild::Object(Object {
                                id: Some(node.clone()),
                                name: None,
                                ty: "Plant/Piping.PipingNode".to_string(),
                                children: Vec::new(),
                            })],
                        }),
                        Child::Data(string_data("PipeConnectorNumber", name)),
                    ],
                })
            })
            .collect();
        // A segment's `Items` is an Object child; the Plant import walks the document, so a
        // `Boundary` collection on the model is enough for it to find them.
        plant_children.push(Child::Components(Components {
            property: "Boundaries".to_string(),
            children: connectors,
        }));
    }

    let plant_model = Object {
        id: Some("PlantModel1".to_string()),
        name: None,
        ty: "Plant/PlantModel".to_string(),
        children: plant_children,
    };
    let engineering_model = Object {
        id: None,
        name: None,
        ty: "Core/EngineeringModel".to_string(),
        children: vec![Child::Components(Components {
            property: "ConceptualModel".to_string(),
            children: vec![ComponentChild::Object(plant_model)],
        })],
    };

    let document = Document {
        model: Model {
            name: "ProteusProjection".to_string(),
            uri: "urn:azoth:dexpi:2.0:plant:proteus-projection".to_string(),
            imports: vec![
                Import {
                    prefix: "Core".to_string(),
                    source: CORE_MODEL.to_string(),
                },
                Import {
                    prefix: "Plant".to_string(),
                    source: PLANT_MODEL.to_string(),
                },
            ],
            children: vec![Child::Object(engineering_model)],
        },
    };
    let drawing = crate::drawing::proteus::extract(&root, &ids).map(|report| {
        diagnostics.extend(report.diagnostics);
        report.drawing
    });
    Ok(Projection { document, drawing })
}

fn collect<'a>(element: &'a RawElement, elements: &mut Vec<&'a RawElement>) {
    elements.push(element);
    for child in element.elements() {
        collect(child, elements);
    }
}

fn node_of(
    by_id: &BTreeMap<String, &RawElement>,
    id: &str,
    index: Option<&str>,
    nodes_of: &impl Fn(&RawElement) -> Vec<String>,
) -> Option<String> {
    let element = by_id.get(id)?;
    let nodes = nodes_of(element);
    let position = index
        .and_then(|value| value.parse::<usize>().ok())
        .and_then(|value| value.checked_sub(1))
        .unwrap_or(0);
    nodes
        .get(position)
        .cloned()
        .or_else(|| nodes.first().cloned())
}

/// Join pass-through components so an equipment-to-equipment edge survives them.
///
/// **The pair is crossed the way the flow reaches it.** A Proteus connection does not say
/// which of a fitting's two nodes is upstream, so an edge arriving at one node leaves by the
/// *other*: incoming to `a` meets outgoing from `b`, and incoming to `b` meets outgoing from
/// `a`. The rewrite runs to a fixpoint, bounded so a ring of fittings cannot spin.
fn contract(
    edges: BTreeMap<(String, String), String>,
    internal: &[(String, String)],
) -> BTreeMap<(String, String), String> {
    let mut edges = edges;
    let bound = internal.len().saturating_mul(4).max(4);
    for _ in 0..bound {
        let mut changed = false;
        for (a, b) in internal {
            let into = |node: &String| -> Vec<((String, String), String)> {
                edges
                    .iter()
                    .filter(|((source, target), _)| target == node && source != a && source != b)
                    .map(|(edge, line)| (edge.clone(), line.clone()))
                    .collect()
            };
            let out = |node: &String| -> Vec<((String, String), String)> {
                edges
                    .iter()
                    .filter(|((source, target), _)| source == node && target != a && target != b)
                    .map(|(edge, line)| (edge.clone(), line.clone()))
                    .collect()
            };
            let into_a = into(a);
            let into_b = into(b);
            let out_a = out(a);
            let out_b = out(b);
            if (into_a.is_empty() && into_b.is_empty()) || (out_a.is_empty() && out_b.is_empty()) {
                continue;
            }
            for ((source, _), line) in &into_a {
                for ((_, target), _) in &out_b {
                    edges
                        .entry((source.clone(), target.clone()))
                        .or_insert_with(|| line.clone());
                    changed = true;
                }
            }
            for ((source, _), line) in &into_b {
                for ((_, target), _) in &out_a {
                    edges
                        .entry((source.clone(), target.clone()))
                        .or_insert_with(|| line.clone());
                    changed = true;
                }
            }
            for ((source, _), _) in &into_a {
                edges.remove(&(source.clone(), a.clone()));
            }
            for ((source, _), _) in &into_b {
                edges.remove(&(source.clone(), b.clone()));
            }
            for ((_, target), _) in &out_a {
                edges.remove(&(a.clone(), target.clone()));
            }
            for ((_, target), _) in &out_b {
                edges.remove(&(b.clone(), target.clone()));
            }
            edges.remove(&(a.clone(), b.clone()));
            edges.remove(&(b.clone(), a.clone()));
        }
        if !changed {
            break;
        }
    }
    edges
}

fn tag_of(element: &RawElement) -> Option<String> {
    if let Some(value) = attribute_value(element, "TagName") {
        return Some(value);
    }
    for label in element.elements().filter(|child| child.name == "Label") {
        if label.attribute("ComponentClass").as_deref() == Some("EquipmentTagNameLabel") {
            for text in label.elements().filter(|child| child.name == "Text") {
                if let Some(value) = text.attribute("String") {
                    return Some(value);
                }
            }
        }
    }
    None
}

fn sub_tag(element: &RawElement) -> Option<String> {
    attribute_value(element, "SubTagName")
}

fn connector_name(element: &RawElement) -> String {
    attribute_value(element, "PipeConnectorNumber")
        .or_else(|| attribute_value(element, "TagName"))
        .or_else(|| element.attribute("ID"))
        .unwrap_or_default()
}

/// A `GenericAttribute` by name, with the `AssignmentClass` suffix the RDL uses.
fn attribute_value(element: &RawElement, name: &str) -> Option<String> {
    for attributes in element
        .elements()
        .filter(|child| child.name == "GenericAttributes")
    {
        for attribute in attributes
            .elements()
            .filter(|child| child.name == "GenericAttribute")
        {
            let attribute_name = attribute.attribute("Name")?;
            if attribute_name == name || attribute_name == format!("{name}AssignmentClass") {
                return attribute.attribute("Value");
            }
        }
    }
    None
}

/// The native DEXPI 2.0 type a Proteus class projects to.
fn native_type(class: &str) -> String {
    match class {
        "GlobeValve"
        | "ControlValve"
        | "GateValve"
        | "BallValve"
        | "ButterflyValve"
        | "CheckValve"
        | "SwingCheckValve"
        | "PressureSafetyValve"
        | "PressureReliefValve"
        | "SpringLoadedGlobeSafetyValve"
        | "RuptureDisk"
        | "Valve" => "Plant/Piping.GlobeValve".to_string(),
        other => format!("Plant/ProcessEquipment.{other}"),
    }
}

fn string_data(property: &str, value: &str) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::String(value.to_string())],
    }
}
