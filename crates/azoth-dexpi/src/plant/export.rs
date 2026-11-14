//! A flowsheet as a DEXPI 2.0 Plant model: tagged items, nozzles, and the piping between them.
//!
//! The Plant model carries no process state, so nothing is run and no quantity is written. What
//! the drawing carries is the topology a PFD does not: the tag of every item, the nozzle each
//! stream leaves or enters by, and the line it belongs to.
//!
//! **One nozzle per palette port.** A `many` inlet is one physical nozzle that several pipes
//! reach; a `many` outlet is one nozzle whose streams leave under one port name and are told
//! apart by the endpoint's index. That is the drawing's own reading, and it round-trips.

use std::collections::BTreeMap;

use azoth_process::flowsheet::{Flowsheet, Instrument, Signal, split_endpoint};
use azoth_process::middleware::graph;
use azoth_process::{Severity, UnitOpSpec, validate};

use crate::diagnostic::DexpiDiagnostic;
use crate::drawing::Drawing;
use crate::drawing::export::{self as drawing_export, DiagramEdge, DiagramNode, NodeShape};
use crate::drawing::native::write_diagram;
use crate::extensions;
use crate::ids::{Mint, model_name};
use crate::mapping::Disposition;
use crate::model::{
    Child, ComponentChild, Components, Data, DataValue, Document, Import, Model, Object, References,
};
use crate::plant::{self, DEFAULT_LINE, PlantMapping};
use crate::profile::{PLANT_SUPPORTED_TYPES, Profile};
use crate::validate::validate as validate_document;
use crate::writer::write;

/// The status an exchange carries when it has not been certified.
pub const APPROVAL_STATUS: &str = "NOT_A_DEXPI_EV_CERTIFICATE";
/// The round-trip status an exchange carries when no named CAE tool has been qualified.
pub const CAE_ROUND_TRIP_STATUS: &str = "QUALIFICATION_REQUIRED";

/// What a Plant export is allowed to do.
#[derive(Debug, Clone, Default)]
pub struct PlantExportOptions {
    /// Include mappings whose Plant class loses something.
    pub allow_lossy: bool,
    /// Per-palette-unit overrides for a refused mapping, as `unit_ops.x=Plant/...`.
    pub overrides: BTreeMap<String, String>,
    /// Write the `Azoth/…` properties that make an azoth round trip exact.
    pub azoth_extensions: bool,
    /// Draw the P&ID as generic `Core/Diagram` primitives.
    pub drawing: bool,
}

/// A completed Plant export.
#[derive(Debug, Clone)]
pub struct PlantExportReport {
    pub xml: String,
    pub diagnostics: Vec<DexpiDiagnostic>,
    /// Tagged items written.
    pub items: usize,
    /// Piping segments written.
    pub segments: usize,
    /// Drawing groups written; zero when the export carried no drawing.
    pub diagram_groups: usize,
    /// Drawing primitives written.
    pub diagram_primitives: usize,
    pub approval_status: &'static str,
    pub cae_round_trip: &'static str,
    pub semantic_valid: bool,
}

/// Export a flowsheet as a DEXPI Plant exchange.
///
/// # Errors
/// The diagnostics that refused the export: a flowsheet the checker rejects, a unit with no
/// reviewed Plant class, or a lossy mapping without [`PlantExportOptions::allow_lossy`].
pub fn export_plant(
    flowsheet: &Flowsheet,
    palette: &[UnitOpSpec],
    options: &PlantExportOptions,
) -> Result<PlantExportReport, Vec<DexpiDiagnostic>> {
    let mut diagnostics = Vec::new();
    let mut fatal = Vec::new();

    for diagnostic in validate(flowsheet, palette) {
        if diagnostic.severity() == Severity::Error {
            fatal.push(DexpiDiagnostic::error(
                "invalid_flowsheet",
                diagnostic.location().path,
                diagnostic.message(),
            ));
        }
    }
    if !fatal.is_empty() {
        return Err(fatal);
    }

    let mapped = resolve_mappings(flowsheet, options, &mut fatal);
    if !fatal.is_empty() {
        return Err(fatal);
    }

    let mut builder = Builder::new(flowsheet, palette, options);
    builder.add_items(&mapped);
    builder.add_segments(&mut fatal);
    if !fatal.is_empty() {
        return Err(fatal);
    }
    let items = builder.items.len();
    let segments = builder.segment_count;
    let document = builder.finish();

    let validation = validate_document(&document);
    diagnostics.extend(validation.diagnostics.clone());
    if !validation.is_valid() {
        return Err(validation.diagnostics);
    }

    let (diagram_groups, diagram_primitives) = crate::drawing::native::counts(&document);
    let xml = write(&document).map_err(|error| {
        vec![DexpiDiagnostic::error(
            "write_failed",
            flowsheet.id.clone(),
            error.to_string(),
        )]
    })?;

    Ok(PlantExportReport {
        xml,
        diagnostics,
        items,
        segments,
        diagram_groups,
        diagram_primitives,
        approval_status: APPROVAL_STATUS,
        cae_round_trip: CAE_ROUND_TRIP_STATUS,
        semantic_valid: true,
    })
}

fn resolve_mappings(
    flowsheet: &Flowsheet,
    options: &PlantExportOptions,
    fatal: &mut Vec<DexpiDiagnostic>,
) -> BTreeMap<String, String> {
    let mut mapped = BTreeMap::new();
    for instance in &flowsheet.instances {
        let Some(row) = plant::for_azoth(&instance.unit) else {
            fatal.push(DexpiDiagnostic::error(
                "unmapped_unit",
                instance.id.clone(),
                format!(
                    "`{}` has no row in specs/dexpi/plant_units.toml",
                    instance.unit
                ),
            ));
            continue;
        };
        if let Some(override_type) = options.overrides.get(&instance.unit) {
            if !PLANT_SUPPORTED_TYPES.contains(&override_type.as_str()) {
                fatal.push(DexpiDiagnostic::error(
                    "override_outside_profile",
                    instance.id.clone(),
                    format!("`{override_type}` is outside the reviewed Plant profile"),
                ));
                continue;
            }
            mapped.insert(instance.id.clone(), override_type.clone());
            continue;
        }
        let allowed = match row.write {
            Disposition::Faithful => true,
            Disposition::Lossy => options.allow_lossy,
            Disposition::Refused => false,
        };
        if !allowed {
            let code = if row.write == Disposition::Refused {
                "refused_mapping"
            } else {
                "lossy_mapping"
            };
            fatal.push(DexpiDiagnostic::error(
                code,
                instance.id.clone(),
                format!(
                    "`{}` maps to {}: {} An export needs a reviewed class or an explicit \
                     `--map {}=<Plant type>`.",
                    instance.unit,
                    describe(row),
                    row.note,
                    instance.unit
                ),
            ));
            continue;
        }
        let Some(plant_type) = row.plant_types.first() else {
            fatal.push(DexpiDiagnostic::error(
                "unmapped_unit",
                instance.id.clone(),
                format!("`{}` names no Plant class", instance.unit),
            ));
            continue;
        };
        mapped.insert(instance.id.clone(), (*plant_type).to_string());
    }
    mapped
}

fn describe(row: &PlantMapping) -> String {
    row.plant_types
        .first()
        .map_or_else(|| "no reviewed class".to_string(), |ty| format!("`{ty}`"))
}

struct Builder<'a> {
    flowsheet: &'a Flowsheet,
    palette: &'a [UnitOpSpec],
    options: &'a PlantExportOptions,
    items: Vec<Object>,
    systems: BTreeMap<String, Vec<Object>>,
    segment_count: usize,
    connector_count: usize,
    mint: Mint,
    port_node: BTreeMap<(String, String), String>,
    /// The instance id to the equipment object id, so an instrument can reference the item it
    /// watches or commands.
    item_ids: BTreeMap<String, String>,
    /// The endpoint name to the off-page connector object id, for the drawing.
    boundaries: BTreeMap<String, String>,
    /// One entry per segment, for the drawing.
    segment_edges: Vec<SegmentEdge>,
    /// The instance id to the Plant type it is written as, including segment-placed units.
    mapped: BTreeMap<String, String>,
    /// A counter for the pipes a segment's valve chain needs.
    pipe_count: usize,
    /// A counter for the valve objects, which are not items.
    valve_count: usize,
}

struct SegmentEdge {
    object_id: String,
    from: String,
    to: String,
    line: String,
    recycle: bool,
}

/// One valve the writer placed inside a segment.
#[derive(Debug, Clone)]
struct ValvePlacement {
    instance: String,
    unit: String,
    ty: String,
}

impl ValvePlacement {
    fn tag(&self) -> &str {
        &self.instance
    }
}

#[derive(Debug, Clone)]
struct Edge {
    from: String,
    to: String,
    line: String,
    recycle: Option<String>,
    /// The valves this edge passes through, in stream order.
    valves: Vec<ValvePlacement>,
}

impl<'a> Builder<'a> {
    fn new(
        flowsheet: &'a Flowsheet,
        palette: &'a [UnitOpSpec],
        options: &'a PlantExportOptions,
    ) -> Self {
        Self {
            flowsheet,
            palette,
            options,
            items: Vec::new(),
            systems: BTreeMap::new(),
            segment_count: 0,
            connector_count: 0,
            mint: Mint::new(),
            port_node: BTreeMap::new(),
            item_ids: BTreeMap::new(),
            boundaries: BTreeMap::new(),
            segment_edges: Vec::new(),
            mapped: BTreeMap::new(),
            pipe_count: 0,
            valve_count: 0,
        }
    }

    fn add_items(&mut self, mapped: &BTreeMap<String, String>) {
        self.mapped = mapped.clone();
        for instance in &self.flowsheet.instances {
            // **A segment-placed unit is not an item.** Its valve object is made inside the
            // segment that passes through it, where the Plant model puts a piping component.
            if plant::is_segment_unit(&instance.unit) {
                continue;
            }
            let Some(ty) = mapped.get(&instance.id) else {
                continue;
            };
            let equipment_id = format!("Equipment{}", self.items.len() + 1);
            let mut nozzles = Vec::new();
            if let Some(spec) = self.palette.iter().find(|spec| spec.id == instance.unit) {
                for port in &spec.ports {
                    let nozzle_id = format!("{equipment_id}Nozzle{}", nozzles.len() + 1);
                    let node = self.mint.next("PipingNode");
                    self.port_node
                        .insert((instance.id.clone(), port.name.clone()), node.clone());
                    let mut nozzle_children = vec![
                        Child::Components(Components {
                            property: "Nodes".to_string(),
                            children: vec![ComponentChild::Object(Object {
                                id: Some(node),
                                name: None,
                                ty: "Plant/Piping.PipingNode".to_string(),
                                children: Vec::new(),
                            })],
                        }),
                        Child::Data(string_data("SubTagName", &port.name)),
                    ];
                    if self.options.azoth_extensions {
                        nozzle_children
                            .push(Child::Data(string_data(extensions::PORT, &port.name)));
                    }
                    nozzles.push(ComponentChild::Object(Object {
                        id: Some(nozzle_id),
                        name: None,
                        ty: "Plant/ProcessEquipment.Nozzle".to_string(),
                        children: nozzle_children,
                    }));
                }
            }
            let mut children = vec![Child::Data(string_data("TagName", &instance.id))];
            if self.options.azoth_extensions {
                children.push(Child::Data(string_data(
                    extensions::UNIT_OP,
                    &instance.unit,
                )));
                children.push(Child::Data(string_data(
                    extensions::INSTANCE_ID,
                    &instance.id,
                )));
            }
            children.push(Child::Components(Components {
                property: "Nozzles".to_string(),
                children: nozzles,
            }));
            self.item_ids
                .insert(instance.id.clone(), equipment_id.clone());
            self.items.push(Object {
                id: Some(equipment_id),
                name: None,
                ty: ty.clone(),
                children,
            });
        }
    }

    fn add_segments(&mut self, fatal: &mut Vec<DexpiDiagnostic>) {
        let mut edges: Vec<Edge> = Vec::new();
        for connection in &self.flowsheet.connections {
            edges.push(Edge {
                from: connection.from.clone(),
                to: connection.to.clone(),
                line: connection
                    .line
                    .clone()
                    .unwrap_or_else(|| DEFAULT_LINE.to_string()),
                recycle: None,
                valves: Vec::new(),
            });
        }
        for recycle in &self.flowsheet.recycles {
            edges.push(Edge {
                from: recycle.from.clone(),
                to: recycle.to.clone(),
                line: DEFAULT_LINE.to_string(),
                recycle: self
                    .options
                    .azoth_extensions
                    .then(|| recycle.stream.clone()),
                valves: Vec::new(),
            });
        }

        // **A valve is inside the segment it sits on.** The two streams that reach its ports
        // become one segment, so the chain source -> valve -> target is one drawing edge with the
        // valve in its `Items`, which is where the Plant model puts a piping component.
        let valves: BTreeMap<String, ValvePlacement> = self
            .flowsheet
            .instances
            .iter()
            .filter(|instance| plant::is_segment_unit(&instance.unit))
            .filter_map(|instance| {
                self.mapped.get(&instance.id).map(|ty| {
                    (
                        instance.id.clone(),
                        ValvePlacement {
                            instance: instance.id.clone(),
                            unit: instance.unit.clone(),
                            ty: ty.clone(),
                        },
                    )
                })
            })
            .collect();
        contract_valves(&mut edges, &valves, fatal);

        for edge in edges {
            let Some(source) = self.endpoint_node(&edge.from, true) else {
                fatal.push(DexpiDiagnostic::error(
                    "unresolved_endpoint",
                    edge.from.clone(),
                    "the stream's source is neither a declared input nor an instance outlet",
                ));
                continue;
            };
            let Some(target) = self.endpoint_node(&edge.to, false) else {
                fatal.push(DexpiDiagnostic::error(
                    "unresolved_endpoint",
                    edge.to.clone(),
                    "the stream's target is neither a declared product nor an instance inlet",
                ));
                continue;
            };
            self.segment_count += 1;
            let number = self.segment_count;
            let mut children = Vec::new();
            let mut boundary_items = Vec::new();
            let source_node = match source {
                Node::Existing(node) => node,
                Node::Boundary { node, item } => {
                    boundary_items.push(item);
                    node
                }
            };
            let target_node = match target {
                Node::Existing(node) => node,
                Node::Boundary { node, item } => {
                    boundary_items.push(item);
                    node
                }
            };
            let mut items: Vec<ComponentChild> = boundary_items
                .into_iter()
                .map(ComponentChild::Object)
                .collect();
            let mut valve_nodes: Vec<(String, String)> = Vec::new();
            for valve in &edge.valves {
                self.valve_count += 1;
                let valve_id = format!("Valve{}", self.valve_count);
                let inlet = self.mint.next("ValveNode");
                let outlet = self.mint.next("ValveNode");
                let mut valve_children = vec![
                    Child::Components(Components {
                        property: "Nodes".to_string(),
                        children: vec![
                            ComponentChild::Object(Object {
                                id: Some(inlet.clone()),
                                name: None,
                                ty: "Plant/Piping.PipingNode".to_string(),
                                children: Vec::new(),
                            }),
                            ComponentChild::Object(Object {
                                id: Some(outlet.clone()),
                                name: None,
                                ty: "Plant/Piping.PipingNode".to_string(),
                                children: Vec::new(),
                            }),
                        ],
                    }),
                    Child::Data(string_data("TagName", valve.tag())),
                ];
                if self.options.azoth_extensions {
                    valve_children.push(Child::Data(string_data(extensions::UNIT_OP, &valve.unit)));
                    valve_children.push(Child::Data(string_data(
                        extensions::INSTANCE_ID,
                        &valve.instance,
                    )));
                    valve_children.push(Child::Data(string_data(extensions::TAG, valve.tag())));
                }
                self.item_ids
                    .insert(valve.instance.clone(), valve_id.clone());
                items.push(ComponentChild::Object(Object {
                    id: Some(valve_id),
                    name: None,
                    ty: valve.ty.clone(),
                    children: valve_children,
                }));
                valve_nodes.push((inlet, outlet));
            }
            if !items.is_empty() {
                children.push(Child::Components(Components {
                    property: "Items".to_string(),
                    children: items,
                }));
            }
            let mut pipes = Vec::new();
            let mut current = source_node.clone();
            for (inlet, outlet) in &valve_nodes {
                self.pipe_count += 1;
                pipes.push(pipe_object(self.pipe_count, &current, inlet));
                current = outlet.clone();
            }
            self.pipe_count += 1;
            pipes.push(pipe_object(self.pipe_count, &current, &target_node));
            let mut segment = vec![
                Child::References(References {
                    property: "SourceNode".to_string(),
                    objects: vec![format!("#{source_node}")],
                }),
                Child::References(References {
                    property: "TargetNode".to_string(),
                    objects: vec![format!("#{target_node}")],
                }),
                Child::Components(Components {
                    property: "Connections".to_string(),
                    children: pipes,
                }),
                Child::Data(string_data("SegmentNumber", &number.to_string())),
            ];
            segment.splice(0..0, children);
            if self.options.azoth_extensions {
                segment.push(Child::Data(string_data(extensions::FROM, &edge.from)));
                segment.push(Child::Data(string_data(extensions::TO, &edge.to)));
                if let Some(recycle) = &edge.recycle {
                    segment.push(Child::Data(string_data(extensions::RECYCLE, recycle)));
                }
            }
            let object_id = format!("PipingNetworkSegment{number}");
            let object = Object {
                id: Some(object_id.clone()),
                name: None,
                ty: "Plant/Piping.PipingNetworkSegment".to_string(),
                children: segment,
            };
            self.segment_edges.push(SegmentEdge {
                object_id,
                from: edge.from.clone(),
                to: edge.to.clone(),
                line: edge.line.clone(),
                recycle: edge.recycle.is_some(),
            });
            self.systems.entry(edge.line).or_default().push(object);
        }
    }

    /// The node an endpoint reaches: an item's nozzle, or a boundary connector minted here.
    fn endpoint_node(&mut self, endpoint: &str, source: bool) -> Option<Node> {
        if let Some((instance, port, _)) = split_endpoint(endpoint) {
            return self
                .port_node
                .get(&(instance.to_string(), port.to_string()))
                .cloned()
                .map(Node::Existing);
        }
        // A bare name is a boundary. It gets an off-page connector inside the segment that
        // uses it, which is where DEXPI puts one.
        let exists = if source {
            self.flowsheet
                .inputs
                .iter()
                .any(|input| input.name == endpoint)
        } else {
            self.flowsheet
                .products
                .iter()
                .any(|product| product == endpoint)
        };
        if !exists {
            return None;
        }
        self.connector_count += 1;
        let number = self.connector_count;
        self.boundaries
            .insert(endpoint.to_string(), format!("BoundaryConnector{number}"));
        let node = self.mint.next("BoundaryPipingNode");
        let item = Object {
            id: Some(format!("BoundaryConnector{number}")),
            name: None,
            ty: if source {
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
                Child::Data(string_data("PipeConnectorNumber", endpoint)),
            ],
        };
        Some(Node::Boundary { node, item })
    }

    /// The P&ID: items, off-page connectors, segments, instruments and signals.
    #[allow(clippy::cast_precision_loss)]
    fn diagram(&self) -> Drawing {
        let positions = graph::positions(self.flowsheet);
        let mut nodes = Vec::new();
        for instance in &self.flowsheet.instances {
            let Some(object_id) = self.item_ids.get(&instance.id) else {
                continue;
            };
            let detail = self
                .items
                .iter()
                .find(|item| item.id.as_deref() == Some(object_id.as_str()))
                .map(|item| item.ty.rsplit('.').next().unwrap_or(&item.ty).to_string())
                .or_else(|| {
                    self.mapped
                        .get(&instance.id)
                        .map(|ty| plant::class_of(ty).to_string())
                });
            let position = positions.get(&instance.id);
            nodes.push(DiagramNode {
                key: instance.id.clone(),
                object_id: object_id.clone(),
                label: instance.id.clone(),
                detail,
                shape: if plant::is_segment_unit(&instance.unit) {
                    NodeShape::Valve
                } else {
                    NodeShape::Unit
                },
                x: position.map_or(0.0, |position| position.x),
                y: position.map_or(0.0, |position| position.y),
            });
        }
        for input in &self.flowsheet.inputs {
            let Some(object_id) = self.boundaries.get(&input.name) else {
                continue;
            };
            let position = positions.get(&input.name);
            nodes.push(DiagramNode {
                key: input.name.clone(),
                object_id: object_id.clone(),
                label: input.name.clone(),
                detail: Some("feed".to_string()),
                shape: NodeShape::Input,
                x: position.map_or(0.0, |position| position.x),
                y: position.map_or(0.0, |position| position.y),
            });
        }
        for product in &self.flowsheet.products {
            let Some(object_id) = self.boundaries.get(product) else {
                continue;
            };
            let position = positions.get(product);
            nodes.push(DiagramNode {
                key: product.clone(),
                object_id: object_id.clone(),
                label: product.clone(),
                detail: Some("product".to_string()),
                shape: NodeShape::Output,
                x: position.map_or(0.0, |position| position.x),
                y: position.map_or(0.0, |position| position.y),
            });
        }

        let max_y = positions
            .values()
            .map(|position| position.y)
            .fold(0.0, f64::max);
        let mut slots: BTreeMap<String, usize> = BTreeMap::new();
        for (index, instrument) in self.flowsheet.instruments.iter().enumerate() {
            let anchor = instrument
                .measures
                .as_deref()
                .or(instrument.actuates.as_deref())
                .map(|target| split_endpoint(target).map_or(target, |(instance, _, _)| instance));
            let (x, y) = match anchor.and_then(|name| positions.get(name).map(|p| (name, p))) {
                Some((name, position)) => {
                    let slot = slots.entry(name.to_string()).or_insert(0);
                    let offset = *slot;
                    *slot += 1;
                    (position.x + 168.0, position.y - 20.0 + 50.0 * offset as f64)
                }
                None => (0.0, max_y + 80.0 + 50.0 * index as f64),
            };
            nodes.push(DiagramNode {
                key: format!("instrument:{}", instrument.id),
                object_id: format!("Instrument{}", index + 1),
                label: instrument.id.clone(),
                detail: Some(instrument.function.clone()),
                shape: NodeShape::Instrument,
                x,
                y,
            });
        }

        let mut edges = Vec::new();
        for segment in &self.segment_edges {
            let label = (segment.line != DEFAULT_LINE).then(|| segment.line.clone());
            edges.push(DiagramEdge {
                object_id: segment.object_id.clone(),
                source: diagram_key(&segment.from),
                target: diagram_key(&segment.to),
                label,
                dashed: segment.recycle,
            });
        }
        let mut signal_number = 0_usize;
        for instrument in &self.flowsheet.instruments {
            for signal in self
                .flowsheet
                .signals
                .iter()
                .filter(|signal| signal.from == instrument.id)
            {
                signal_number += 1;
                edges.push(DiagramEdge {
                    object_id: format!("Signal{signal_number}"),
                    source: format!("instrument:{}", signal.from),
                    target: format!("instrument:{}", signal.to),
                    label: signal.kind.clone(),
                    dashed: true,
                });
            }
        }
        drawing_export::build(&self.flowsheet.name, &nodes, &edges)
    }

    fn finish(mut self) -> Document {
        let diagram = self.options.drawing.then(|| self.diagram());
        let mut plant_children = vec![Child::Components(Components {
            property: "TaggedPlantItems".to_string(),
            children: self.items.into_iter().map(ComponentChild::Object).collect(),
        })];
        if !self.systems.is_empty() {
            let mut systems = Vec::new();
            for (index, (line, segments)) in self.systems.into_iter().enumerate() {
                systems.push(ComponentChild::Object(Object {
                    id: Some(format!("PipingNetworkSystem{}", index + 1)),
                    name: None,
                    ty: "Plant/Piping.PipingNetworkSystem".to_string(),
                    children: vec![
                        Child::Data(string_data("LineNumber", &line)),
                        Child::Components(Components {
                            property: "Segments".to_string(),
                            children: segments.into_iter().map(ComponentChild::Object).collect(),
                        }),
                    ],
                }));
            }
            plant_children.push(Child::Components(Components {
                property: "PipingNetworkSystems".to_string(),
                children: systems,
            }));
        }

        if !self.flowsheet.instruments.is_empty() {
            // The instrument objects first, so a signal can reference either end by object id.
            let by_tag: BTreeMap<&str, String> = self
                .flowsheet
                .instruments
                .iter()
                .enumerate()
                .map(|(index, instrument)| {
                    (instrument.id.as_str(), format!("Instrument{}", index + 1))
                })
                .collect();
            let mut built = Vec::new();
            let mut signal_number = 0usize;
            for instrument in &self.flowsheet.instruments {
                let id = by_tag[instrument.id.as_str()].clone();
                let mut children =
                    instrument_children(instrument, &self.item_ids, self.options.azoth_extensions);
                let mut signals: Vec<ComponentChild> = Vec::new();
                for signal in self
                    .flowsheet
                    .signals
                    .iter()
                    .filter(|signal| signal.from == instrument.id)
                {
                    signal_number += 1;
                    signals.push(ComponentChild::Object(signal_object(
                        signal,
                        &by_tag,
                        signal_number,
                        self.options.azoth_extensions,
                    )));
                }
                if !signals.is_empty() {
                    children.push(Child::Components(Components {
                        property: "SignalConveyingFunctions".to_string(),
                        children: signals,
                    }));
                }
                built.push(ComponentChild::Object(Object {
                    id: Some(id),
                    name: None,
                    ty: "Plant/Instrumentation.ProcessInstrumentationFunction".to_string(),
                    children,
                }));
            }
            plant_children.push(Child::Components(Components {
                property: "ProcessInstrumentationFunctions".to_string(),
                children: built,
            }));
        }

        let plant_model = Object {
            id: Some("PlantModel1".to_string()),
            name: None,
            ty: "Plant/PlantModel".to_string(),
            children: plant_children,
        };

        let mut engineering_children = Vec::new();
        if self.options.azoth_extensions {
            engineering_children.push(Child::Data(string_data(
                extensions::FLOWSHEET_ID,
                &self.flowsheet.id,
            )));
            engineering_children.push(Child::Data(string_data(
                extensions::FLOWSHEET_NAME,
                &self.flowsheet.name,
            )));
        }
        engineering_children.push(Child::Components(Components {
            property: "ConceptualModel".to_string(),
            children: vec![ComponentChild::Object(plant_model)],
        }));
        if let Some(drawing) = diagram {
            let diagram = write_diagram(&drawing, "PlantModel1", &mut self.mint);
            engineering_children.push(Child::Components(Components {
                property: "Diagram".to_string(),
                children: vec![ComponentChild::Object(diagram)],
            }));
        }
        let engineering_model = Object {
            id: None,
            name: None,
            ty: "Core/EngineeringModel".to_string(),
            children: engineering_children,
        };

        let name = model_name(&self.flowsheet.id);
        Document {
            model: Model {
                uri: format!("urn:azoth:dexpi:2.0:plant:{name}"),
                name,
                imports: Profile::Plant
                    .imports()
                    .into_iter()
                    .map(|(prefix, source)| Import {
                        prefix: prefix.to_string(),
                        source: source.to_string(),
                    })
                    .collect(),
                children: vec![Child::Object(engineering_model)],
            },
        }
    }
}

/// Join the two streams that reach a valve's ports into one segment.
///
/// A valve's inlet and outlet are two connections in a flowsheet and one segment in a P&ID.
/// The contraction is repeated so a chain of valves becomes one segment with several items;
/// a valve reached on only one side is refused, because a stream cannot cross it.
fn contract_valves(
    edges: &mut Vec<Edge>,
    valves: &BTreeMap<String, ValvePlacement>,
    fatal: &mut Vec<DexpiDiagnostic>,
) {
    loop {
        let mut merged = false;
        for (instance, valve) in valves {
            let incoming = edges
                .iter()
                .position(|edge| endpoint_instance(&edge.to) == Some(instance.as_str()));
            let outgoing = edges
                .iter()
                .position(|edge| endpoint_instance(&edge.from) == Some(instance.as_str()));
            let (Some(incoming), Some(outgoing)) = (incoming, outgoing) else {
                continue;
            };
            if incoming == outgoing {
                continue;
            }
            let into = edges[incoming].clone();
            let out = edges[outgoing].clone();
            let (high, low) = if incoming > outgoing {
                (incoming, outgoing)
            } else {
                (outgoing, incoming)
            };
            edges.remove(high);
            edges.remove(low);
            let mut through = into.valves;
            through.push(valve.clone());
            through.extend(out.valves);
            edges.push(Edge {
                from: into.from,
                to: out.to,
                line: merged_line(&into.line, &out.line),
                recycle: into.recycle.or(out.recycle),
                valves: through,
            });
            merged = true;
            break;
        }
        if !merged {
            break;
        }
    }
    for edge in edges.iter() {
        for endpoint in [&edge.from, &edge.to] {
            if let Some(instance) = endpoint_instance(endpoint) {
                if valves.contains_key(instance) {
                    fatal.push(DexpiDiagnostic::error(
                        "unpaired_segment_unit",
                        instance,
                        "a valve needs one inlet and one outlet stream to sit in a segment",
                    ));
                }
            }
        }
    }
}

fn endpoint_instance(endpoint: &str) -> Option<&str> {
    split_endpoint(endpoint).map(|(instance, _, _)| instance)
}

/// A segment carries one line; where a valve joins two lines, the incoming one is kept.
fn merged_line(incoming: &str, outgoing: &str) -> String {
    if incoming == DEFAULT_LINE {
        outgoing.to_string()
    } else {
        incoming.to_string()
    }
}

enum Node {
    Existing(String),
    Boundary { node: String, item: Object },
}

/// One pipe between two nodes, as the segment's `Connections` holds it.
fn pipe_object(number: usize, source: &str, target: &str) -> ComponentChild {
    ComponentChild::Object(Object {
        id: Some(format!("Pipe{number}")),
        name: None,
        ty: "Plant/Piping.Pipe".to_string(),
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
    })
}

/// The drawing node an endpoint names: the instance, or the boundary name.
fn diagram_key(endpoint: &str) -> String {
    split_endpoint(endpoint)
        .map_or(endpoint, |(instance, _, _)| instance)
        .to_string()
}

fn string_data(property: &str, value: &str) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::String(value.to_string())],
    }
}

fn number_data(property: &str, value: f64) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::Double(value)],
    }
}

/// The DEXPI `ProcessInstrumentationFunction`'s properties, with the parts the model has no
/// field for carried in `Azoth/…` extensions.
fn instrument_children(
    instrument: &Instrument,
    item_ids: &BTreeMap<String, String>,
    extensions: bool,
) -> Vec<Child> {
    let mut children = Vec::new();
    let mut letters = instrument.function.chars();
    if let Some(category) = letters.next() {
        children.push(Child::Data(string_data(
            "ProcessInstrumentationFunctionCategory",
            &category.to_string(),
        )));
        let modifier: String = letters.collect();
        if !modifier.is_empty() {
            children.push(Child::Data(string_data(
                "ProcessInstrumentationFunctionModifier",
                &modifier,
            )));
        }
    }
    let number = instrument
        .loop_number
        .clone()
        .or_else(|| trailing_digits(&instrument.id));
    if let Some(number) = number.filter(|number| !number.is_empty()) {
        children.push(Child::Data(string_data(
            "ProcessInstrumentationFunctionNumber",
            &number,
        )));
    }
    if let Some(reference) = instrument.location.as_deref().and_then(location_reference) {
        children.push(Child::Data(Data {
            property: "Location".to_string(),
            values: vec![DataValue::Reference(reference.to_string())],
        }));
    }
    if let Some(safety) = &instrument.safety {
        children.push(Child::Data(string_data(
            "SafetyRelevanceClass",
            &format!("SIL {}", safety.sil),
        )));
        if let Some(voting) = &safety.voting {
            children.push(Child::Data(string_data(
                "VotingSystemRepresentation",
                voting,
            )));
        }
        if extensions {
            children.push(Child::Data(string_data(
                extensions::SAFETY_FUNCTION,
                &safety.function,
            )));
            if let Some(state) = &safety.safe_state {
                children.push(Child::Data(string_data(extensions::SAFE_STATE, state)));
            }
            if let Some(hours) = safety.proof_test_hours {
                children.push(Child::Data(number_data(
                    extensions::PROOF_TEST_HOURS,
                    hours,
                )));
            }
            if let Some(seconds) = safety.closure_seconds {
                children.push(Child::Data(number_data(
                    extensions::CLOSURE_SECONDS,
                    seconds,
                )));
            }
        }
    }
    if extensions {
        children.push(Child::Data(string_data(extensions::TAG, &instrument.id)));
        if let Some(service) = &instrument.service {
            children.push(Child::Data(string_data(extensions::SERVICE, service)));
        }
        if let Some(measures) = &instrument.measures {
            children.push(Child::Data(string_data(extensions::MEASURES, measures)));
        }
        if let Some(actuates) = &instrument.actuates {
            children.push(Child::Data(string_data(extensions::ACTUATES, actuates)));
        }
        if let Some(range) = &instrument.range {
            children.push(Child::Data(number_data(extensions::RANGE_MIN, range.min)));
            children.push(Child::Data(number_data(extensions::RANGE_MAX, range.max)));
            children.push(Child::Data(string_data(
                extensions::RANGE_UNIT,
                &range.unit,
            )));
        }
    }
    for (property, target) in [
        ("SensingLocation", &instrument.measures),
        ("ActuatingLocation", &instrument.actuates),
    ] {
        if let Some(target) = target {
            if let Some(item) = item_for(target, item_ids) {
                children.push(Child::References(References {
                    property: property.to_string(),
                    objects: vec![format!("#{item}")],
                }));
            }
        }
    }
    children
}

fn signal_object(
    signal: &Signal,
    by_tag: &BTreeMap<&str, String>,
    number: usize,
    _extensions: bool,
) -> Object {
    let mut children = Vec::new();
    if let Some(reference) = signal.kind.as_deref().and_then(signal_conveying_reference) {
        children.push(Child::Data(Data {
            property: "SignalConveyingType".to_string(),
            values: vec![DataValue::Reference(reference.to_string())],
        }));
    }
    for (property, tag) in [("Source", &signal.from), ("Target", &signal.to)] {
        if let Some(id) = by_tag.get(tag.as_str()) {
            children.push(Child::References(References {
                property: property.to_string(),
                objects: vec![format!("#{id}")],
            }));
        }
    }
    Object {
        id: Some(format!("Signal{number}")),
        name: None,
        ty: "Plant/Instrumentation.SignalConveyingFunction".to_string(),
        children,
    }
}

/// The item an endpoint's instance is, where the document has one.
fn item_for<'a>(target: &str, item_ids: &'a BTreeMap<String, String>) -> Option<&'a str> {
    let instance = split_endpoint(target).map_or(target, |(instance, _, _)| instance);
    item_ids.get(instance).map(String::as_str)
}

fn location_reference(location: &str) -> Option<&'static str> {
    match location {
        "field" => Some("Plant/Enumerations.LocationClassification.Field"),
        "panel" => Some("Plant/Enumerations.LocationClassification.ControlPanel"),
        "central" => Some("Plant/Enumerations.LocationClassification.CentralLocation"),
        _ => None,
    }
}

fn signal_conveying_reference(kind: &str) -> Option<&'static str> {
    match kind {
        "electric" => {
            Some("Plant/Enumerations.SignalConveyingTypeClassification.ElectricalSignalConveying")
        }
        "pneumatic" => {
            Some("Plant/Enumerations.SignalConveyingTypeClassification.PneumaticSignalConveying")
        }
        "hydraulic" => {
            Some("Plant/Enumerations.SignalConveyingTypeClassification.HydraulicSignalConveying")
        }
        "capillary" => {
            Some("Plant/Enumerations.SignalConveyingTypeClassification.CapillarySignalConveying")
        }
        "radiation" => Some(
            "Plant/Enumerations.SignalConveyingTypeClassification.ConductedRadiationSignalConveying",
        ),
        _ => None,
    }
}

fn trailing_digits(tag: &str) -> Option<String> {
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
