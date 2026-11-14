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

use azoth_process::flowsheet::{Flowsheet, split_endpoint};
use azoth_process::{Severity, UnitOpSpec, validate};

use crate::diagnostic::DexpiDiagnostic;
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
}

struct Edge {
    from: String,
    to: String,
    line: String,
    recycle: Option<String>,
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
        }
    }

    fn add_items(&mut self, mapped: &BTreeMap<String, String>) {
        for instance in &self.flowsheet.instances {
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
            });
        }

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
            if !boundary_items.is_empty() {
                children.push(Child::Components(Components {
                    property: "Items".to_string(),
                    children: boundary_items
                        .into_iter()
                        .map(ComponentChild::Object)
                        .collect(),
                }));
            }
            let pipe_id = format!("Pipe{number}");
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
                    children: vec![ComponentChild::Object(Object {
                        id: Some(pipe_id),
                        name: None,
                        ty: "Plant/Piping.Pipe".to_string(),
                        children: vec![
                            Child::References(References {
                                property: "SourceNode".to_string(),
                                objects: vec![format!("#{source_node}")],
                            }),
                            Child::References(References {
                                property: "TargetNode".to_string(),
                                objects: vec![format!("#{target_node}")],
                            }),
                        ],
                    })],
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
            let object = Object {
                id: Some(format!("PipingNetworkSegment{number}")),
                name: None,
                ty: "Plant/Piping.PipingNetworkSegment".to_string(),
                children: segment,
            };
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

    fn finish(self) -> Document {
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

enum Node {
    Existing(String),
    Boundary { node: String, item: Object },
}

fn string_data(property: &str, value: &str) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::String(value.to_string())],
    }
}
