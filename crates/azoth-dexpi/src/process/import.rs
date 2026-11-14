//! A DEXPI 2.0 Process model as a flowsheet scaffold.
//!
//! **Import cannot produce a runnable flowsheet, and does not pretend to.** DEXPI carries no
//! composition, so mass flow cannot become molar flow and a feed record cannot be written. The
//! importer therefore emits the instances, connections and products it can prove, leaves every
//! feed *undeclared* (the checker reports `unknown_feed`), and gives each step no parameters
//! (the checker reports `missing_parameter`). The document's own `azoth check` is the statement
//! that it is not runnable yet, and the report names every site.
//!
//! Port names are the lossy part. DEXPI's material ports are connection-scoped, not named the
//! way a palette's ports are. Azoth's own extension properties recover them exactly; without
//! them the importer assigns by direction and document order against the palette, and says so.

use std::collections::{BTreeMap, BTreeSet};

use azoth_process::flowsheet::{Connection, Flowsheet, Instance, Recycle};
use azoth_process::{Direction, Multiplicity, UnitOpSpec};

use crate::diagnostic::{DexpiDiagnostic, Severity};
use crate::ids::{model_name, sanitize_instance_id, uniquify};
use crate::mapping;
use crate::model::{Child, ComponentChild, DataValue, Document, Object};
use crate::process::extensions;
use crate::reader;
use crate::validate::validate;

/// How strict an import is.
#[derive(Debug, Clone, Copy, Default)]
pub struct ImportOptions {
    /// Turn every warning into a refusal.
    pub strict: bool,
}

/// The flowsheet an exchange became, and what it cost.
#[derive(Debug, Clone)]
pub struct ImportReport {
    pub flowsheet: Flowsheet,
    pub diagnostics: Vec<DexpiDiagnostic>,
    pub document: Document,
}

/// Read a DEXPI Process exchange as a flowsheet scaffold.
///
/// `palette` is what names the ports: a DEXPI material port is not named the way an azoth port
/// is, so an exchange without the azoth extension properties is mapped by direction and order
/// against the palette's own port declarations.
///
/// # Errors
/// The diagnostics that refused the import: a document that will not read, a document that fails
/// the semantic rules, or - under [`ImportOptions::strict`] - any warning.
pub fn import_process(
    xml: &str,
    palette: &[UnitOpSpec],
    options: &ImportOptions,
) -> Result<ImportReport, Vec<DexpiDiagnostic>> {
    let read = reader::read(xml).map_err(|diagnostic| vec![diagnostic])?;
    let validation = validate(&read.document);
    if !validation.is_valid() {
        return Err(validation.diagnostics);
    }
    let process_model = find_process_model(&read.document).ok_or_else(|| {
        vec![DexpiDiagnostic::error(
            "missing_process_model",
            "",
            "the document holds no `Process/ProcessModel`",
        )]
    })?;

    let mut diagnostics = Vec::new();
    let mut taken = BTreeSet::new();
    let steps = read_steps(process_model, palette, &mut diagnostics, &mut taken);
    let step_index: BTreeMap<&str, usize> = steps
        .iter()
        .enumerate()
        .filter_map(|(index, step)| step.object_id.as_deref().map(|id| (id, index)))
        .collect();
    let ports = read_ports(process_model, &step_index);
    let (streams, mut names) = read_streams(process_model, &ports, &mut diagnostics);
    name_ports(&steps, &streams, palette, &mut names, &mut diagnostics);

    let mut connections = Vec::new();
    let mut touched_by_recycle: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (index, stream) in streams.iter().enumerate() {
        match (
            &steps[stream.source.step].kind,
            &steps[stream.target.step].kind,
        ) {
            (StepKind::Recycle, _) => {
                touched_by_recycle
                    .entry(stream.source.step)
                    .or_default()
                    .push(index);
                continue;
            }
            (_, StepKind::Recycle) => {
                touched_by_recycle
                    .entry(stream.target.step)
                    .or_default()
                    .push(index);
                continue;
            }
            _ => {}
        }
        let from = endpoint(
            &steps[stream.source.step],
            names.get(&(index, Role::Source)),
        );
        let to = endpoint(
            &steps[stream.target.step],
            names.get(&(index, Role::Target)),
        );
        connections.push(Connection {
            from,
            to,
            line: None,
        });
    }

    let mut recycles = Vec::new();
    for (step, indices) in &touched_by_recycle {
        let incoming: Vec<usize> = indices
            .iter()
            .copied()
            .filter(|index| streams[*index].target.step == *step)
            .collect();
        let outgoing: Vec<usize> = indices
            .iter()
            .copied()
            .filter(|index| streams[*index].source.step == *step)
            .collect();
        let (Some(incoming), Some(outgoing)) = (incoming.first(), outgoing.first()) else {
            diagnostics.push(DexpiDiagnostic::warning(
                "recycle_incomplete",
                steps[*step].identifier.clone(),
                "a recycle marker needs one incoming and one outgoing stream; the loop is not \
                 reconstructed",
            ));
            continue;
        };
        let from = endpoint(
            &steps[streams[*incoming].source.step],
            names.get(&(*incoming, Role::Source)),
        );
        let to = endpoint(
            &steps[streams[*outgoing].target.step],
            names.get(&(*outgoing, Role::Target)),
        );
        recycles.push(Recycle::new(&steps[*step].identifier, &from, &to));
    }
    if !recycles.is_empty() {
        diagnostics.push(DexpiDiagnostic::info(
            "recycle_detected",
            "",
            format!(
                "{} loop(s) were read through an `Azoth.Recycle` marker",
                recycles.len()
            ),
        ));
    }

    for step in steps.iter().filter(|step| step.kind == StepKind::Source) {
        diagnostics.push(DexpiDiagnostic::warning(
            "input_composition_missing",
            step.identifier.clone(),
            "DEXPI carries no composition, so this feed is left undeclared: add an [[inputs]] \
             record with components, z, n, P and T, and `azoth check` stops reporting \
             `unknown_feed`",
        ));
        diagnostics.push(DexpiDiagnostic::warning(
            "input_flow_unresolved",
            step.identifier.clone(),
            "mass flow cannot become molar flow without a composition, so the feed's `n` is not \
             recovered",
        ));
    }

    let instances = steps
        .iter()
        .filter_map(|step| match &step.kind {
            StepKind::Unit(unit) => Some(Instance {
                id: step.instance_id.clone()?,
                unit: unit.clone(),
                parameters: BTreeMap::new(),
            }),
            _ => None,
        })
        .collect();
    let products = steps
        .iter()
        .filter(|step| step.kind == StepKind::Sink)
        .map(|step| step.identifier.clone())
        .collect();

    let root = find_engineering_model(&read.document);
    let flowsheet = Flowsheet {
        id: root
            .and_then(|root| root.data_string(extensions::FLOWSHEET_ID))
            .map_or_else(
                || model_name(&read.document.model.name),
                |id| id.to_string(),
            ),
        name: root
            .and_then(|root| root.data_string(extensions::FLOWSHEET_NAME))
            .map_or_else(|| read.document.model.name.clone(), |name| name.to_string()),
        products,
        inputs: Vec::new(),
        instances,
        connections,
        recycles,
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

    Ok(ImportReport {
        flowsheet,
        diagnostics,
        document: read.document,
    })
}

#[derive(Debug, Clone, PartialEq, Eq)]
enum StepKind {
    Source,
    Sink,
    Recycle,
    Unit(String),
}

struct Step {
    object_id: Option<String>,
    identifier: String,
    kind: StepKind,
    instance_id: Option<String>,
}

struct PortRow {
    step: usize,
    explicit: Option<String>,
    explicit_index: Option<usize>,
}

struct StreamRec {
    source: PortRef,
    target: PortRef,
}

struct PortRef {
    step: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
enum Role {
    Source,
    Target,
}

type Names = BTreeMap<(usize, Role), (String, Option<usize>)>;

fn read_steps(
    process_model: &Object,
    _palette: &[UnitOpSpec],
    diagnostics: &mut Vec<DexpiDiagnostic>,
    taken: &mut BTreeSet<String>,
) -> Vec<Step> {
    let mut steps = Vec::new();
    let Some(component) = process_model.components("ProcessSteps") else {
        return steps;
    };
    for child in &component.children {
        let ComponentChild::Object(step) = child else {
            continue;
        };
        let identifier = step
            .data_string("Identifier")
            .or(step.id.as_deref())
            .unwrap_or_default()
            .to_string();
        let kind = if step.data_string(extensions::RECYCLE).is_some() {
            StepKind::Recycle
        } else {
            match step.ty.as_str() {
                "Process/Process.Source" => StepKind::Source,
                "Process/Process.Sink" => StepKind::Sink,
                ty => StepKind::Unit(unit_for(step, ty, &identifier, diagnostics)),
            }
        };
        let instance_id = match &kind {
            StepKind::Unit(_) => {
                let sanitized = sanitize_instance_id(&identifier);
                let id = uniquify(&sanitized, taken);
                if id != sanitized {
                    diagnostics.push(DexpiDiagnostic::warning(
                        "renamed_on_import",
                        identifier.clone(),
                        format!(
                            "`{identifier}` collides with another name; this instance is `{id}`"
                        ),
                    ));
                }
                Some(id)
            }
            _ => None,
        };
        steps.push(Step {
            object_id: step.id.clone(),
            identifier,
            kind,
            instance_id,
        });
    }
    steps
}

fn unit_for(
    step: &Object,
    dexpi_type: &str,
    identifier: &str,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) -> String {
    if let Some(unit) = step.data_string(extensions::UNIT_OP) {
        if mapping::for_azoth(unit).is_some() {
            return unit.to_string();
        }
        diagnostics.push(DexpiDiagnostic::warning(
            "unknown_unit_op",
            identifier.to_string(),
            format!(
                "`Azoth.UnitOp` names `{unit}`, which has no mapping row; the type is used instead"
            ),
        ));
    }
    let candidates = mapping::candidates_for(dexpi_type);
    if candidates.is_empty() {
        diagnostics.push(DexpiDiagnostic::warning(
            "unsupported_step_type",
            identifier.to_string(),
            format!(
                "`{dexpi_type}` has no reviewed azoth unit; the step is kept as `dexpi.unmapped`"
            ),
        ));
        return "dexpi.unmapped".to_string();
    }
    let chosen = preferred(dexpi_type, candidates);
    if candidates.len() > 1 {
        diagnostics.push(DexpiDiagnostic::warning(
            "ambiguous_step_type",
            identifier.to_string(),
            format!(
                "`{dexpi_type}` may stand for {candidates:?}; `{chosen}` is used, and an export \
                 with --with-azoth-extensions preserves the exact unit"
            ),
        ));
    }
    chosen.to_string()
}

fn preferred(dexpi_type: &str, candidates: &[&'static str]) -> &'static str {
    let preferred = match dexpi_type {
        "Process/Process.Distilling" => "unit_ops.distillation_column",
        "Process/Process.SeparatingByGravity" => "unit_ops.separator",
        "Process/Process.TransportingFluids" => "unit_ops.pipe",
        "Process/Process.SplittingMaterial" => "unit_ops.splitter",
        _ => candidates[0],
    };
    candidates
        .iter()
        .copied()
        .find(|candidate| *candidate == preferred)
        .unwrap_or(candidates[0])
}

fn read_ports(
    process_model: &Object,
    step_index: &BTreeMap<&str, usize>,
) -> BTreeMap<String, PortRow> {
    let mut ports = BTreeMap::new();
    let Some(component) = process_model.components("ProcessSteps") else {
        return ports;
    };
    for child in &component.children {
        let ComponentChild::Object(step) = child else {
            continue;
        };
        let Some(object_id) = step.id.as_deref() else {
            continue;
        };
        let Some(&index) = step_index.get(object_id) else {
            continue;
        };
        let Some(port_component) = step.components("Ports") else {
            continue;
        };
        for port in port_component
            .children
            .iter()
            .filter_map(|child| match child {
                ComponentChild::Object(object) => Some(object),
                _ => None,
            })
        {
            let Some(id) = port.id.clone() else {
                continue;
            };
            ports.insert(
                id,
                PortRow {
                    step: index,
                    explicit: port.data_string(extensions::PORT).map(ToString::to_string),
                    explicit_index: match port.data_value(extensions::INDEX) {
                        Some(DataValue::Integer(value)) => usize::try_from(*value).ok(),
                        _ => None,
                    },
                },
            );
        }
    }
    ports
}

fn read_streams(
    process_model: &Object,
    ports: &BTreeMap<String, PortRow>,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) -> (Vec<StreamRec>, Names) {
    let mut streams = Vec::new();
    let mut names = Names::new();
    let Some(component) = process_model.components("ProcessConnections") else {
        return (streams, names);
    };
    for child in &component.children {
        let ComponentChild::Object(stream) = child else {
            continue;
        };
        let Some(source) = stream
            .references("Source")
            .and_then(|references| references.objects.first())
        else {
            continue;
        };
        let Some(target) = stream
            .references("Target")
            .and_then(|references| references.objects.first())
        else {
            continue;
        };
        let source = source.trim_start_matches('#');
        let target = target.trim_start_matches('#');
        let (Some(source), Some(target)) = (ports.get(source), ports.get(target)) else {
            let label = stream
                .data_string("Label")
                .or_else(|| stream.data_string("Identifier"))
                .unwrap_or_default();
            diagnostics.push(DexpiDiagnostic::warning(
                "connection_dropped",
                label.to_string(),
                "the stream names a port the document does not hold",
            ));
            continue;
        };
        let index = streams.len();
        if let Some(explicit) = &source.explicit {
            names.insert(
                (index, Role::Source),
                (explicit.clone(), source.explicit_index),
            );
        }
        if let Some(explicit) = &target.explicit {
            names.insert(
                (index, Role::Target),
                (explicit.clone(), target.explicit_index),
            );
        }
        streams.push(StreamRec {
            source: PortRef { step: source.step },
            target: PortRef { step: target.step },
        });
    }
    (streams, names)
}

fn name_ports(
    steps: &[Step],
    streams: &[StreamRec],
    palette: &[UnitOpSpec],
    names: &mut Names,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) {
    let mut outgoing: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    let mut incoming: BTreeMap<usize, Vec<usize>> = BTreeMap::new();
    for (index, stream) in streams.iter().enumerate() {
        outgoing.entry(stream.source.step).or_default().push(index);
        incoming.entry(stream.target.step).or_default().push(index);
    }
    for (step_index, step) in steps.iter().enumerate() {
        let StepKind::Unit(unit) = &step.kind else {
            continue;
        };
        for (role, catalog, is_outgoing) in [
            (Role::Source, &outgoing, true),
            (Role::Target, &incoming, false),
        ] {
            let Some(indices) = catalog.get(&step_index) else {
                continue;
            };
            if indices
                .iter()
                .all(|index| names.contains_key(&(*index, role)))
            {
                // Every stream on this side named itself through the extension properties.
                continue;
            }
            let assigned = assign(
                unit,
                is_outgoing,
                indices.len(),
                palette,
                &step.identifier,
                diagnostics,
            );
            for (ordinal, stream_index) in indices.iter().enumerate() {
                if names.contains_key(&(*stream_index, role)) {
                    continue;
                }
                if let Some((name, index)) = assigned.get(ordinal) {
                    names.insert((*stream_index, role), (name.clone(), *index));
                }
            }
        }
    }
}

fn assign(
    unit: &str,
    outgoing: bool,
    count: usize,
    palette: &[UnitOpSpec],
    identifier: &str,
    diagnostics: &mut Vec<DexpiDiagnostic>,
) -> Vec<(String, Option<usize>)> {
    if count == 0 {
        return Vec::new();
    }
    let Some(spec) = palette.iter().find(|spec| spec.id == unit) else {
        diagnostics.push(DexpiDiagnostic::warning(
            "port_assignment_assumed",
            identifier.to_string(),
            format!("`{unit}` is not in the palette, so its ports are named `port`"),
        ));
        return vec![("port".to_string(), None); count];
    };
    let direction = if outgoing {
        Direction::Out
    } else {
        Direction::In
    };
    let ports: Vec<_> = spec
        .ports
        .iter()
        .filter(|port| port.direction == direction)
        .collect();
    match ports.as_slice() {
        [] => vec![("port".to_string(), None); count],
        // **The index names a position on either side of a `many` port.** A two-way split
        // returns `products[0]` and `products[1]`; a mixer's feed is a row of connection points
        // too, and `mix1.feed[0]` and `mix1.feed[1]` are two of them. An index on a `one` port
        // is the checker's refusal, not this function's.
        [port] if port.multiplicity == Multiplicity::Many => (0..count)
            .map(|ordinal| (port.name.clone(), Some(ordinal)))
            .collect(),
        [port] => vec![(port.name.clone(), None); count],
        many => {
            diagnostics.push(DexpiDiagnostic::warning(
                "port_assignment_assumed",
                identifier.to_string(),
                format!(
                    "`{unit}` declares {} {} ports and the exchange does not name which is \
                     which; assignment is by document order and unverified (a heat exchanger's \
                     hot and cold sides are the dangerous case)",
                    many.len(),
                    if outgoing { "outlet" } else { "inlet" }
                ),
            ));
            (0..count)
                .map(|ordinal| {
                    let port = many.get(ordinal).copied().unwrap_or(many[0]);
                    (port.name.clone(), None)
                })
                .collect()
        }
    }
}

fn endpoint(step: &Step, name: Option<&(String, Option<usize>)>) -> String {
    match &step.kind {
        StepKind::Source | StepKind::Sink | StepKind::Recycle => step.identifier.clone(),
        StepKind::Unit(_) => {
            let instance = step
                .instance_id
                .clone()
                .unwrap_or_else(|| sanitize_instance_id(&step.identifier));
            let (port, index) = name.cloned().unwrap_or_else(|| ("port".to_string(), None));
            match index {
                Some(index) => format!("{instance}.{port}[{index}]"),
                None => format!("{instance}.{port}"),
            }
        }
    }
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

fn find_process_model(document: &Document) -> Option<&Object> {
    for child in &document.model.children {
        let Child::Object(root) = child else {
            continue;
        };
        let Some(conceptual) = root.components("ConceptualModel") else {
            continue;
        };
        for model in conceptual.children.iter().filter_map(|child| match child {
            ComponentChild::Object(object) => Some(object),
            _ => None,
        }) {
            if model.ty == "Process/ProcessModel" {
                return Some(model);
            }
        }
    }
    None
}
