//! A flowsheet, and its run, as a DEXPI 2.0 Process model.
//!
//! The projection is deliberately plain: one DEXPI step per boundary or instance, one
//! `MaterialPort` per stream endpoint, one `Process.Stream` per material edge, and the run's
//! numbers as the quantities NeqSim's writer emits - kg/h, bar absolute and °C.
//!
//! **The run is the normal path.** A DEXPI Process exchange exists to carry process state, so a
//! caller who has not asked for topology only gets the flowsheet run first; a run that does not
//! converge is refused unless the caller says the untrusted numbers are wanted, and a stream the
//! run did not produce simply carries no quantity. Nothing is defaulted to zero.

use std::collections::BTreeMap;

use azoth_process::executor::session::{RunReport, run};
use azoth_process::flowsheet::{Flowsheet, split_endpoint};
use azoth_process::{ExecutionOrder, Severity, UnitOpSpec, validate};

use crate::diagnostic::DexpiDiagnostic;
use crate::ids::{Mint, model_name};
use crate::mapping::{self, Disposition};
use crate::model::{
    Child, ComponentChild, Components, Data, DataValue, Document, Import, Model, Object, References,
};
use crate::process::extensions;
use crate::profile::Profile;
use crate::validate::validate as validate_document;
use crate::writer::write;

/// The status an exchange carries when it has not been certified.
pub const APPROVAL_STATUS: &str = "NOT_A_DEXPI_EV_CERTIFICATE";
/// The round-trip status an exchange carries when no named CAE tool has been qualified.
pub const CAE_ROUND_TRIP_STATUS: &str = "QUALIFICATION_REQUIRED";

/// What an export is allowed to do.
///
/// The default is the strict one: every flag that relaxes a refusal is off, and a caller that
/// wants an unconverged run, a lossy mapping or a projected recycle says so.
#[derive(Debug, Clone, Default)]
pub struct ExportOptions {
    /// Include lossy mappings; refused ones still refuse.
    pub allow_lossy: bool,
    /// Per-palette-unit overrides for a refused mapping, as `unit_ops.x=<DEXPI type>`.
    pub overrides: BTreeMap<String, String>,
    /// Turn each `[[recycles]]` entry into a `TransportingFluids` step.
    pub project_recycles: bool,
    /// Write the `Azoth.*` properties that make an azoth round trip exact.
    pub azoth_extensions: bool,
    /// Write topology only, without running the flowsheet.
    pub topology_only: bool,
    /// Publish an unconverged run's numbers rather than refusing.
    pub allow_unconverged: bool,
}

/// A completed export.
#[derive(Debug, Clone)]
pub struct ExportReport {
    pub xml: String,
    pub diagnostics: Vec<DexpiDiagnostic>,
    /// Whether the run converged; `false` for a topology-only exchange.
    pub converged: bool,
    /// Always [`APPROVAL_STATUS`]: schema conformance is not a DEXPI EV certificate.
    pub approval_status: &'static str,
    /// Always [`CAE_ROUND_TRIP_STATUS`]: no named CAE tool has been imported into.
    pub cae_round_trip: &'static str,
    /// Whether the built document passed the semantic validator. The XSD is a separate gate.
    pub semantic_valid: bool,
    /// Process steps written.
    pub steps: usize,
    /// Material streams written.
    pub streams: usize,
}

/// Export a flowsheet as a DEXPI Process exchange.
///
/// # Errors
/// The diagnostics that refused the export: a flowsheet the checker rejects, a unit with no
/// reviewed mapping, a lossy mapping without [`ExportOptions::allow_lossy`], a recycle without
/// [`ExportOptions::project_recycles`], or a run that did not converge.
pub fn export_process(
    flowsheet: &Flowsheet,
    palette: &[UnitOpSpec],
    options: &ExportOptions,
) -> Result<ExportReport, Vec<DexpiDiagnostic>> {
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

    let run = if options.topology_only {
        None
    } else {
        let report = run(
            flowsheet,
            palette,
            &BTreeMap::new(),
            ExecutionOrder::Topological,
        )
        .map_err(|error| {
            vec![DexpiDiagnostic::error(
                "run_failed",
                flowsheet.id.clone(),
                error.to_string(),
            )]
        })?;
        if !report.converged && !options.allow_unconverged {
            return Err(vec![DexpiDiagnostic::error(
                "unconverged_run",
                flowsheet.id.clone(),
                "the run did not converge, so its numbers are not evidence; export with \
                 --allow-unconverged to publish them with a warning anyway",
            )]);
        }
        if !report.converged {
            diagnostics.push(DexpiDiagnostic::warning(
                "unconverged_run",
                flowsheet.id.clone(),
                "the run did not converge and its numbers are published anyway",
            ));
        }
        Some(report)
    };

    let mut builder = Builder::new(flowsheet, options);
    builder.add_boundaries();
    builder.add_instances(&mapped);
    builder.add_edges(&mut fatal);
    if !fatal.is_empty() {
        return Err(fatal);
    }
    if let Some(run) = &run {
        builder.add_quantities(run, &mut diagnostics);
    }
    let document = builder.finish();
    let counts = crate::inventory::inventory(&document);
    let steps = counts
        .get("step.count")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0);
    let streams = counts
        .get("stream.count")
        .and_then(|value| value.parse().ok())
        .unwrap_or(0);

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

    Ok(ExportReport {
        xml,
        diagnostics,
        converged: run.as_ref().is_some_and(|report| report.converged),
        approval_status: APPROVAL_STATUS,
        cae_round_trip: CAE_ROUND_TRIP_STATUS,
        semantic_valid: true,
        steps,
        streams,
    })
}

/// The DEXPI type each instance is exchanged as.
fn resolve_mappings(
    flowsheet: &Flowsheet,
    options: &ExportOptions,
    fatal: &mut Vec<DexpiDiagnostic>,
) -> BTreeMap<String, String> {
    let mut mapped = BTreeMap::new();
    for instance in &flowsheet.instances {
        let Some(row) = mapping::for_azoth(&instance.unit) else {
            fatal.push(DexpiDiagnostic::error(
                "unmapped_unit",
                instance.id.clone(),
                format!(
                    "`{}` has no row in specs/dexpi/process_units.toml",
                    instance.unit
                ),
            ));
            continue;
        };
        if let Some(override_type) = options.overrides.get(&instance.unit) {
            if !crate::profile::supports(override_type) {
                fatal.push(DexpiDiagnostic::error(
                    "override_outside_profile",
                    instance.id.clone(),
                    format!("`{override_type}` is outside the reviewed Process profile"),
                ));
                continue;
            }
            mapped.insert(instance.id.clone(), override_type.clone());
            continue;
        }
        match row.disposition {
            Disposition::Faithful => {}
            Disposition::Lossy if options.allow_lossy => {}
            Disposition::Lossy => {
                fatal.push(DexpiDiagnostic::error(
                    "lossy_mapping",
                    instance.id.clone(),
                    format!(
                        "`{}` maps to `{}` but loses something: {} Export with \
                         --allow-lossy-mappings, or override the mapping with --map.",
                        instance.unit,
                        row.dexpi_types.first().copied().unwrap_or_default(),
                        row.note
                    ),
                ));
                continue;
            }
            Disposition::Refused => {
                fatal.push(DexpiDiagnostic::error(
                    "refused_mapping",
                    instance.id.clone(),
                    format!(
                        "`{}` has no reviewed DEXPI Process type: {} An export needs an \
                         explicit `--map {}=<type>` and that override is recorded in the report.",
                        instance.unit, row.note, instance.unit
                    ),
                ));
                continue;
            }
        }
        let Some(dexpi_type) = row.dexpi_types.first() else {
            fatal.push(DexpiDiagnostic::error(
                "unmapped_unit",
                instance.id.clone(),
                format!("`{}` names no DEXPI type", instance.unit),
            ));
            continue;
        };
        mapped.insert(instance.id.clone(), (*dexpi_type).to_string());
    }
    mapped
}

struct Builder<'a> {
    flowsheet: &'a Flowsheet,
    options: &'a ExportOptions,
    steps: BTreeMap<String, StepBuild>,
    order: Vec<String>,
    streams: Vec<Object>,
    mint: Mint,
    edges: Vec<Edge>,
}

struct StepBuild {
    identifier: String,
    ty: String,
    leaf_unit: Option<String>,
    azoth_instance: Option<String>,
    recycle: bool,
    ports: Vec<ComponentChild>,
}

struct Edge {
    /// The run key the stream's quantities are read from.
    quantity_key: String,
    stream_id: String,
    label: String,
}

impl<'a> Builder<'a> {
    fn new(flowsheet: &'a Flowsheet, options: &'a ExportOptions) -> Self {
        Self {
            flowsheet,
            options,
            steps: BTreeMap::new(),
            order: Vec::new(),
            streams: Vec::new(),
            mint: Mint::new(),
            edges: Vec::new(),
        }
    }

    fn add_step(&mut self, key: String, identifier: String, ty: String) {
        if self.steps.contains_key(&key) {
            return;
        }
        self.order.push(key.clone());
        self.steps.insert(
            key,
            StepBuild {
                identifier,
                ty,
                leaf_unit: None,
                azoth_instance: None,
                recycle: false,
                ports: Vec::new(),
            },
        );
    }

    fn add_boundaries(&mut self) {
        for input in &self.flowsheet.inputs {
            self.add_step(
                format!("input:{}", input.name),
                input.name.clone(),
                "Process/Process.Source".to_string(),
            );
        }
        for product in &self.flowsheet.products {
            self.add_step(
                format!("product:{product}"),
                product.clone(),
                "Process/Process.Sink".to_string(),
            );
        }
    }

    fn add_instances(&mut self, mapped: &BTreeMap<String, String>) {
        for instance in &self.flowsheet.instances {
            let Some(ty) = mapped.get(&instance.id) else {
                continue;
            };
            self.add_step(
                format!("instance:{}", instance.id),
                instance.id.clone(),
                ty.clone(),
            );
            if let Some(step) = self.steps.get_mut(&format!("instance:{}", instance.id)) {
                step.leaf_unit = Some(instance.unit.clone());
                step.azoth_instance = Some(instance.id.clone());
            }
        }
    }

    fn add_edges(&mut self, fatal: &mut Vec<DexpiDiagnostic>) {
        let connections: Vec<(String, String)> = self
            .flowsheet
            .connections
            .iter()
            .map(|connection| (connection.from.clone(), connection.to.clone()))
            .collect();
        for (from, to) in connections {
            self.push_edge(&from, &to, &from, &from, fatal);
        }
        for recycle in &self.flowsheet.recycles {
            if !self.options.project_recycles {
                fatal.push(DexpiDiagnostic::error(
                    "recycle_not_projected",
                    recycle.stream.clone(),
                    "DEXPI Process has no recycle; export with --project-recycles to write the \
                     tear as a TransportingFluids step, or remove the recycle",
                ));
                continue;
            }
            let key = format!("recycle:{}", recycle.stream);
            self.add_step(
                key.clone(),
                recycle.stream.clone(),
                "Process/Process.TransportingFluids".to_string(),
            );
            if let Some(step) = self.steps.get_mut(&key) {
                step.azoth_instance = Some(recycle.stream.clone());
                step.recycle = true;
            }
            let inlet = format!("recycle:{}.in", recycle.stream);
            let outlet = format!("recycle:{}.out", recycle.stream);
            self.push_edge(&recycle.from, &inlet, &recycle.stream, &recycle.from, fatal);
            self.push_edge(&outlet, &recycle.to, &recycle.stream, &recycle.from, fatal);
        }
    }

    fn push_edge(
        &mut self,
        from: &str,
        to: &str,
        label: &str,
        quantity_key: &str,
        fatal: &mut Vec<DexpiDiagnostic>,
    ) {
        let Some(source_key) = self.step_key(from, true) else {
            fatal.push(DexpiDiagnostic::error(
                "unresolved_endpoint",
                from,
                "the stream's source is neither a declared input nor an instance outlet",
            ));
            return;
        };
        let Some(target_key) = self.step_key(to, false) else {
            fatal.push(DexpiDiagnostic::error(
                "unresolved_endpoint",
                to,
                "the stream's target is neither a declared product nor an instance inlet",
            ));
            return;
        };
        let stream_id = self.mint.next("Stream");
        let source_port = format!("{stream_id}SourcePort");
        let target_port = format!("{stream_id}TargetPort");

        let source_endpoint = endpoint_parts(from);
        let target_endpoint = endpoint_parts(to);
        self.add_port(
            &source_key,
            &source_port,
            &stream_id,
            "Outlet",
            source_endpoint.as_ref(),
        );
        self.add_port(
            &target_key,
            &target_port,
            &stream_id,
            "Inlet",
            target_endpoint.as_ref(),
        );

        let mut children = vec![
            Child::Data(string_data("Identifier", &stream_id)),
            Child::Data(string_data("Label", label)),
            Child::References(References {
                property: "Source".to_string(),
                objects: vec![format!("#{source_port}")],
            }),
            Child::References(References {
                property: "Target".to_string(),
                objects: vec![format!("#{target_port}")],
            }),
        ];
        if self.options.azoth_extensions {
            children.push(Child::Data(string_data(extensions::FROM, from)));
            children.push(Child::Data(string_data(extensions::TO, to)));
        }
        self.streams.push(Object {
            id: Some(stream_id.clone()),
            name: None,
            ty: "Process/Process.Stream".to_string(),
            children,
        });
        self.edges.push(Edge {
            quantity_key: quantity_key.to_string(),
            stream_id,
            label: label.to_string(),
        });
    }

    fn add_port(
        &mut self,
        step_key: &str,
        port_id: &str,
        stream_id: &str,
        direction: &str,
        endpoint: Option<&Endpoint>,
    ) {
        let Some(step) = self.steps.get_mut(step_key) else {
            return;
        };
        let mut children = vec![
            Child::Data(string_data("Identifier", port_id)),
            Child::Data(Data {
                property: "NominalDirection".to_string(),
                values: vec![DataValue::Reference(format!(
                    "Process/Enumerations.PortDirection.{direction}"
                ))],
            }),
            Child::References(References {
                property: "ConnectorReference".to_string(),
                objects: vec![format!("#{stream_id}")],
            }),
        ];
        if self.options.azoth_extensions {
            if let Some(endpoint) = endpoint {
                if let Some(port) = &endpoint.port {
                    children.push(Child::Data(string_data(extensions::PORT, port)));
                }
                if let Some(index) = endpoint.index {
                    children.push(Child::Data(Data {
                        property: extensions::INDEX.to_string(),
                        values: vec![DataValue::Integer(index as i64)],
                    }));
                }
            }
        }
        step.ports.push(ComponentChild::Object(Object {
            id: Some(port_id.to_string()),
            name: None,
            ty: "Process/Process.MaterialPort".to_string(),
            children,
        }));
    }

    fn step_key(&self, endpoint: &str, source: bool) -> Option<String> {
        if let Some(rest) = endpoint.strip_prefix("recycle:") {
            let tear = rest.split('.').next().unwrap_or(rest);
            let key = format!("recycle:{tear}");
            return self.steps.contains_key(&key).then_some(key);
        }
        match split_endpoint(endpoint) {
            Some((instance, _, _)) => {
                let key = format!("instance:{instance}");
                self.steps.contains_key(&key).then_some(key)
            }
            None if source => {
                let key = format!("input:{endpoint}");
                self.steps.contains_key(&key).then_some(key)
            }
            None => {
                let key = format!("product:{endpoint}");
                self.steps.contains_key(&key).then_some(key)
            }
        }
    }

    fn add_quantities(&mut self, run: &RunReport, diagnostics: &mut Vec<DexpiDiagnostic>) {
        for index in 0..self.streams.len() {
            let key = self.edges[index].quantity_key.clone();
            let stream_id = self.edges[index].stream_id.clone();
            let label = self.edges[index].label.clone();
            let Some(stream) = run.streams.get(&key) else {
                diagnostics.push(DexpiDiagnostic::warning(
                    "operating_value_missing",
                    stream_id,
                    format!("the run produced no stream at `{key}` (label `{label}`), so the exchange carries none of its quantities"),
                ));
                continue;
            };
            let quantities = [
                (
                    "MassFlow",
                    stream.mass_flow().map(|value| value * 3600.0),
                    "Core/PhysicalQuantities.MassFlowRateUnit.KilogramPerHour",
                ),
                (
                    "Pressure",
                    Ok(stream.p.value / 1.0e5),
                    "Core/PhysicalQuantities.PressureAbsoluteUnit.Bar",
                ),
                (
                    "Temperature",
                    Ok(stream.t.value - 273.15),
                    "Core/PhysicalQuantities.TemperatureUnit.DegreeCelsius",
                ),
            ];
            for (property, value, unit) in quantities {
                match value {
                    Ok(value) if value.is_finite() => {
                        self.streams[index]
                            .children
                            .push(Child::Components(quantity_components(
                                property, value, unit,
                            )));
                    }
                    Ok(value) => diagnostics.push(DexpiDiagnostic::warning(
                        "non_finite_quantity",
                        stream_id.clone(),
                        format!("`{property}` is `{value}`, which is not a finite number"),
                    )),
                    Err(error) => diagnostics.push(DexpiDiagnostic::warning(
                        "operating_value_missing",
                        stream_id.clone(),
                        format!("`{property}`: {error}"),
                    )),
                }
            }
        }
    }

    fn finish(self) -> Document {
        let mut process_steps = Vec::new();
        for (ordinal, key) in self.order.iter().enumerate() {
            let Some(step) = self.steps.get(key) else {
                continue;
            };
            let mut children = vec![
                Child::Data(string_data("Identifier", &step.identifier)),
                Child::Data(string_data("Label", &step.identifier)),
            ];
            if self.options.azoth_extensions {
                if let Some(unit) = &step.leaf_unit {
                    children.push(Child::Data(string_data(extensions::UNIT_OP, unit)));
                }
                if let Some(instance) = &step.azoth_instance {
                    children.push(Child::Data(string_data(extensions::INSTANCE_ID, instance)));
                }
                if step.recycle {
                    children.push(Child::Data(string_data(
                        extensions::RECYCLE,
                        &step.identifier,
                    )));
                }
            }
            if !step.ports.is_empty() {
                children.push(Child::Components(Components {
                    property: "Ports".to_string(),
                    children: step.ports.clone(),
                }));
            }
            process_steps.push(ComponentChild::Object(Object {
                id: Some(format!("ProcessStep{}", ordinal + 1)),
                name: None,
                ty: step.ty.clone(),
                children,
            }));
        }

        let mut process_model_children = vec![Child::Components(Components {
            property: "ProcessSteps".to_string(),
            children: process_steps,
        })];
        if !self.streams.is_empty() {
            process_model_children.push(Child::Components(Components {
                property: "ProcessConnections".to_string(),
                children: self
                    .streams
                    .into_iter()
                    .map(ComponentChild::Object)
                    .collect(),
            }));
        }

        let process_model = Object {
            id: Some("ProcessModel1".to_string()),
            name: None,
            ty: "Process/ProcessModel".to_string(),
            children: process_model_children,
        };

        let mut engineering_children = Vec::new();
        if self.options.azoth_extensions {
            // The `Model` name is a DEXPI `name` and may carry no dot, so the flowsheet's own
            // id travels here; without the extension the name is the sanitised one.
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
            children: vec![ComponentChild::Object(process_model)],
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
                uri: format!("urn:azoth:dexpi:2.0:process:{name}"),
                name,
                imports: Profile::Process
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

struct Endpoint {
    port: Option<String>,
    index: Option<usize>,
}

fn endpoint_parts(endpoint: &str) -> Option<Endpoint> {
    let (_, port, index) = split_endpoint(endpoint)?;
    Some(Endpoint {
        port: Some(port.to_string()),
        index,
    })
}

fn string_data(property: &str, value: &str) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::String(value.to_string())],
    }
}

fn quantity_components(property: &str, magnitude: f64, unit: &str) -> Components {
    let quantity = Object {
        id: None,
        name: None,
        ty: "Core/QualifiedValue".to_string(),
        children: vec![Child::Data(Data {
            property: "Value".to_string(),
            values: vec![DataValue::Aggregated {
                ty: "Core/PhysicalQuantities.PhysicalQuantity".to_string(),
                data: vec![
                    Data {
                        property: "Unit".to_string(),
                        values: vec![DataValue::Reference(unit.to_string())],
                    },
                    Data {
                        property: "Value".to_string(),
                        values: vec![DataValue::Double(magnitude)],
                    },
                ],
            }],
        })],
    };
    Components {
        property: property.to_string(),
        children: vec![ComponentChild::Object(quantity)],
    }
}

/// The endpoint key a stream is published under, as the executor names it.
#[must_use]
pub fn endpoint_key(instance: &str, port: &str, index: Option<usize>) -> String {
    azoth_process::flowsheet::stream_path(instance, port, index)
}
