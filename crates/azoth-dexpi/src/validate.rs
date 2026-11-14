//! The portable half of DEXPI's semantic validation, for the Process profile.
//!
//! This is what NeqSim's `Dexpi20SemanticValidator` checks that does not need a Java DOM: the
//! exact information-model imports, unique identities, resolvable references, the supported
//! type set, a Process model with at least one step and one stream, the port/stream
//! reciprocity, and finite, unit-qualified quantities.
//!
//! **It is not schema validation.** XML well-formedness and the XSD are `xmllint`'s job, run by
//! the test gate and by `azoth dexpi validate --xsd`; a document can pass this and fail the
//! schema, and both results are reported separately rather than folded into one "valid".

use std::collections::{BTreeMap, BTreeSet};

use crate::diagnostic::DexpiDiagnostic;
use crate::model::{Child, ComponentChild, DataValue, Document, Object};
use crate::profile::{self, CORE_MODEL, PROCESS_MODEL};
use crate::units;

/// Every finding the semantic rules produced.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct ValidationReport {
    pub diagnostics: Vec<DexpiDiagnostic>,
}

impl ValidationReport {
    /// Whether no finding is an error.
    #[must_use]
    pub fn is_valid(&self) -> bool {
        !self.diagnostics.iter().any(DexpiDiagnostic::is_error)
    }

    fn error(
        &mut self,
        code: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) {
        self.diagnostics
            .push(DexpiDiagnostic::error(code, subject, message));
    }

    fn warning(
        &mut self,
        code: &'static str,
        subject: impl Into<String>,
        message: impl Into<String>,
    ) {
        self.diagnostics
            .push(DexpiDiagnostic::warning(code, subject, message));
    }
}

/// Validate a Process-profile document against the semantic rules.
#[must_use]
pub fn validate(document: &Document) -> ValidationReport {
    let mut report = ValidationReport::default();
    let model = &document.model;

    if model.imports.is_empty() {
        report.error(
            "missing_import",
            "Model",
            "a DEXPI 2.0 exchange imports the Core and Process information models",
        );
    }
    let mut prefixes = BTreeSet::new();
    for import in &model.imports {
        if !prefixes.insert(import.prefix.clone()) {
            report.error(
                "duplicate_import",
                &import.prefix,
                "this import prefix is declared twice",
            );
        }
        let expected = match import.prefix.as_str() {
            "Core" => Some(CORE_MODEL),
            "Process" => Some(PROCESS_MODEL),
            _ => None,
        };
        match expected {
            Some(uri) if import.source == uri => {}
            Some(uri) => report.error(
                "wrong_import",
                &import.prefix,
                format!(
                    "`{}` must resolve to `{uri}` at this release",
                    import.source
                ),
            ),
            None => report.error(
                "unknown_import",
                &import.prefix,
                "the Process profile imports exactly `Core` and `Process`",
            ),
        }
    }
    for (prefix, uri) in profile::Profile::Process.imports() {
        if !model
            .imports
            .iter()
            .any(|import| import.prefix == prefix && import.source == uri)
        {
            report.error(
                "missing_import",
                prefix,
                format!("the Process profile requires the import `{prefix}` at `{uri}`"),
            );
        }
    }

    let mut objects: BTreeMap<String, &Object> = BTreeMap::new();
    let mut tag_names: BTreeSet<String> = BTreeSet::new();
    collect_objects(
        document,
        &model.children,
        &mut objects,
        &prefixes,
        &mut report,
        &mut tag_names,
    );

    // References resolve, against ids and names.
    for object in objects.values() {
        for child in &object.children {
            let Child::References(references) = child else {
                continue;
            };
            if references.objects.is_empty() {
                report.error(
                    "empty_reference",
                    object.id.clone().unwrap_or_default(),
                    format!("`{}` names no object", references.property),
                );
            }
            for target in &references.objects {
                if resolve(target, &objects).is_none() {
                    report.error(
                        "unresolved_reference",
                        object.id.clone().unwrap_or_default(),
                        format!(
                            "`{}` points at `{target}`, which the document does not hold",
                            references.property
                        ),
                    );
                }
            }
        }
    }

    // One Process model, with content.
    let models: Vec<&Object> = objects
        .values()
        .copied()
        .filter(|object| object.ty == "Process/ProcessModel")
        .collect();
    match models.len() {
        0 => report.error(
            "missing_process_model",
            "Model",
            "the document has no `Process/ProcessModel`",
        ),
        1 => {}
        count => report.error(
            "duplicate_process_model",
            "Model",
            format!(
                "the document has {count} `Process/ProcessModel` objects; one exchange holds one"
            ),
        ),
    }

    let steps: Vec<&Object> = objects
        .values()
        .copied()
        .filter(|object| is_step(object))
        .collect();
    let streams: Vec<&Object> = objects
        .values()
        .copied()
        .filter(|object| object.ty == "Process/Process.Stream")
        .collect();
    if steps.is_empty() {
        report.error(
            "missing_step",
            "ProcessModel",
            "a Process exchange requires at least one process step",
        );
    }
    if streams.is_empty() {
        report.error(
            "missing_stream",
            "ProcessModel",
            "a Process exchange requires at least one material stream",
        );
    }

    for step in &steps {
        if step.data_string("Identifier").is_none_or(str::is_empty) {
            report.error(
                "missing_identifier",
                step.id.clone().unwrap_or_default(),
                "a process step requires a non-empty `Identifier`",
            );
        }
        match step.components("Ports") {
            Some(ports) if !ports.children.is_empty() => {}
            _ => report.warning(
                "missing_ports",
                step.id.clone().unwrap_or_default(),
                "a tagged process step carries no material port",
            ),
        }
    }

    // Port <-> stream reciprocity.
    for stream in &streams {
        let subject = stream.id.clone().unwrap_or_default();
        if stream.data_string("Identifier").is_none_or(str::is_empty) {
            report.error(
                "missing_identifier",
                &subject,
                "a stream requires a non-empty `Identifier`",
            );
        }
        let sources = reference_targets(stream, "Source");
        let targets = reference_targets(stream, "Target");
        if sources.len() != 1 {
            report.error(
                "bad_stream_reference",
                &subject,
                format!(
                    "a stream has exactly one `Source`; this one has {}",
                    sources.len()
                ),
            );
        }
        if targets.len() != 1 {
            report.error(
                "bad_stream_reference",
                &subject,
                format!(
                    "a stream has exactly one `Target`; this one has {}",
                    targets.len()
                ),
            );
        }
        if sources.first() == targets.first() && !sources.is_empty() {
            report.error(
                "bad_stream_reference",
                &subject,
                "a stream's source and target are the same port",
            );
        }
        for (property, expected_direction) in [
            ("Source", profile::DIRECTION_OUTLET),
            ("Target", profile::DIRECTION_INLET),
        ] {
            for target in reference_targets(stream, property) {
                let Some(port) = resolve(target, &objects) else {
                    continue; // reported by the reference pass
                };
                if port.ty != "Process/Process.MaterialPort" {
                    report.error(
                        "bad_stream_reference",
                        &subject,
                        format!("`{property}` points at `{target}`, which is not a material port"),
                    );
                    continue;
                }
                let direction = port
                    .data_value("NominalDirection")
                    .and_then(|value| match value {
                        DataValue::Reference(reference) => Some(reference.as_str()),
                        _ => None,
                    });
                if direction != Some(expected_direction) {
                    report.error(
                        "bad_port_direction",
                        port.id.clone().unwrap_or_default(),
                        format!("a `{property}` port must be `{expected_direction}`"),
                    );
                }
                if !reference_targets(port, "ConnectorReference")
                    .iter()
                    .any(|target| target.trim_start_matches('#') == subject)
                {
                    report.error(
                        "missing_connector_reference",
                        port.id.clone().unwrap_or_default(),
                        format!(
                            "`{}` does not point back at stream `{subject}`",
                            port.id.clone().unwrap_or_default()
                        ),
                    );
                }
            }
        }
        for quantity in ["MassFlow", "Pressure", "Temperature"] {
            validate_quantity(stream, quantity, &subject, &mut report);
        }
    }

    for port in objects
        .values()
        .filter(|object| object.ty == "Process/Process.MaterialPort")
    {
        let subject = port.id.clone().unwrap_or_default();
        if port.data_string("Identifier").is_none_or(str::is_empty) {
            report.error(
                "missing_identifier",
                &subject,
                "a material port requires an `Identifier`",
            );
        }
        match port.data_value("NominalDirection") {
            Some(DataValue::Reference(reference))
                if reference == profile::DIRECTION_INLET
                    || reference == profile::DIRECTION_OUTLET => {}
            _ => report.error(
                "bad_port_direction",
                &subject,
                "a material port's `NominalDirection` is `Inlet` or `Outlet`",
            ),
        }
        if reference_targets(port, "ConnectorReference").is_empty() {
            report.error(
                "missing_connector_reference",
                &subject,
                "a material port names the stream it connects",
            );
        }
    }

    report
}

fn validate_quantity(
    stream: &Object,
    property: &str,
    subject: &str,
    report: &mut ValidationReport,
) {
    let Some(components) = stream.components(property) else {
        return;
    };
    for child in &components.children {
        let ComponentChild::Object(qualified) = child else {
            report.error(
                "bad_quantity",
                subject,
                format!(
                    "`{property}` holds an object reference; a quantity is a `Core/QualifiedValue`"
                ),
            );
            continue;
        };
        let Some(value) = qualified.data_value("Value") else {
            report.error(
                "bad_quantity",
                subject,
                format!("`{property}` has no `Value`"),
            );
            continue;
        };
        if value.quantity_unit().is_none() {
            report.error(
                "bad_quantity",
                subject,
                format!("`{property}` carries no `Core/PhysicalQuantities.PhysicalQuantity`"),
            );
            continue;
        }
        let unit = value.quantity_unit().unwrap_or_default();
        let Some(expected) = quantity_kind(property) else {
            continue;
        };
        match units::kind(unit) {
            Some(kind) if kind == expected => {}
            Some(_) => report.error(
                "wrong_quantity_unit",
                subject,
                format!("`{property}` carries `{unit}`, which measures something else"),
            ),
            None => {
                let message = match units::from_si(unit, 1.0) {
                    Err(error) => error.message(),
                    Ok(_) => String::new(),
                };
                report.error("unknown_quantity_unit", subject, message);
            }
        }
        if value.quantity_magnitude().is_none() {
            report.error(
                "non_finite_quantity",
                subject,
                format!("`{property}` carries a value that is not a finite number"),
            );
        }
    }
}

fn quantity_kind(property: &str) -> Option<units::QuantityKind> {
    match property {
        "MassFlow" => Some(units::QuantityKind::MassFlow),
        "Pressure" => Some(units::QuantityKind::PressureAbsolute),
        "Temperature" => Some(units::QuantityKind::Temperature),
        _ => None,
    }
}

fn is_step(object: &Object) -> bool {
    object.ty.starts_with("Process/Process.")
        && !matches!(
            object.ty.as_str(),
            "Process/Process.MaterialPort" | "Process/Process.Stream"
        )
}

fn reference_targets<'a>(object: &'a Object, property: &str) -> &'a [String] {
    object
        .references(property)
        .map_or(&[], |references| references.objects.as_slice())
}

fn resolve<'a>(target: &str, objects: &'a BTreeMap<String, &'a Object>) -> Option<&'a Object> {
    let local = target.trim_start_matches('#');
    if let Some(object) = objects.get(local) {
        return Some(object);
    }
    objects
        .values()
        .copied()
        .find(|object| object.name.as_deref() == Some(local))
}

fn collect_objects<'a>(
    _document: &'a Document,
    children: &'a [Child],
    objects: &mut BTreeMap<String, &'a Object>,
    prefixes: &BTreeSet<String>,
    report: &mut ValidationReport,
    tag_names: &mut BTreeSet<String>,
) {
    for child in children {
        if let Child::Object(object) = child {
            collect_object(object, objects, prefixes, report, tag_names);
        }
    }
}

fn collect_object<'a>(
    object: &'a Object,
    objects: &mut BTreeMap<String, &'a Object>,
    prefixes: &BTreeSet<String>,
    report: &mut ValidationReport,
    tag_names: &mut BTreeSet<String>,
) {
    let subject = object.id.clone().unwrap_or_default();
    match object.ty.split_once('/') {
        Some((prefix, _)) => {
            if !prefixes.contains(prefix) {
                report.error(
                    "undeclared_prefix",
                    &subject,
                    format!(
                        "type `{}` uses prefix `{prefix}`, which no import declares",
                        object.ty
                    ),
                );
            }
        }
        None => report.error(
            "bad_type",
            &subject,
            format!("`{}` is not a `prefix/name` reference", object.ty),
        ),
    }
    if !profile::supports(&object.ty) {
        report.error(
            "unsupported_type",
            &subject,
            format!("`{}` is outside the reviewed Process profile", object.ty),
        );
    }
    if let Some(id) = &object.id {
        if objects.insert(id.clone(), object).is_some() {
            report.error("duplicate_id", id, "this object id is used twice");
        }
    }
    if let Some(tag) = object.data_string("TagName") {
        if !tag_names.insert(tag.to_string()) {
            report.error("duplicate_tag_name", tag, "this TagName is used twice");
        }
    }
    for child in &object.children {
        if let Child::Components(components) = child {
            for nested in &components.children {
                if let ComponentChild::Object(nested) = nested {
                    collect_object(nested, objects, prefixes, report, tag_names);
                }
            }
        }
    }
}
