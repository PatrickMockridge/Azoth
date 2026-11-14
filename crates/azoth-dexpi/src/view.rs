//! A drawing and a coverage report, for a viewer.
//!
//! **The question this answers is "what would azoth model?"** A DEXPI file carries more than
//! a simulation can: tagged items with no unit, instruments, junctions, labels and the
//! drawing's furniture. [`inspect`] reads the file and its graphics, runs the same import the
//! CLI does, and marks every conceptual object with what it became - modelled, lossy,
//! ambiguous, unmapped, topology, control - and every drawn group with nothing behind it as
//! drawing-only. The statuses come from the same mapping tables the importers read; this
//! module never holds a second mapping.

use std::collections::{BTreeMap, BTreeSet};

use azoth_process::UnitOpSpec;
use azoth_process::flowsheet::Flowsheet;
use serde::Serialize;

use crate::diagnostic::DexpiDiagnostic;
use crate::drawing::{Drawing, DrawingGroup};
use crate::extensions;
use crate::mapping;
use crate::model::{Child, ComponentChild, Document, Object};
use crate::plant;
use crate::profile::Profile;

/// The serialisation a view came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum ExchangeFormat {
    Dexpi20,
    Proteus,
}

/// What an object in the drawing would become in azoth.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum CoverageStatus {
    Modelled,
    Lossy,
    Ambiguous,
    Unmapped,
    Topology,
    Control,
    DrawingOnly,
}

/// One conceptual object, or one drawn group with nothing behind it.
#[derive(Debug, Clone, Serialize)]
pub struct ObjectCoverage {
    pub id: String,
    pub ty: String,
    pub label: String,
    pub unit: Option<String>,
    pub azoth: Option<String>,
    pub status: CoverageStatus,
    pub reason: String,
    pub alternatives: Vec<String>,
    pub imported: bool,
    pub exportable: bool,
    pub drawn: bool,
    pub primitives: usize,
}

/// How many objects fall in each category.
#[derive(Debug, Clone, Default, Serialize)]
pub struct CoverageCounts {
    pub semantic: usize,
    pub modelled: usize,
    pub lossy: usize,
    pub ambiguous: usize,
    pub unmapped: usize,
    pub topology: usize,
    pub control: usize,
    pub drawing_only: usize,
    pub drawn: usize,
    pub not_drawn: usize,
    pub groups: usize,
    pub primitives: usize,
    pub layout_positions: usize,
}

/// The scaffold an import produced, and what it still owes.
#[derive(Debug, Clone, Serialize)]
pub struct Scaffold {
    pub toml: String,
    pub flowsheet_id: String,
    pub instances: usize,
    pub connections: usize,
    pub recycles: usize,
    pub instruments: usize,
    pub signals: usize,
    pub layout_positions: usize,
    pub check_errors: usize,
    pub check_warnings: usize,
}

/// A drawing, its coverage, and the scaffold the file would become.
#[derive(Debug, Clone, Serialize)]
pub struct ViewReport {
    pub format: ExchangeFormat,
    pub profile: Option<Profile>,
    pub model_name: String,
    pub model_uri: String,
    pub drawing: Option<Drawing>,
    pub objects: Vec<ObjectCoverage>,
    pub counts: CoverageCounts,
    pub diagnostics: Vec<DexpiDiagnostic>,
    pub scaffold: Option<Scaffold>,
}

#[derive(Default)]
struct Imported {
    flowsheet: Option<Flowsheet>,
    bindings: BTreeMap<String, String>,
    drawing: Option<Drawing>,
    layout_positions: usize,
    diagnostics: Vec<DexpiDiagnostic>,
}

struct Coverage {
    status: CoverageStatus,
    unit: Option<String>,
    alternatives: Vec<String>,
    reason: String,
    exportable: bool,
}

/// Read a DEXPI 2.0 or Proteus drawing as a view.
///
/// # Errors
/// Only a file that cannot be parsed at all. A document the semantic rules refuse still
/// returns a view: the drawing is what a person opens it to see, and the diagnostics say why
/// the scaffold is missing.
pub fn inspect(xml: &str, palette: &[UnitOpSpec]) -> Result<ViewReport, Vec<DexpiDiagnostic>> {
    if crate::proteus::is_proteus(xml) {
        inspect_proteus(xml, palette)
    } else {
        inspect_native(xml, palette)
    }
}

fn inspect_native(xml: &str, palette: &[UnitOpSpec]) -> Result<ViewReport, Vec<DexpiDiagnostic>> {
    let read = crate::reader::read(xml).map_err(|diagnostic| vec![diagnostic])?;
    let document = read.document;
    let profile = Profile::of(&document);
    let mut imported = Imported::default();
    match profile {
        Some(Profile::Process) => {
            match crate::process::import::import_process(
                xml,
                palette,
                &crate::process::import::ImportOptions::default(),
            ) {
                Ok(report) => {
                    imported.flowsheet = Some(report.flowsheet);
                    imported.bindings = report.bindings;
                    imported.drawing = report.drawing;
                    imported.layout_positions = report.layout_positions;
                    imported.diagnostics = report.diagnostics;
                }
                Err(diagnostics) => imported.diagnostics = diagnostics,
            }
        }
        Some(Profile::Plant) => {
            match crate::plant::import::import_plant_document(
                document.clone(),
                palette,
                &crate::plant::import::PlantImportOptions::default(),
            ) {
                Ok(report) => {
                    imported.flowsheet = Some(report.flowsheet);
                    imported.bindings = report.bindings;
                    imported.drawing = report.drawing;
                    imported.layout_positions = report.layout_positions;
                    imported.diagnostics = report.diagnostics;
                }
                Err(diagnostics) => imported.diagnostics = diagnostics,
            }
        }
        None => imported.diagnostics.push(DexpiDiagnostic::warning(
            "unknown_profile",
            document.model.name.clone(),
            "the document imports neither the Process nor the Plant model",
        )),
    }
    if imported.drawing.is_none() {
        if let Some(report) = crate::drawing::native::read_document(&document) {
            imported.diagnostics.extend(report.diagnostics);
            imported.drawing = Some(report.drawing);
        }
    }
    let drawing = imported.drawing.clone().or_else(|| {
        imported
            .flowsheet
            .as_ref()
            .map(crate::drawing::export::schematic)
    });
    Ok(build_view(
        ExchangeFormat::Dexpi20,
        profile,
        document.model.name.clone(),
        document.model.uri.clone(),
        &document,
        drawing,
        imported,
        palette,
    ))
}

fn inspect_proteus(xml: &str, palette: &[UnitOpSpec]) -> Result<ViewReport, Vec<DexpiDiagnostic>> {
    let mut diagnostics = Vec::new();
    let projection = match crate::proteus::project(xml, &mut diagnostics) {
        Ok(projection) => projection,
        Err(diagnostic) => {
            diagnostics.push(diagnostic);
            return Err(diagnostics);
        }
    };
    let document = projection.document.clone();
    let drawing = projection.drawing;
    let mut imported = Imported {
        diagnostics,
        ..Imported::default()
    };
    match crate::plant::import::import_plant_document(
        projection.document,
        palette,
        &crate::plant::import::PlantImportOptions::default(),
    ) {
        Ok(report) => {
            imported.flowsheet = Some(report.flowsheet);
            imported.bindings = report.bindings;
            imported.layout_positions = report.layout_positions;
            imported.diagnostics.extend(report.diagnostics);
        }
        Err(report) => {
            imported.diagnostics.extend(report);
        }
    }
    if let Some(drawing) = drawing.clone() {
        if let Some(flowsheet) = &mut imported.flowsheet {
            imported.layout_positions =
                crate::drawing::layout::apply(flowsheet, &drawing, &imported.bindings);
        }
    }
    let name = drawing.as_ref().map_or_else(
        || "Proteus drawing".to_string(),
        |drawing| drawing.name.clone(),
    );
    let view = build_view(
        ExchangeFormat::Proteus,
        Some(Profile::Plant),
        name,
        "urn:azoth:dexpi:1:proteus".to_string(),
        &document,
        drawing,
        imported,
        palette,
    );
    Ok(view)
}

#[allow(clippy::too_many_arguments)]
fn build_view(
    format: ExchangeFormat,
    profile: Option<Profile>,
    model_name: String,
    model_uri: String,
    document: &Document,
    drawing: Option<Drawing>,
    imported: Imported,
    palette: &[UnitOpSpec],
) -> ViewReport {
    let imported_ok = imported.flowsheet.is_some();
    let mut objects = Vec::new();
    let mut seen = BTreeSet::new();
    for child in &document.model.children {
        if let Child::Object(object) = child {
            walk(
                object,
                &mut objects,
                &mut seen,
                &imported.bindings,
                imported_ok,
            );
        }
    }

    let groups = drawing.as_ref().map_or(&[][..], |drawing| &drawing.groups);
    for object in &mut objects {
        let (drawn, primitives) = drawn(groups, &object.id, object.azoth.as_deref());
        object.drawn = drawn;
        object.primitives = primitives;
    }
    let represented: BTreeSet<String> = objects
        .iter()
        .flat_map(|object| std::iter::once(object.id.clone()).chain(object.azoth.clone()))
        .collect();
    for group in groups {
        let Some(represents) = group.represents.as_deref() else {
            objects.push(drawing_only(group));
            continue;
        };
        if !represented.contains(represents) {
            objects.push(drawing_only(group));
        }
    }

    let mut counts = CoverageCounts {
        groups: groups.len(),
        primitives: drawing.as_ref().map_or(0, Drawing::primitive_count),
        layout_positions: imported.layout_positions,
        ..CoverageCounts::default()
    };
    for object in &objects {
        match object.status {
            CoverageStatus::Modelled => counts.modelled += 1,
            CoverageStatus::Lossy => counts.lossy += 1,
            CoverageStatus::Ambiguous => counts.ambiguous += 1,
            CoverageStatus::Unmapped => counts.unmapped += 1,
            CoverageStatus::Topology => counts.topology += 1,
            CoverageStatus::Control => counts.control += 1,
            CoverageStatus::DrawingOnly => counts.drawing_only += 1,
        }
        if object.status != CoverageStatus::DrawingOnly {
            counts.semantic += 1;
            if object.drawn {
                counts.drawn += 1;
            }
        }
    }
    counts.not_drawn = counts.semantic.saturating_sub(counts.drawn);

    let has_flowsheet = imported.flowsheet.is_some();
    let scaffold = imported
        .flowsheet
        .as_ref()
        .and_then(|flowsheet| scaffold(flowsheet, palette, imported.layout_positions));
    let mut diagnostics = imported.diagnostics;
    if scaffold.is_none() && has_flowsheet {
        diagnostics.push(DexpiDiagnostic::warning(
            "scaffold_unreadable",
            "",
            "the imported flowsheet could not be written as TOML",
        ));
    }
    ViewReport {
        format,
        profile,
        model_name,
        model_uri,
        drawing,
        objects,
        counts,
        diagnostics,
        scaffold,
    }
}

fn scaffold(
    flowsheet: &Flowsheet,
    palette: &[UnitOpSpec],
    layout_positions: usize,
) -> Option<Scaffold> {
    let toml = flowsheet.to_toml().ok()?;
    let check = azoth_process::validate(flowsheet, palette);
    let check_errors = check
        .iter()
        .filter(|diagnostic| diagnostic.severity() == azoth_process::Severity::Error)
        .count();
    Some(Scaffold {
        toml,
        flowsheet_id: flowsheet.id.clone(),
        instances: flowsheet.instances.len(),
        connections: flowsheet.connections.len(),
        recycles: flowsheet.recycles.len(),
        instruments: flowsheet.instruments.len(),
        signals: flowsheet.signals.len(),
        layout_positions,
        check_errors,
        check_warnings: check.len().saturating_sub(check_errors),
    })
}

fn walk(
    object: &Object,
    objects: &mut Vec<ObjectCoverage>,
    seen: &mut BTreeSet<String>,
    bindings: &BTreeMap<String, String>,
    imported_ok: bool,
) {
    if crate::profile::is_graphics(&object.ty) {
        return;
    }
    if let Some(coverage) = coverage_of(object) {
        let id = object.id.clone().unwrap_or_else(|| object.ty.clone());
        if seen.insert(id.clone()) {
            let azoth = bindings.get(&id).cloned();
            objects.push(ObjectCoverage {
                id,
                ty: object.ty.clone(),
                label: label(object),
                unit: coverage.unit,
                azoth,
                status: coverage.status,
                reason: coverage.reason,
                alternatives: coverage.alternatives,
                imported: imported_ok,
                exportable: coverage.exportable,
                drawn: false,
                primitives: 0,
            });
        }
    }
    for child in &object.children {
        match child {
            Child::Object(nested) => walk(nested, objects, seen, bindings, imported_ok),
            Child::Components(components) => {
                for nested in &components.children {
                    if let ComponentChild::Object(nested) = nested {
                        walk(nested, objects, seen, bindings, imported_ok);
                    }
                }
            }
            _ => {}
        }
    }
}

fn coverage_of(object: &Object) -> Option<Coverage> {
    let ty = object.ty.as_str();
    if ty == "Process/Process.Stream" {
        return Some(topology("a material connection"));
    }
    if ty == "Process/Process.Source" || ty == "Process/Process.Sink" {
        return Some(topology("a boundary step"));
    }
    if ty.starts_with("Process/Process.") && ty != "Process/Process.MaterialPort" {
        return Some(process_step(object, ty));
    }
    if ty.starts_with("Plant/ProcessEquipment.") && ty != "Plant/ProcessEquipment.Nozzle" {
        return Some(plant_item(ty));
    }
    if ty == "Plant/Instrumentation.ProcessInstrumentationFunction"
        || ty == "Plant/Instrumentation.SignalConveyingFunction"
    {
        return Some(Coverage {
            status: CoverageStatus::Control,
            unit: None,
            alternatives: Vec::new(),
            reason: "carried as control evidence, not simulated".to_string(),
            exportable: true,
        });
    }
    if ty.starts_with("Plant/Piping.") {
        let class = plant::class_of(ty);
        if class == "PipingNetworkSegment" || class.ends_with("PipeOffPageConnector") {
            return Some(topology("piping topology"));
        }
        if plant::for_class(class).is_some() {
            return Some(plant_item(ty));
        }
    }
    None
}

fn process_step(object: &Object, ty: &str) -> Coverage {
    if object.data_string(extensions::RECYCLE).is_some() {
        return topology("a projected recycle tear");
    }
    if let Some(unit) = object.data_string(extensions::UNIT_OP) {
        return from_unit(unit, "the exchange states `Azoth/UnitOp`");
    }
    let candidates = mapping::candidates_for(ty);
    match candidates {
        [] => Coverage {
            status: CoverageStatus::Unmapped,
            unit: None,
            alternatives: Vec::new(),
            reason: format!("no reviewed azoth unit for `{ty}`"),
            exportable: false,
        },
        [unit] => from_unit(unit, "one azoth unit maps to this type"),
        _ => {
            let chosen = crate::process::import::preferred(ty, candidates);
            Coverage {
                status: CoverageStatus::Ambiguous,
                unit: None,
                alternatives: candidates.iter().map(|unit| (*unit).to_string()).collect(),
                reason: format!(
                    "several azoth units share `{ty}`; the importer chooses `{chosen}`"
                ),
                exportable: true,
            }
        }
    }
}

fn from_unit(unit: &str, reason: &str) -> Coverage {
    match mapping::for_azoth(unit) {
        Some(row) => Coverage {
            status: if row.disposition == mapping::Disposition::Lossy {
                CoverageStatus::Lossy
            } else {
                CoverageStatus::Modelled
            },
            unit: Some(unit.to_string()),
            alternatives: Vec::new(),
            reason: if row.disposition == mapping::Disposition::Lossy {
                format!("lossy: {}", row.note)
            } else {
                reason.to_string()
            },
            exportable: row.disposition != mapping::Disposition::Refused,
        },
        None => Coverage {
            status: CoverageStatus::Unmapped,
            unit: None,
            alternatives: Vec::new(),
            reason: format!("`{unit}` has no row in the Process mapping"),
            exportable: false,
        },
    }
}

fn plant_item(ty: &str) -> Coverage {
    let class = plant::class_of(ty);
    let Some(unit) = plant::for_class(class) else {
        return Coverage {
            status: CoverageStatus::Unmapped,
            unit: None,
            alternatives: Vec::new(),
            reason: format!("no reviewed azoth unit for `{class}`; kept as `dexpi.unmapped`"),
            exportable: false,
        };
    };
    let Some(row) = plant::for_azoth(unit) else {
        return Coverage {
            status: CoverageStatus::Unmapped,
            unit: None,
            alternatives: Vec::new(),
            reason: format!("`{unit}` has no row in the Plant mapping"),
            exportable: false,
        };
    };
    let (status, reason, exportable) = match row.write {
        crate::mapping::Disposition::Lossy => {
            (CoverageStatus::Lossy, format!("lossy: {}", row.note), true)
        }
        crate::mapping::Disposition::Refused => (
            CoverageStatus::Modelled,
            format!("readable, not writable: {}", row.note),
            false,
        ),
        crate::mapping::Disposition::Faithful => {
            (CoverageStatus::Modelled, "mapped".to_string(), true)
        }
    };
    Coverage {
        status,
        unit: Some(unit.to_string()),
        alternatives: Vec::new(),
        reason,
        exportable,
    }
}

fn topology(reason: &str) -> Coverage {
    Coverage {
        status: CoverageStatus::Topology,
        unit: None,
        alternatives: Vec::new(),
        reason: reason.to_string(),
        exportable: true,
    }
}

fn drawing_only(group: &DrawingGroup) -> ObjectCoverage {
    ObjectCoverage {
        id: group.id.clone().unwrap_or_else(|| "drawing".to_string()),
        ty: "drawing".to_string(),
        label: group
            .id
            .clone()
            .unwrap_or_else(|| "unrepresented".to_string()),
        unit: None,
        azoth: None,
        status: CoverageStatus::DrawingOnly,
        reason: "graphics with no conceptual object behind them".to_string(),
        alternatives: Vec::new(),
        imported: false,
        exportable: false,
        drawn: true,
        primitives: group.primitives.len(),
    }
}

fn drawn(groups: &[DrawingGroup], id: &str, azoth: Option<&str>) -> (bool, usize) {
    let mut count = 0;
    for group in groups {
        let represents = group.represents.as_deref();
        if represents == Some(id) || azoth.is_some() && represents == azoth {
            count += group.primitives.len();
        }
    }
    (count > 0, count)
}

fn label(object: &Object) -> String {
    if let Some(tag) = object.data_string("TagName") {
        return tag.to_string();
    }
    if object.ty.contains("ProcessInstrumentationFunction") {
        if let Some(tag) = instrument_tag(object) {
            return tag;
        }
    }
    object
        .data_string("Identifier")
        .or_else(|| object.data_string("Label"))
        .or_else(|| object.data_string(extensions::TAG))
        .or(object.id.as_deref())
        .unwrap_or("(unnamed)")
        .to_string()
}

fn instrument_tag(object: &Object) -> Option<String> {
    let category = object.data_string("ProcessInstrumentationFunctionCategory")?;
    let number = object.data_string("ProcessInstrumentationFunctionNumber")?;
    let modifier = object
        .data_string("ProcessInstrumentationFunctionModifier")
        .unwrap_or("");
    Some(format!("{category}{modifier}-{number}"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::Object;

    fn object(ty: &str) -> Object {
        Object {
            id: Some("O1".to_string()),
            name: None,
            ty: ty.to_string(),
            children: Vec::new(),
        }
    }

    #[test]
    fn a_type_several_units_share_is_ambiguous() {
        let coverage = coverage_of(&object("Process/Process.Distilling")).expect("a step");
        assert_eq!(coverage.status, CoverageStatus::Ambiguous);
        assert!(coverage.unit.is_none());
        assert!(
            !coverage.alternatives.is_empty(),
            "{:?}",
            coverage.alternatives
        );
    }

    #[test]
    fn a_boundary_step_is_topology() {
        let coverage = coverage_of(&object("Process/Process.Source")).expect("a step");
        assert_eq!(coverage.status, CoverageStatus::Topology);
    }

    #[test]
    fn a_lossy_plant_item_says_what_it_loses() {
        let coverage =
            coverage_of(&object("Plant/ProcessEquipment.AirCoolingSystem")).expect("an item");
        assert_eq!(coverage.status, CoverageStatus::Lossy);
        assert_eq!(coverage.unit.as_deref(), Some("unit_ops.cooler"));
        assert!(!coverage.reason.is_empty());
    }

    #[test]
    fn a_plant_class_with_no_unit_is_unmapped() {
        let coverage = coverage_of(&object("Plant/ProcessEquipment.Dryer")).expect("an item");
        assert_eq!(coverage.status, CoverageStatus::Unmapped);
        assert!(coverage.unit.is_none());
    }
}
