//! DEXPI exchange, exposed to Python.
//!
//! Thin by design: the projection and the XML live in `azoth-dexpi`, and Python supplies a
//! flowsheet as its TOML text and receives back the exchange text plus the report as JSON. The
//! split is the CLI's, so a notebook and a shell read the same two documents.

use std::collections::BTreeMap;
use std::path::Path;

use azoth_dexpi::plant::export::PlantExportOptions;
use azoth_dexpi::plant::import::PlantImportOptions;
use azoth_dexpi::process::export::ExportOptions;
use azoth_dexpi::process::import::ImportOptions;
use azoth_dexpi::{
    export_plant, export_process, import_plant, import_process, import_proteus, writer,
};
use azoth_process::{load_palette, parse_flowsheet};
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Write a flowsheet as a DEXPI 2.0 Process exchange.
///
/// Returns `(xml, report_json)`. A refused exchange raises `ValueError` whose message is the
/// report JSON, so the diagnostics survive the boundary.
#[pyfunction]
#[pyo3(signature = (flowsheet, palette_dir="specs/unit_ops", allow_lossy_mappings=false, project_recycles=false, with_azoth_extensions=false, topology_only=false, allow_unconverged=false, drawing=false, maps=None))]
#[allow(clippy::too_many_arguments)]
pub fn dexpi_export(
    flowsheet: &str,
    palette_dir: &str,
    allow_lossy_mappings: bool,
    project_recycles: bool,
    with_azoth_extensions: bool,
    topology_only: bool,
    allow_unconverged: bool,
    drawing: bool,
    maps: Option<BTreeMap<String, String>>,
) -> PyResult<(String, String)> {
    let palette = load_palette(Path::new(palette_dir))
        .map_err(|error| PyValueError::new_err(format!("palette: {error}")))?;
    let flowsheet = parse_flowsheet(flowsheet)
        .map_err(|error| PyValueError::new_err(format!("flowsheet: {error}")))?;
    let options = ExportOptions {
        allow_lossy: allow_lossy_mappings,
        overrides: maps.unwrap_or_default(),
        project_recycles,
        azoth_extensions: with_azoth_extensions,
        topology_only,
        allow_unconverged,
        drawing,
    };
    match export_process(&flowsheet, &palette, &options) {
        Ok(report) => {
            let json = serde_json::json!({
                "steps": report.steps,
                "streams": report.streams,
                "diagram_groups": report.diagram_groups,
                "diagram_primitives": report.diagram_primitives,
                "converged": report.converged,
                "semantic_valid": report.semantic_valid,
                "approval_status": report.approval_status,
                "cae_round_trip": report.cae_round_trip,
                "diagnostics": report.diagnostics,
            });
            Ok((report.xml, json.to_string()))
        }
        Err(diagnostics) => Err(PyValueError::new_err(
            serde_json::json!({ "diagnostics": diagnostics }).to_string(),
        )),
    }
}

/// Read a DEXPI 2.0 Process exchange as a flowsheet scaffold.
///
/// Returns `(flowsheet_toml, report_json)`. A refused import raises `ValueError` whose message is
/// the report JSON.
#[pyfunction]
#[pyo3(signature = (xml, palette_dir="specs/unit_ops", strict=false, layout=true))]
pub fn dexpi_import(
    xml: &str,
    palette_dir: &str,
    strict: bool,
    layout: bool,
) -> PyResult<(String, String)> {
    let palette = load_palette(Path::new(palette_dir))
        .map_err(|error| PyValueError::new_err(format!("palette: {error}")))?;
    match import_process(xml, &palette, &ImportOptions { strict }) {
        Ok(mut report) => {
            if !layout {
                report.flowsheet.layout = None;
            }
            let toml = report
                .flowsheet
                .to_toml()
                .map_err(|error| PyValueError::new_err(error.to_string()))?;
            let json = serde_json::json!({
                "flowsheet": report.flowsheet,
                "layout_positions": report.layout_positions,
                "diagnostics": report.diagnostics,
            });
            Ok((toml, json.to_string()))
        }
        Err(diagnostics) => Err(PyValueError::new_err(
            serde_json::json!({ "diagnostics": diagnostics }).to_string(),
        )),
    }
}

/// Write a flowsheet as a DEXPI 2.0 Plant exchange - a P&ID.
///
/// Returns `(xml, report_json)`. A refused exchange raises `ValueError` whose message is the
/// report JSON.
#[pyfunction]
#[pyo3(signature = (flowsheet, palette_dir="specs/unit_ops", allow_lossy_mappings=false, with_azoth_extensions=false, drawing=false, maps=None))]
pub fn dexpi_export_plant(
    flowsheet: &str,
    palette_dir: &str,
    allow_lossy_mappings: bool,
    with_azoth_extensions: bool,
    drawing: bool,
    maps: Option<BTreeMap<String, String>>,
) -> PyResult<(String, String)> {
    let palette = load_palette(Path::new(palette_dir))
        .map_err(|error| PyValueError::new_err(format!("palette: {error}")))?;
    let flowsheet = parse_flowsheet(flowsheet)
        .map_err(|error| PyValueError::new_err(format!("flowsheet: {error}")))?;
    let options = PlantExportOptions {
        allow_lossy: allow_lossy_mappings,
        overrides: maps.unwrap_or_default(),
        azoth_extensions: with_azoth_extensions,
        drawing,
    };
    match export_plant(&flowsheet, &palette, &options) {
        Ok(report) => {
            let json = serde_json::json!({
                "items": report.items,
                "segments": report.segments,
                "diagram_groups": report.diagram_groups,
                "diagram_primitives": report.diagram_primitives,
                "semantic_valid": report.semantic_valid,
                "approval_status": report.approval_status,
                "cae_round_trip": report.cae_round_trip,
                "diagnostics": report.diagnostics,
            });
            Ok((report.xml, json.to_string()))
        }
        Err(diagnostics) => Err(PyValueError::new_err(
            serde_json::json!({ "diagnostics": diagnostics }).to_string(),
        )),
    }
}

/// Read a DEXPI 2.0 Plant exchange as a flowsheet scaffold.
///
/// Returns `(flowsheet_toml, report_json)`; a refusal raises `ValueError`.
#[pyfunction]
#[pyo3(signature = (xml, palette_dir="specs/unit_ops", strict=false, layout=true))]
pub fn dexpi_import_plant(
    xml: &str,
    palette_dir: &str,
    strict: bool,
    layout: bool,
) -> PyResult<(String, String)> {
    let palette = load_palette(Path::new(palette_dir))
        .map_err(|error| PyValueError::new_err(format!("palette: {error}")))?;
    match import_plant(xml, &palette, &PlantImportOptions { strict }) {
        Ok(report) => plant_report(report, layout),
        Err(diagnostics) => Err(PyValueError::new_err(
            serde_json::json!({ "diagnostics": diagnostics }).to_string(),
        )),
    }
}

/// Read a Proteus `PlantModel` as a flowsheet scaffold.
///
/// Returns `(flowsheet_toml, report_json)`; a refusal raises `ValueError`.
#[pyfunction]
#[pyo3(signature = (xml, palette_dir="specs/unit_ops", strict=false, layout=true))]
pub fn dexpi_import_proteus(
    xml: &str,
    palette_dir: &str,
    strict: bool,
    layout: bool,
) -> PyResult<(String, String)> {
    let palette = load_palette(Path::new(palette_dir))
        .map_err(|error| PyValueError::new_err(format!("palette: {error}")))?;
    match import_proteus(xml, &palette, &PlantImportOptions { strict }) {
        Ok(report) => plant_report(report, layout),
        Err(diagnostics) => Err(PyValueError::new_err(
            serde_json::json!({ "diagnostics": diagnostics }).to_string(),
        )),
    }
}

fn plant_report(
    mut report: azoth_dexpi::PlantImportReport,
    layout: bool,
) -> PyResult<(String, String)> {
    if !layout {
        report.flowsheet.layout = None;
    }
    let toml = report
        .flowsheet
        .to_toml()
        .map_err(|error| PyValueError::new_err(error.to_string()))?;
    let json = serde_json::json!({
        "flowsheet": report.flowsheet,
        "layout_positions": report.layout_positions,
        "diagnostics": report.diagnostics,
    });
    Ok((toml, json.to_string()))
}

/// Check a DEXPI 2.0 document against the semantic rules.
///
/// Returns the report as JSON; a document with an error is still reported rather than raised,
/// because the report is what a caller reads.
#[pyfunction]
pub fn dexpi_validate(xml: &str) -> PyResult<String> {
    match azoth_dexpi::reader::read(xml) {
        Ok(read) => {
            let report = azoth_dexpi::validate::validate(&read.document);
            Ok(serde_json::json!({
                "valid": report.is_valid(),
                "diagnostics": report.diagnostics,
            })
            .to_string())
        }
        Err(diagnostic) => Ok(serde_json::json!({
            "valid": false,
            "diagnostics": [diagnostic],
        })
        .to_string()),
    }
}

/// The vendored schema's digest, so a caller can pin what it validated against.
#[pyfunction]
pub fn dexpi_schema_sha256() -> &'static str {
    azoth_dexpi::xsd::SHA256
}

/// `writer::format_double` is re-exported for the round-trip tests a caller writes.
#[pyfunction]
pub fn dexpi_format_double(value: f64) -> String {
    writer::format_double(value)
}
