//! One facade an editor, a server and a binding can all call.
//!
//! **A refused export is a report, not an exception.** The doors below this one differ in
//! transport, not in answer: a flowsheet in TOML and a profile name out of the UI become a
//! DEXPI file or the diagnostics that refused it. Only a document that does not parse is an
//! error, because there is no report to give about bytes that are not a flowsheet.

use std::collections::BTreeMap;

use azoth_process::{UnitOpSpec, parse_flowsheet};
use serde::{Deserialize, Serialize};

use crate::diagnostic::DexpiDiagnostic;
use crate::plant::export::{PlantExportOptions, export_plant};
use crate::process::export::{ExportOptions, export_process};

/// What a caller may ask an export to do, as JSON.
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct ExchangeOptions {
    #[serde(default)]
    pub allow_lossy: bool,
    #[serde(default)]
    pub project_recycles: bool,
    #[serde(default)]
    pub azoth_extensions: bool,
    #[serde(default)]
    pub topology_only: bool,
    #[serde(default)]
    pub allow_unconverged: bool,
    #[serde(default)]
    pub drawing: bool,
    #[serde(default)]
    pub overrides: BTreeMap<String, String>,
}

/// The file, or the diagnostics that refused it.
#[derive(Debug, Clone, Serialize)]
pub struct ExchangeExport {
    pub ok: bool,
    pub xml: Option<String>,
    pub diagnostics: Vec<DexpiDiagnostic>,
    pub steps: usize,
    pub streams: usize,
    pub items: usize,
    pub segments: usize,
    pub diagram_groups: usize,
    pub diagram_primitives: usize,
}

impl ExchangeExport {
    fn refused(diagnostics: Vec<DexpiDiagnostic>) -> Self {
        Self {
            ok: false,
            xml: None,
            diagnostics,
            steps: 0,
            streams: 0,
            items: 0,
            segments: 0,
            diagram_groups: 0,
            diagram_primitives: 0,
        }
    }
}

/// Export a flowsheet document as the named profile.
///
/// # Errors
/// A TOML document that does not read, or a profile name this port does not carry. A refused
/// export is `Ok` with `ok == false` and the diagnostics that refused it.
pub fn export(
    document: &str,
    palette: &[UnitOpSpec],
    profile: &str,
    options: &ExchangeOptions,
) -> Result<ExchangeExport, Vec<DexpiDiagnostic>> {
    let flowsheet = parse_flowsheet(document).map_err(|error| {
        vec![DexpiDiagnostic::error(
            "flowsheet_unreadable",
            "",
            error.to_string(),
        )]
    })?;
    match profile {
        "process" => {
            let options = ExportOptions {
                allow_lossy: options.allow_lossy,
                overrides: options.overrides.clone(),
                project_recycles: options.project_recycles,
                azoth_extensions: options.azoth_extensions,
                topology_only: options.topology_only,
                allow_unconverged: options.allow_unconverged,
                drawing: options.drawing,
            };
            match export_process(&flowsheet, palette, &options) {
                Ok(report) => Ok(ExchangeExport {
                    ok: true,
                    xml: Some(report.xml),
                    diagnostics: report.diagnostics,
                    steps: report.steps,
                    streams: report.streams,
                    items: 0,
                    segments: 0,
                    diagram_groups: report.diagram_groups,
                    diagram_primitives: report.diagram_primitives,
                }),
                Err(diagnostics) => Ok(ExchangeExport::refused(diagnostics)),
            }
        }
        "plant" => {
            if options.project_recycles || options.topology_only || options.allow_unconverged {
                return Ok(ExchangeExport::refused(vec![DexpiDiagnostic::error(
                    "plant_flag_mismatch",
                    "",
                    "a Plant exchange carries no process state: project_recycles, topology_only \
                     and allow_unconverged are Process options",
                )]));
            }
            let options = PlantExportOptions {
                allow_lossy: options.allow_lossy,
                overrides: options.overrides.clone(),
                azoth_extensions: options.azoth_extensions,
                drawing: options.drawing,
            };
            match export_plant(&flowsheet, palette, &options) {
                Ok(report) => Ok(ExchangeExport {
                    ok: true,
                    xml: Some(report.xml),
                    diagnostics: report.diagnostics,
                    steps: 0,
                    streams: 0,
                    items: report.items,
                    segments: report.segments,
                    diagram_groups: report.diagram_groups,
                    diagram_primitives: report.diagram_primitives,
                }),
                Err(diagnostics) => Ok(ExchangeExport::refused(diagnostics)),
            }
        }
        other => Err(vec![DexpiDiagnostic::error(
            "unknown_profile",
            other,
            "the profile is `process` or `plant`",
        )]),
    }
}
