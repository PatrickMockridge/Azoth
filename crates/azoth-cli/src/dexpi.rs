//! `azoth dexpi` - the file exchange with plant-engineering tools.
//!
//! Three verbs. `export` writes a flowsheet and its run as a DEXPI 2.0 Process model; `import`
//! reads one back as a scaffold; `validate` holds a document to the semantic rules and, with
//! `--xsd`, to the official schema through `xmllint`.
//!
//! **The report is the product as much as the file is.** A lossy mapping, a refusal, a port
//! assignment the importer had to guess, and a feed whose composition DEXPI cannot carry are all
//! printed - a file that validates while its author does not know what was lost is the failure
//! this surface exists to prevent.

use std::path::Path;

use azoth_dexpi::plant::export::PlantExportOptions;
use azoth_dexpi::plant::import::PlantImportOptions;
use azoth_dexpi::process::export::ExportOptions;
use azoth_dexpi::process::import::ImportOptions;
use azoth_dexpi::{
    DexpiDiagnostic, ExportReport, ImportReport, PlantExportReport, PlantImportReport, Profile,
    Severity, export_plant, export_process, import_plant, import_process, xsd,
};
use azoth_process::{load_palette, parse_flowsheet};

use crate::cli::{DexpiExportArgs, DexpiImportArgs, DexpiValidateArgs};

/// What a verb produced: whether it succeeded, and the text to print.
pub struct Outcome {
    pub ok: bool,
    pub text: String,
}

/// Write a DEXPI exchange - a PFD (`--profile process`) or a P&ID (`--profile plant`).
///
/// # Errors
/// A message naming what could not be read or written. A document the exchange refuses is an
/// [`Outcome`] with `ok == false`, because the refusal is the report.
pub fn export(args: &DexpiExportArgs) -> Result<Outcome, String> {
    let palette = load_palette(&args.palette)?;
    let text = std::fs::read_to_string(&args.flowsheet)
        .map_err(|error| format!("{}: {error}", args.flowsheet.display()))?;
    let flowsheet =
        parse_flowsheet(&text).map_err(|error| format!("{}: {error}", args.flowsheet.display()))?;

    let mut overrides = std::collections::BTreeMap::new();
    for entry in &args.maps {
        let Some((unit, dexpi_type)) = entry.split_once('=') else {
            return Err(format!(
                "--map `{entry}` is not `unit_ops.<name>=<prefix>/<class>`"
            ));
        };
        if !unit.starts_with("unit_ops.") {
            return Err(format!("--map `{entry}` must name a palette id"));
        }
        overrides.insert(unit.to_string(), dexpi_type.to_string());
    }

    match args.profile.as_str() {
        "process" => {
            let options = ExportOptions {
                allow_lossy: args.allow_lossy_mappings,
                overrides,
                project_recycles: args.project_recycles,
                azoth_extensions: args.with_azoth_extensions,
                topology_only: args.topology_only,
                allow_unconverged: args.allow_unconverged,
                drawing: args.with_drawing,
            };
            match export_process(&flowsheet, &palette, &options) {
                Ok(report) => {
                    write_file(&args.out, &report.xml)?;
                    let text = if args.json {
                        render_export_json(&report, &args.out)
                    } else {
                        render_export(&report, &args.out)
                    };
                    Ok(Outcome { ok: true, text })
                }
                Err(diagnostics) => Ok(failure(&diagnostics, args.json)),
            }
        }
        "plant" => {
            if args.topology_only || args.project_recycles || args.allow_unconverged {
                return Err(
                    "a Plant exchange carries no process state: --topology-only, \
                     --project-recycles and --allow-unconverged are Process flags"
                        .to_string(),
                );
            }
            let options = PlantExportOptions {
                allow_lossy: args.allow_lossy_mappings,
                overrides,
                azoth_extensions: args.with_azoth_extensions,
                drawing: args.with_drawing,
            };
            match export_plant(&flowsheet, &palette, &options) {
                Ok(report) => {
                    write_file(&args.out, &report.xml)?;
                    let text = if args.json {
                        render_plant_export_json(&report, &args.out)
                    } else {
                        render_plant_export(&report, &args.out)
                    };
                    Ok(Outcome { ok: true, text })
                }
                Err(diagnostics) => Ok(failure(&diagnostics, args.json)),
            }
        }
        other => Err(format!("--profile `{other}` is not `process` or `plant`")),
    }
}

/// Read a DEXPI exchange as a flowsheet scaffold, detecting its profile from its imports.
///
/// # Errors
/// A message naming what could not be read or written. A refused document is an [`Outcome`] with
/// `ok == false`.
pub fn import(args: &DexpiImportArgs) -> Result<Outcome, String> {
    let palette = load_palette(&args.palette)?;
    let xml = std::fs::read_to_string(&args.input)
        .map_err(|error| format!("{}: {error}", args.input.display()))?;
    // **A Proteus file is a different serialisation, not a broken DEXPI 2.0 one.** It is
    // projected into the Plant metamodel and read through the same importer.
    let document = match azoth_dexpi::reader::read(&xml) {
        Ok(read) => Some(read.document),
        Err(diagnostic)
            if diagnostic.code == "unsupported_root" && azoth_dexpi::proteus::is_proteus(&xml) =>
        {
            None
        }
        Err(diagnostic) => return Ok(failure(&[diagnostic], args.json)),
    };
    let Some(document) = document else {
        return match args.profile.as_str() {
            "process" => {
                Err("a Proteus PlantModel is a P&ID; read it with `--profile plant`".to_string())
            }
            "auto" | "plant" => match azoth_dexpi::import_proteus(
                &xml,
                &palette,
                &PlantImportOptions {
                    strict: args.strict,
                },
            ) {
                Ok(report) => write_plant_import(report, args),
                Err(diagnostics) => Ok(failure(&diagnostics, args.json)),
            },
            other => Err(format!(
                "--profile `{other}` is not `auto`, `process` or `plant`"
            )),
        };
    };
    let detected = Profile::of(&document);
    let selected = match args.profile.as_str() {
        "auto" => detected.ok_or_else(|| {
            "the document imports neither the Plant nor the Process model, so its profile \
             cannot be taken from its imports"
                .to_string()
        })?,
        "process" => Profile::Process,
        "plant" => Profile::Plant,
        other => {
            return Err(format!(
                "--profile `{other}` is not `auto`, `process` or `plant`"
            ));
        }
    };
    if args.profile != "auto" && detected != Some(selected) {
        return Err(format!(
            "the document imports the {detected:?} model, not the {selected:?} one"
        ));
    }

    match selected {
        Profile::Process => match import_process(
            &xml,
            &palette,
            &ImportOptions {
                strict: args.strict,
            },
        ) {
            Ok(report) => write_import(report, args),
            Err(diagnostics) => Ok(failure(&diagnostics, args.json)),
        },
        Profile::Plant => match import_plant(
            &xml,
            &palette,
            &PlantImportOptions {
                strict: args.strict,
            },
        ) {
            Ok(report) => write_plant_import(report, args),
            Err(diagnostics) => Ok(failure(&diagnostics, args.json)),
        },
    }
}

fn write_file(path: &Path, contents: &str) -> Result<(), String> {
    std::fs::write(path, contents).map_err(|error| format!("{}: {error}", path.display()))
}

fn failure(diagnostics: &[DexpiDiagnostic], json: bool) -> Outcome {
    let text = if json {
        render_failure_json(diagnostics)
    } else {
        render_failure(diagnostics)
    };
    Outcome { ok: false, text }
}

fn write_import(mut report: ImportReport, args: &DexpiImportArgs) -> Result<Outcome, String> {
    if args.no_layout {
        report.flowsheet.layout = None;
    }
    let toml = report
        .flowsheet
        .to_toml()
        .map_err(|error| format!("the flowsheet could not be written: {error}"))?;
    write_file(&args.out, &toml)?;
    let text = if args.json {
        render_import_json(&report, &args.out)
    } else {
        render_import(&report, &args.out)
    };
    Ok(Outcome { ok: true, text })
}

fn write_plant_import(
    mut report: PlantImportReport,
    args: &DexpiImportArgs,
) -> Result<Outcome, String> {
    if args.no_layout {
        report.flowsheet.layout = None;
    }
    let toml = report
        .flowsheet
        .to_toml()
        .map_err(|error| format!("the flowsheet could not be written: {error}"))?;
    write_file(&args.out, &toml)?;
    let text = if args.json {
        render_plant_import_json(&report, &args.out)
    } else {
        render_plant_import(&report, &args.out)
    };
    Ok(Outcome { ok: true, text })
}

/// Check a document against the semantic rules, and optionally the schema.
///
/// # Errors
/// A message naming what could not be read, or that `xmllint` is absent when `--xsd` asked for
/// the schema check. A document that fails is an [`Outcome`] with `ok == false`.
pub fn validate(args: &DexpiValidateArgs) -> Result<Outcome, String> {
    let xml = std::fs::read_to_string(&args.input)
        .map_err(|error| format!("{}: {error}", args.input.display()))?;
    let read = match azoth_dexpi::reader::read(&xml) {
        Ok(read) => read,
        Err(diagnostic) => {
            let text = if args.json {
                render_failure_json(&[diagnostic])
            } else {
                render_failure(&[diagnostic])
            };
            return Ok(Outcome { ok: false, text });
        }
    };
    let report = azoth_dexpi::validate::validate(&read.document);
    let mut lines = Vec::new();
    for diagnostic in &report.diagnostics {
        lines.push(render_diagnostic(diagnostic));
    }
    let mut ok = report.is_valid();
    let mut schema = None;
    if args.xsd {
        let path = args.input.clone();
        schema = Some(xmllint(&path)?);
        ok &= schema.as_deref() == Some("valid");
    }
    if args.json {
        let value = serde_json::json!({
            "input": args.input.display().to_string(),
            "semantic": if report.is_valid() { "valid" } else { "invalid" },
            "xsd": schema,
            "diagnostics": report.diagnostics,
        });
        return Ok(Outcome {
            ok,
            text: serde_json::to_string_pretty(&value).unwrap_or_default(),
        });
    }
    let mut text = format!(
        "{}: semantic rules {}\n",
        args.input.display(),
        if report.is_valid() { "pass" } else { "fail" }
    );
    if let Some(schema) = &schema {
        text.push_str(&format!("  official XSD: {schema}\n"));
    }
    text.push_str(&lines.join("\n"));
    if !lines.is_empty() {
        text.push('\n');
    }
    Ok(Outcome { ok, text })
}

/// Run `xmllint` against the vendored schema, through a temporary file.
fn xmllint(input: &Path) -> Result<String, String> {
    use std::io::Write;
    let mut schema_path = std::env::temp_dir();
    schema_path.push(format!("azoth-dexpi-schema-{}.xsd", std::process::id()));
    {
        let mut file = std::fs::File::create(&schema_path).map_err(|error| {
            format!(
                "could not write the schema to {}: {error}",
                schema_path.display()
            )
        })?;
        file.write_all(xsd::SCHEMA.as_bytes())
            .map_err(|error| format!("could not write the schema: {error}"))?;
    }
    let output = std::process::Command::new("xmllint")
        .args(["--nonet", "--noout", "--schema"])
        .arg(&schema_path)
        .arg(input)
        .output()
        .map_err(|error| {
            format!(
                "xmllint is required for --xsd and could not be run ({error}); install \
                 libxml2-utils. The schema check is not reported as passed when it did not run"
            )
        })?;
    let _ = std::fs::remove_file(&schema_path);
    if output.status.success() {
        Ok("valid".to_string())
    } else {
        Ok(format!(
            "invalid:\n{}",
            String::from_utf8_lossy(&output.stderr).trim()
        ))
    }
}

fn render_export(report: &ExportReport, out: &Path) -> String {
    let mut text = format!(
        "{}: {} steps, {} streams, run {}\n",
        out.display(),
        report.steps,
        report.streams,
        if report.converged {
            "converged"
        } else {
            "not run"
        }
    );
    text.push_str("  semantic rules: pass\n");
    if report.diagram_groups > 0 {
        text.push_str(&format!(
            "  drawing: {} groups, {} primitives\n",
            report.diagram_groups, report.diagram_primitives
        ));
    }
    text.push_str(&format!("  approval: {}\n", report.approval_status));
    text.push_str(&format!(
        "  named-CAE round trip: {}\n",
        report.cae_round_trip
    ));
    for diagnostic in &report.diagnostics {
        text.push_str(&render_diagnostic(diagnostic));
        text.push('\n');
    }
    text
}

fn render_export_json(report: &ExportReport, out: &Path) -> String {
    let value = serde_json::json!({
        "out": out.display().to_string(),
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
    serde_json::to_string_pretty(&value).unwrap_or_default()
}

fn render_import(report: &ImportReport, out: &Path) -> String {
    let flowsheet = &report.flowsheet;
    let mut text = format!(
        "{}: {} instances, {} connections, {} recycle(s)\n",
        out.display(),
        flowsheet.instances.len(),
        flowsheet.connections.len(),
        flowsheet.recycles.len()
    );
    text.push_str(
        "  this is a scaffold: feeds carry no composition and no step carries parameters\n",
    );
    if report.layout_positions > 0 {
        text.push_str(&format!(
            "  layout: {} position(s) read from the drawing\n",
            report.layout_positions
        ));
    }
    for diagnostic in &report.diagnostics {
        text.push_str(&render_diagnostic(diagnostic));
        text.push('\n');
    }
    text
}

fn render_plant_export(report: &PlantExportReport, out: &Path) -> String {
    let mut text = format!(
        "{}: {} tagged items, {} piping segments\n",
        out.display(),
        report.items,
        report.segments
    );
    text.push_str("  semantic rules: pass\n");
    if report.diagram_groups > 0 {
        text.push_str(&format!(
            "  drawing: {} groups, {} primitives\n",
            report.diagram_groups, report.diagram_primitives
        ));
    }
    text.push_str(&format!("  approval: {}\n", report.approval_status));
    text.push_str(&format!(
        "  named-CAE round trip: {}\n",
        report.cae_round_trip
    ));
    for diagnostic in &report.diagnostics {
        text.push_str(&render_diagnostic(diagnostic));
        text.push('\n');
    }
    text
}

fn render_plant_export_json(report: &PlantExportReport, out: &Path) -> String {
    let value = serde_json::json!({
        "out": out.display().to_string(),
        "items": report.items,
        "segments": report.segments,
        "diagram_groups": report.diagram_groups,
        "diagram_primitives": report.diagram_primitives,
        "semantic_valid": report.semantic_valid,
        "approval_status": report.approval_status,
        "cae_round_trip": report.cae_round_trip,
        "diagnostics": report.diagnostics,
    });
    serde_json::to_string_pretty(&value).unwrap_or_default()
}

fn render_plant_import(report: &PlantImportReport, out: &Path) -> String {
    let flowsheet = &report.flowsheet;
    let mut text = format!(
        "{}: {} instances, {} connections, {} recycle(s)\n",
        out.display(),
        flowsheet.instances.len(),
        flowsheet.connections.len(),
        flowsheet.recycles.len()
    );
    text.push_str(
        "  this is a scaffold: feeds carry no composition and no item carries parameters\n",
    );
    if report.layout_positions > 0 {
        text.push_str(&format!(
            "  layout: {} position(s) read from the drawing\n",
            report.layout_positions
        ));
    }
    for diagnostic in &report.diagnostics {
        text.push_str(&render_diagnostic(diagnostic));
        text.push('\n');
    }
    text
}

fn render_plant_import_json(report: &PlantImportReport, out: &Path) -> String {
    let value = serde_json::json!({
        "out": out.display().to_string(),
        "flowsheet": report.flowsheet,
        "layout_positions": report.layout_positions,
        "diagnostics": report.diagnostics,
    });
    serde_json::to_string_pretty(&value).unwrap_or_default()
}

fn render_import_json(report: &ImportReport, out: &Path) -> String {
    let value = serde_json::json!({
        "out": out.display().to_string(),
        "flowsheet": report.flowsheet,
        "layout_positions": report.layout_positions,
        "diagnostics": report.diagnostics,
    });
    serde_json::to_string_pretty(&value).unwrap_or_default()
}

fn render_failure(diagnostics: &[DexpiDiagnostic]) -> String {
    let mut text: Vec<String> = diagnostics.iter().map(render_diagnostic).collect();
    text.push(format!(
        "{} finding(s) - the exchange was refused",
        diagnostics.len()
    ));
    text.join("\n")
}

fn render_failure_json(diagnostics: &[DexpiDiagnostic]) -> String {
    serde_json::to_string_pretty(&serde_json::json!({ "diagnostics": diagnostics }))
        .unwrap_or_default()
}

fn render_diagnostic(diagnostic: &DexpiDiagnostic) -> String {
    let level = match diagnostic.severity {
        Severity::Error => "error",
        Severity::Warning => "warning",
        Severity::Info => "note",
    };
    if diagnostic.subject.is_empty() {
        format!("  {level}: {}: {}", diagnostic.code, diagnostic.message)
    } else {
        format!(
            "  {level}: {}: {}: {}",
            diagnostic.code, diagnostic.subject, diagnostic.message
        )
    }
}
