//! `azoth edit` — one command against a flowsheet, and what it left.
//!
//! **The editing half of the wire, from a shell.** `run` answers what a document computes; this
//! answers what an edit makes of it, which is the check and, on request, the run. The document on
//! disk is never written: a command line has nowhere to keep a session, so the answer is the
//! document the command produced and saving it is the caller's business.

use std::path::Path;

use azoth_process::load_palette;
use azoth_process::middleware::command::Command;
use azoth_process::middleware::envelope;
use azoth_process::middleware::session::Workspace;

/// The rendered answer, and whether the document it made can run.
pub struct EditReport {
    pub text: String,
    pub ok: bool,
}

/// Apply one command to a flowsheet, and render what it left.
///
/// # Errors
/// A palette that does not load, a document that does not read, a command that is not a command,
/// or a run that was asked for and failed. **A document the edit broke is not an error**: it comes
/// back with `ok` false and the diagnostics that say so, which is the whole point of an editor.
pub fn edit(
    flowsheet_path: &Path,
    palette_dir: &Path,
    command: &str,
    run: bool,
    json: bool,
) -> Result<EditReport, String> {
    let palette = load_palette(palette_dir)?;
    let text = std::fs::read_to_string(flowsheet_path)
        .map_err(|error| format!("{}: {error}", flowsheet_path.display()))?;
    let mut workspace = Workspace::open(&text, palette)
        .map_err(|error| format!("{}: {error}", flowsheet_path.display()))?;

    let command: Command = serde_json::from_str(command)
        .map_err(|error| format!("the command is not one: {error}"))?;
    workspace
        .apply(&command)
        .map_err(|error| error.to_string())?;
    if run {
        // A run the edit made impossible is reported rather than raised: what a caller wants is
        // the document and why, not an exit code.
        let _ = workspace.run();
    }

    let text = if json {
        envelope::to_json(&workspace).map_err(|error| error.to_string())?
    } else {
        render(&workspace)
    };
    Ok(EditReport {
        text,
        ok: workspace.ok() && workspace.run_error().is_none(),
    })
}

/// The human reading: what the edit left, and what the checker says about it.
fn render(workspace: &Workspace) -> String {
    let mut lines = vec![format!(
        "{}: {}",
        workspace.flowsheet().id,
        if workspace.ok() { "OK" } else { "refused" }
    )];
    for diagnostic in workspace.diagnostics() {
        lines.push(format!("  {diagnostic}"));
    }
    if let Some(error) = workspace.run_error() {
        lines.push(format!("  the run refused it: {error}"));
    }
    if !workspace.dirty() {
        lines.push(format!(
            "  {} values, {} pass(es)",
            workspace.paths().len(),
            workspace
                .report()
                .and_then(Result::ok)
                .map_or(0, |report| report.iterations)
        ));
    }
    lines.join("\n")
}
