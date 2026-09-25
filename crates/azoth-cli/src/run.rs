//! `azoth run` - run a flowsheet to its steady state.
//!
//! `azoth check` says a flowsheet is wired correctly; this says what it *reaches*. The document
//! is self-contained - its `[[inputs]]` declare the fluid and the state, so nothing is passed on
//! the command line but the two paths - and the report is every named stream plus the tear's
//! convergence.
//!
//! **It validates before it runs.** A document with an error diagnostic is one the executor would
//! refuse or, worse, misread, so the checker's verdict is printed and the run does not start.
//! Warnings do not stop it: `UnusedFeed` is a document that runs and is probably not what was
//! meant, which is the line the two severities are drawn on.
//!
//! **`--json` prints the codec's document rather than a second rendering of it.** The same string
//! is what `azoth.process.run_flowsheet(...).document` hands a Python caller, so a script that
//! reads one reads the other.

use std::path::Path;

use azoth_process::executor::{Session, to_json};
use azoth_process::{ExecutionOrder, load_palette, parse_flowsheet, validate};

/// The rendered report, and whether the flowsheet ran.
pub struct RunReport {
    pub text: String,
    pub ran: bool,
}

/// Run one flowsheet file against a palette directory.
///
/// `as_json` selects the codec's document instead of the rendered report.
pub fn run(flowsheet_path: &Path, palette_dir: &Path, as_json: bool) -> Result<RunReport, String> {
    let palette = load_palette(palette_dir)?;
    let text = std::fs::read_to_string(flowsheet_path)
        .map_err(|e| format!("{}: {e}", flowsheet_path.display()))?;
    let flowsheet =
        parse_flowsheet(&text).map_err(|e| format!("{}: {e}", flowsheet_path.display()))?;

    let diags = validate(&flowsheet, &palette);
    let mut lines: Vec<String> = diags
        .iter()
        .filter(|d| d.severity() == azoth_process::Severity::Error)
        .map(|d| format!("{}: {d:?}", flowsheet.id))
        .collect();
    if !lines.is_empty() {
        lines.push(format!(
            "{}: {} error(s) - a flowsheet that does not validate is not one to run",
            flowsheet.id,
            lines.len()
        ));
        return Ok(RunReport {
            text: lines.join("\n"),
            ran: false,
        });
    }

    // The document's own inputs; nothing is overridden from a command line, because a boundary
    // stated in two places is a boundary that goes stale in one of them.
    let session = Session::run(
        flowsheet,
        &palette,
        &Default::default(),
        ExecutionOrder::Insertion,
    )
    .map_err(|e| format!("{}: {e}", flowsheet_path.display()))?;

    if as_json {
        return Ok(RunReport {
            text: to_json(&session).map_err(|e| e.to_string())?,
            ran: true,
        });
    }

    let mut lines = vec![session.flowsheet().id.clone()];
    for (endpoint, stream) in &session.report().streams {
        lines.push(format!(
            "  {endpoint}: n={} mol/s  z=[{}]  P={} Pa  T={} K  h={} J/mol",
            stream.n,
            stream
                .z
                .iter()
                .map(|z| format!("{z:.6}"))
                .collect::<Vec<_>>()
                .join(", "),
            stream.p.value,
            stream.t.value,
            stream.h.value,
        ));
    }
    for tear in &session.report().tears {
        // **`active` is printed beside `solved` rather than folded into it.** A tear the low-flow
        // cutoff switched off reports `solved` by being absent - its residuals are declared zero,
        // not measured - and a reader who saw only `solved=true` would take a loop that carries
        // nothing for one that closed.
        let state = if tear.active { "active" } else { "deactivated" };
        match tear.residuals {
            Some(residuals) => lines.push(format!(
                "  {}: iterations={} solved={} {state} flow={} composition={} temperature={} \
                 pressure={}",
                tear.stream,
                tear.iterations,
                tear.solved,
                residuals.flow,
                residuals.composition,
                residuals.temperature,
                residuals.pressure,
            )),
            None => lines.push(format!(
                "  {}: iterations={} solved={} {state} residuals=unmeasured",
                tear.stream, tear.iterations, tear.solved
            )),
        }
    }
    lines.push(format!(
        "  converged={} iterations={}",
        session.converged(),
        session.iterations()
    ));

    Ok(RunReport {
        text: lines.join("\n"),
        ran: true,
    })
}
