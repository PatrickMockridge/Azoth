//! `azoth check` - validate a flowsheet against the palette.
//!
//! The same checks the build runs, exposed as a command: load the palette, load a
//! flowsheet, and report every way the flowsheet breaks the calculus's rules. A
//! flowsheet that prints `OK` is one a GUI or an executor could consume.

use std::path::Path;

use azoth_process::{load_palette, parse_flowsheet, validate, validate_palette};

/// The rendered report, and how many defects it found.
pub struct CheckReport {
    pub text: String,
    pub defects: usize,
}

/// Validate one flowsheet file against a palette directory.
pub fn check(flowsheet_path: &Path, palette_dir: &Path) -> Result<CheckReport, String> {
    let palette = load_palette(palette_dir)?;
    let palette_diags = validate_palette(&palette);

    let text = std::fs::read_to_string(flowsheet_path)
        .map_err(|e| format!("{}: {e}", flowsheet_path.display()))?;
    let flowsheet =
        parse_flowsheet(&text).map_err(|e| format!("{}: {e}", flowsheet_path.display()))?;

    let mut lines = Vec::new();
    for d in &palette_diags {
        lines.push(format!("palette: {d:?}"));
    }
    let diags = validate(&flowsheet, &palette);
    if diags.is_empty() {
        lines.push(format!("{}: OK", flowsheet.id));
    } else {
        for d in &diags {
            lines.push(format!("{}: {d:?}", flowsheet.id));
        }
    }

    Ok(CheckReport {
        text: lines.join("\n"),
        defects: palette_diags.len() + diags.len(),
    })
}
