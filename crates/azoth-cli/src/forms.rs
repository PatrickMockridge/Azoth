//! `azoth forms` — the palette as one form per unit operation.
//!
//! **The output *is* the schema, so there is no `--json`.** `run` prints a report a person reads
//! and writes JSON on request; this prints the form a front-end renders, and a human-readable
//! form of it would be a second shape for the same declaration. What it is for is the thing on
//! the other side of the pipe: a widget generator, a test fixture, an agent's tool list.
//!
//! The document itself is assembled in `middleware::catalogue`, which is also where the wasm
//! module gets it - so a command line and a browser cannot describe the palette differently.

use std::path::Path;

use azoth_process::load_palette;
use azoth_process::middleware::catalogue;

/// The palette's forms, as JSON.
///
/// `with_tools` adds the agent's tool schema, which is the same declaration read as the commands
/// an agent may send rather than as the fields a person fills.
pub fn forms_json(palette_dir: &Path, with_tools: bool) -> Result<String, String> {
    let palette = load_palette(palette_dir)?;
    catalogue(&palette, with_tools).map_err(|error| error.to_string())
}
