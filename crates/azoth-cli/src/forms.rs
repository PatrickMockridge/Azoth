//! `azoth forms` — the palette as one form per unit operation.
//!
//! **The output *is* the schema, so there is no `--json`.** `run` prints a report a person reads
//! and writes JSON on request; this prints the form a front-end renders, and a human-readable
//! form of it would be a second shape for the same declaration. What it is for is the thing on
//! the other side of the pipe: a widget generator, a test fixture, an agent's tool list.

use std::path::Path;

use azoth_process::load_palette;
use azoth_process::middleware::form::forms;
use azoth_process::middleware::tools::tools;
use serde_json::json;

/// The palette's forms, as JSON.
///
/// `tools` asks for the agent's tool schema as well, which is the same declaration read as the
/// commands an agent may send rather than as the fields a person fills.
pub fn forms_json(palette_dir: &Path, with_tools: bool) -> Result<String, String> {
    let palette = load_palette(palette_dir)?;
    let mut document = json!({ "unit_ops": forms(&palette) });
    if with_tools {
        let document = document
            .as_object_mut()
            .expect("the object this function just built");
        document.insert("tools".to_string(), json!(tools(&palette)));
    }
    serde_json::to_string_pretty(&document).map_err(|error| format!("the forms: {error}"))
}
