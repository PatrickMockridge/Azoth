//! Reading the process layer's format: the palette and the flowsheet, as TOML.
//!
//! The Rust types are the schema, and this is how a document becomes one of them.
//! Both the `azoth check` command and the Python binding read through here, so
//! there is one loader rather than two that can disagree about the format.

use std::path::{Path, PathBuf};

use crate::flowsheet::Flowsheet;
use crate::unit_op::UnitOpSpec;

/// Every palette entry under a directory, recursively.
pub fn load_palette(dir: &Path) -> Result<Vec<UnitOpSpec>, String> {
    let mut paths = Vec::new();
    collect_toml(dir, &mut paths).map_err(|e| format!("{}: {e}", dir.display()))?;
    paths.sort();

    let mut specs = Vec::new();
    for path in paths {
        let text =
            std::fs::read_to_string(&path).map_err(|e| format!("{}: {e}", path.display()))?;
        let spec: UnitOpSpec =
            toml::from_str(&text).map_err(|e| format!("{}: {e}", path.display()))?;
        specs.push(spec);
    }
    Ok(specs)
}

/// Every palette entry in an embedded bundle, in the bundle's own order.
///
/// **The other source of the same loader.** `load_palette` walks a directory, which a browser
/// cannot do; a wasm module carries `palette_gen::PALETTE` and parses it through here, so there is
/// one place that knows what a palette entry is. The bundle names its entries by their path
/// relative to `specs/unit_ops`, which is what makes the two orders the same order.
///
/// # Errors
/// A bundle entry that is not a palette entry, named with the file it came from.
pub fn load_palette_text(entries: &[(&str, &str)]) -> Result<Vec<UnitOpSpec>, String> {
    entries
        .iter()
        .map(|(name, text)| {
            toml::from_str(text).map_err(|error| format!("specs/unit_ops/{name}: {error}"))
        })
        .collect()
}

/// A flowsheet, from its TOML text.
pub fn parse_flowsheet(text: &str) -> Result<Flowsheet, String> {
    // One parse rather than two spellings of it: `Flowsheet::from_toml` is where the document's
    // reader lives, and a second `toml::from_str` here would be a second place a schema change
    // could be made in one of them.
    Flowsheet::from_toml(text).map_err(|error| error.to_string())
}

fn collect_toml(dir: &Path, out: &mut Vec<PathBuf>) -> std::io::Result<()> {
    for entry in std::fs::read_dir(dir)? {
        let path = entry?.path();
        if path.is_dir() {
            collect_toml(&path, out)?;
        } else if path.extension().and_then(|e| e.to_str()) == Some("toml") {
            out.push(path);
        }
    }
    Ok(())
}
