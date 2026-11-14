//! Shared fixtures: the shipped palette and flowsheet, read from the workspace.
//!
//! Each integration test compiles this module separately, so a helper only one of them uses is
//! dead code in the others.
#![allow(dead_code)]

use std::path::{Path, PathBuf};

use azoth_process::{UnitOpSpec, load_palette};

pub fn workspace_root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("..")
        .join("..")
        .canonicalize()
        .expect("the workspace root resolves")
}

pub fn palette() -> Vec<UnitOpSpec> {
    load_palette(&workspace_root().join("specs").join("unit_ops")).expect("the palette loads")
}

/// The topology `validation/neqsim/Dexpi20ProcessProbe.java` builds.
pub fn oracle_simple() -> azoth_process::Flowsheet {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("oracle_simple.toml"),
    )
    .expect("the oracle fixture reads");
    azoth_process::parse_flowsheet(&text).expect("the oracle fixture parses")
}

/// The Plant export fixture: only units the Plant writer reviewed.
pub fn plant_simple() -> azoth_process::Flowsheet {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("plant_simple.toml"),
    )
    .expect("the Plant fixture reads");
    azoth_process::parse_flowsheet(&text).expect("the Plant fixture parses")
}

/// A pump, a valve and a separator: the valve is a segment item, not a tagged item.
pub fn plant_valve() -> azoth_process::Flowsheet {
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("plant_valve.toml"),
    )
    .expect("the valve fixture reads");
    azoth_process::parse_flowsheet(&text).expect("the valve fixture parses")
}

pub fn demo() -> azoth_process::Flowsheet {
    let text = std::fs::read_to_string(
        workspace_root()
            .join("specs")
            .join("flowsheets")
            .join("demo.toml"),
    )
    .expect("the demo flowsheet reads");
    azoth_process::parse_flowsheet(&text).expect("the demo flowsheet parses")
}

/// A minimal valid flowsheet around one unit, for the refusal tests.
pub fn single_unit(unit_id: &str) -> azoth_process::Flowsheet {
    let text = format!(
        r#"
id = "flowsheets.refusal"
name = "One unit"

products = ["vent"]

[[inputs]]
name = "feed_1"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e5
T = 300.0

[[instances]]
id = "u1"
unit = "{unit_id}"

[[connections]]
from = "feed_1"
to = "u1.inlet"

[[connections]]
from = "u1.product"
to = "vent"
"#
    );
    azoth_process::parse_flowsheet(&text).expect("the fixture parses")
}

/// A unique path under the system temporary directory.
pub fn temp_path(name: &str) -> PathBuf {
    let mut path = std::env::temp_dir();
    path.push(format!("azoth-dexpi-{}-{name}", std::process::id()));
    path
}

/// Validate a file with `xmllint`, failing (not skipping) when the tool is absent.
pub fn xmllint(xml_path: &Path) {
    let schema = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("schema")
        .join("DEXPI_XML_Schema.xsd");
    let output = std::process::Command::new("xmllint")
        .args(["--nonet", "--noout", "--schema"])
        .arg(&schema)
        .arg(xml_path)
        .output()
        .expect(
            "xmllint is required: install libxml2-utils. A schema gate that silently skips when \
             the tool is missing is a check nobody knows did not run",
        );
    assert!(
        output.status.success(),
        "the document failed the official DEXPI schema:\n{}",
        String::from_utf8_lossy(&output.stderr)
    );
}
