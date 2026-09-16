//! The shipped flowsheets, loaded from `specs/flowsheets/`, must validate against
//! the shipped palette.

use std::path::Path;

use azoth_process::{Flowsheet, UnitOpSpec, validate};

fn repo_root() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../.."))
}

fn load_palette() -> Vec<UnitOpSpec> {
    let root = repo_root().join("specs/unit_ops");
    let mut specs = Vec::new();
    let mut families: Vec<_> = std::fs::read_dir(&root)
        .expect("specs/unit_ops must exist")
        .map(|e| e.expect("readable entry").path())
        .filter(|p| p.is_dir())
        .collect();
    families.sort();
    for family in families {
        for entry in std::fs::read_dir(&family).expect("readable family dir") {
            let path = entry.expect("readable entry").path();
            if path.extension().and_then(|e| e.to_str()) != Some("toml") {
                continue;
            }
            let text = std::fs::read_to_string(&path).unwrap();
            specs.push(toml::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display())));
        }
    }
    specs
}

fn load_flowsheets() -> Vec<Flowsheet> {
    let root = repo_root().join("specs/flowsheets");
    let mut sheets = Vec::new();
    for entry in std::fs::read_dir(&root).expect("specs/flowsheets must exist") {
        let path = entry.expect("readable entry").path();
        if path.extension().and_then(|e| e.to_str()) != Some("toml") {
            continue;
        }
        let text = std::fs::read_to_string(&path).unwrap();
        sheets.push(toml::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display())));
    }
    sheets
}

#[test]
fn every_shipped_flowsheet_validates_against_the_palette() {
    let palette = load_palette();
    let flowsheets = load_flowsheets();
    assert!(
        !flowsheets.is_empty(),
        "no flowsheets found under specs/flowsheets"
    );

    for flowsheet in &flowsheets {
        let diags = validate(flowsheet, &palette);
        assert!(
            diags.is_empty(),
            "{} does not validate: {diags:?}",
            flowsheet.id
        );
    }
}
