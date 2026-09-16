//! The shipped palette, loaded from `specs/unit_ops/`, must be well-formed.
//!
//! This is the CI gate for the palette: every entry parses as a `UnitOpSpec`,
//! names only dimensions the vocabulary knows, carries a provenance line, and has
//! ports with fields. It reads the files from disk rather than an embedded copy,
//! so it validates what a GUI would read, not a snapshot.

use std::collections::HashSet;
use std::path::Path;

use azoth_process::{UnitOpSpec, validate_palette};

fn unit_op_root() -> &'static Path {
    Path::new(concat!(env!("CARGO_MANIFEST_DIR"), "/../../specs/unit_ops"))
}

fn load_palette() -> Vec<UnitOpSpec> {
    let mut specs = Vec::new();

    let mut families: Vec<_> = std::fs::read_dir(unit_op_root())
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
            let text = std::fs::read_to_string(&path)
                .unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            let spec: UnitOpSpec =
                toml::from_str(&text).unwrap_or_else(|e| panic!("{}: {e}", path.display()));
            specs.push(spec);
        }
    }

    specs
}

#[test]
fn every_palette_entry_is_well_formed() {
    let specs = load_palette();
    assert!(
        !specs.is_empty(),
        "no palette entries found under specs/unit_ops"
    );

    let diags = validate_palette(&specs);
    assert!(diags.is_empty(), "palette defects: {diags:?}");
}

#[test]
fn every_palette_entry_is_named_sourced_and_typed() {
    let mut ids: HashSet<String> = HashSet::new();
    for spec in load_palette() {
        assert!(
            ids.insert(spec.id.clone()),
            "duplicate unit-op id: {}",
            spec.id
        );
        assert!(
            spec.source.is_some(),
            "{} carries no `source`, so its port shape has no provenance",
            spec.id
        );
        assert!(!spec.ports.is_empty(), "{} declares no ports", spec.id);
        for port in &spec.ports {
            assert!(
                !port.fields.is_empty(),
                "{}.{} declares no fields",
                spec.id,
                port.name
            );
        }
    }
}
