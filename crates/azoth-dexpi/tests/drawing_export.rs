//! A flowsheet drawn as a PFD and a simplified P&ID, and the arrangement read back.

mod common;

use std::collections::BTreeSet;

use azoth_dexpi::plant::export::{PlantExportOptions, export_plant};
use azoth_dexpi::plant::import::{PlantImportOptions, import_plant};
use azoth_dexpi::process::export::{ExportOptions, export_process};
use azoth_dexpi::process::import::{ImportOptions, import_process};
use azoth_dexpi::reader::read;

fn pfd() -> azoth_dexpi::ExportReport {
    export_process(
        &common::demo(),
        &common::palette(),
        &ExportOptions {
            allow_lossy: true,
            project_recycles: true,
            azoth_extensions: true,
            drawing: true,
            ..ExportOptions::default()
        },
    )
    .expect("the demo exports a PFD")
}

fn pid() -> azoth_dexpi::PlantExportReport {
    export_plant(
        &common::plant_simple(),
        &common::palette(),
        &PlantExportOptions {
            allow_lossy: true,
            azoth_extensions: true,
            drawing: true,
            ..PlantExportOptions::default()
        },
    )
    .expect("the plant fixture exports a P&ID")
}

#[test]
fn the_demo_pfd_draws_every_step_and_stream() {
    let report = pfd();
    assert_eq!(
        report.diagram_groups, 14,
        "seven nodes (input, four units, product, recycle) and seven streams"
    );
    assert!(report.diagram_primitives >= 14);
    assert!(report.xml.contains("Core/Diagram.Diagram"));
    assert!(report.xml.contains("Core/Diagram.PolyLine"));
    assert!(report.xml.contains("kg/h"));
    let path = common::temp_path("pfd_drawing.xml");
    std::fs::write(&path, &report.xml).expect("the PFD writes");
    common::xmllint(&path);
    assert!(unresolved(&report.xml).is_empty());
}

#[test]
fn the_pfd_is_deterministic() {
    assert_eq!(pfd().xml, pfd().xml);
}

#[test]
fn a_topology_only_pfd_carries_the_drawing_but_no_operating_values() {
    let report = export_process(
        &common::demo(),
        &common::palette(),
        &ExportOptions {
            allow_lossy: true,
            project_recycles: true,
            drawing: true,
            topology_only: true,
            ..ExportOptions::default()
        },
    )
    .expect("the topology-only PFD exports");
    assert!(report.diagram_groups > 0);
    assert!(!report.xml.contains("kg/h"));
}

#[test]
fn a_pfd_import_places_the_drawing() {
    let xml = pfd().xml;
    let imported = import_process(&xml, &common::palette(), &ImportOptions::default())
        .expect("the PFD imports");
    let layout = imported
        .flowsheet
        .layout
        .as_ref()
        .expect("the drawing places the scaffold");
    assert_eq!(layout.instances.len(), 4);
    assert_eq!(layout.products.len(), 1);
    assert_eq!(imported.layout_positions, 5);
    let x = |name: &str| layout.instances.get(name).expect("placed")[0];
    assert!(x("mix1") < x("p1"));
    assert!(x("p1") < x("hx1"));
    assert!(x("hx1") < x("sep1"));
}

#[test]
fn the_plant_pid_draws_items_segments_and_instruments() {
    let report = pid();
    assert_eq!(
        report.diagram_groups, 14,
        "three items, three boundaries, five segments, two instruments and one signal"
    );
    assert!(report.xml.contains("PSHH-101"));
    assert!(report.xml.contains("L-1001"));
    let path = common::temp_path("pid_drawing.xml");
    std::fs::write(&path, &report.xml).expect("the P&ID writes");
    common::xmllint(&path);
    assert!(unresolved(&report.xml).is_empty());
}

#[test]
fn a_pid_import_places_the_drawing() {
    let imported = import_plant(
        &pid().xml,
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the P&ID imports");
    let layout = imported
        .flowsheet
        .layout
        .as_ref()
        .expect("the drawing places the scaffold");
    assert_eq!(layout.instances.len(), 3);
    assert_eq!(layout.products.len(), 2);
    assert_eq!(imported.layout_positions, 5);
    let x = |name: &str| layout.instances.get(name).expect("placed")[0];
    assert!(x("V-101") < x("K-101"));
    assert!(x("K-101") < x("E-101"));
}

#[test]
fn a_drawing_export_still_passes_the_semantic_rules() {
    let pfd = pfd();
    assert!(pfd.semantic_valid);
    let pid = pid();
    assert!(pid.semantic_valid);
}

/// Every `Represents` target in a drawing resolves to an object in the document.
fn unresolved(xml: &str) -> Vec<String> {
    let document = read(xml).expect("the export reads").document;
    let mut ids = BTreeSet::new();
    collect_ids(&document.model.children, &mut ids);
    let mut missing = Vec::new();
    collect_references(&document.model.children, &ids, &mut missing);
    missing
}

fn collect_ids(children: &[azoth_dexpi::Child], ids: &mut BTreeSet<String>) {
    for child in children {
        if let azoth_dexpi::Child::Object(object) = child {
            if let Some(id) = &object.id {
                ids.insert(id.clone());
            }
            collect_ids(&object.children, ids);
        }
        if let azoth_dexpi::Child::Components(components) = child {
            for nested in &components.children {
                if let azoth_dexpi::ComponentChild::Object(object) = nested {
                    if let Some(id) = &object.id {
                        ids.insert(id.clone());
                    }
                    collect_ids(&object.children, ids);
                }
            }
        }
    }
}

fn collect_references(
    children: &[azoth_dexpi::Child],
    ids: &BTreeSet<String>,
    missing: &mut Vec<String>,
) {
    for child in children {
        match child {
            azoth_dexpi::Child::References(references) if references.property == "Represents" => {
                for target in &references.objects {
                    let target = target.trim_start_matches('#');
                    if !ids.contains(target) {
                        missing.push(target.to_string());
                    }
                }
            }
            azoth_dexpi::Child::Object(object) => {
                collect_references(&object.children, ids, missing)
            }
            azoth_dexpi::Child::Components(components) => {
                for nested in &components.children {
                    if let azoth_dexpi::ComponentChild::Object(object) = nested {
                        collect_references(&object.children, ids, missing);
                    }
                }
            }
            _ => {}
        }
    }
}
