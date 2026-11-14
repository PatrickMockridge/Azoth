//! A valve in a segment: written inside its `Items`, split around its nodes, and read back.

mod common;

use azoth_dexpi::plant;
use azoth_dexpi::plant::export::{PlantExportOptions, export_plant};
use azoth_dexpi::plant::import::{PlantImportOptions, import_plant};
use azoth_dexpi::reader::read;
use azoth_dexpi::view::{CoverageStatus, inspect};

fn pid() -> azoth_dexpi::PlantExportReport {
    export_plant(
        &common::plant_valve(),
        &common::palette(),
        &PlantExportOptions {
            allow_lossy: true,
            azoth_extensions: true,
            drawing: true,
            ..PlantExportOptions::default()
        },
    )
    .expect("the valve fixture exports a P&ID")
}

#[test]
fn the_valve_is_written_inside_a_segment_not_as_an_item() {
    let report = pid();
    assert!(report.semantic_valid);
    assert_eq!(
        report.items, 2,
        "the pump and the separator are the tagged items"
    );
    assert_eq!(
        report.segments, 4,
        "the two streams through the valve are one segment"
    );
    assert!(report.xml.contains("Plant/Piping.GlobeValve"));
    assert!(report.xml.contains("<String>V-101</String>"));
    let path = common::temp_path("plant_valve.xml");
    std::fs::write(&path, &report.xml).expect("the P&ID writes");
    common::xmllint(&path);
}

#[test]
fn the_segment_carries_the_valve_and_two_pipes() {
    let report = pid();
    let document = read(&report.xml).expect("the export reads").document;
    let model = find_plant(&document).expect("the plant model");
    let systems = model.components("PipingNetworkSystems").expect("systems");
    let valve_segment = systems
        .children
        .iter()
        .filter_map(|child| match child {
            azoth_dexpi::ComponentChild::Object(object) => Some(object),
            _ => None,
        })
        .flat_map(|system| {
            system
                .components("Segments")
                .map(|segments| segments.children.iter())
                .into_iter()
                .flatten()
        })
        .filter_map(|child| match child {
            azoth_dexpi::ComponentChild::Object(object) => Some(object),
            _ => None,
        })
        .find(|segment| {
            segment.components("Items").is_some_and(|items| {
                items.children.iter().any(|child| {
                    matches!(child, azoth_dexpi::ComponentChild::Object(object)
                        if object.ty == "Plant/Piping.GlobeValve")
                })
            })
        })
        .expect("a segment carries the valve");
    let items = valve_segment.components("Items").expect("items");
    let valve = items
        .children
        .iter()
        .find_map(|child| match child {
            azoth_dexpi::ComponentChild::Object(object)
                if object.ty == "Plant/Piping.GlobeValve" =>
            {
                Some(object)
            }
            _ => None,
        })
        .expect("the valve item");
    let nodes = valve.components("Nodes").expect("the valve's nodes");
    assert_eq!(nodes.children.len(), 2);
    let pipes = valve_segment.components("Connections").expect("pipes");
    assert_eq!(
        pipes.children.len(),
        2,
        "one pipe to the valve, one from it"
    );
}

#[test]
fn the_valve_survives_the_round_trip() {
    let imported = import_plant(
        &pid().xml,
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the P&ID imports");
    let valve = imported
        .flowsheet
        .instances
        .iter()
        .find(|instance| instance.id == "V-101")
        .expect("the valve is an instance again");
    assert_eq!(valve.unit, "unit_ops.throttling_valve");
    let around: Vec<&str> = imported
        .flowsheet
        .connections
        .iter()
        .filter(|connection| connection.from.contains("V-101") || connection.to.contains("V-101"))
        .map(|connection| {
            if connection.from.contains("V-101") {
                connection.to.as_str()
            } else {
                connection.from.as_str()
            }
        })
        .collect();
    assert_eq!(around.len(), 2, "{:?}", imported.flowsheet.connections);
    assert!(around.contains(&"P-101.outlet"));
    assert!(around.contains(&"SEP-101.feed"));
}

#[test]
fn the_valve_is_modelled_and_writable_in_the_view() {
    let view = inspect(&pid().xml, &common::palette()).expect("the P&ID is viewable");
    let valve = view
        .objects
        .iter()
        .find(|object| object.label == "V-101")
        .expect("the valve is in the coverage");
    assert_eq!(valve.status, CoverageStatus::Modelled);
    assert!(valve.exportable);
}

#[test]
fn the_mapping_counts_the_valve_as_writable() {
    let writable = plant::units()
        .iter()
        .filter(|row| row.write != azoth_dexpi::Disposition::Refused)
        .count();
    assert_eq!(writable, 10);
    assert!(plant::is_segment_unit("unit_ops.throttling_valve"));
    assert!(!plant::is_segment_unit("unit_ops.pump"));
}

fn find_plant(document: &azoth_dexpi::Document) -> Option<&azoth_dexpi::Object> {
    fn walk(object: &azoth_dexpi::Object) -> Option<&azoth_dexpi::Object> {
        if object.ty == "Plant/PlantModel" {
            return Some(object);
        }
        for child in &object.children {
            if let azoth_dexpi::Child::Components(components) = child {
                for nested in &components.children {
                    if let azoth_dexpi::ComponentChild::Object(nested) = nested {
                        if let Some(found) = walk(nested) {
                            return Some(found);
                        }
                    }
                }
            }
        }
        None
    }
    document
        .model
        .children
        .iter()
        .find_map(|child| match child {
            azoth_dexpi::Child::Object(object) => walk(object),
            _ => None,
        })
}
