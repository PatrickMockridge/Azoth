//! A Proteus drawing's geometry: shape catalogue, placement, segment centre lines and the
//! join from a drawn group to the projected object.

mod common;

use azoth_dexpi::drawing::{GroupKind, Primitive};
use azoth_dexpi::import_proteus;
use azoth_dexpi::plant::import::PlantImportOptions;

fn fixture() -> String {
    std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("proteus_simple.xml"),
    )
    .expect("the Proteus fixture reads")
}

#[test]
fn the_proteus_drawing_resolves_its_shapes_and_positions() {
    let report = import_proteus(
        &fixture(),
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the Proteus fixture imports");
    let drawing = report.drawing.as_ref().expect("the file carries graphics");
    assert!(drawing.groups.len() >= 8, "{:?}", drawing.groups.len());

    let pump = report
        .bindings
        .iter()
        .find(|(_, name)| name.as_str() == "P-100")
        .map(|(id, _)| id.clone())
        .expect("the pump is bound");
    let group = drawing
        .groups
        .iter()
        .find(|group| group.represents.as_deref() == Some(pump.as_str()))
        .expect("the pump is drawn");
    assert!(
        group
            .primitives
            .iter()
            .any(|primitive| matches!(primitive, Primitive::Ellipse { .. })),
        "the pump's circle became an ellipse"
    );
    assert!(
        group
            .primitives
            .iter()
            .any(|primitive| matches!(primitive, Primitive::Text { text, .. } if text == "P-100")),
        "the pump's label is attached to its group"
    );

    let edges = drawing
        .groups
        .iter()
        .filter(|group| group.kind == GroupKind::Edge)
        .count();
    assert!(edges >= 4, "each segment's centre line is an edge");
    assert!(
        drawing
            .groups
            .iter()
            .any(|group| group.represents.is_none() && group.kind == GroupKind::Node),
        "the reducer has no azoth unit and is drawing-only"
    );

    let switch = report
        .bindings
        .iter()
        .find(|(_, name)| name.as_str() == "PSHH-101")
        .map(|(id, _)| id.clone())
        .expect("the switch is bound");
    let bubble = drawing
        .groups
        .iter()
        .find(|group| group.represents.as_deref() == Some(switch.as_str()))
        .expect("the instrument bubble is drawn");
    assert!(
        bubble
            .primitives
            .iter()
            .any(|primitive| matches!(primitive, Primitive::Ellipse { .. })),
        "the bubble's circle became an ellipse"
    );
    assert!(
        bubble
            .primitives
            .iter()
            .any(|primitive| matches!(primitive, Primitive::Text { text, .. } if text == "PSHH")),
        "the instrument's label is attached to its group"
    );
}

#[test]
fn the_proteus_drawing_places_the_instances() {
    let report = import_proteus(
        &fixture(),
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the Proteus fixture imports");
    let layout = report
        .flowsheet
        .layout
        .as_ref()
        .expect("the drawing places the scaffold");
    assert!(layout.instances.contains_key("P-100"));
    assert!(layout.instances.contains_key("H-100"));
    assert!(layout.instances.contains_key("T-100"));
    assert!(layout.instances.contains_key("GlobeValve-1"));
    let x = |name: &str| layout.instances.get(name).expect("placed")[0];
    assert!(x("P-100") < x("H-100"));
    assert!(x("H-100") < x("T-100"));
    assert!(layout.products.contains_key("PRODUCT"));
    assert!(!layout.inputs.contains_key("FEED"));
}

#[test]
fn a_proteus_import_with_graphics_still_reads_the_topology() {
    let report = import_proteus(
        &fixture(),
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the Proteus fixture imports");
    let units: Vec<&str> = report
        .flowsheet
        .instances
        .iter()
        .map(|instance| instance.unit.as_str())
        .collect();
    assert!(units.contains(&"unit_ops.pump"));
    assert!(units.contains(&"unit_ops.tank"));
    assert!(units.contains(&"unit_ops.throttling_valve"));
}
