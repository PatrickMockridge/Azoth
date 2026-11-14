//! The Proteus projection: a DEXPI 1.x `PlantModel` read as the Plant metamodel.

mod common;

use std::path::Path;

use azoth_dexpi::plant::import::PlantImportOptions;
use azoth_dexpi::proteus::import_proteus;

fn fixture() -> String {
    std::fs::read_to_string(
        Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("proteus_simple.xml"),
    )
    .expect("the Proteus fixture reads")
}

#[test]
fn the_projection_recovers_the_items_and_their_units() {
    let report = import_proteus(
        &fixture(),
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the Proteus fixture imports");
    let instances: Vec<(String, String)> = report
        .flowsheet
        .instances
        .iter()
        .map(|instance| (instance.id.clone(), instance.unit.clone()))
        .collect();
    for expected in [
        ("P-100".to_string(), "unit_ops.pump".to_string()),
        ("H-100".to_string(), "unit_ops.heat_exchanger".to_string()),
        ("T-100".to_string(), "unit_ops.tank".to_string()),
        (
            "GlobeValve-1".to_string(),
            "unit_ops.throttling_valve".to_string(),
        ),
    ] {
        assert!(
            instances.contains(&expected),
            "{expected:?} in {instances:?}"
        );
    }
}

#[test]
fn the_projection_recovers_the_chain_through_the_valve_and_the_pass_through() {
    let report = import_proteus(
        &fixture(),
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the Proteus fixture imports");
    let edges: Vec<(String, String)> = report
        .flowsheet
        .connections
        .iter()
        .map(|connection| (connection.from.clone(), connection.to.clone()))
        .collect();
    for expected in [
        ("FEED".to_string(), "P-100.inlet".to_string()),
        ("P-100.outlet".to_string(), "H-100.hot_in".to_string()),
        (
            "H-100.hot_out".to_string(),
            "GlobeValve-1.inlet".to_string(),
        ),
        ("GlobeValve-1.outlet".to_string(), "T-100.inlet".to_string()),
        ("T-100.gas".to_string(), "PRODUCT".to_string()),
    ] {
        assert!(edges.contains(&expected), "{expected:?} in {edges:?}");
    }
    assert_eq!(report.flowsheet.products, ["PRODUCT"]);
    assert!(
        report
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code == "inline_component_skipped"
                && diagnostic.subject == "PipeReducer-1"),
        "{:?}",
        report.diagnostics
    );
}

#[test]
fn the_line_number_survives_every_connection() {
    let report = import_proteus(
        &fixture(),
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the Proteus fixture imports");
    assert!(!report.flowsheet.connections.is_empty());
    for connection in &report.flowsheet.connections {
        assert_eq!(connection.line.as_deref(), Some("L-1001"), "{connection:?}");
    }
}

#[test]
fn a_dexpi_2_document_is_not_proteus() {
    let failure = import_proteus(
        "<Model name=\"m\" uri=\"u\"/>",
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect_err("a DEXPI 2.0 root is not a Proteus file");
    assert_eq!(failure[0].code, "unsupported_root");
}
