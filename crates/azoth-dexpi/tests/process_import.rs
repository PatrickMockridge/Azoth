//! The Process import: the scaffold it produces, and what it refuses.

mod common;

use azoth_dexpi::export_process;
use azoth_dexpi::process::export::ExportOptions;
use azoth_dexpi::process::import::{ImportOptions, import_process};
use azoth_process::validate;

fn exported_demo() -> String {
    let palette = common::palette();
    let flowsheet = common::demo();
    let options = ExportOptions {
        allow_lossy: true,
        project_recycles: true,
        azoth_extensions: true,
        ..ExportOptions::default()
    };
    export_process(&flowsheet, &palette, &options)
        .expect("the demo exports")
        .xml
}

#[test]
fn azoths_own_extension_exchange_round_trips_its_topology() {
    let palette = common::palette();
    let original = common::demo();
    let xml = exported_demo();
    let imported = import_process(&xml, &palette, &ImportOptions::default())
        .expect("the extension exchange imports");
    let flowsheet = &imported.flowsheet;

    let expected: Vec<(String, String)> = original
        .instances
        .iter()
        .map(|instance| (instance.id.clone(), instance.unit.clone()))
        .collect();
    let actual: Vec<(String, String)> = flowsheet
        .instances
        .iter()
        .map(|instance| (instance.id.clone(), instance.unit.clone()))
        .collect();
    assert_eq!(actual, expected);

    let mut expected_edges: Vec<(String, String)> = original
        .connections
        .iter()
        .map(|connection| (connection.from.clone(), connection.to.clone()))
        .collect();
    expected_edges.sort();
    let mut actual_edges: Vec<(String, String)> = flowsheet
        .connections
        .iter()
        .map(|connection| (connection.from.clone(), connection.to.clone()))
        .collect();
    actual_edges.sort();
    assert_eq!(actual_edges, expected_edges);

    assert_eq!(flowsheet.products, original.products);
    assert_eq!(flowsheet.recycles.len(), 1);
    assert_eq!(flowsheet.recycles[0].stream, original.recycles[0].stream);
    assert_eq!(flowsheet.recycles[0].from, original.recycles[0].from);
    assert_eq!(flowsheet.recycles[0].to, original.recycles[0].to);
}

#[test]
fn an_imported_scaffold_is_not_runnable_and_the_checker_says_so() {
    let palette = common::palette();
    let xml = exported_demo();
    let imported =
        import_process(&xml, &palette, &ImportOptions::default()).expect("the exchange imports");
    let diagnostics = validate(&imported.flowsheet, &palette);
    let codes: Vec<&str> = diagnostics
        .iter()
        .map(|diagnostic| diagnostic.code())
        .collect();
    assert!(
        codes.contains(&"unknown_feed"),
        "a feed is left undeclared: {codes:?}"
    );
    assert!(
        codes.contains(&"missing_parameter"),
        "no parameter is invented: {codes:?}"
    );
}

#[test]
fn the_report_names_the_gap_composition() {
    let palette = common::palette();
    let xml = exported_demo();
    let imported =
        import_process(&xml, &palette, &ImportOptions::default()).expect("the exchange imports");
    assert!(
        imported
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code == "input_composition_missing"
                && diagnostic.subject == "feed_1"),
        "{:?}",
        imported.diagnostics
    );
}

#[test]
fn strict_import_refuses_the_scaffold() {
    let palette = common::palette();
    let xml = exported_demo();
    let failure = import_process(&xml, &palette, &ImportOptions { strict: true })
        .expect_err("a strict import refuses warnings");
    assert!(failure.iter().any(|diagnostic| diagnostic.is_error()));
}

#[test]
fn a_plant_document_is_refused_as_the_wrong_serialisation() {
    let palette = common::palette();
    let failure = import_process(
        "<PlantModel><PlantInformation/></PlantModel>",
        &palette,
        &ImportOptions::default(),
    )
    .expect_err("a Proteus document is not a Process exchange");
    assert_eq!(failure[0].code, "unsupported_root");
}

#[test]
fn a_doctype_is_refused() {
    let palette = common::palette();
    let xml =
        "<?xml version=\"1.0\"?><!DOCTYPE Model [<!ENTITY x \"y\">]><Model name=\"m\" uri=\"u\"/>";
    let failure =
        import_process(xml, &palette, &ImportOptions::default()).expect_err("a DTD is refused");
    assert_eq!(failure[0].code, "doctype");
}
