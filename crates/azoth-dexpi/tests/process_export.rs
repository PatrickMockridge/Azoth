//! The Process export: what it writes, what it refuses, and the schema it must pass.

mod common;

use azoth_dexpi::export_process;
use azoth_dexpi::process::export::ExportOptions;

fn export_demo() -> azoth_dexpi::ExportReport {
    let palette = common::palette();
    let flowsheet = common::demo();
    let options = ExportOptions {
        allow_lossy: true,
        project_recycles: true,
        azoth_extensions: true,
        ..ExportOptions::default()
    };
    export_process(&flowsheet, &palette, &options).expect("the demo exports")
}

#[test]
fn the_demo_flowsheet_exports_the_steps_its_units_map_to() {
    let report = export_demo();
    assert!(report.semantic_valid);
    assert!(report.converged);
    assert_eq!(report.approval_status, "NOT_A_DEXPI_EV_CERTIFICATE");
    assert_eq!(report.cae_round_trip, "QUALIFICATION_REQUIRED");

    let inventory = azoth_dexpi::inventory::inventory(
        &azoth_dexpi::reader::read(&report.xml)
            .expect("the export reads back")
            .document,
    );
    assert_eq!(
        inventory.get("step.feed_1.type").map(String::as_str),
        Some("Process/Process.Source")
    );
    assert_eq!(
        inventory.get("step.mix1.type").map(String::as_str),
        Some("Process/Process.MixingSimple")
    );
    assert_eq!(
        inventory.get("step.p1.type").map(String::as_str),
        Some("Process/Process.Pumping")
    );
    assert_eq!(
        inventory.get("step.hx1.type").map(String::as_str),
        Some("Process/Process.HeatingInFurnace")
    );
    assert_eq!(
        inventory.get("step.sep1.type").map(String::as_str),
        Some("Process/Process.SeparatingByGravity")
    );
    assert_eq!(
        inventory
            .get("step.vapour_product.type")
            .map(String::as_str),
        Some("Process/Process.Sink")
    );
}

#[test]
fn the_demo_export_passes_the_official_schema() {
    let report = export_demo();
    let path = common::temp_path("demo.dexpi.xml");
    std::fs::write(&path, &report.xml).expect("the fixture writes");
    common::xmllint(&path);
}

#[test]
fn a_lossy_mapping_needs_the_flag() {
    let palette = common::palette();
    let flowsheet = common::demo();
    let options = ExportOptions {
        project_recycles: true,
        ..ExportOptions::default()
    };
    let failure = export_process(&flowsheet, &palette, &options).expect_err("a lossy unit refuses");
    assert!(
        failure
            .iter()
            .any(|diagnostic| diagnostic.code == "lossy_mapping" && diagnostic.subject == "hx1"),
        "the refusal names the heater: {failure:?}"
    );
}

#[test]
fn a_refused_mapping_is_refused_by_name() {
    let palette = common::palette();
    let flowsheet = common::single_unit("unit_ops.flare");
    let options = ExportOptions::default();
    let failure = export_process(&flowsheet, &palette, &options).expect_err("a flare refuses");
    assert_eq!(failure[0].code, "refused_mapping");
    assert_eq!(failure[0].subject, "u1");
    assert!(failure[0].message.contains("unit_ops.flare"));
}

#[test]
fn a_recycle_needs_the_projection_flag() {
    let palette = common::palette();
    let flowsheet = common::demo();
    let options = ExportOptions {
        allow_lossy: true,
        ..ExportOptions::default()
    };
    let failure = export_process(&flowsheet, &palette, &options).expect_err("the tear refuses");
    assert!(
        failure
            .iter()
            .any(|diagnostic| diagnostic.code == "recycle_not_projected"),
        "{failure:?}"
    );
}

#[test]
fn the_export_is_deterministic() {
    let first = export_demo();
    let second = export_demo();
    assert_eq!(first.xml, second.xml);
}
