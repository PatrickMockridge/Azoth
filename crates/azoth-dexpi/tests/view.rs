//! The view: what a drawing would become, and what it would not.

mod common;

use azoth_dexpi::Profile;
use azoth_dexpi::drawing::DrawingOrigin;
use azoth_dexpi::plant::export::{PlantExportOptions, export_plant};
use azoth_dexpi::process::export::{ExportOptions, export_process};
use azoth_dexpi::view::{CoverageStatus, ExchangeFormat, inspect};

fn pfd() -> String {
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
    .expect("the demo exports")
    .xml
}

fn pid() -> String {
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
    .expect("the plant fixture exports")
    .xml
}

fn fixture(name: &str) -> String {
    std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join(name),
    )
    .expect("the fixture reads")
}

#[test]
fn a_pfd_view_marks_what_the_import_models() {
    let view = inspect(&pfd(), &common::palette()).expect("the PFD is viewable");
    assert_eq!(view.format, ExchangeFormat::Dexpi20);
    assert_eq!(view.profile, Some(Profile::Process));
    assert!(view.drawing.is_some());
    let scaffold = view.scaffold.as_ref().expect("the PFD imports");
    assert_eq!(scaffold.instances, 4);
    assert_eq!(scaffold.layout_positions, 5);
    assert!(scaffold.check_errors > 0, "the feeds are still undeclared");
    assert!(view.counts.modelled >= 3, "{:?}", view.counts);
    assert_eq!(view.counts.lossy, 1, "the heater is lossy");
    assert_eq!(view.counts.unmapped, 0);
    assert!(view.counts.topology > 0);
    assert!(view.counts.drawn > 0);
    let heater = view
        .objects
        .iter()
        .find(|object| object.label == "hx1")
        .expect("the heater is in the report");
    assert_eq!(heater.status, CoverageStatus::Lossy);
    assert_eq!(heater.unit.as_deref(), Some("unit_ops.heater"));
}

#[test]
fn a_pid_view_marks_control_and_loss() {
    let view = inspect(&pid(), &common::palette()).expect("the P&ID is viewable");
    assert_eq!(view.profile, Some(Profile::Plant));
    let scaffold = view.scaffold.as_ref().expect("the P&ID imports");
    assert_eq!(scaffold.instruments, 2);
    assert_eq!(scaffold.signals, 1);
    assert_eq!(scaffold.layout_positions, 5);
    assert_eq!(view.counts.control, 3, "{:?}", view.counts);
    assert!(
        view.counts.lossy >= 1,
        "the cooler is lossy: {:?}",
        view.counts
    );
    assert!(view.counts.modelled >= 2);
    let switch = view
        .objects
        .iter()
        .find(|object| object.label == "PSHH-101")
        .expect("the switch is in the report");
    assert_eq!(switch.status, CoverageStatus::Control);
}

#[test]
fn a_proteus_view_reads_the_graphics() {
    let view = inspect(&fixture("proteus_simple.xml"), &common::palette())
        .expect("the Proteus file is viewable");
    assert_eq!(view.format, ExchangeFormat::Proteus);
    assert_eq!(view.profile, Some(Profile::Plant));
    assert!(view.drawing.is_some());
    assert!(view.counts.modelled >= 3, "{:?}", view.counts);
    assert!(view.counts.lossy >= 1, "{:?}", view.counts);
    assert!(view.counts.topology >= 5);
    assert!(view.counts.drawing_only >= 1, "{:?}", view.counts);
    let scaffold = view.scaffold.as_ref().expect("the file imports");
    assert!(scaffold.layout_positions > 0);
}

#[test]
fn a_file_without_graphics_gets_a_derived_schematic() {
    let xml = export_process(
        &common::demo(),
        &common::palette(),
        &ExportOptions {
            allow_lossy: true,
            project_recycles: true,
            azoth_extensions: true,
            ..ExportOptions::default()
        },
    )
    .expect("the demo exports without a drawing")
    .xml;
    let view = inspect(&xml, &common::palette()).expect("the file is viewable");
    let drawing = view.drawing.as_ref().expect("a schematic is derived");
    assert_eq!(drawing.origin, DrawingOrigin::Derived);
    assert!(view.counts.groups > 0);
}

#[test]
fn a_refused_model_still_shows_its_drawing() {
    let xml = fixture("native_projection_shape.dexpi.xml")
        .replace("id=\"Equipment1\"", "id=\"Equipment2\"");
    let view = inspect(&xml, &common::palette()).expect("the drawing is still viewable");
    assert!(view.drawing.is_some());
    assert!(view.scaffold.is_none());
    assert!(!view.diagnostics.is_empty());
}
