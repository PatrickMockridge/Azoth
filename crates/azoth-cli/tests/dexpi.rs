//! `azoth dexpi`, end to end: the three verbs against the shipped flowsheet.
//!
//! The round trip is the one assertion that spans them: export the demo, validate the file,
//! import it back, and read the topology that comes out. What is *not* asserted here is the
//! run's numbers - those are `azoth-dexpi`'s tests - but that the CLI writes the file, reads
//! it, and exits with the status a script switches on.

use std::path::{Path, PathBuf};
use std::process::Command;

use azoth_process::parse_flowsheet;

fn root() -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR")).join("../..")
}

fn demo() -> String {
    root()
        .join("specs/flowsheets/demo.toml")
        .display()
        .to_string()
}

fn palette() -> String {
    root().join("specs/unit_ops").display().to_string()
}

fn temp(name: &str) -> String {
    std::env::temp_dir()
        .join(format!("azoth-cli-dexpi-{}-{name}", std::process::id()))
        .display()
        .to_string()
}

fn azoth(args: &[&str]) -> std::process::Output {
    Command::new(env!("CARGO_BIN_EXE_azoth"))
        .args(args)
        .output()
        .expect("the built binary should run")
}

#[test]
fn the_demo_exports_validates_and_imports_again() {
    let xml = temp("demo.dexpi.xml");
    let toml = temp("demo.imported.toml");

    let export = azoth(&[
        "dexpi",
        "export",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--out",
        &xml,
        "--allow-lossy-mappings",
        "--project-recycles",
        "--with-azoth-extensions",
    ]);
    assert!(
        export.status.success(),
        "export failed: {}",
        String::from_utf8_lossy(&export.stderr)
    );
    assert!(Path::new(&xml).is_file(), "the exchange was written");

    let validate = azoth(&["dexpi", "validate", "--input", &xml, "--xsd"]);
    assert!(
        validate.status.success(),
        "validation failed: {}",
        String::from_utf8_lossy(&validate.stdout)
    );

    let import = azoth(&[
        "dexpi",
        "import",
        "--input",
        &xml,
        "--palette",
        &palette(),
        "--out",
        &toml,
    ]);
    assert!(
        import.status.success(),
        "import failed: {}",
        String::from_utf8_lossy(&import.stdout)
    );
    let text = std::fs::read_to_string(&toml).expect("the imported flowsheet is written");
    let flowsheet = parse_flowsheet(&text).expect("the imported document parses");
    let ids: Vec<&str> = flowsheet
        .instances
        .iter()
        .map(|instance| instance.id.as_str())
        .collect();
    assert_eq!(ids, ["mix1", "p1", "hx1", "sep1"]);
    assert_eq!(flowsheet.recycles.len(), 1);
    assert_eq!(flowsheet.recycles[0].stream, "recycle_1");
}

#[test]
fn the_drawing_flag_writes_a_pfd_and_the_import_places_it() {
    let xml = temp("demo.drawing.xml");
    let placed = temp("demo.placed.toml");
    let unplaced = temp("demo.unplaced.toml");

    let export = azoth(&[
        "dexpi",
        "export",
        "--flowsheet",
        &demo(),
        "--palette",
        &palette(),
        "--out",
        &xml,
        "--allow-lossy-mappings",
        "--project-recycles",
        "--with-azoth-extensions",
        "--with-drawing",
    ]);
    assert!(
        export.status.success(),
        "export failed: {}",
        String::from_utf8_lossy(&export.stderr)
    );
    let text = std::fs::read_to_string(&xml).expect("the exchange is written");
    assert!(
        text.contains("Core/Diagram.Diagram"),
        "the PFD carries a drawing"
    );

    let import = azoth(&[
        "dexpi",
        "import",
        "--input",
        &xml,
        "--palette",
        &palette(),
        "--out",
        &placed,
    ]);
    assert!(import.status.success());
    let text = std::fs::read_to_string(&placed).expect("the import is written");
    assert!(text.contains("[layout.instances]"), "{text}");

    let import = azoth(&[
        "dexpi",
        "import",
        "--input",
        &xml,
        "--palette",
        &palette(),
        "--out",
        &unplaced,
        "--no-layout",
    ]);
    assert!(import.status.success());
    let text = std::fs::read_to_string(&unplaced).expect("the import is written");
    assert!(!text.contains("[layout.instances]"), "{text}");
}

#[test]
fn the_plant_profile_exports_validates_and_imports_again() {
    let fixture = root()
        .join("crates/azoth-dexpi/tests/fixtures/plant_simple.toml")
        .display()
        .to_string();
    let xml = temp("plant.dexpi.xml");
    let toml = temp("plant.imported.toml");

    let export = azoth(&[
        "dexpi",
        "export",
        "--profile",
        "plant",
        "--flowsheet",
        &fixture,
        "--palette",
        &palette(),
        "--out",
        &xml,
        "--allow-lossy-mappings",
        "--with-azoth-extensions",
    ]);
    assert!(
        export.status.success(),
        "export failed: {}",
        String::from_utf8_lossy(&export.stderr)
    );

    let validate = azoth(&["dexpi", "validate", "--input", &xml, "--xsd"]);
    assert!(
        validate.status.success(),
        "validation failed: {}",
        String::from_utf8_lossy(&validate.stdout)
    );

    let import = azoth(&[
        "dexpi",
        "import",
        "--input",
        &xml,
        "--palette",
        &palette(),
        "--out",
        &toml,
    ]);
    assert!(
        import.status.success(),
        "import failed: {}",
        String::from_utf8_lossy(&import.stdout)
    );
    let text = std::fs::read_to_string(&toml).expect("the imported flowsheet is written");
    let flowsheet = parse_flowsheet(&text).expect("the imported document parses");
    assert_eq!(flowsheet.instances.len(), 3);
    assert_eq!(flowsheet.connections.len(), 5);
    assert_eq!(flowsheet.instruments.len(), 2);
    assert_eq!(flowsheet.signals.len(), 1);
    assert_eq!(flowsheet.connections[0].line.as_deref(), Some("L-1001"));
}

#[test]
fn a_valve_exports_inside_a_segment_and_imports_again() {
    let fixture = root()
        .join("crates/azoth-dexpi/tests/fixtures/plant_valve.toml")
        .display()
        .to_string();
    let xml = temp("valve.dexpi.xml");
    let toml = temp("valve.imported.toml");
    let export = azoth(&[
        "dexpi",
        "export",
        "--profile",
        "plant",
        "--flowsheet",
        &fixture,
        "--palette",
        &palette(),
        "--out",
        &xml,
        "--with-azoth-extensions",
        "--with-drawing",
    ]);
    assert!(
        export.status.success(),
        "export failed: {}",
        String::from_utf8_lossy(&export.stderr)
    );
    let text = std::fs::read_to_string(&xml).expect("the P&ID is written");
    assert!(text.contains("Plant/Piping.GlobeValve"), "{text}");

    let import = azoth(&[
        "dexpi",
        "import",
        "--input",
        &xml,
        "--palette",
        &palette(),
        "--out",
        &toml,
    ]);
    assert!(import.status.success());
    let text = std::fs::read_to_string(&toml).expect("the import is written");
    let flowsheet = parse_flowsheet(&text).expect("the scaffold parses");
    assert!(
        flowsheet
            .instances
            .iter()
            .any(|instance| instance.id == "V-101" && instance.unit == "unit_ops.throttling_valve"),
        "{:?}",
        flowsheet.instances
    );
}

#[test]
fn a_proteus_plant_model_imports_from_the_command_line() {
    let fixture = root()
        .join("crates/azoth-dexpi/tests/fixtures/proteus_simple.xml")
        .display()
        .to_string();
    let toml = temp("proteus.imported.toml");
    let import = azoth(&[
        "dexpi",
        "import",
        "--input",
        &fixture,
        "--palette",
        &palette(),
        "--out",
        &toml,
    ]);
    assert!(
        import.status.success(),
        "import failed: {}",
        String::from_utf8_lossy(&import.stdout)
    );
    let text = std::fs::read_to_string(&toml).expect("the imported flowsheet is written");
    let flowsheet = parse_flowsheet(&text).expect("the imported document parses");
    let units: Vec<&str> = flowsheet
        .instances
        .iter()
        .map(|instance| instance.unit.as_str())
        .collect();
    assert!(units.contains(&"unit_ops.pump"));
    assert!(units.contains(&"unit_ops.tank"));
}

#[test]
fn the_view_reports_coverage_layout_and_the_scaffold() {
    let fixture = root()
        .join("crates/azoth-dexpi/tests/fixtures/proteus_simple.xml")
        .display()
        .to_string();
    let view = azoth(&[
        "dexpi",
        "view",
        "--input",
        &fixture,
        "--palette",
        &palette(),
        "--json",
    ]);
    assert!(
        view.status.success(),
        "view failed: {}",
        String::from_utf8_lossy(&view.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&view.stdout).expect("the view is JSON");
    assert_eq!(report["format"], serde_json::json!("proteus"));
    assert_eq!(report["profile"], serde_json::json!("plant"));
    assert!(
        report["drawing"]["groups"]
            .as_array()
            .is_some_and(|g| !g.is_empty())
    );
    assert!(report["counts"]["modelled"].as_u64().is_some_and(|n| n > 0));
    assert!(
        report["scaffold"]["layout_positions"]
            .as_u64()
            .is_some_and(|n| n > 0)
    );
}

#[test]
fn a_refused_export_exits_non_zero_and_prints_the_refusal() {
    let fixture = temp("flare.toml");
    std::fs::write(
        &fixture,
        r#"
id = "flowsheets.flare"
name = "One flare"

products = ["vent"]

[[inputs]]
name = "feed_1"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e5
T = 300.0

[[instances]]
id = "f1"
unit = "unit_ops.flare"

[[connections]]
from = "feed_1"
to = "f1.inlet"

[[connections]]
from = "f1.product"
to = "vent"
"#,
    )
    .expect("the fixture writes");

    let out = temp("flare.dexpi.xml");
    let export = azoth(&[
        "dexpi",
        "export",
        "--flowsheet",
        &fixture,
        "--palette",
        &palette(),
        "--out",
        &out,
    ]);
    assert_eq!(
        export.status.code(),
        Some(2),
        "a refusal is not a success status"
    );
    let text = String::from_utf8_lossy(&export.stdout);
    assert!(text.contains("refused_mapping"), "{text}");
    assert!(
        !Path::new(&out).exists(),
        "nothing is written for a refused exchange"
    );
}
