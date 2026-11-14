//! The export against NeqSim's pinned Process writer.
//!
//! NeqSim has no native DEXPI 2.0 Process **reader**, so this is the half that has an oracle:
//! the same process, built by NeqSim's `Dexpi20ProcessProbe`, exported by its reviewed writer,
//! and exported by this crate. What must agree is what the documents say - the steps and their
//! types, the edges between them, and the three quantities - not their bytes: NeqSim's DOM
//! transformer and this writer indent differently and mint their own ids.
//!
//! The capture is the committed file `validation/neqsim/captures/dexpi_process_simple.tsv`, the
//! same shape `tools/oracle_sweep.py` refreshes for every other probe.

mod common;

use std::collections::BTreeMap;

use azoth_dexpi::export_process;
use azoth_dexpi::process::export::ExportOptions;

fn capture() -> BTreeMap<String, String> {
    let path = common::workspace_root()
        .join("validation")
        .join("neqsim")
        .join("captures")
        .join("dexpi_process_simple.tsv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .filter_map(|line| line.split_once('='))
        .map(|(key, value)| (key.trim().to_string(), value.trim().to_string()))
        .collect()
}

struct Edge {
    source: String,
    target: String,
    quantity: BTreeMap<String, f64>,
    unit: BTreeMap<String, String>,
}

fn edges(inventory: &BTreeMap<String, String>) -> Vec<Edge> {
    let mut edges: BTreeMap<String, Edge> = BTreeMap::new();
    for (key, value) in inventory {
        let Some(rest) = key.strip_prefix("stream.") else {
            continue;
        };
        let Some((left, field)) = rest.rsplit_once('.') else {
            continue;
        };
        // `<label>.<property>` and `<label>.<property>.unit` are the two shapes.
        let (label, property) = if field == "unit" {
            match left.rsplit_once('.') {
                Some((label, property)) => (label, Some(property)),
                None => continue,
            }
        } else {
            (left, None)
        };
        let edge = edges.entry(label.to_string()).or_insert_with(|| Edge {
            source: String::new(),
            target: String::new(),
            quantity: BTreeMap::new(),
            unit: BTreeMap::new(),
        });
        match (field, property) {
            ("source", _) => edge.source = value.clone(),
            ("target", _) => edge.target = value.clone(),
            ("unit", Some(property)) => {
                edge.unit.insert(property.to_string(), value.clone());
            }
            (other, None) => {
                if let Ok(number) = value.parse::<f64>() {
                    edge.quantity.insert(other.to_string(), number);
                }
            }
            _ => {}
        }
    }
    edges.into_values().collect()
}

#[test]
fn the_same_process_exports_the_same_steps_and_streams() {
    let expected = capture();
    let flowsheet = common::oracle_simple();
    let report = export_process(
        &flowsheet,
        &common::palette(),
        &ExportOptions {
            azoth_extensions: true,
            ..ExportOptions::default()
        },
    )
    .expect("the oracle fixture exports");
    let actual = azoth_dexpi::inventory::inventory(
        &azoth_dexpi::reader::read(&report.xml)
            .expect("the export reads")
            .document,
    );

    assert_eq!(actual.get("step.count"), expected.get("step.count"));
    assert_eq!(actual.get("stream.count"), expected.get("stream.count"));

    for (key, value) in &expected {
        if key.starts_with("step.") && key.ends_with(".type") {
            assert_eq!(actual.get(key), Some(value), "{key}");
        }
    }

    let expected_edges: Vec<Edge> = edges(&expected);
    let actual_edges: Vec<Edge> = edges(&actual);
    assert_eq!(actual_edges.len(), expected_edges.len());

    for expected in &expected_edges {
        let actual_edge = actual_edges
            .iter()
            .find(|edge| edge.source == expected.source && edge.target == expected.target)
            .unwrap_or_else(|| {
                panic!(
                    "no exported stream {} -> {}; exported {:?}",
                    expected.source,
                    expected.target,
                    actual_edges
                        .iter()
                        .map(|edge| format!("{} -> {}", edge.source, edge.target))
                        .collect::<Vec<_>>()
                )
            });
        for (property, unit) in &expected.unit {
            assert_eq!(
                actual_edge.unit.get(property),
                Some(unit),
                "{} -> {} {property} unit",
                expected.source,
                expected.target
            );
        }
        for (property, value) in &expected.quantity {
            let actual_value = actual_edge.quantity.get(property).unwrap_or_else(|| {
                panic!(
                    "{} -> {} carries no {property}",
                    expected.source, expected.target
                )
            });
            // **Only the values the exchange itself controls are compared, and that is a
            // finding, not a convenience.** The two runs are two property models: the
            // separator's split on this feed is 723.5/276.5 kg/h here against NeqSim's
            // 715.2/284.8, a 1-3% difference that belongs to the flash. The separator's own
            // oracle is `process_separator.tsv`; a DEXPI probe is not where that disagreement
            // should be either hidden under a wide tolerance or mistaken for an export defect.
            //
            // What the exchange controls is compared exactly: every quantity and unit is
            // present (above), every pressure is the equipment's set point, the feed's mass
            // flow is its stated composition, and the cooler's outlet temperature is the one
            // the flowsheet states. Everything else is the run's, and `inventory.rs` proves it
            // was carried without proving the physics.
            let agrees = match property.as_str() {
                "Pressure" => (actual_value - value).abs() <= 1e-9 * value.abs().max(1.0),
                "MassFlow" if expected.source == "FEED" => {
                    (actual_value - value).abs() <= 1e-4 * value.abs().max(1.0)
                }
                "Temperature" if expected.source == "E-101" => (actual_value - value).abs() <= 1e-9,
                _ => true,
            };
            assert!(
                agrees,
                "{} -> {} {property}: {actual_value} against {value}",
                expected.source, expected.target
            );
        }
    }
}

#[test]
fn the_neqsim_fixture_passes_this_crates_semantic_rules() {
    let path = common::workspace_root()
        .join("validation")
        .join("dexpi")
        .join("neqsim")
        .join("process_simple.dexpi.xml");
    let xml = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    let document = azoth_dexpi::reader::read(&xml)
        .expect("this reader reads NeqSim's export")
        .document;
    let report = azoth_dexpi::validate::validate(&document);
    assert!(report.is_valid(), "{:?}", report.diagnostics);
}
