//! The Plant export against NeqSim's pinned Plant writer.
//!
//! The reference is `validation/neqsim/captures/dexpi_plant_simple.tsv`: the same process,
//! built and written by NeqSim, read back and printed. What must agree is the tagged items -
//! tag and type - the nozzle count, and the internal edges between items. NeqSim's default
//! writer draws no boundary segments, so the boundary connections this port adds are not part
//! of the comparison.

mod common;

use std::collections::{BTreeMap, BTreeSet};

use azoth_dexpi::plant::export::{PlantExportOptions, export_plant};
use azoth_dexpi::plant::import::{PlantImportOptions, import_plant};

fn capture() -> BTreeMap<String, String> {
    let path = common::workspace_root()
        .join("validation")
        .join("neqsim")
        .join("captures")
        .join("dexpi_plant_simple.tsv");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|error| panic!("{}: {error}", path.display()));
    text.lines()
        .filter(|line| !line.starts_with('#') && !line.trim().is_empty())
        .map(|line| match line.split_once('=') {
            Some((key, value)) => (key.trim().to_string(), value.trim().to_string()),
            None => (line.trim().to_string(), String::new()),
        })
        .collect()
}

#[test]
fn the_same_process_exports_the_same_items_and_internal_edges() {
    let expected = capture();
    let flowsheet = common::plant_simple();
    let report = export_plant(
        &flowsheet,
        &common::palette(),
        &PlantExportOptions {
            allow_lossy: true,
            azoth_extensions: true,
            ..PlantExportOptions::default()
        },
    )
    .expect("the fixture exports");

    let document = azoth_dexpi::reader::read(&report.xml)
        .expect("the export reads")
        .document;
    let items: BTreeMap<String, String> = document
        .model
        .children
        .iter()
        .filter_map(|child| match child {
            azoth_dexpi::Child::Object(root) => Some(root),
            _ => None,
        })
        .flat_map(plant_items)
        .collect();
    for (key, value) in &expected {
        if key.starts_with("item.") && key.ends_with(".type") {
            let tag = key
                .trim_start_matches("item.")
                .trim_end_matches(".type")
                .to_string();
            assert_eq!(items.get(&tag), Some(value), "{key}");
        }
    }
    assert_eq!(
        items.len().to_string(),
        expected["item.count"],
        "item count"
    );

    // The internal edges: both endpoints are instance endpoints, so a boundary stream is not
    // one of them.
    let imported = import_plant(
        &report.xml,
        &common::palette(),
        &PlantImportOptions::default(),
    )
    .expect("the export imports");
    let edges: BTreeSet<(String, String)> = imported
        .flowsheet
        .connections
        .iter()
        .filter(|connection| connection.from.contains('.') && connection.to.contains('.'))
        .map(|connection| {
            (
                connection
                    .from
                    .split('.')
                    .next()
                    .unwrap_or_default()
                    .to_string(),
                connection
                    .to
                    .split('.')
                    .next()
                    .unwrap_or_default()
                    .to_string(),
            )
        })
        .collect();
    let expected_edges: BTreeSet<(String, String)> = expected
        .keys()
        .filter_map(|key| {
            key.strip_prefix("segment.")
                .and_then(|rest| rest.split_once("->"))
                .map(|(source, target)| (source.to_string(), target.to_string()))
        })
        .collect();
    assert_eq!(edges, expected_edges);
}

fn plant_items(root: &azoth_dexpi::Object) -> Vec<(String, String)> {
    let Some(conceptual) = root.components("ConceptualModel") else {
        return Vec::new();
    };
    let mut items = Vec::new();
    for plant in conceptual.children.iter().filter_map(|child| match child {
        azoth_dexpi::ComponentChild::Object(object) => Some(object),
        _ => None,
    }) {
        let Some(tagged) = plant.components("TaggedPlantItems") else {
            continue;
        };
        for item in tagged.children.iter().filter_map(|child| match child {
            azoth_dexpi::ComponentChild::Object(object) => Some(object),
            _ => None,
        }) {
            if let Some(tag) = item.data_string("TagName") {
                items.push((tag.to_string(), item.ty.clone()));
            }
        }
    }
    items
}
