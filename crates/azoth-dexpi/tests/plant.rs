//! The Plant profile: tagged items, nozzles, lines, and the topology an import recovers.

mod common;

use azoth_dexpi::plant::export::{PlantExportOptions, export_plant};
use azoth_dexpi::plant::import::{PlantImportOptions, import_plant};
use azoth_dexpi::reader::read;
use azoth_dexpi::validate::validate;
use azoth_process::{parse_flowsheet, validate as validate_flowsheet};

fn export_oracle(extensions: bool) -> azoth_dexpi::PlantExportReport {
    let palette = common::palette();
    let flowsheet = common::plant_simple();
    let options = PlantExportOptions {
        allow_lossy: true,
        azoth_extensions: extensions,
        ..PlantExportOptions::default()
    };
    export_plant(&flowsheet, &palette, &options).expect("the oracle fixture exports as a Plant")
}

#[test]
fn the_oracle_fixture_exports_and_passes_the_official_schema() {
    let report = export_oracle(true);
    assert!(report.semantic_valid);
    assert_eq!(report.items, 3);
    assert_eq!(report.segments, 5);
    let path = common::temp_path("plant.dexpi.xml");
    std::fs::write(&path, &report.xml).expect("the fixture writes");
    common::xmllint(&path);
}

#[test]
fn the_export_names_the_items_and_their_nozzles() {
    let report = export_oracle(true);
    let inventory =
        azoth_dexpi::inventory::inventory(&read(&report.xml).expect("the export reads").document);
    let _ = inventory;
    assert!(report.xml.contains("<String>V-101</String>"));
    assert!(report.xml.contains("Plant/ProcessEquipment.Separator"));
    assert!(report.xml.contains("<String>vapour</String>"));
    assert!(report.xml.contains("<String>FEED</String>"));
}

#[test]
fn the_round_trip_recovers_the_topology() {
    let palette = common::palette();
    let original = common::plant_simple();
    let xml = export_oracle(true).xml;
    let imported = import_plant(&xml, &palette, &PlantImportOptions::default())
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
}

#[test]
fn line_numbers_group_the_systems_and_survive_the_round_trip() {
    let text = r#"
id = "flowsheets.lines"
name = "Two lines"

products = ["GAS"]

[[inputs]]
name = "FEED"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e6
T = 300.0

[[instances]]
id = "P-101"
unit = "unit_ops.pump"
[instances.parameters]
outlet_pressure = 6.0e6
isentropic_efficiency = 0.75

[[instances]]
id = "K-101"
unit = "unit_ops.compressor"
[instances.parameters]
outlet_pressure = 8.0e6
isentropic_efficiency = 0.75

[[connections]]
from = "FEED"
to = "P-101.inlet"
line = "L-1001"

[[connections]]
from = "P-101.outlet"
to = "K-101.inlet"
line = "L-1002"

[[connections]]
from = "K-101.outlet"
to = "GAS"
line = "L-1003"
"#;
    let flowsheet = parse_flowsheet(text).expect("the fixture parses");
    let palette = common::palette();
    let report = export_plant(
        &flowsheet,
        &palette,
        &PlantExportOptions {
            allow_lossy: true,
            azoth_extensions: true,
            ..PlantExportOptions::default()
        },
    )
    .expect("the fixture exports");
    let document = read(&report.xml).expect("the export reads").document;
    let model = document
        .model
        .children
        .iter()
        .find_map(|child| match child {
            azoth_dexpi::Child::Object(object) => Some(object),
            _ => None,
        })
        .expect("the engineering model");
    let conceptual = model
        .components("ConceptualModel")
        .expect("conceptual model");
    let plant = conceptual
        .children
        .iter()
        .find_map(|child| match child {
            azoth_dexpi::ComponentChild::Object(object) => Some(object),
            _ => None,
        })
        .expect("the plant model");
    let systems = plant.components("PipingNetworkSystems").expect("systems");
    assert_eq!(systems.children.len(), 3, "one system per line");

    let imported = import_plant(&report.xml, &palette, &PlantImportOptions::default())
        .expect("the exchange imports");
    let lines: Vec<Option<String>> = imported
        .flowsheet
        .connections
        .iter()
        .map(|connection| connection.line.clone())
        .collect();
    for line in ["L-1001", "L-1002", "L-1003"] {
        assert!(
            lines.contains(&Some(line.to_string())),
            "{line} survived: {lines:?}"
        );
    }
}

#[test]
fn a_mixer_is_refused_because_it_is_a_junction() {
    let text = r#"
id = "flowsheets.mix"
name = "A mixer"

products = ["OUT"]

[[inputs]]
name = "A"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e5
T = 300.0

[[inputs]]
name = "B"
components = ["methane"]
n = 1.0
z = [1.0]
P = 5.0e5
T = 300.0

[[instances]]
id = "M-1"
unit = "unit_ops.mixer"

[[connections]]
from = "A"
to = "M-1.feed"

[[connections]]
from = "B"
to = "M-1.feed"

[[connections]]
from = "M-1.product"
to = "OUT"
"#;
    let flowsheet = parse_flowsheet(text).expect("the fixture parses");
    assert!(validate_flowsheet(&flowsheet, &common::palette()).is_empty());
    let failure = export_plant(
        &flowsheet,
        &common::palette(),
        &PlantExportOptions::default(),
    )
    .expect_err("a mixer refuses");
    assert_eq!(failure[0].code, "refused_mapping");
    assert_eq!(failure[0].subject, "M-1");
}

#[test]
fn an_unmarked_loop_is_inferred_and_reported() {
    let palette = common::palette();
    let xml = r##"<?xml version="1.0" encoding="UTF-8"?>
<Model name="looped" uri="urn:test:looped">
    <Import prefix="Core" source="https://data.dexpi.org/models/2.0.0/Core.xml"/>
    <Import prefix="Plant" source="https://data.dexpi.org/models/2.0.0/Plant.xml"/>
    <Object type="Core/EngineeringModel">
        <Components property="ConceptualModel">
            <Object id="PlantModel1" type="Plant/PlantModel">
                <Components property="TaggedPlantItems">
                    <Object id="P1" type="Plant/ProcessEquipment.CentrifugalPump">
                        <Data property="TagName"><String>P-100</String></Data>
                        <Components property="Nozzles">
                            <Object id="P1In" type="Plant/ProcessEquipment.Nozzle">
                                <Components property="Nodes"><Object id="N1" type="Plant/Piping.PipingNode"/></Components>
                                <Data property="SubTagName"><String>inlet</String></Data>
                            </Object>
                            <Object id="P1Out" type="Plant/ProcessEquipment.Nozzle">
                                <Components property="Nodes"><Object id="N2" type="Plant/Piping.PipingNode"/></Components>
                                <Data property="SubTagName"><String>outlet</String></Data>
                            </Object>
                        </Components>
                    </Object>
                    <Object id="K1" type="Plant/ProcessEquipment.CentrifugalCompressor">
                        <Data property="TagName"><String>K-100</String></Data>
                        <Components property="Nozzles">
                            <Object id="K1In" type="Plant/ProcessEquipment.Nozzle">
                                <Components property="Nodes"><Object id="N3" type="Plant/Piping.PipingNode"/></Components>
                                <Data property="SubTagName"><String>inlet</String></Data>
                            </Object>
                            <Object id="K1Out" type="Plant/ProcessEquipment.Nozzle">
                                <Components property="Nodes"><Object id="N4" type="Plant/Piping.PipingNode"/></Components>
                                <Data property="SubTagName"><String>outlet</String></Data>
                            </Object>
                        </Components>
                    </Object>
                </Components>
                <Components property="PipingNetworkSystems">
                    <Object id="S1" type="Plant/Piping.PipingNetworkSystem">
                        <Data property="LineNumber"><String>L-1</String></Data>
                        <Components property="Segments">
                            <Object id="Seg1" type="Plant/Piping.PipingNetworkSegment">
                                <References property="SourceNode" objects="#N2"/>
                                <References property="TargetNode" objects="#N3"/>
                                <Data property="SegmentNumber"><String>1</String></Data>
                            </Object>
                            <Object id="Seg2" type="Plant/Piping.PipingNetworkSegment">
                                <References property="SourceNode" objects="#N4"/>
                                <References property="TargetNode" objects="#N1"/>
                                <Data property="SegmentNumber"><String>2</String></Data>
                            </Object>
                        </Components>
                    </Object>
                </Components>
            </Object>
        </Components>
    </Object>
</Model>
"##;
    let document = read(xml).expect("the document reads").document;
    assert!(
        validate(&document).is_valid(),
        "{:?}",
        validate(&document).diagnostics
    );
    let imported = import_plant(xml, &palette, &PlantImportOptions::default())
        .expect("the looped drawing imports");
    assert_eq!(imported.flowsheet.recycles.len(), 1);
    assert!(
        imported
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code == "recycle_inferred"),
        "{:?}",
        imported.diagnostics
    );
}
