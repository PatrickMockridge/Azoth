//! The drawing model: DEXPI 2.0 `Core/Diagram` read, write, and the warnings that keep a
//! real drawing viewable when a primitive or a shape is not one this port renders.

mod common;

use azoth_dexpi::Profile;
use azoth_dexpi::drawing::model::{
    Color, Drawing, DrawingGroup, DrawingOrigin, GroupKind, Point, Primitive, Stroke,
};
use azoth_dexpi::drawing::native::{read_document, write_diagram};
use azoth_dexpi::model::{
    Child, ComponentChild, Components, Data, DataValue, Document, Import, Model, Object,
};
use azoth_dexpi::reader::read;
use azoth_dexpi::validate::validate;

fn fixture() -> String {
    std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("tests")
            .join("fixtures")
            .join("native_projection_shape.dexpi.xml"),
    )
    .expect("the drawing fixture reads")
}

fn fixture_document() -> Document {
    read(&fixture()).expect("the fixture reads").document
}

#[test]
fn the_fixture_reads_its_groups_and_flattens_the_shape_usage() {
    let report = read_document(&fixture_document()).expect("the fixture carries a diagram");
    let drawing = report.drawing;
    assert_eq!(drawing.name, "Drawing for this test");
    assert_eq!(drawing.origin, DrawingOrigin::File);
    assert_eq!(drawing.groups.len(), 2, "{:?}", drawing.groups);
    assert!(report.diagnostics.is_empty(), "{:?}", report.diagnostics);

    let first = &drawing.groups[0];
    assert_eq!(first.represents.as_deref(), Some("Equipment1"));
    assert_eq!(first.kind, GroupKind::Node);
    assert_eq!(first.primitives.len(), 3);
    assert!(
        first
            .primitives
            .iter()
            .any(|primitive| matches!(primitive, Primitive::Text { text, .. } if text == "V-101"))
    );

    let second = &drawing.groups[1];
    assert_eq!(second.represents.as_deref(), Some("Equipment2"));
    let polygon = second
        .primitives
        .iter()
        .find_map(|primitive| match primitive {
            Primitive::Polygon { points, .. } => Some(points.clone()),
            _ => None,
        })
        .expect("the shape usage became a polygon");
    assert_eq!(polygon.len(), 3);
    assert_eq!(polygon[1], Point { x: 96.0, y: 20.0 });
}

#[test]
fn the_fixture_passes_the_official_schema_and_the_semantic_rules() {
    let document = fixture_document();
    let report = validate(&document);
    assert!(report.is_valid(), "{:?}", report.diagnostics);
    let path = common::temp_path("native_projection_shape.xml");
    std::fs::write(&path, fixture()).expect("the fixture writes");
    common::xmllint(&path);
}

#[test]
fn write_then_read_round_trips_the_drawing() {
    let drawing = Drawing::from_groups(
        "round trip".to_string(),
        Some("mm".to_string()),
        Color::WHITE,
        DrawingOrigin::Derived,
        vec![DrawingGroup {
            kind: GroupKind::Node,
            id: None,
            represents: Some("Equipment1".to_string()),
            primitives: vec![
                Primitive::Polygon {
                    id: Some("Primitive1".to_string()),
                    points: vec![
                        Point { x: 10.0, y: 10.0 },
                        Point { x: 30.0, y: 10.0 },
                        Point { x: 30.0, y: 20.0 },
                        Point { x: 10.0, y: 20.0 },
                    ],
                    fill: "Solid".to_string(),
                    stroke: Stroke {
                        color: Color { r: 255, g: 0, b: 0 },
                        width: 0.5,
                        dash: "Solid".to_string(),
                    },
                },
                Primitive::Polyline {
                    id: Some("Primitive2".to_string()),
                    points: vec![Point { x: 30.0, y: 15.0 }, Point { x: 40.0, y: 15.0 }],
                    stroke: Stroke {
                        color: Color::BLACK,
                        width: 1.0,
                        dash: "Dash".to_string(),
                    },
                },
                Primitive::Text {
                    id: Some("Primitive3".to_string()),
                    position: Point { x: 20.0, y: 5.0 },
                    text: "V-101".to_string(),
                    size: 2.5,
                    rotation: 0.0,
                    alignment: "CenterCenter".to_string(),
                    color: Color::BLACK,
                },
            ],
        }],
    );
    let mut mint = azoth_dexpi::ids::Mint::new();
    let diagram = write_diagram(&drawing, "Equipment1", &mut mint);
    let document = document_with(diagram);
    let xml = azoth_dexpi::writer::write(&document).expect("the round trip writes");
    let path = common::temp_path("drawing_round_trip.xml");
    std::fs::write(&path, &xml).expect("the round trip writes");
    common::xmllint(&path);

    let read_back = read(&xml).expect("the round trip reads").document;
    assert!(validate(&read_back).is_valid());
    let report = read_document(&read_back).expect("the round trip carries a diagram");
    let recovered = report.drawing;
    assert_eq!(recovered.groups.len(), 1);
    assert_eq!(
        recovered.groups[0].represents.as_deref(),
        Some("Equipment1")
    );
    assert_eq!(recovered.groups[0].primitives.len(), 3);
    assert_eq!(recovered.bounds, drawing.bounds);
    assert_eq!(recovered.groups[0].primitives, drawing.groups[0].primitives);
}

#[test]
fn a_group_with_no_represented_object_is_not_written() {
    let drawing = Drawing::from_groups(
        "unrepresented".to_string(),
        None,
        Color::WHITE,
        DrawingOrigin::Derived,
        vec![DrawingGroup {
            kind: GroupKind::Annotation,
            id: None,
            represents: None,
            primitives: vec![Primitive::Polyline {
                id: None,
                points: vec![Point { x: 0.0, y: 0.0 }, Point { x: 1.0, y: 0.0 }],
                stroke: Stroke::default(),
            }],
        }],
    );
    let mut mint = azoth_dexpi::ids::Mint::new();
    let diagram = write_diagram(&drawing, "Equipment1", &mut mint);
    let document = document_with(diagram);
    let xml = azoth_dexpi::writer::write(&document).expect("the document writes");
    let read_back = read(&xml).expect("the document reads").document;
    let report = read_document(&read_back).expect("the diagram reads");
    assert!(report.drawing.groups.is_empty());
}

#[test]
fn an_unsupported_primitive_is_a_warning_not_a_refusal() {
    let mut document = fixture_document();
    let diagram = find_type(&mut document, "Core/Diagram.Diagram").expect("the diagram");
    let groups = diagram
        .children
        .iter_mut()
        .find_map(|child| match child {
            Child::Components(components) if components.property == "Groups" => Some(components),
            _ => None,
        })
        .expect("the groups");
    let ComponentChild::Object(group) = &mut groups.children[0] else {
        panic!("the first group is an object");
    };
    let static_group = group
        .children
        .iter_mut()
        .find_map(|child| match child {
            Child::Components(components) if components.property == "Groups" => {
                components.children.first_mut()
            }
            _ => None,
        })
        .expect("the static group");
    let ComponentChild::Object(static_group) = static_group else {
        panic!("the static group is an object");
    };
    let elements = static_group
        .children
        .iter_mut()
        .find_map(|child| match child {
            Child::Components(components) if components.property == "Elements" => Some(components),
            _ => None,
        })
        .expect("the elements");
    elements.children.push(ComponentChild::Object(Object {
        id: Some("Mystery1".to_string()),
        name: None,
        ty: "Core/Diagram.Mystery".to_string(),
        children: Vec::new(),
    }));

    let report = read_document(&document).expect("the drawing still reads");
    assert_eq!(report.drawing.groups[0].primitives.len(), 3);
    assert!(
        report
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code == "DEXPI_GRAPHICS_UNSUPPORTED"),
        "{:?}",
        report.diagnostics
    );
    let validation = validate(&document);
    assert!(validation.is_valid(), "{:?}", validation.diagnostics);
}

#[test]
fn a_shape_usage_with_no_catalogue_entry_is_a_warning() {
    let mut document = fixture_document();
    let usage = find_type(&mut document, "Core/Diagram.ShapeUsage").expect("the usage");
    for child in &mut usage.children {
        if let Child::References(references) = child {
            if references.property == "Shape" {
                references.objects = vec!["#Missing".to_string()];
            }
        }
    }
    let report = read_document(&document).expect("the drawing still reads");
    assert!(
        report
            .diagnostics
            .iter()
            .any(|diagnostic| diagnostic.code == "DEXPI_GRAPHICS_UNSUPPORTED"),
        "{:?}",
        report.diagnostics
    );
}

#[test]
fn a_three_point_polyline_reads_three_points() {
    let mut document = fixture_document();
    let line = find_type(&mut document, "Core/Diagram.PolyLine").expect("the line");
    let points = line
        .children
        .iter_mut()
        .find_map(|child| match child {
            Child::Data(data) if data.property == "Points" => Some(data),
            _ => None,
        })
        .expect("the points");
    let first = points.values[0].clone();
    points.values.push(first);
    let report = read_document(&document).expect("the drawing still reads");
    let read = report.drawing.groups[0]
        .primitives
        .iter()
        .find_map(|primitive| match primitive {
            Primitive::Polyline { points, .. } => Some(points.len()),
            _ => None,
        })
        .expect("the polyline");
    assert_eq!(read, 3);
}

#[test]
fn the_writer_is_deterministic() {
    let drawing = Drawing::from_groups(
        "deterministic".to_string(),
        None,
        Color::WHITE,
        DrawingOrigin::Derived,
        vec![DrawingGroup {
            kind: GroupKind::Node,
            id: None,
            represents: Some("Equipment1".to_string()),
            primitives: vec![Primitive::Polyline {
                id: None,
                points: vec![Point { x: 0.0, y: 0.0 }, Point { x: 10.0, y: 10.0 }],
                stroke: Stroke::default(),
            }],
        }],
    );
    let first = {
        let mut mint = azoth_dexpi::ids::Mint::new();
        azoth_dexpi::writer::write(&document_with(write_diagram(
            &drawing,
            "Equipment1",
            &mut mint,
        )))
        .expect("the first write")
    };
    let second = {
        let mut mint = azoth_dexpi::ids::Mint::new();
        azoth_dexpi::writer::write(&document_with(write_diagram(
            &drawing,
            "Equipment1",
            &mut mint,
        )))
        .expect("the second write")
    };
    assert_eq!(first, second);
}

fn document_with(diagram: Object) -> Document {
    let equipment = Object {
        id: Some("Equipment1".to_string()),
        name: None,
        ty: "Plant/ProcessEquipment.Separator".to_string(),
        children: vec![Child::Data(Data {
            property: "TagName".to_string(),
            values: vec![DataValue::String("V-101".to_string())],
        })],
    };
    let plant = Object {
        id: Some("PlantModel1".to_string()),
        name: None,
        ty: "Plant/PlantModel".to_string(),
        children: vec![Child::Components(Components {
            property: "TaggedPlantItems".to_string(),
            children: vec![ComponentChild::Object(equipment)],
        })],
    };
    let engineering = Object {
        id: None,
        name: None,
        ty: "Core/EngineeringModel".to_string(),
        children: vec![
            Child::Components(Components {
                property: "ConceptualModel".to_string(),
                children: vec![ComponentChild::Object(plant)],
            }),
            Child::Components(Components {
                property: "Diagram".to_string(),
                children: vec![ComponentChild::Object(diagram)],
            }),
        ],
    };
    Document {
        model: Model {
            name: "drawing_round_trip".to_string(),
            uri: "urn:test:drawing_round_trip".to_string(),
            imports: Profile::Plant
                .imports()
                .into_iter()
                .map(|(prefix, source)| Import {
                    prefix: prefix.to_string(),
                    source: source.to_string(),
                })
                .collect(),
            children: vec![Child::Object(engineering)],
        },
    }
}

fn find_type<'a>(document: &'a mut Document, ty: &str) -> Option<&'a mut Object> {
    fn walk<'a>(object: &'a mut Object, ty: &str) -> Option<&'a mut Object> {
        if object.ty == ty {
            return Some(object);
        }
        for child in &mut object.children {
            if let Child::Components(components) = child {
                for nested in &mut components.children {
                    if let ComponentChild::Object(nested) = nested {
                        if let Some(found) = walk(nested, ty) {
                            return Some(found);
                        }
                    }
                }
            }
        }
        None
    }
    for child in &mut document.model.children {
        if let Child::Object(object) = child {
            if let Some(found) = walk(object, ty) {
                return Some(found);
            }
        }
    }
    None
}
