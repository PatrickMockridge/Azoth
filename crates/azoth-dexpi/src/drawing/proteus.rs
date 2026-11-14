//! Proteus 1.x graphics, read into the neutral drawing.
//!
//! A vendor P&ID in the field is usually Proteus: a `ShapeCatalogue` of symbol primitives,
//! each component placed by `Position`/`Scale`, segments drawn as `CenterLine` coordinates
//! and labels as `Text`. The 2.0 projection keeps the topology and drops all of this; this
//! module reads it from the same parsed tree and joins each group to the projection's
//! object ids, so a viewer can show the drawing and a coverage report can find the object.

use std::collections::{BTreeMap, BTreeSet};

use crate::diagnostic::DexpiDiagnostic;
use crate::reader::RawElement;

use super::model::{
    Color, Drawing, DrawingGroup, DrawingOrigin, DrawingReport, GroupKind, Point, Primitive, Stroke,
};

/// Read the graphics a Proteus `PlantModel` carries.
///
/// `ids` maps a source element id to the projected DEXPI object id; a group whose source id
/// has no entry is drawing furniture (an inline fitting, a segment) and carries no
/// representation.
pub(crate) fn extract(root: &RawElement, ids: &BTreeMap<String, String>) -> Option<DrawingReport> {
    let shapes = shape_catalogue(root);
    let mut draft = Draft::default();
    walk(root, None, &shapes, &mut draft);
    if draft.groups.is_empty() && draft.annotations.is_empty() {
        return None;
    }
    for (owner, primitive) in draft.labels {
        if let Some(index) = owner.as_ref().and_then(|id| draft.index.get(id)).copied() {
            draft.groups[index].primitives.push(primitive);
        } else {
            draft.annotations.push(primitive);
        }
    }

    let mut groups = Vec::new();
    for group in draft.groups {
        groups.push(DrawingGroup {
            kind: group.kind,
            id: group.id,
            represents: ids.get(&group.source).cloned(),
            primitives: group.primitives,
        });
    }
    if !draft.annotations.is_empty() {
        groups.push(DrawingGroup {
            kind: GroupKind::Annotation,
            id: None,
            represents: None,
            primitives: draft.annotations,
        });
    }

    let max_y = extent_y(root)
        .or_else(|| {
            groups
                .iter()
                .filter_map(DrawingGroup::geometry_bounds)
                .reduce(super::model::Bounds::union)
                .map(|bounds| bounds.max_y)
        })
        .unwrap_or(0.0);
    for group in &mut groups {
        for primitive in &mut group.primitives {
            primitive.flip_y(max_y);
        }
    }
    let name = drawing_name(root).unwrap_or_else(|| "Proteus drawing".to_string());
    let units = plant_information(root).and_then(|info| info.attribute("Units"));
    let drawing = Drawing::from_groups(name, units, Color::WHITE, DrawingOrigin::File, groups);
    Some(DrawingReport {
        drawing,
        diagnostics: draft.diagnostics,
    })
}

#[derive(Default)]
struct Draft {
    groups: Vec<DraftGroup>,
    index: BTreeMap<String, usize>,
    labels: Vec<(Option<String>, Primitive)>,
    annotations: Vec<Primitive>,
    missing: BTreeSet<String>,
    diagnostics: Vec<DexpiDiagnostic>,
}

struct DraftGroup {
    source: String,
    id: Option<String>,
    kind: GroupKind,
    primitives: Vec<Primitive>,
}

fn walk(
    element: &RawElement,
    owner: Option<&str>,
    shapes: &BTreeMap<String, Vec<Primitive>>,
    draft: &mut Draft,
) {
    if element.name == "ShapeCatalogue" {
        return;
    }
    let mut owner = owner.map(ToString::to_string);
    if element.name == "PipingNetworkSegment" {
        if let Some(id) = element.attribute("ID") {
            owner = Some(id);
        }
    }
    if let Some(id) = element.attribute("ID") {
        if let Some(component) = element.attribute("ComponentName") {
            if let Some(position) = location(element) {
                owner = Some(id.clone());
                if let Some(primitives) = shapes.get(&component) {
                    let scale = scale_of(element);
                    let mut placed = Vec::new();
                    for mut primitive in primitives.clone() {
                        place(&mut primitive, position, scale);
                        placed.push(primitive);
                    }
                    push_group(draft, id, GroupKind::Node, placed);
                } else if draft.missing.insert(component.clone()) {
                    draft.diagnostics.push(DexpiDiagnostic::warning(
                        "DEXPI_GRAPHICS_UNSUPPORTED",
                        component.clone(),
                        format!("`{component}` is not in the file's `ShapeCatalogue`"),
                    ));
                }
            }
        }
    }
    if element.name == "CenterLine" {
        let points: Vec<Point> = element
            .elements()
            .filter(|child| child.name == "Coordinate")
            .filter_map(point_attrs)
            .collect();
        if points.len() >= 2 {
            let stroke = stroke_of(element);
            let id = owner.clone().unwrap_or_default();
            push_group(
                draft,
                id,
                GroupKind::Edge,
                vec![Primitive::Polyline {
                    id: None,
                    points,
                    stroke,
                }],
            );
        }
    }
    if element.name == "Label" {
        if let Some(primitive) = label_primitive(element) {
            draft.labels.push((owner.clone(), primitive));
        }
    }
    for child in element.elements() {
        walk(child, owner.as_deref(), shapes, draft);
    }
}

fn push_group(draft: &mut Draft, source: String, kind: GroupKind, primitives: Vec<Primitive>) {
    if primitives.is_empty() {
        return;
    }
    if source.is_empty() {
        draft.annotations.extend(primitives);
        return;
    }
    if let Some(index) = draft.index.get(&source).copied() {
        draft.groups[index].primitives.extend(primitives);
        return;
    }
    draft.index.insert(source.clone(), draft.groups.len());
    draft.groups.push(DraftGroup {
        id: Some(source.clone()),
        source,
        kind,
        primitives,
    });
}

fn shape_catalogue(root: &RawElement) -> BTreeMap<String, Vec<Primitive>> {
    let mut shapes = BTreeMap::new();
    for catalogue in descendants(root, "ShapeCatalogue") {
        for definition in catalogue.elements() {
            let Some(name) = definition.attribute("ComponentName") else {
                continue;
            };
            let mut primitives = Vec::new();
            shape_primitives(definition, &mut primitives);
            if !primitives.is_empty() {
                shapes.insert(name, primitives);
            }
        }
    }
    shapes
}

fn shape_primitives(element: &RawElement, out: &mut Vec<Primitive>) {
    for child in element.elements() {
        match child.name.as_str() {
            "PolyLine" => {
                let points: Vec<Point> = child
                    .elements()
                    .filter(|node| node.name == "Coordinate")
                    .filter_map(point_attrs)
                    .collect();
                if points.len() >= 2 {
                    out.push(Primitive::Polyline {
                        id: None,
                        points,
                        stroke: stroke_of(child),
                    });
                }
            }
            "Polygon" => {
                let points: Vec<Point> = child
                    .elements()
                    .filter(|node| node.name == "Coordinate")
                    .filter_map(point_attrs)
                    .collect();
                if points.len() >= 3 {
                    out.push(Primitive::Polygon {
                        id: None,
                        points,
                        fill: "Transparent".to_string(),
                        stroke: stroke_of(child),
                    });
                }
            }
            "Circle" => {
                if let Some(primitive) = circle(child) {
                    out.push(primitive);
                }
            }
            "TrimmedCurve" => {
                if let Some(primitive) = arc(child) {
                    out.push(primitive);
                }
            }
            "Text" => {
                if let Some(primitive) = text_primitive(child) {
                    out.push(primitive);
                }
            }
            "Position" | "Scale" | "Presentation" | "GenericAttributes" | "ConnectionPoints" => {}
            _ => shape_primitives(child, out),
        }
    }
}

fn circle(element: &RawElement) -> Option<Primitive> {
    let radius = number(element, "Radius")?;
    let center = location(element).unwrap_or(Point { x: 0.0, y: 0.0 });
    Some(Primitive::Ellipse {
        id: None,
        center,
        radius_x: radius,
        radius_y: radius,
        rotation: 0.0,
        fill: "Transparent".to_string(),
        stroke: stroke_of(element),
    })
}

fn arc(element: &RawElement) -> Option<Primitive> {
    let start = number(element, "StartAngle")?;
    let end = number(element, "EndAngle")?;
    let circle = element.elements().find(|child| child.name == "Circle")?;
    let radius = number(circle, "Radius")?;
    let center = location(circle).unwrap_or(Point { x: 0.0, y: 0.0 });
    Some(Primitive::EllipseArc {
        id: None,
        center,
        radius_x: radius,
        radius_y: radius,
        start_angle: start,
        end_angle: end,
        rotation: 0.0,
        stroke: stroke_of(element),
    })
}

fn label_primitive(label: &RawElement) -> Option<Primitive> {
    let text = label.elements().find(|child| child.name == "Text")?;
    text_primitive(text)
}

fn text_primitive(text: &RawElement) -> Option<Primitive> {
    let value = text.attribute("String")?;
    let position = location(text)?;
    let size = number(text, "Height").unwrap_or(3.0);
    let rotation = number(text, "TextAngle").unwrap_or(0.0);
    let alignment = text
        .attribute("Justification")
        .unwrap_or_else(|| "LeftCenter".to_string());
    let (color, _, _) = presentation(text);
    Some(Primitive::Text {
        id: None,
        position,
        text: value,
        size,
        rotation,
        alignment,
        color,
    })
}

fn place(primitive: &mut Primitive, position: Point, scale: Point) {
    let apply = |point: Point| Point {
        x: position.x + point.x * scale.x,
        y: position.y + point.y * scale.y,
    };
    match primitive {
        Primitive::Polyline { points, .. } | Primitive::Polygon { points, .. } => {
            for point in points {
                *point = apply(*point);
            }
        }
        Primitive::Ellipse {
            center,
            radius_x,
            radius_y,
            ..
        }
        | Primitive::EllipseArc {
            center,
            radius_x,
            radius_y,
            ..
        } => {
            *center = apply(*center);
            *radius_x *= scale.x.abs();
            *radius_y *= scale.y.abs();
        }
        Primitive::Text {
            position: anchor,
            size,
            ..
        } => {
            *anchor = apply(*anchor);
            *size *= scale.y.abs();
        }
    }
}

fn stroke_of(element: &RawElement) -> Stroke {
    let (color, width, dash) = presentation(element);
    Stroke { color, width, dash }
}

fn presentation(element: &RawElement) -> (Color, f64, String) {
    let Some(presentation) = element
        .elements()
        .find(|child| child.name == "Presentation")
    else {
        return (Color::BLACK, 1.0, "Solid".to_string());
    };
    let color = Color {
        r: channel(number(presentation, "R").unwrap_or(0.0)),
        g: channel(number(presentation, "G").unwrap_or(0.0)),
        b: channel(number(presentation, "B").unwrap_or(0.0)),
    };
    let width = number(presentation, "LineWeight").unwrap_or(1.0);
    let dash = match presentation.attribute("LineType").as_deref() {
        Some("0") | None => "Solid",
        Some("1") => "Dash",
        Some("2") => "Dot",
        Some(_) => "Dash",
    };
    (color, width, dash.to_string())
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn channel(value: f64) -> u8 {
    let scaled = if value > 1.0 { value } else { value * 255.0 };
    scaled.round().clamp(0.0, 255.0) as u8
}

fn scale_of(element: &RawElement) -> Point {
    match element.elements().find(|child| child.name == "Scale") {
        Some(scale) => Point {
            x: number(scale, "X").unwrap_or(1.0),
            y: number(scale, "Y").unwrap_or(1.0),
        },
        None => Point { x: 1.0, y: 1.0 },
    }
}

fn location(element: &RawElement) -> Option<Point> {
    let position = element.elements().find(|child| child.name == "Position")?;
    let location = position.elements().find(|child| child.name == "Location")?;
    point_attrs(location)
}

fn point_attrs(element: &RawElement) -> Option<Point> {
    Some(Point {
        x: number(element, "X")?,
        y: number(element, "Y")?,
    })
}

fn number(element: &RawElement, name: &str) -> Option<f64> {
    element
        .attribute(name)?
        .trim()
        .parse::<f64>()
        .ok()
        .filter(|value: &f64| value.is_finite())
}

fn extent_y(root: &RawElement) -> Option<f64> {
    let drawing = descendants(root, "Drawing").into_iter().next()?;
    let extent = drawing.elements().find(|child| child.name == "Extent")?;
    let max = extent.elements().find(|child| child.name == "Max")?;
    number(max, "Y")
}

fn drawing_name(root: &RawElement) -> Option<String> {
    descendants(root, "Drawing")
        .into_iter()
        .next()
        .and_then(|drawing| drawing.attribute("Name"))
}

fn plant_information(root: &RawElement) -> Option<&RawElement> {
    root.elements()
        .find(|child| child.name == "PlantInformation")
}

fn descendants<'a>(element: &'a RawElement, name: &str) -> Vec<&'a RawElement> {
    let mut found = Vec::new();
    if element.name == name {
        found.push(element);
    }
    for child in element.elements() {
        found.extend(descendants(child, name));
    }
    found
}
