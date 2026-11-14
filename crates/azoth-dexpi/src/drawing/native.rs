//! DEXPI 2.0's `Core/Diagram` package, read into the neutral drawing and written back.
//!
//! **The writer emits the reviewed generic vocabulary and no profile symbol.** One
//! `Diagram`, one `RepresentationGroup` per represented object, a `Static` under it, and
//! `PolyLine`/`Polygon`/`Text` under that - the shape NeqSim's
//! `Dexpi20GraphicalProjectionWriter` writes, which pyDEXPI's 2.0 classes load. The reader
//! is wider: it understands ellipses, arcs, connectors, labels and shape usages, and turns
//! anything it cannot draw into a diagnostic rather than refusing the file.

use std::collections::BTreeMap;

use crate::diagnostic::DexpiDiagnostic;
use crate::ids::Mint;
use crate::model::{
    Child, ComponentChild, Components, Data, DataValue, Document, Object, References,
};

use super::model::{
    Bounds, Color, Drawing, DrawingGroup, DrawingOrigin, DrawingReport, GroupKind, Point,
    Primitive, Stroke,
};

/// Whether a diagram type is one this reader renders.
#[must_use]
pub fn renders(ty: &str) -> bool {
    matches!(
        ty,
        "Core/Diagram.PolyLine"
            | "Core/Diagram.Polygon"
            | "Core/Diagram.Ellipse"
            | "Core/Diagram.EllipseArc"
            | "Core/Diagram.Text"
            | "Core/Diagram.ConnectorLine"
            | "Core/Diagram.ShapeUsage"
            | "Core/Diagram.Static"
            | "Core/Diagram.Symbol"
            | "Core/Diagram.CustomSymbol"
            | "Core/Diagram.Label"
            | "Core/Diagram.Border"
            | "Core/Diagram.RepresentationGroup"
            | "Core/Diagram.Diagram"
    ) || ty.ends_with(".Label")
        || ty.ends_with(".Static")
        || ty.ends_with(".Symbol")
        || ty.ends_with(".RepresentationGroup")
        || ty.ends_with(".NodePosition")
        || matches!(
            ty,
            "Core/Diagram.MetaData"
                | "Core/Diagram.ShapeCatalogue"
                | "Core/Diagram.Shape"
                | "Core/Diagram.TextTemplate"
                | "Plant/Diagram.PlantMetaData"
        )
}

/// Read the first diagram a document carries.
///
/// `None` means the document states no `Core/Diagram.Diagram`; a diagram whose primitives
/// are all unsupported still returns a report, with its diagnostics.
#[must_use]
pub fn read_document(document: &Document) -> Option<DrawingReport> {
    let diagram = find_diagram(document)?;
    let mut reader = Reader::new(document);
    let groups = reader.read_top(diagram);
    let name = text(diagram, "Name")
        .unwrap_or(&document.model.name)
        .to_string();
    let background = color_of(diagram, "BackgroundColor").unwrap_or(Color::WHITE);
    let max_y = number(diagram, "MaxY")
        .filter(|value| value.is_finite())
        .or_else(|| {
            groups
                .iter()
                .filter_map(DrawingGroup::geometry_bounds)
                .reduce(Bounds::union)
                .map(|bounds| bounds.max_y)
        })
        .unwrap_or(0.0);
    let mut groups = groups;
    for group in &mut groups {
        for primitive in &mut group.primitives {
            primitive.flip_y(max_y);
        }
    }
    let drawing = Drawing::from_groups(
        name,
        Some("mm".to_string()),
        background,
        DrawingOrigin::File,
        groups,
    );
    Some(DrawingReport {
        drawing,
        diagnostics: reader.diagnostics,
    })
}

/// Build a `Core/Diagram.Diagram` from a drawing.
///
/// **This is the only drawing vocabulary azoth writes**: polylines, polygons and text. A
/// group with no represented object is not written, because the model requires `Represents`.
#[must_use]
pub fn write_diagram(drawing: &Drawing, represents: &str, mint: &mut Mint) -> Object {
    let height = drawing.bounds.height();
    let mut groups = Vec::new();
    for group in &drawing.groups {
        let Some(represented) = group.represents.as_deref() else {
            continue;
        };
        let mut elements = Vec::new();
        for primitive in &group.primitives {
            if let Some(written) = primitive_object(primitive, height, mint) {
                elements.push(ComponentChild::Object(written));
            }
        }
        if elements.is_empty() {
            continue;
        }
        let static_group = Object {
            id: Some(mint.next("Static")),
            name: None,
            ty: "Core/Diagram.Static".to_string(),
            children: vec![Child::Components(Components {
                property: "Elements".to_string(),
                children: elements,
            })],
        };
        groups.push(ComponentChild::Object(Object {
            id: Some(mint.next("RepresentationGroup")),
            name: None,
            ty: "Core/Diagram.RepresentationGroup".to_string(),
            children: vec![
                Child::References(References {
                    property: "Represents".to_string(),
                    objects: vec![format!("#{represented}")],
                }),
                Child::Components(Components {
                    property: "Groups".to_string(),
                    children: vec![ComponentChild::Object(static_group)],
                }),
            ],
        }));
    }

    let mut children = vec![
        Child::Data(string_data("Name", &drawing.name)),
        Child::Data(color_data("BackgroundColor", drawing.background)),
        Child::Data(double_data("MinX", 0.0)),
        Child::Data(double_data("MinY", 0.0)),
        Child::Data(double_data("MaxX", drawing.bounds.width())),
        Child::Data(double_data("MaxY", height)),
        Child::References(References {
            property: "Represents".to_string(),
            objects: vec![format!("#{represents}")],
        }),
    ];
    if !groups.is_empty() {
        children.push(Child::Components(Components {
            property: "Groups".to_string(),
            children: groups,
        }));
    }
    Object {
        id: Some(mint.next("Diagram")),
        name: None,
        ty: "Core/Diagram.Diagram".to_string(),
        children,
    }
}

fn primitive_object(primitive: &Primitive, height: f64, mint: &mut Mint) -> Option<Object> {
    let mut primitive = primitive.clone();
    primitive.flip_y(height);
    let id = Some(mint.next("Primitive"));
    match primitive {
        Primitive::Polyline { points, stroke, .. } if points.len() >= 2 => Some(Object {
            id,
            name: None,
            ty: "Core/Diagram.PolyLine".to_string(),
            children: vec![
                Child::Data(points_data("Points", &points)),
                Child::Data(stroke_data(&stroke)),
            ],
        }),
        Primitive::Polygon {
            points,
            fill,
            stroke,
            ..
        } if points.len() >= 3 => Some(Object {
            id,
            name: None,
            ty: "Core/Diagram.Polygon".to_string(),
            children: vec![
                Child::Data(reference_data("FillStyle", &fill_style(&fill))),
                Child::Data(points_data("Points", &points)),
                Child::Data(stroke_data(&stroke)),
            ],
        }),
        Primitive::Text {
            position,
            text: value,
            size,
            rotation,
            alignment,
            color,
            ..
        } => Some(Object {
            id,
            name: None,
            ty: "Core/Diagram.Text".to_string(),
            children: vec![
                Child::Data(reference_data(
                    "Alignment",
                    &format!("Core/Diagram.TextAlignment.{}", text_alignment(&alignment)),
                )),
                Child::Data(color_data("Color", color)),
                Child::Data(string_data("Font", "Arial")),
                Child::Data(point_data("Position", position)),
                Child::Data(double_data("Rotation", rotation)),
                Child::Data(double_data("Size", size)),
                Child::Data(string_data("Text", &value)),
            ],
        }),
        // A reader can hold these, but the writer's reviewed vocabulary does not emit them.
        Primitive::Ellipse { .. } | Primitive::EllipseArc { .. } => None,
        Primitive::Polyline { .. } | Primitive::Polygon { .. } => None,
    }
}

fn fill_style(fill: &str) -> String {
    match fill {
        "Hatch" | "Solid" | "Transparent" => format!("Core/Diagram.FillStyle.{fill}"),
        _ => "Core/Diagram.FillStyle.Transparent".to_string(),
    }
}

fn text_alignment(alignment: &str) -> String {
    match alignment {
        "CenterBottom" | "CenterCenter" | "CenterTop" | "LeftBottom" | "LeftCenter" | "LeftTop"
        | "RightBottom" | "RightCenter" | "RightTop" => alignment.to_string(),
        _ => "LeftCenter".to_string(),
    }
}

// --- reading -----------------------------------------------------------------------------

struct Reader {
    types: BTreeMap<String, String>,
    nodes: BTreeMap<String, Point>,
    shapes: BTreeMap<String, Vec<Primitive>>,
    diagnostics: Vec<DexpiDiagnostic>,
}

impl Reader {
    fn new(document: &Document) -> Self {
        let mut types = BTreeMap::new();
        let mut nodes = BTreeMap::new();
        let mut shapes = Vec::new();
        for child in &document.model.children {
            if let Child::Object(object) = child {
                collect(object, &mut types, &mut nodes, &mut shapes);
            }
        }
        let mut reader = Self {
            types,
            nodes,
            shapes: BTreeMap::new(),
            diagnostics: Vec::new(),
        };
        for shape in shapes {
            let mut primitives = Vec::new();
            reader.read_shape(shape, &mut primitives);
            if primitives.is_empty() {
                continue;
            }
            if let Some(id) = &shape.id {
                reader.shapes.insert(id.clone(), primitives.clone());
            }
            if let Some(name) = text(shape, "Name") {
                reader.shapes.insert(name.to_string(), primitives);
            }
        }
        reader
    }

    fn read_top(&mut self, diagram: &Object) -> Vec<DrawingGroup> {
        let mut groups = Vec::new();
        self.read_group(diagram, &mut groups);
        groups
    }

    fn read_group(&mut self, object: &Object, groups: &mut Vec<DrawingGroup>) {
        let represents = object
            .references("Represents")
            .and_then(|references| references.objects.first())
            .map(|target| target.trim_start_matches('#').to_string());
        let mut primitives = Vec::new();
        if let Some(elements) = object.components("Elements") {
            for child in &elements.children {
                if let ComponentChild::Object(child) = child {
                    self.read_container(child, &mut primitives, groups);
                }
            }
        }
        if let Some(children) = object.components("Groups") {
            for child in &children.children {
                let ComponentChild::Object(child) = child else {
                    continue;
                };
                if is_representation_group(&child.ty) && child.references("Represents").is_some() {
                    self.read_group(child, groups);
                } else {
                    self.read_container(child, &mut primitives, groups);
                }
            }
        }
        if primitives.is_empty() {
            return;
        }
        groups.push(DrawingGroup {
            kind: self.kind_of(represents.as_deref()),
            id: object.id.clone(),
            represents,
            primitives,
        });
    }

    fn read_container(
        &mut self,
        object: &Object,
        out: &mut Vec<Primitive>,
        groups: &mut Vec<DrawingGroup>,
    ) {
        if is_container(&object.ty) {
            if let Some(elements) = object.components("Elements") {
                for child in &elements.children {
                    if let ComponentChild::Object(child) = child {
                        self.read_container(child, out, groups);
                    }
                }
            }
            if let Some(children) = object.components("Groups") {
                for child in &children.children {
                    let ComponentChild::Object(child) = child else {
                        continue;
                    };
                    if is_representation_group(&child.ty)
                        && child.references("Represents").is_some()
                    {
                        self.read_group(child, groups);
                    } else {
                        self.read_container(child, out, groups);
                    }
                }
            }
            return;
        }
        self.read_primitive(object, out);
    }

    fn read_shape(&mut self, shape: &Object, out: &mut Vec<Primitive>) {
        for property in ["Primitives", "Elements"] {
            if let Some(children) = shape.components(property) {
                for child in &children.children {
                    if let ComponentChild::Object(child) = child {
                        self.read_primitive(child, out);
                    }
                }
            }
        }
        if let Some(children) = shape.components("Groups") {
            for child in &children.children {
                if let ComponentChild::Object(child) = child {
                    self.read_container(child, out, &mut Vec::new());
                }
            }
        }
    }

    fn read_primitive(&mut self, object: &Object, out: &mut Vec<Primitive>) {
        let id = object.id.clone();
        let stroke = stroke_of(object);
        let primitive = match object.ty.as_str() {
            "Core/Diagram.PolyLine" => {
                let points = self.points_of(object, "Points");
                if points.len() < 2 {
                    self.unsupported(object, "a polyline needs at least two points");
                    return;
                }
                Primitive::Polyline { id, points, stroke }
            }
            "Core/Diagram.Polygon" => {
                let points = self.points_of(object, "Points");
                if points.len() < 3 {
                    self.unsupported(object, "a polygon needs at least three points");
                    return;
                }
                Primitive::Polygon {
                    id,
                    points,
                    fill: fill_of(object),
                    stroke,
                }
            }
            "Core/Diagram.Ellipse" => {
                let Some(center) = point_property(object, "Center") else {
                    self.unsupported(object, "an ellipse needs a center");
                    return;
                };
                let (Some(radius_x), Some(radius_y)) = (
                    number(object, "HorizontalSemiAxis"),
                    number(object, "VerticalSemiAxis"),
                ) else {
                    self.unsupported(object, "an ellipse needs two semi-axes");
                    return;
                };
                Primitive::Ellipse {
                    id,
                    center,
                    radius_x,
                    radius_y,
                    rotation: number(object, "Rotation").unwrap_or(0.0),
                    fill: fill_of(object),
                    stroke,
                }
            }
            "Core/Diagram.EllipseArc" => {
                let Some(center) = point_property(object, "Center") else {
                    self.unsupported(object, "an ellipse arc needs a center");
                    return;
                };
                let (Some(radius_x), Some(radius_y), Some(start_angle), Some(end_angle)) = (
                    number(object, "HorizontalSemiAxis"),
                    number(object, "VerticalSemiAxis"),
                    number(object, "StartAngle"),
                    number(object, "EndAngle"),
                ) else {
                    self.unsupported(object, "an ellipse arc needs its axes and angles");
                    return;
                };
                Primitive::EllipseArc {
                    id,
                    center,
                    radius_x,
                    radius_y,
                    start_angle,
                    end_angle,
                    rotation: number(object, "Rotation").unwrap_or(0.0),
                    stroke,
                }
            }
            "Core/Diagram.Text" => {
                let Some(position) = point_property(object, "Position") else {
                    self.unsupported(object, "a text needs a position");
                    return;
                };
                let Some(value) = text(object, "Text") else {
                    self.unsupported(object, "a text needs its string");
                    return;
                };
                Primitive::Text {
                    id,
                    position,
                    text: value.to_string(),
                    size: number(object, "Size").unwrap_or(3.0),
                    rotation: number(object, "Rotation").unwrap_or(0.0),
                    alignment: literal_of(object, "Alignment")
                        .unwrap_or("LeftCenter")
                        .to_string(),
                    color: color_of(object, "Color").unwrap_or(Color::BLACK),
                }
            }
            "Core/Diagram.ConnectorLine" => {
                let mut points = Vec::new();
                if let Some(source) = self.node_point(object, "Source") {
                    points.push(source);
                }
                points.extend(self.points_of(object, "InnerPoints"));
                if let Some(target) = self.node_point(object, "Target") {
                    points.push(target);
                }
                if points.len() < 2 {
                    self.unsupported(object, "a connector line needs its two ends");
                    return;
                }
                Primitive::Polyline { id, points, stroke }
            }
            "Core/Diagram.ShapeUsage" => {
                self.read_shape_usage(object, out);
                return;
            }
            other => {
                self.unsupported(object, other);
                return;
            }
        };
        if primitive.is_finite() {
            out.push(primitive);
        } else {
            self.unsupported(object, "a coordinate is not finite");
        }
    }

    fn read_shape_usage(&mut self, object: &Object, out: &mut Vec<Primitive>) {
        let Some(shape) = object
            .references("Shape")
            .and_then(|references| references.objects.first())
            .map(|target| target.trim_start_matches('#').to_string())
        else {
            self.unsupported(object, "a shape usage names no shape");
            return;
        };
        let Some(primitives) = self.shapes.get(&shape).cloned() else {
            self.unsupported(object, &format!("shape `{shape}` is not in the catalogue"));
            return;
        };
        let position = point_property(object, "Position").unwrap_or(Point { x: 0.0, y: 0.0 });
        let rotation = number(object, "Rotation").unwrap_or(0.0);
        let scale_x = number(object, "ScaleX").unwrap_or(1.0);
        let scale_y = number(object, "ScaleY").unwrap_or(1.0);
        let mirrored = boolean(object, "IsMirrored").unwrap_or(false);
        for mut primitive in primitives {
            transform(
                &mut primitive,
                position,
                rotation,
                scale_x,
                scale_y,
                mirrored,
            );
            if primitive.is_finite() {
                out.push(primitive);
            }
        }
    }

    fn points_of(&self, object: &Object, property: &str) -> Vec<Point> {
        object
            .children
            .iter()
            .filter_map(|child| match child {
                Child::Data(data) if data.property == property => Some(&data.values),
                _ => None,
            })
            .flatten()
            .filter_map(point_of)
            .collect()
    }

    fn node_point(&self, object: &Object, property: &str) -> Option<Point> {
        let target = object
            .references(property)?
            .objects
            .first()?
            .trim_start_matches('#');
        self.nodes.get(target).copied()
    }

    fn kind_of(&self, represents: Option<&str>) -> GroupKind {
        let Some(ty) = represents.and_then(|id| self.types.get(id)) else {
            return GroupKind::Annotation;
        };
        if ty.ends_with(".Stream")
            || ty.ends_with(".PipingNetworkSegment")
            || ty.ends_with(".SignalConveyingFunction")
        {
            GroupKind::Edge
        } else if ty.contains("ProcessEquipment.")
            || ty.contains(".ProcessInstrumentationFunction")
            || ty.contains(".Piping.GlobeValve")
            || ty.ends_with(".Valve")
            || ty.contains("Process/Process.")
        {
            GroupKind::Node
        } else {
            GroupKind::Annotation
        }
    }

    fn unsupported(&mut self, object: &Object, what: &str) {
        self.diagnostics.push(DexpiDiagnostic::warning(
            "DEXPI_GRAPHICS_UNSUPPORTED",
            object.id.clone().unwrap_or_default(),
            format!("`{}` is not drawn: {what}", object.ty),
        ));
    }
}

fn collect<'a>(
    object: &'a Object,
    types: &mut BTreeMap<String, String>,
    nodes: &mut BTreeMap<String, Point>,
    shapes: &mut Vec<&'a Object>,
) {
    if let Some(id) = &object.id {
        types.insert(id.clone(), object.ty.clone());
        if object.ty.ends_with(".NodePosition") {
            if let Some(point) = point_property(object, "Position") {
                nodes.insert(id.clone(), point);
            }
        }
    }
    if object.ty.ends_with(".Shape") {
        shapes.push(object);
    }
    for child in &object.children {
        if let Child::Components(components) = child {
            for nested in &components.children {
                if let ComponentChild::Object(nested) = nested {
                    collect(nested, types, nodes, shapes);
                }
            }
        }
    }
}

fn find_diagram(document: &Document) -> Option<&Object> {
    fn walk(object: &Object) -> Option<&Object> {
        if object.ty == "Core/Diagram.Diagram" {
            return Some(object);
        }
        for child in &object.children {
            if let Child::Components(components) = child {
                for nested in &components.children {
                    if let ComponentChild::Object(nested) = nested {
                        if let Some(found) = walk(nested) {
                            return Some(found);
                        }
                    }
                }
            }
        }
        None
    }
    for child in &document.model.children {
        if let Child::Object(object) = child {
            if let Some(found) = walk(object) {
                return Some(found);
            }
        }
    }
    None
}

fn is_representation_group(ty: &str) -> bool {
    ty.ends_with(".RepresentationGroup") || ty.ends_with(".Diagram")
}

fn is_container(ty: &str) -> bool {
    ty.ends_with(".Static")
        || ty.ends_with(".Symbol")
        || ty.ends_with(".CustomSymbol")
        || ty.ends_with(".Label")
        || ty.ends_with(".Border")
}

fn transform(
    primitive: &mut Primitive,
    position: Point,
    shape_rotation: f64,
    scale_x: f64,
    scale_y: f64,
    mirrored: bool,
) {
    let (sin, cos) = shape_rotation.to_radians().sin_cos();
    let apply = |point: Point| -> Point {
        let x = if mirrored { -point.x } else { point.x } * scale_x;
        let y = point.y * scale_y;
        Point {
            x: position.x + x * cos - y * sin,
            y: position.y + x * sin + y * cos,
        }
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
            rotation,
            ..
        }
        | Primitive::EllipseArc {
            center,
            radius_x,
            radius_y,
            rotation,
            ..
        } => {
            *center = apply(*center);
            *radius_x *= scale_x.abs();
            *radius_y *= scale_y.abs();
            *rotation += shape_rotation;
        }
        Primitive::Text {
            position,
            rotation,
            size,
            ..
        } => {
            *position = apply(*position);
            *size *= scale_y.abs();
            *rotation += shape_rotation;
        }
    }
}

// --- values ------------------------------------------------------------------------------

fn text<'a>(object: &'a Object, property: &str) -> Option<&'a str> {
    match object.data_value(property)? {
        DataValue::String(value) => Some(value),
        _ => None,
    }
}

fn number(object: &Object, property: &str) -> Option<f64> {
    object.data_value(property).and_then(DataValue::as_f64)
}

fn boolean(object: &Object, property: &str) -> Option<bool> {
    match object.data_value(property)? {
        DataValue::Boolean(value) => Some(*value),
        _ => None,
    }
}

fn literal_of<'a>(object: &'a Object, property: &str) -> Option<&'a str> {
    match object.data_value(property)? {
        DataValue::Reference(value) => Some(literal(value)),
        DataValue::String(value) => Some(literal(value)),
        _ => None,
    }
}

fn literal(reference: &str) -> &str {
    reference.rsplit('.').next().unwrap_or(reference)
}

fn point_property(object: &Object, property: &str) -> Option<Point> {
    object.data_values(property).find_map(point_of)
}

fn point_of(value: &DataValue) -> Option<Point> {
    let DataValue::Aggregated { ty, data } = value else {
        return None;
    };
    if !ty.ends_with(".Point") {
        return None;
    }
    let x = number_in(data, "X")?;
    let y = number_in(data, "Y")?;
    if x.is_finite() && y.is_finite() {
        Some(Point { x, y })
    } else {
        None
    }
}

fn number_in(data: &[Data], property: &str) -> Option<f64> {
    data.iter()
        .find(|entry| entry.property == property)
        .and_then(|entry| entry.values.first())
        .and_then(DataValue::as_f64)
}

fn color_of(object: &Object, property: &str) -> Option<Color> {
    object.data_values(property).find_map(|value| match value {
        DataValue::Aggregated { ty, data } if ty.ends_with(".Color") => {
            let r = number_in(data, "R")?;
            let g = number_in(data, "G")?;
            let b = number_in(data, "B")?;
            Some(Color {
                r: channel(r),
                g: channel(g),
                b: channel(b),
            })
        }
        DataValue::String(text) => parse_hex(text),
        _ => None,
    })
}

fn parse_hex(text: &str) -> Option<Color> {
    let digits = text.trim().strip_prefix('#')?;
    if digits.len() != 6 {
        return None;
    }
    let r = u8::from_str_radix(digits.get(0..2)?, 16).ok()?;
    let g = u8::from_str_radix(digits.get(2..4)?, 16).ok()?;
    let b = u8::from_str_radix(digits.get(4..6)?, 16).ok()?;
    Some(Color { r, g, b })
}

#[allow(clippy::cast_possible_truncation, clippy::cast_sign_loss)]
fn channel(value: f64) -> u8 {
    value.round().clamp(0.0, 255.0) as u8
}

fn stroke_of(object: &Object) -> Stroke {
    let Some(DataValue::Aggregated { ty, data }) = object.data_value("Stroke") else {
        return Stroke::default();
    };
    if !ty.ends_with(".Stroke") {
        return Stroke::default();
    }
    let color = data
        .iter()
        .find(|entry| entry.property == "Color")
        .and_then(|entry| entry.values.first())
        .and_then(|value| match value {
            DataValue::Aggregated { ty, data } if ty.ends_with(".Color") => {
                let r = number_in(data, "R")?;
                let g = number_in(data, "G")?;
                let b = number_in(data, "B")?;
                Some(Color {
                    r: channel(r),
                    g: channel(g),
                    b: channel(b),
                })
            }
            _ => None,
        })
        .unwrap_or(Color::BLACK);
    let dash = data
        .iter()
        .find(|entry| entry.property == "DashStyle")
        .and_then(|entry| entry.values.first())
        .and_then(|value| match value {
            DataValue::Reference(reference) => Some(literal(reference).to_string()),
            _ => None,
        })
        .unwrap_or_else(|| "Solid".to_string());
    let width = number_in(data, "Width").unwrap_or(1.0);
    Stroke {
        color,
        width: if width.is_finite() { width } else { 1.0 },
        dash,
    }
}

fn fill_of(object: &Object) -> String {
    literal_of(object, "FillStyle")
        .unwrap_or("Transparent")
        .to_string()
}

// --- writing -----------------------------------------------------------------------------

fn string_data(property: &str, value: &str) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::String(value.to_string())],
    }
}

fn double_data(property: &str, value: f64) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::Double(value)],
    }
}

fn reference_data(property: &str, value: &str) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::Reference(value.to_string())],
    }
}

fn point_data(property: &str, point: Point) -> Data {
    Data {
        property: property.to_string(),
        values: vec![point_value(point)],
    }
}

fn points_data(property: &str, points: &[Point]) -> Data {
    Data {
        property: property.to_string(),
        values: points.iter().copied().map(point_value).collect(),
    }
}

fn point_value(point: Point) -> DataValue {
    DataValue::Aggregated {
        ty: "Core/Diagram.Point".to_string(),
        data: vec![double_data("X", point.x), double_data("Y", point.y)],
    }
}

fn color_data(property: &str, color: Color) -> Data {
    Data {
        property: property.to_string(),
        values: vec![DataValue::Aggregated {
            ty: "Core/Diagram.Color".to_string(),
            data: vec![
                Data {
                    property: "B".to_string(),
                    values: vec![DataValue::Integer(i64::from(color.b))],
                },
                Data {
                    property: "G".to_string(),
                    values: vec![DataValue::Integer(i64::from(color.g))],
                },
                Data {
                    property: "R".to_string(),
                    values: vec![DataValue::Integer(i64::from(color.r))],
                },
            ],
        }],
    }
}

fn stroke_data(stroke: &Stroke) -> Data {
    let color = color_data("Color", stroke.color);
    Data {
        property: "Stroke".to_string(),
        values: vec![DataValue::Aggregated {
            ty: "Core/Diagram.Stroke".to_string(),
            data: vec![
                color,
                reference_data(
                    "DashStyle",
                    &format!("Core/Diagram.DashStyle.{}", dash_style(&stroke.dash)),
                ),
                double_data("Width", stroke.width),
            ],
        }],
    }
}

fn dash_style(dash: &str) -> String {
    match dash {
        "Dash"
        | "DashShortDash"
        | "Dot"
        | "LongDash"
        | "LongDashShortDash"
        | "LongDashShortDashShortDash"
        | "ShortDash"
        | "Solid" => dash.to_string(),
        _ => "Solid".to_string(),
    }
}
