//! A resolved, profile-neutral drawing.
//!
//! The exchange's two graphic dialects - DEXPI 2.0 `Core/Diagram` and Proteus
//! `ShapeCatalogue` + `Position` - both reduce to this: groups of polylines, polygons,
//! ellipses and text. **The model is presentation, not process state**: nothing here is
//! validated as a quantity, and a viewer draws it rather than a simulation reading it.
//!
//! Coordinates are y-down and normalised to a zero origin. DEXPI and CAD are y-up, so the
//! readers flip once and translate; the writer flips back. Position survives an exchange up
//! to a translation, which is what the flowsheet's `[layout]` view needs and all it claims.

use serde::{Deserialize, Serialize};

/// Where a drawing came from.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum DrawingOrigin {
    /// The file carries graphics.
    File,
    /// The file carries none; this is a schematic derived from its topology.
    Derived,
}

/// What a group draws, which decides the layer order and what a selection means.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum GroupKind {
    Node,
    Edge,
    Annotation,
}

/// A point in drawing units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Point {
    pub x: f64,
    pub y: f64,
}

/// An RGB colour, one byte per channel.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Color {
    pub r: u8,
    pub g: u8,
    pub b: u8,
}

impl Color {
    pub const BLACK: Self = Self { r: 0, g: 0, b: 0 };
    pub const WHITE: Self = Self {
        r: 255,
        g: 255,
        b: 255,
    };
}

/// How a line is drawn.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Stroke {
    pub color: Color,
    pub width: f64,
    /// The DEXPI `DashStyle` literal, kept as read so an unreviewed value survives.
    pub dash: String,
}

impl Default for Stroke {
    fn default() -> Self {
        Self {
            color: Color::BLACK,
            width: 1.0,
            dash: "Solid".to_string(),
        }
    }
}

/// A rectangular extent in drawing units.
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub struct Bounds {
    pub min_x: f64,
    pub min_y: f64,
    pub max_x: f64,
    pub max_y: f64,
}

impl Bounds {
    #[must_use]
    pub fn new(min_x: f64, min_y: f64, max_x: f64, max_y: f64) -> Self {
        Self {
            min_x,
            min_y,
            max_x,
            max_y,
        }
    }

    #[must_use]
    pub fn width(self) -> f64 {
        self.max_x - self.min_x
    }

    #[must_use]
    pub fn height(self) -> f64 {
        self.max_y - self.min_y
    }

    /// The smallest box holding both. A non-finite box is ignored.
    #[must_use]
    pub fn union(self, other: Self) -> Self {
        Self {
            min_x: self.min_x.min(other.min_x),
            min_y: self.min_y.min(other.min_y),
            max_x: self.max_x.max(other.max_x),
            max_y: self.max_y.max(other.max_y),
        }
    }

    /// Whether every corner is finite.
    #[must_use]
    pub fn is_finite(self) -> bool {
        self.min_x.is_finite()
            && self.min_y.is_finite()
            && self.max_x.is_finite()
            && self.max_y.is_finite()
    }
}

/// One primitive. `kind` is the serde tag the wire and the viewer read.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
pub enum Primitive {
    Polyline {
        id: Option<String>,
        points: Vec<Point>,
        stroke: Stroke,
    },
    Polygon {
        id: Option<String>,
        points: Vec<Point>,
        fill: String,
        stroke: Stroke,
    },
    Ellipse {
        id: Option<String>,
        center: Point,
        radius_x: f64,
        radius_y: f64,
        rotation: f64,
        fill: String,
        stroke: Stroke,
    },
    EllipseArc {
        id: Option<String>,
        center: Point,
        radius_x: f64,
        radius_y: f64,
        start_angle: f64,
        end_angle: f64,
        rotation: f64,
        stroke: Stroke,
    },
    Text {
        id: Option<String>,
        position: Point,
        text: String,
        size: f64,
        rotation: f64,
        alignment: String,
        color: Color,
    },
}

impl Primitive {
    /// The primitive's own id, when the file gave one.
    #[must_use]
    pub fn id(&self) -> Option<&str> {
        match self {
            Self::Polyline { id, .. }
            | Self::Polygon { id, .. }
            | Self::Ellipse { id, .. }
            | Self::EllipseArc { id, .. }
            | Self::Text { id, .. } => id.as_deref(),
        }
    }

    /// Whether every coordinate is finite. A primitive that is not is dropped and reported.
    #[must_use]
    pub fn is_finite(&self) -> bool {
        match self {
            Self::Polyline { points, .. } | Self::Polygon { points, .. } => {
                !points.is_empty() && points.iter().all(|point| finite(*point))
            }
            Self::Ellipse {
                center,
                radius_x,
                radius_y,
                rotation,
                ..
            }
            | Self::EllipseArc {
                center,
                radius_x,
                radius_y,
                rotation,
                ..
            } => {
                finite(*center)
                    && radius_x.is_finite()
                    && radius_y.is_finite()
                    && rotation.is_finite()
            }
            Self::Text {
                position,
                size,
                rotation,
                ..
            } => finite(*position) && size.is_finite() && rotation.is_finite(),
        }
    }

    /// Translate every point.
    pub fn translate(&mut self, dx: f64, dy: f64) {
        match self {
            Self::Polyline { points, .. } | Self::Polygon { points, .. } => {
                for point in points {
                    point.x += dx;
                    point.y += dy;
                }
            }
            Self::Ellipse { center, .. }
            | Self::EllipseArc { center, .. }
            | Self::Text {
                position: center, ..
            } => {
                center.x += dx;
                center.y += dy;
            }
        }
    }

    /// Flip the y axis about `axis`.
    pub fn flip_y(&mut self, axis: f64) {
        match self {
            Self::Polyline { points, .. } | Self::Polygon { points, .. } => {
                for point in points {
                    point.y = axis - point.y;
                }
            }
            Self::Ellipse { center, .. }
            | Self::EllipseArc { center, .. }
            | Self::Text {
                position: center, ..
            } => {
                center.y = axis - center.y;
            }
        }
    }

    /// The primitive's geometry box. A text carries only its anchor, which is why a group
    /// ignores text when it has anything else to measure.
    #[must_use]
    pub fn bounds(&self) -> Option<Bounds> {
        match self {
            Self::Polyline { points, .. } | Self::Polygon { points, .. } => points_box(points),
            Self::Ellipse {
                center,
                radius_x,
                radius_y,
                ..
            }
            | Self::EllipseArc {
                center,
                radius_x,
                radius_y,
                ..
            } => Some(Bounds::new(
                center.x - radius_x,
                center.y - radius_y,
                center.x + radius_x,
                center.y + radius_y,
            )),
            Self::Text { position, .. } => {
                Some(Bounds::new(position.x, position.y, position.x, position.y))
            }
        }
    }

    /// Whether this is text, which a group excludes from its placement box when it can.
    #[must_use]
    pub fn is_text(&self) -> bool {
        matches!(self, Self::Text { .. })
    }
}

fn finite(point: Point) -> bool {
    point.x.is_finite() && point.y.is_finite()
}

fn points_box(points: &[Point]) -> Option<Bounds> {
    let first = points.first()?;
    let mut bounds = Bounds::new(first.x, first.y, first.x, first.y);
    for point in points.iter().skip(1) {
        bounds = bounds.union(Bounds::new(point.x, point.y, point.x, point.y));
    }
    Some(bounds)
}

/// One represented object's graphics, or drawing furniture with no represented object.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct DrawingGroup {
    pub kind: GroupKind,
    pub id: Option<String>,
    /// The DEXPI object id this group represents, or an azoth key for a derived drawing.
    pub represents: Option<String>,
    pub primitives: Vec<Primitive>,
}

impl DrawingGroup {
    /// The placement box: the geometry of every non-text primitive, or the text when that
    /// is all the group has.
    #[must_use]
    pub fn geometry_bounds(&self) -> Option<Bounds> {
        let geometry = self
            .primitives
            .iter()
            .filter(|primitive| !primitive.is_text())
            .filter_map(Primitive::bounds)
            .reduce(Bounds::union);
        geometry.or_else(|| {
            self.primitives
                .iter()
                .filter_map(Primitive::bounds)
                .reduce(Bounds::union)
        })
    }
}

/// A drawing, resolved and normalised.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Drawing {
    pub name: String,
    /// The source's stated unit (`mm` for Proteus), never converted.
    pub units: Option<String>,
    pub background: Color,
    pub bounds: Bounds,
    pub origin: DrawingOrigin,
    pub groups: Vec<DrawingGroup>,
}

impl Drawing {
    /// Assemble a drawing from groups, translating every primitive to a zero origin.
    #[must_use]
    pub fn from_groups(
        name: String,
        units: Option<String>,
        background: Color,
        origin: DrawingOrigin,
        mut groups: Vec<DrawingGroup>,
    ) -> Self {
        let bounds = groups
            .iter()
            .filter_map(DrawingGroup::geometry_bounds)
            .reduce(Bounds::union)
            .unwrap_or(Bounds::new(0.0, 0.0, 0.0, 0.0));
        let bounds = if bounds.is_finite() {
            bounds
        } else {
            Bounds::new(0.0, 0.0, 0.0, 0.0)
        };
        if bounds.min_x != 0.0 || bounds.min_y != 0.0 {
            for group in &mut groups {
                for primitive in &mut group.primitives {
                    primitive.translate(-bounds.min_x, -bounds.min_y);
                }
            }
        }
        let bounds = Bounds::new(0.0, 0.0, bounds.width(), bounds.height());
        Self {
            name,
            units,
            background,
            bounds,
            origin,
            groups,
        }
    }

    /// Every primitive, in group order.
    pub fn primitives(&self) -> impl Iterator<Item = &Primitive> {
        self.groups.iter().flat_map(|group| group.primitives.iter())
    }

    /// Number of primitives, for the report.
    #[must_use]
    pub fn primitive_count(&self) -> usize {
        self.primitives().count()
    }
}

/// A drawing and the losses reading it incurred.
#[derive(Debug, Clone, PartialEq)]
pub struct DrawingReport {
    pub drawing: Drawing,
    pub diagnostics: Vec<crate::diagnostic::DexpiDiagnostic>,
}
