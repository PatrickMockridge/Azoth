//! A flowsheet drawn as generic Core primitives.
//!
//! **This is a schematic, not a symbol library.** A unit is a box carrying its name, a
//! stream is an orthogonal polyline with an arrow, a boundary is a pentagon and a recycle
//! is a diamond. A CAD tool replaces the symbols; azoth's job is to carry the arrangement,
//! the labels and the values so that replacement has something to start from.
//!
//! Routing is one orthogonal step, with a return lane for edges that go backwards. There is
//! no obstacle avoidance and no `DexpiLayoutEngine`: the drawing places what the canvas
//! places, which is the document's own `[layout]` or the deterministic fallback both read.

use std::collections::BTreeMap;

use super::model::{
    Color, Drawing, DrawingGroup, DrawingOrigin, GroupKind, Point, Primitive, Stroke,
};

/// The box a unit or boundary occupies in drawing units.
const NODE_WIDTH: f64 = 128.0;
const NODE_HEIGHT: f64 = 54.0;
/// The instrument bubble's edge.
const INSTRUMENT: f64 = 40.0;
/// A return lane's clearance below the nodes it joins.
const RETURN_LANE: f64 = 60.0;

/// What a node is drawn as.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(crate) enum NodeShape {
    Unit,
    Input,
    Output,
    Recycle,
    Instrument,
}

/// One node the drawing places.
#[derive(Debug, Clone)]
pub(crate) struct DiagramNode {
    /// The flowsheet key the node is addressed by.
    pub(crate) key: String,
    /// The DEXPI object id the group represents.
    pub(crate) object_id: String,
    pub(crate) label: String,
    pub(crate) detail: Option<String>,
    pub(crate) shape: NodeShape,
    pub(crate) x: f64,
    pub(crate) y: f64,
}

impl DiagramNode {
    fn width(&self) -> f64 {
        if self.shape == NodeShape::Instrument {
            INSTRUMENT
        } else {
            NODE_WIDTH
        }
    }

    fn height(&self) -> f64 {
        if self.shape == NodeShape::Instrument {
            INSTRUMENT
        } else {
            NODE_HEIGHT
        }
    }

    fn center(&self) -> Point {
        Point {
            x: self.x + self.width() / 2.0,
            y: self.y + self.height() / 2.0,
        }
    }

    fn right(&self) -> Point {
        Point {
            x: self.x + self.width(),
            y: self.y + self.height() / 2.0,
        }
    }

    fn left(&self) -> Point {
        Point {
            x: self.x,
            y: self.y + self.height() / 2.0,
        }
    }
}

/// One edge the drawing routes.
#[derive(Debug, Clone)]
pub(crate) struct DiagramEdge {
    /// The DEXPI object id the group represents.
    pub(crate) object_id: String,
    /// The source node's flowsheet key.
    pub(crate) source: String,
    /// The target node's flowsheet key.
    pub(crate) target: String,
    pub(crate) label: Option<String>,
    pub(crate) dashed: bool,
}

/// Build a drawing from placed nodes and edges.
pub(crate) fn build(name: &str, nodes: &[DiagramNode], edges: &[DiagramEdge]) -> Drawing {
    let index: BTreeMap<&str, &DiagramNode> =
        nodes.iter().map(|node| (node.key.as_str(), node)).collect();
    let mut groups = Vec::new();
    for edge in edges {
        let (Some(source), Some(target)) = (
            index.get(edge.source.as_str()),
            index.get(edge.target.as_str()),
        ) else {
            continue;
        };
        groups.push(edge_group(edge, source, target));
    }
    for node in nodes {
        groups.push(node_group(node));
    }
    Drawing::from_groups(
        name.to_string(),
        None,
        Color::WHITE,
        DrawingOrigin::Derived,
        groups,
    )
}

fn node_group(node: &DiagramNode) -> DrawingGroup {
    let mut primitives = vec![Primitive::Polygon {
        id: None,
        points: outline(node),
        fill: "Transparent".to_string(),
        stroke: Stroke::default(),
    }];
    let label_size = if node.shape == NodeShape::Instrument {
        8.0
    } else {
        10.0
    };
    primitives.push(Primitive::Text {
        id: None,
        position: Point {
            x: node.center().x,
            y: node.center().y - 2.0,
        },
        text: node.label.clone(),
        size: label_size,
        rotation: 0.0,
        alignment: "CenterCenter".to_string(),
        color: Color::BLACK,
    });
    if let Some(detail) = &node.detail {
        primitives.push(Primitive::Text {
            id: None,
            position: Point {
                x: node.center().x,
                y: node.center().y + 11.0,
            },
            text: detail.clone(),
            size: 7.0,
            rotation: 0.0,
            alignment: "CenterCenter".to_string(),
            color: Color::BLACK,
        });
    }
    DrawingGroup {
        kind: GroupKind::Node,
        id: Some(node.object_id.clone()),
        represents: Some(node.object_id.clone()),
        primitives,
    }
}

fn outline(node: &DiagramNode) -> Vec<Point> {
    let (x, y, width, height) = (node.x, node.y, node.width(), node.height());
    match node.shape {
        NodeShape::Unit | NodeShape::Instrument => vec![
            Point { x, y },
            Point { x: x + width, y },
            Point {
                x: x + width,
                y: y + height,
            },
            Point { x, y: y + height },
            Point { x, y },
        ],
        NodeShape::Recycle => vec![
            Point {
                x: x + width / 2.0,
                y,
            },
            Point {
                x: x + width,
                y: y + height / 2.0,
            },
            Point {
                x: x + width / 2.0,
                y: y + height,
            },
            Point {
                x,
                y: y + height / 2.0,
            },
            Point {
                x: x + width / 2.0,
                y,
            },
        ],
        NodeShape::Input => vec![
            Point { x, y },
            Point {
                x: x + width - 12.0,
                y,
            },
            Point {
                x: x + width,
                y: y + height / 2.0,
            },
            Point {
                x: x + width - 12.0,
                y: y + height,
            },
            Point { x, y: y + height },
            Point { x, y },
        ],
        NodeShape::Output => vec![
            Point { x: x + 12.0, y },
            Point { x: x + width, y },
            Point {
                x: x + width,
                y: y + height,
            },
            Point {
                x: x + 12.0,
                y: y + height,
            },
            Point {
                x,
                y: y + height / 2.0,
            },
            Point { x: x + 12.0, y },
        ],
    }
}

fn edge_group(edge: &DiagramEdge, source: &DiagramNode, target: &DiagramNode) -> DrawingGroup {
    let (start, end) = if source.shape == NodeShape::Instrument {
        (source.center(), target.center())
    } else if target.shape == NodeShape::Instrument {
        (source.right(), target.center())
    } else {
        (source.right(), target.left())
    };
    let points = route(start, end);
    let stroke = Stroke {
        color: Color::BLACK,
        width: 1.0,
        dash: if edge.dashed {
            "Dash".to_string()
        } else {
            "Solid".to_string()
        },
    };
    let mut primitives = vec![Primitive::Polyline {
        id: None,
        points: points.clone(),
        stroke,
    }];
    if let Some(arrow) = arrow(&points) {
        primitives.push(Primitive::Polygon {
            id: None,
            points: arrow,
            fill: "Solid".to_string(),
            stroke: Stroke::default(),
        });
    }
    if let Some(label) = &edge.label {
        if !label.is_empty() {
            primitives.push(Primitive::Text {
                id: None,
                position: Point {
                    x: (start.x + end.x) / 2.0,
                    y: (start.y + end.y) / 2.0 - 6.0,
                },
                text: label.clone(),
                size: 7.0,
                rotation: 0.0,
                alignment: "CenterCenter".to_string(),
                color: Color::BLACK,
            });
        }
    }
    DrawingGroup {
        kind: GroupKind::Edge,
        id: Some(edge.object_id.clone()),
        represents: Some(edge.object_id.clone()),
        primitives,
    }
}

fn route(start: Point, end: Point) -> Vec<Point> {
    if end.x >= start.x {
        let middle = (start.x + end.x) / 2.0;
        vec![
            start,
            Point {
                x: middle,
                y: start.y,
            },
            Point {
                x: middle,
                y: end.y,
            },
            end,
        ]
    } else {
        // A return lane, below both nodes.
        let lane = start.y.max(end.y) + RETURN_LANE;
        vec![
            start,
            Point {
                x: start.x + 24.0,
                y: start.y,
            },
            Point {
                x: start.x + 24.0,
                y: lane,
            },
            Point {
                x: end.x - 24.0,
                y: lane,
            },
            Point {
                x: end.x - 24.0,
                y: end.y,
            },
            end,
        ]
    }
}

fn arrow(points: &[Point]) -> Option<Vec<Point>> {
    let end = *points.last()?;
    let before = *points.get(points.len().checked_sub(2)?)?;
    let (dx, dy) = (end.x - before.x, end.y - before.y);
    let length = (dx * dx + dy * dy).sqrt();
    if length == 0.0 {
        return None;
    }
    let (ux, uy) = (dx / length, dy / length);
    let size = 7.0;
    let base = Point {
        x: end.x - ux * size,
        y: end.y - uy * size,
    };
    Some(vec![
        end,
        Point {
            x: base.x - uy * size / 2.0,
            y: base.y + ux * size / 2.0,
        },
        Point {
            x: base.x + uy * size / 2.0,
            y: base.y - ux * size / 2.0,
        },
        end,
    ])
}
