//! A drawing's arrangement, read into the flowsheet's `[layout]`.
//!
//! **The arrangement is the part of a CAD drawing azoth keeps.** Shape and symbol geometry
//! belong to the drawing, but where a unit sits is the editor's placement, so a P&ID opened
//! as a simulation draws the way it was drawn. The mapping is relative: the drawing is
//! normalised to a zero origin, so only positions relative to one another survive, which is
//! all `[layout]` claims.

use std::collections::BTreeMap;

use azoth_process::flowsheet::{Flowsheet, Layout};

use super::model::Drawing;

/// Place every represented group on the flowsheet node its binding names.
///
/// Returns how many positions were written. A group with no binding is drawing furniture and
/// is skipped; an instrument has no `[layout]` entry and is skipped too.
pub(crate) fn apply(
    flowsheet: &mut Flowsheet,
    drawing: &Drawing,
    bindings: &BTreeMap<String, String>,
) -> usize {
    let mut layout = Layout::default();
    let mut placed = 0_usize;
    for group in &drawing.groups {
        let Some(represented) = group.represents.as_deref() else {
            continue;
        };
        let Some(name) = bindings.get(represented) else {
            continue;
        };
        let Some(bounds) = group.geometry_bounds() else {
            continue;
        };
        let position = [bounds.min_x, bounds.min_y];
        let target = if flowsheet.instances.iter().any(|item| &item.id == name) {
            Some(&mut layout.instances)
        } else if flowsheet.inputs.iter().any(|input| &input.name == name) {
            Some(&mut layout.inputs)
        } else if flowsheet.products.iter().any(|product| product == name) {
            Some(&mut layout.products)
        } else {
            None
        };
        let Some(target) = target else {
            continue;
        };
        target
            .entry(name.clone())
            .and_modify(|current| {
                current[0] = current[0].min(position[0]);
                current[1] = current[1].min(position[1]);
            })
            .or_insert(position);
        placed += 1;
    }
    if placed > 0 {
        flowsheet.layout = Some(layout);
    }
    placed
}
