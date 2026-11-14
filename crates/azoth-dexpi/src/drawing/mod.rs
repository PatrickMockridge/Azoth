//! The drawing half of the exchange: a neutral, resolved figure and its DEXPI 2.0 form.
//!
//! A drawing is **presentation, not process state**. It links each group to the conceptual
//! object it represents, so a viewer can colour what a model would take, but no simulation
//! reads a coordinate. The neutral model is in [`model`]; [`native`] reads and writes
//! `Core/Diagram`.

pub mod model;
pub mod native;

pub use model::{
    Bounds, Color, Drawing, DrawingGroup, DrawingOrigin, DrawingReport, GroupKind, Point,
    Primitive, Stroke,
};
