//! The process layer: unit operations and flowsheets, on typed, directional
//! channels.
//!
//! A unit operation is a process with inlets and outlets; a channel carries a
//! field record (names → dimensions) with a polarity; conservation is linearity.
//! This crate states the schema — the palette of unit operations and the
//! flowsheet that wires them — and checks that a flowsheet honours the calculus's
//! rules. The executor that turns a flowsheet into something that runs is tranche
//! P12, and is not built.

pub mod channel;
pub mod check;
pub mod flowsheet;
pub mod kernels;
pub mod load;
pub mod model_gen;
pub mod models;
pub mod stream;
pub mod unit_op;

pub use channel::{Direction, FieldType, Multiplicity, Port, Shape};
pub use check::{Diagnostic, validate, validate_palette};
pub use flowsheet::{Connection, Flowsheet, Instance, Recycle};
pub use load::{load_palette, parse_flowsheet};
// **The kernels stay under `kernels` and the ids take the flat names**, which is the split
// the Python package makes too. Two call shapes under one name is not a naming problem to
// work around: a kernel takes and returns a `Stream`, which is what a flowsheet's
// connection carries, while a model takes the record field by field so that a case, a
// cross-impl test and a NeqSim capture can address it.
pub use models::{
    MixerResult, PumpResult, SeparatorResult, SplitterResult, mixer, pump, separator, splitter,
};
pub use stream::Stream;
pub use unit_op::{Param, Source, UnitOpSpec};
