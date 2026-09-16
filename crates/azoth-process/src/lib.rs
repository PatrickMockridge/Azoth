//! The process layer: unit operations and flowsheets, on typed, directional
//! channels.
//!
//! A unit operation is a process with inlets and outlets; a channel carries a
//! field record (names → dimensions) with a polarity; conservation is linearity.
//! This crate states the schema — the palette of unit operations and the
//! flowsheet that wires them — and checks that a flowsheet honours the calculus's
//! rules. The executor, and the kernels that give each unit op its arithmetic,
//! are a later tranche; what lives here now is the declaration and the check.

pub mod channel;
pub mod check;
pub mod flowsheet;
pub mod kernels;
pub mod stream;
pub mod unit_op;

pub use channel::{Direction, FieldType, Multiplicity, Port, Shape};
pub use check::{Diagnostic, validate, validate_palette};
pub use flowsheet::{Connection, Flowsheet, Instance, Recycle};
pub use kernels::{mixer, separator, splitter};
pub use stream::Stream;
pub use unit_op::{Param, Source, UnitOpSpec};
