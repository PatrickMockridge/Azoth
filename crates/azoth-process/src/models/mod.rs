//! The registry surface: one module per unit operation, wrapping its kernel.
//!
//! A kernel ([`crate::kernels`]) is the flowsheet's arithmetic — it reads and writes
//! [`crate::stream::Stream`] values, which is what a connection carries. A model is the
//! same unit operation as an id the rest of the library can address: a declared input
//! list, a [`CalcResult`](azoth_core::CalcResult) impl, a cross-implementation test and a
//! NeqSim oracle.
//!
//! **Two surfaces, one arithmetic**, and the split is deliberate: the executor P12 will
//! build calls kernels, while a case, a cross-impl test and an oracle capture call models.
//! Neither wraps the other's types, so a change to the flowsheet's stream record cannot
//! silently move a published id.
//!
//! # How a port crosses
//!
//! A port carries the record `n, z, P, T, h`. The fluid's `components` is shared by every
//! port of a unit operation, so it is declared once; each port then contributes
//! `<port>_n`, `<port>_z`, `<port>_p`, `<port>_t`, and — on an outlet — `<port>_h`.
//!
//! **An inlet takes four of the record's five fields, and that is stated rather than
//! assumed.** `h` is a state function of `(T, P, z)`: the same flash that `Stream::from_pt`
//! runs decides it. Accepting an `h` as well would let a case hand over a state that does
//! not exist, so the inlet's enthalpy is computed and an outlet's is reported.
//!
//! **A `many` port crosses as a matrix and vectors.** `z` becomes one row per stream and
//! each scalar field one entry per stream, which is the port's multiplicity written out —
//! the calculus makes the multiplicity part of the declaration an implementation must
//! agree with, so folding a mixer to a fixed two inlets would answer a smaller question.
//! [`splitter`] is the first id with such a port, and its result is the shape in full:
//! `products_n`, `products_p`, `products_t` and `products_h` are vectors with one entry
//! per outlet and `products_z` a matrix with one row per outlet.

pub mod heat_exchanger;
pub mod mixer;
pub mod pump;
pub mod separator;
pub mod splitter;
pub mod throttling_valve;

pub use heat_exchanger::{HeatExchangerResult, heat_exchanger};
pub use mixer::{MixerResult, mixer};
pub use pump::{PumpResult, pump};
pub use separator::{SeparatorResult, separator};
pub use splitter::{SplitterResult, splitter};
pub use throttling_valve::{ThrottlingValveResult, throttling_valve};
