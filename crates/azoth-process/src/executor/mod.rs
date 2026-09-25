//! The executor: turning a flowsheet into something that runs.
//!
//! **It is a port, and `ProcessSystem.runSequential` + `Recycle` are the source.**
//! `ProcessSystem.java` at the pin is 10,522 lines and most of it is not the executor - the
//! DEXPI/JSON/graphviz exporters, `RunStatus`, threading, the field-development hooks. What is
//! ported is the **sequential path the class calls legacy**, plus `Recycle` and the slice of
//! `RecycleController` that path uses; the rest is named in `ROADMAP.md` as not ported, each with
//! the class that would close it.
//!
//! The layers, in the order the middleware's gap table names them:
//!
//! * [`dispatch`] - the `unit_ops.*` id to kernel table, with typed parameter coercion.
//! * `order` - the execution order, which is `getTopologicalOrder()` when
//!   `useGraphBasedExecution` and insertion order otherwise.
//! * `recycle` - the tear's fixed point: `Recycle`'s four tolerances, its `solved()` and its
//!   three acceleration methods, driven by the outer loop's hundred iterations.
//! * `session` - the named results a widget and an agent both point at.
//!
//! **A recycle stays a declared tear** (`[[recycles]]` in the flowsheet), which is the calculus's
//! "feedback is restriction" and the checker's `UnrecycledLoop` rule. The declaration carries
//! `Recycle`'s own convergence parameters so that the class's numbers are the tear's rather than
//! invented.

pub mod dispatch;
pub mod json;
pub mod session;

pub use dispatch::{DISPATCH, Kernel, Parameters, UNRUNNABLE, dispatch, kernel_for};
pub use json::{Quantity, SessionReport, StreamRecord, TearReport, to_json};
pub use session::{MAX_PASSES, RunReport, STREAM_FIELDS, Session, TearRecord, run};
