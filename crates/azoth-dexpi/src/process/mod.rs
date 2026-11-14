//! The projection between a flowsheet and a DEXPI 2.0 Process document.
//!
//! [`export`] writes a flowsheet - and, unless the exchange is topology only, its run - as a
//! Process model. [`import`] reads one back as a scaffold, naming everything it could not
//! recover. Neither direction invents a value: a quantity the run did not compute is omitted,
//! a fluid DEXPI does not carry is left undeclared, and a unit with no reviewed mapping is
//! refused by name.

pub mod export;
pub mod import;

/// The extension property names the Process projection writes; they live at the crate root
/// because the Plant projection writes the same set.
pub use crate::extensions;
