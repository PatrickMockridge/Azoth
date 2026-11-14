//! DEXPI 2.0 exchange for azoth: a flowsheet as a DEXPI Process model, and back.
//!
//! The crate is two layers. [`model`], [`reader`] and [`writer`] are the DEXPI XML metamodel -
//! `Model`, `Import`, `Object`, `Components`, `Data` and `References` - and nothing about
//! process engineering. [`process`] is the projection between that tree and a
//! [`azoth_process::Flowsheet`], and [`mapping`] is the one table that says which azoth unit
//! operation is which DEXPI Process step.
//!
//! **What this crate does not claim.** A DEXPI file that validates is not a DEXPI EV
//! certificate, and schema conformance is not a named-tool round trip. [`ExportReport`]
//! carries both statements the way NeqSim's own conformance report does, so a caller cannot
//! read "conformant" as either.

pub mod diagnostic;
pub mod ids;
pub mod inventory;
pub mod mapping;
pub mod model;
pub mod process;
pub mod profile;
pub mod reader;
pub mod units;
pub mod validate;
pub mod writer;
pub mod xsd;

mod process_mapping_gen;

pub use diagnostic::{DexpiDiagnostic, Severity};
pub use mapping::{Disposition, UnitMapping};
pub use model::{Child, Components, Data, DataValue, Document, Import, Model, Object, References};
pub use process::export::{ExportOptions, ExportReport, export_process};
pub use process::import::{ImportOptions, ImportReport, import_process};
pub use profile::Profile;
