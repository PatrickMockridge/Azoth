//! What a run rests on, as a front end receives it.
//!
//! `azoth_core::provenance` is the value and does not serialise, for the reason
//! `azoth_core::check`'s `Diagnostic` does not: the JSON form of a thing is the front end's
//! business, and this is that form. Every key is always present - `null` rather than a
//! missing key - because a reader switching on `verification` should not also have to test
//! for the field's existence.
//!
//! **On the envelope rather than inside the session report.** `SessionReport` is the one
//! stream encoder and its published bytes are pinned by tests; what code a document rests on
//! is a fact about the *run* and not about any stream, so putting it there would move a
//! document that has nothing to do with it.

use azoth_core::provenance_gen;
use serde::Serialize;

use crate::flowsheet::Instance;
use crate::middleware::form::form;
use crate::middleware::session::Workspace;
use crate::unit_op::UnitOpSpec;

/// One instantiated unit operation, and what backs it.
#[derive(Debug, Clone, Serialize)]
pub struct UnitProvenance {
    /// The document's own name for this instance, as `Envelope::paths` spells it.
    pub instance: String,
    /// The palette entry it was added as, e.g. `unit_ops.pump`.
    pub unit: String,
    /// The registered model id, where the palette entry has one.
    pub model: Option<&'static str>,
    /// The provenance line the palette carries, e.g. the NeqSim class.
    pub source: Option<String>,
    /// Whether the executor has a kernel for it, and why not where it does not.
    pub runnable: bool,
    pub refusal: Option<&'static str>,
    /// The model's spec and Rust kernel hashes, and how far its answer is checked.
    ///
    /// Absent together, and absent rather than defaulted: a palette entry with no model has
    /// no calculation to describe, and `null` says so where an empty string would look like a
    /// hash nobody could compare.
    pub spec_sha256: Option<&'static str>,
    pub rust_sha256: Option<&'static str>,
    pub verification: Option<&'static str>,
}

/// What a document as a whole rests on.
#[derive(Debug, Clone, Serialize)]
pub struct ProvenanceRecord {
    /// The library that produced the answer, and the version this binary reports.
    pub library: &'static str,
    pub version: &'static str,
    /// One entry per instantiated unit operation, in the document's own order.
    ///
    /// The document's order rather than a sort: the envelope is what a canvas diffs between
    /// edits, and a list that reordered itself would read as a change that did not happen.
    pub units: Vec<UnitProvenance>,
}

impl ProvenanceRecord {
    /// Whether every unit the document instantiates can say what backs it.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.units
            .iter()
            .all(|unit| unit.model.is_none() || unit.spec_sha256.is_some())
    }
}

/// One instance's block, or its own account of why it has none.
fn unit_provenance(instance: &Instance, palette: &[UnitOpSpec]) -> UnitProvenance {
    let entry = palette.iter().find(|spec| spec.id == instance.unit);
    let Some(spec) = entry else {
        // An instance naming a palette entry this build does not carry. The checker reports
        // it as a diagnostic; here it is simply a unit with nothing behind it.
        return UnitProvenance {
            instance: instance.id.clone(),
            unit: instance.unit.clone(),
            model: None,
            source: None,
            runnable: false,
            refusal: Some("not a palette entry this build carries"),
            spec_sha256: None,
            rust_sha256: None,
            verification: None,
        };
    };

    let shape = form(spec);
    let block = shape.model.and_then(provenance_gen::provenance);
    UnitProvenance {
        instance: instance.id.clone(),
        unit: instance.unit.clone(),
        model: shape.model,
        source: shape.source,
        runnable: shape.runnable,
        refusal: shape.refusal,
        spec_sha256: block.map(|entry| entry.spec_sha256),
        rust_sha256: block.map(|entry| entry.rust_sha256),
        verification: block.map(|entry| entry.verification.as_str()),
    }
}

/// The run's provenance: the library, and every unit the document instantiates.
#[must_use]
pub fn session_provenance(workspace: &Workspace) -> ProvenanceRecord {
    ProvenanceRecord {
        library: "azoth",
        version: env!("CARGO_PKG_VERSION"),
        units: workspace
            .flowsheet()
            .instances
            .iter()
            .map(|instance| unit_provenance(instance, workspace.palette()))
            .collect(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn a_document_with_no_instances_is_complete_and_empty() {
        let record = ProvenanceRecord {
            library: "azoth",
            version: "0.1.0",
            units: Vec::new(),
        };
        assert!(record.is_complete());
    }

    #[test]
    fn a_unit_with_no_model_does_not_make_a_record_incomplete() {
        // A boundary entry - a feed or a product - has no calculation behind it. Treating
        // that as incompleteness would make every document incomplete and the flag useless.
        let record = ProvenanceRecord {
            library: "azoth",
            version: "0.1.0",
            units: vec![UnitProvenance {
                instance: "feed".to_string(),
                unit: "unit_ops.feed".to_string(),
                model: None,
                source: None,
                runnable: true,
                refusal: None,
                spec_sha256: None,
                rust_sha256: None,
                verification: None,
            }],
        };
        assert!(record.is_complete());
    }
}
