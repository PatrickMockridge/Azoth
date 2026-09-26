//! The interoperation surface: what a front-end reads and writes.
//!
//! `docs/src/architecture/middleware.md` states what sits between the schema and the things
//! that drive it - a flowsheet editor, a notebook, an agent. This module is that surface, and
//! every one of those is a **binding** of it rather than a second implementation: the session
//! is [`session::Workspace`], the graph is [`graph`], a unit op's declaration is [`form`], an
//! edit is a [`command::Command`], and one call's answer is the [`envelope::Envelope`].
//!
//! **Nothing here is a second rule set.** The checker decides what a document may say, the
//! executor decides what it computes, and `executor::json` writes the values; this layer
//! projects those and adds no verdict of its own.

pub mod command;
pub mod diagnostic;
pub mod envelope;
pub mod form;
pub mod graph;
pub mod session;
pub mod tools;

use azoth_core::Result;
use azoth_core::unit_vocab_gen::{DISPLAY_UNITS, UNIT_NAMES, UNIT_SETS, si_factor};

use crate::middleware::form::dimension_id;
use crate::unit_op::UnitOpSpec;

/// The palette as a front-end reads it: one form per entry, and optionally the agent's tools.
///
/// **One assembly, so a binding cannot spell the document differently.** The CLI and the wasm
/// module differ in where the palette comes from and in nothing else: both call this with the
/// specs they loaded, and `tools` is `false` for a widget generator and `true` for an agent's
/// tool list - the same declaration read two ways rather than two documents.
///
/// The units and the unit sets ride here too, because they are the same kind of fact as the
/// palette: something about the *library* rather than about a document, which a front end has
/// nowhere else to ask for. A document says what a stream is worth in the unit the codec chose;
/// this says what that unit is, and which unit a reader would rather see it in.
///
/// # Errors
/// A write that fails, which for this shape it cannot.
pub fn catalogue(palette: &[UnitOpSpec], with_tools: bool) -> Result<String> {
    let mut document = serde_json::json!({
        "unit_ops": form::forms(palette),
        "units": units(),
        "unit_sets": unit_sets(),
    });
    if with_tools {
        let object = document
            .as_object_mut()
            .expect("the object this function just built");
        object.insert(
            "tools".to_string(),
            serde_json::json!(tools::tools(palette)),
        );
    }
    serde_json::to_string_pretty(&document).map_err(|error| {
        azoth_core::AzothError::invalid_input(
            "catalogue",
            format!("it could not be written: {error}"),
        )
    })
}

/// Every canonical unit **and every display unit**, keyed by the string a spec and a stream
/// record spell: its dimension, the SI base magnitude one of it is worth, and the affine
/// constant where it has one.
///
/// **The factor is computed and not written down.** `si_factor` runs the same conversion a
/// calculation runs, so a front end converting for display and a calculation converting on the
/// way in reach one function - and the number a reader is shown is the number `pint` was used
/// to check, not a copy of it. A **display unit** is the one case where the table carries a
/// number, because an offset is a definition neither library exposes as data; the conversion
/// itself is `affine_si`, held to `pint` by the same test.
///
/// Both lists land in one map because that is what a reader looks a unit up in: a set names a
/// unit's id and a quantity carries one, and the front end has no way to know - and no reason
/// to care - which table the name came from.
fn units() -> serde_json::Value {
    let mut out = serde_json::Map::new();
    for name in UNIT_NAMES {
        out.insert(
            (*name).to_string(),
            serde_json::json!({
                "dimension": dimension_id(name),
                "factor": si_factor(name),
                // A scale is the same expression with no shift, so the front end has one
                // formula rather than two branches.
                "offset": 0.0,
            }),
        );
    }
    for unit in DISPLAY_UNITS {
        out.insert(
            unit.id.to_string(),
            serde_json::json!({
                "dimension": unit.dimension,
                "factor": unit.factor,
                "offset": unit.offset,
            }),
        );
    }
    serde_json::Value::Object(out)
}

/// The named unit sets, as the vocabulary declares them: a unit per dimension, per set.
///
/// `serde_json::json!` cannot key an object by a run-time string, so these two maps are built
/// rather than written - which is also the shape that keeps the wire a projection: nothing here
/// decides what a set contains.
fn unit_sets() -> serde_json::Value {
    serde_json::Value::Array(
        UNIT_SETS
            .iter()
            .map(|set| {
                serde_json::json!({
                    "id": set.id,
                    "name": set.name,
                    "units": set
                        .units
                        .iter()
                        .map(|(dimension, unit)| ((*dimension).to_string(), serde_json::json!(unit)))
                        .collect::<serde_json::Map<String, serde_json::Value>>(),
                })
            })
            .collect(),
    )
}
