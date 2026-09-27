//! The provenance block, as the Rust side computes it, exposed for comparison.

use azoth_core::provenance::{Provenance, ProvenanceStatic};
use azoth_core::warning::{Warning, WarningCode};
use pyo3::PyResult;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::pyfunction;

/// One call's block, merged by this side, for the cross-language test.
///
/// **The comparison the feature's claim rests on.** Both languages decide independently what
/// a block says - which checks were skipped and in what order, which warning codes were
/// raised, whether the call was clean - and until this existed nothing compared them. The
/// *tables* they read are held equal by id, but the arithmetic over them was not, and
/// arithmetic written twice is arithmetic that can disagree.
///
/// Returns the four fields the merge computes rather than the whole block: the static half
/// comes from the generated table, which the same test compares separately, and re-sending it
/// would be comparing the same bytes twice.
///
/// The warnings cross as `(code, field)` pairs and the message is dropped, because the merge
/// reads no message - a message is for a reader, and the block carries codes.
///
/// # Errors
/// [`PyValueError`] for a code this build does not know. `None` rather than an error for an
/// id it does not ship, which is what the by-id lookup means.
#[pyfunction]
#[pyo3(name = "provenance_block")]
pub fn provenance_block(
    calc_id: &str,
    warnings: Vec<(String, Option<String>)>,
) -> PyResult<Option<(String, Vec<String>, Vec<String>, bool)>> {
    let Some(static_half): Option<&'static ProvenanceStatic> =
        azoth_core::provenance_gen::provenance(calc_id)
    else {
        return Ok(None);
    };

    let mut parsed = Vec::with_capacity(warnings.len());
    for (code, field) in warnings {
        let Some(code) = WarningCode::parse(&code) else {
            return Err(PyValueError::new_err(format!(
                "{code:?} is not a warning code this build knows"
            )));
        };
        parsed.push(Warning {
            code,
            message: String::new(),
            field,
        });
    }

    let block = Provenance::of(*static_half, &parsed);
    Ok(Some((
        block.static_half.calc_id.to_owned(),
        block
            .skipped_checks
            .iter()
            .map(|s| (*s).to_owned())
            .collect(),
        block
            .warning_codes
            .iter()
            .map(|s| (*s).to_owned())
            .collect(),
        block.clean,
    )))
}
