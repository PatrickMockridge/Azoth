//! A keycard overlay, as a value crossing into the extension.
//!
//! `azoth_eos::databank::Overlay` is the Rust type; this is how Python builds one and
//! reads back what it resolves to. Nothing here decides anything - the merge rule lives
//! in `azoth-eos`, and this module is a window onto it.
//!
//! **Why a window rather than a parameter.** `python/src/azoth/_rust_bridge.py` already
//! resolves names to numbers before crossing, so the model functions have no need of an
//! overlay and their signatures are untouched. What the overlay *does* need to be is
//! **readable from Python**, because the merge rule - a card naming only `omega` keeps
//! the shipped `Tc` and `Pc` - is implemented twice, once in
//! `python/src/azoth/eos/components.py` and once in `azoth-eos::databank`, and two
//! implementations of one rule are only worth having if something compares them.
//! `python/tests/test_data_agreement.py` is that comparison.

use pyo3::prelude::*;

use azoth_eos::databank::{ComponentOverride, Overlay};

use crate::data::{PyComponentRow, row_of};

/// One substance as a card states it: `(name, Tc, Pc, omega)`, each parameter optional
/// because a card naming only `omega` keeps the shipped `Tc` and `Pc`.
pub type ComponentArguments = (String, Option<f64>, Option<f64>, Option<f64>);

/// A keycard's sections, as this crate can read them.
///
/// Built either entry by entry, from the values `azoth.keycard` resolved to, or by
/// reading a card's text with [`card_overlay`].
/// Opaque on purpose: a card crosses as a value and *what it resolves to* is read
/// through the functions below, so there is nothing about it Python can inspect
/// directly and nothing here for the stub to describe.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "Overlay")]
pub struct PyOverlay {
    inner: Overlay,
}

impl PyOverlay {
    /// The Rust overlay this wraps.
    pub(crate) fn as_overlay(&self) -> &Overlay {
        &self.inner
    }
}

/// A keycard, built from the values `azoth.keycard` resolved to.
///
/// A function rather than a builder object, so the card is a value the caller holds
/// from the moment it exists: there is no half-built overlay, and no way to keep one
/// and add to it after it has been used. A keycard is immutable once read, and this is
/// the shape that says so.
///
/// `components` is `(name, Tc, Pc, omega)` per substance, each parameter optional
/// because a card naming only `omega` keeps the shipped `Tc` and `Pc`. `kij` is
/// `(first, second, value)` per pair, and a value of exactly zero is a caller stating
/// ideal mixing rather than a card saying nothing.
///
/// # Errors
/// * `InvalidInputError` if a `kij` pair names one substance twice.
#[pyfunction]
#[pyo3(signature = (components, kij))]
pub fn overlay(
    py: Python<'_>,
    components: Vec<ComponentArguments>,
    kij: Vec<(String, String, f64)>,
) -> PyResult<PyOverlay> {
    let mut inner = Overlay::new();
    for (name, tc, pc, omega) in components {
        inner.set_component(&name, ComponentOverride { tc, pc, omega });
    }
    for (first, second, value) in kij {
        inner
            .set_kij(&first, &second, value)
            .map_err(|error| crate::errors::to_pyerr(py, error))?;
    }
    Ok(PyOverlay { inner })
}

/// One substance as a card resolves it, in the table's own row shape.
///
/// # Errors
/// * `PropertyUnavailableError` if neither the databank nor the card carries the name,
///   or if the card adds one without every parameter a cubic reads.
#[pyfunction]
pub fn overlay_entry_row(
    py: Python<'_>,
    name: &str,
    overlay: &PyOverlay,
) -> PyResult<PyComponentRow> {
    let entry = azoth_eos::databank::entry(name, Some(overlay.as_overlay()))
        .map_err(|error| crate::errors::to_pyerr(py, error))?;
    Ok(row_of(&entry))
}

/// Every substance a card *names*, resolved, in name order.
///
/// Only the card's own names, because those are the ones the baseline does not answer
/// for: comparing the whole table would compare the file against itself. A card naming a
/// substance that is also shipped resolves it exactly as `overlay_entry_row` does.
///
/// # Errors
/// * `PropertyUnavailableError` as [`overlay_entry_row`].
#[pyfunction]
pub fn overlay_component_rows(
    py: Python<'_>,
    overlay: &PyOverlay,
) -> PyResult<Vec<PyComponentRow>> {
    let mut names = overlay.as_overlay().component_names();
    names.sort_unstable();
    names
        .into_iter()
        .map(|name| overlay_entry_row(py, name, overlay))
        .collect()
}

/// Every interaction pair a card states, resolved, each pair once, in name order.
///
/// Resolved rather than as written, so a comparison reads the number a calculation
/// would actually use. That is how a card's zero winning over a fitted value is visible
/// rather than merely stored.
#[pyfunction]
pub fn overlay_kij_rows(overlay: &PyOverlay) -> Vec<(String, String, f64)> {
    overlay
        .as_overlay()
        .kij_pairs()
        .into_iter()
        .map(|(first, second)| {
            // The *effective* value, not the stated one: a card's zero overriding a
            // fitted parameter is the rule worth comparing, and a row carrying only what
            // the card wrote could not show it.
            let value = azoth_eos::databank::kij(&first, &second, Some(overlay.as_overlay()));
            (first, second, value)
        })
        .collect()
}

/// A card read by **Rust**, as the overlay it resolves to.
///
/// The second reader of one document. `python/tests/test_card_agreement.py` hands one
/// card's text to both this and `azoth.keycard`, and compares what each says - one case
/// per rule, on the shipped template. That is the arrangement
/// `python/tests/test_data_agreement.py` already uses for the merge rule, and the reason
/// is the same: two implementations of one document are only worth having if something
/// holds them to each other on real data.
///
/// # Errors
/// * `InvalidInputError` if the text is not a card this build reads, naming the section,
///   parameter, unit or model choice that was refused.
#[pyfunction]
pub fn card_overlay(py: Python<'_>, text: &str) -> PyResult<PyOverlay> {
    let card = azoth_eos::card::Card::from_toml(text)
        .map_err(|error| crate::errors::to_pyerr(py, error))?;
    Ok(PyOverlay {
        inner: card.overlay().clone(),
    })
}
