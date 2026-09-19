//! A keycard overlay, as a value crossing into the extension.
//!
//! `azoth_eos::databank::Overlay` is the Rust type; this is how Python builds one and
//! reads back what it resolves to. Nothing here decides anything - the merge rule lives
//! in `azoth-eos`, and this module is a window onto it. The rule is implemented twice,
//! here and in `python/src/azoth/eos/components.py`, and
//! `python/tests/test_card_agreement.py` compares them.

use pyo3::prelude::*;

use azoth_eos::Cubic;
use azoth_eos::card::CoefficientValue;
use azoth_eos::databank::{ComponentOverride, Overlay};

use crate::data::{PyComponentRow, row_of};

/// One substance as a card states it: `(name, Tc, Pc, omega)`, each parameter optional
/// because a card naming only `omega` keeps the shipped `Tc` and `Pc`.
pub type ComponentArguments = (String, Option<f64>, Option<f64>, Option<f64>);

/// A keycard's sections, as this crate can read them. Opaque: a card crosses as a value,
/// and *what it resolves to* is read through the functions below.
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

/// A keycard, built from the values `azoth.keycard` resolved to. A function rather than a
/// builder, so the card is a value from the moment it exists.
///
/// `components` is `(name, Tc, Pc, omega)` per substance, each optional because a card
/// naming only `omega` keeps the shipped `Tc` and `Pc`. `kij` is
/// `(first, second, value)` per pair, where exactly zero is a caller stating ideal mixing.
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
        inner.set_component(
            &name,
            ComponentOverride {
                tc,
                pc,
                omega,
                ..Default::default()
            },
        );
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

/// Every substance a card names, resolved, in name order.
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
            let value =
                azoth_eos::databank::kij(&first, &second, Cubic::Pr, Some(overlay.as_overlay()));
            (first, second, value)
        })
        .collect()
}

/// Every *associating* interaction pair a card states, resolved, each pair once.
///
/// The sibling of [`overlay_kij_rows`] for the other column, and resolved for the same
/// reason: what is compared is the value an associating mixing rule would read rather
/// than what the card wrote, so a card's override - including a zero - is visible.
///
/// # Errors
/// * `InvalidInputError` if `family` is not `"srk"` or `"pr"`.
#[pyfunction]
pub fn overlay_cpa_kij_rows(
    py: Python<'_>,
    overlay: &PyOverlay,
    family: &str,
) -> PyResult<Vec<(String, String, f64)>> {
    let cubic = match family {
        "srk" => azoth_eos::association::AssociationCubic::Srk,
        "pr" => azoth_eos::association::AssociationCubic::Pr,
        other => {
            return Err(crate::errors::to_pyerr(
                py,
                azoth_core::AzothError::invalid_input(
                    "family",
                    format!(
                        "{other:?} is not an associating cubic family; expected \"srk\" or \"pr\""
                    ),
                ),
            ));
        }
    };
    Ok(overlay
        .as_overlay()
        .kij_pairs()
        .into_iter()
        .map(|(first, second)| {
            // A pair resolved through the pair-wise entry point, which is what a
            // `Mixture`'s own reduced parameters call; the two names are the whole list.
            let names = [first.as_str(), second.as_str()];
            let value = azoth_eos::databank::cpa_kij(&names, cubic, Some(overlay.as_overlay()))[1];
            (first, second, value)
        })
        .collect())
}

/// One coefficient as it crosses: `(calc_id, name, unit, rows, cols, values)`.
///
/// A tuple rather than a class because nothing on this side reads it - the Python side
/// unpacks it into a comparison - and `values` is row-major beside the shape it belongs
/// to rather than nested, because the boundary carries numbers.
type CoefficientRow = (String, String, String, usize, usize, Vec<f64>);

/// Every coefficient a card states, in a shape a comparison can read.
///
/// **This exists because nothing compared the two readers on coefficients.** The card
/// format's own test compares `components` and `kij`; a coefficient was parsed by both
/// and checked by neither, which is how a value one reader accepts and the other
/// refuses stays invisible. It was a scalar until a card had to carry a matrix, so the
/// gap was harmless until it was not.
///
/// The value keeps its shape: `rows` and `cols` are `1` and `1` for a number, `1` and
/// `n` for a vector, and the matrix's own size otherwise, with `values` flattened
/// row-major beside them. Flattened because the boundary carries numbers; shaped because
/// a matrix read as a vector is the defect this is here to catch.
///
/// The value is the **SI** one, `Coefficient::si_value`'s, not the number as written -
/// the same choice `overlay_kij_rows` makes, and for the same reason: what two readers
/// have to agree about is the number a calculation would use.
///
/// # Errors
/// * `InvalidInputError` if the text is not a card this build reads.
#[pyfunction]
pub fn card_coefficients(py: Python<'_>, text: &str) -> PyResult<Vec<CoefficientRow>> {
    let card = azoth_eos::card::Card::from_toml(text)
        .map_err(|error| crate::errors::to_pyerr(py, error))?;
    let mut out = Vec::new();
    for (calc_id, arguments) in &card.coefficients {
        for (name, body) in arguments {
            let (rows, cols, values) = flatten(&body.si_value());
            out.push((
                calc_id.clone(),
                name.clone(),
                body.unit.clone(),
                rows,
                cols,
                values,
            ));
        }
    }
    Ok(out)
}

/// A coefficient value as `(rows, cols, values row-major)`, which is how it crosses.
fn flatten(value: &CoefficientValue) -> (usize, usize, Vec<f64>) {
    match value {
        CoefficientValue::Scalar(number) => (1, 1, vec![*number]),
        CoefficientValue::List(entries) => match entries.first() {
            Some(CoefficientValue::List(row)) => {
                let cols = row.len();
                let values = entries
                    .iter()
                    .flat_map(|entry| match entry {
                        CoefficientValue::List(row) => {
                            row.iter().filter_map(CoefficientValue::scalar).collect()
                        }
                        CoefficientValue::Scalar(_) => Vec::new(),
                    })
                    .collect();
                (entries.len(), cols, values)
            }
            _ => (
                1,
                entries.len(),
                entries
                    .iter()
                    .filter_map(CoefficientValue::scalar)
                    .collect(),
            ),
        },
    }
}

/// A card read by **Rust**, as the overlay it resolves to. The second reader of one
/// document; `python/tests/test_card_agreement.py` holds it to `azoth.keycard`.
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
