//! The result objects the extension hands back to Python.
//!
//! These are a *transport* representation, not the Python API. The adapter in
//! ``azoth._rust_bridge`` converts each one into the result dataclass the pure
//! Python implementation returns, so a caller cannot tell which backend
//! answered - and, more to the point, cannot accidentally depend on the
//! difference.
//!
//! Field names here must match both the Rust struct's `CalcResult::FIELDS` and
//! the Python dataclass. Three-way agreement is asserted by
//! `python/tests/test_cross_impl.py::test_result_shapes_agree_across_languages`,
//! because a mismatch is invisible to every numerical test: the numbers agree
//! perfectly and only the attribute name differs.

use azoth_core::solver::SolverKind;
use azoth_core::units::UNIT_NAMES;
use azoth_core::warning::{Warning, WarningCode};
use pyo3::prelude::*;

/// A physical quantity crossing the FFI boundary.
///
/// `magnitude_si` is always the SI base value, and `unit` is the canonical
/// display unit. Keeping both means the Python side can rebuild a genuine pint
/// quantity rather than being handed a bare float and having to remember what it
/// was in.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "Qty")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyQty {
    /// The magnitude, in SI base units.
    #[pyo3(get)]
    pub magnitude_si: f64,
    /// The canonical unit string, e.g. `Pa`.
    #[pyo3(get)]
    pub unit: String,
}

#[pymethods]
impl PyQty {
    #[new]
    fn new(magnitude_si: f64, unit: String) -> Self {
        Self { magnitude_si, unit }
    }

    fn __repr__(&self) -> String {
        format!("Qty({} {})", self.magnitude_si, self.unit)
    }

    /// Always raises.
    ///
    /// A dimensioned value must never be silently coerced to a bare number: that
    /// is exactly how a pressure in bar becomes a pressure in pascal with
    /// nothing to show for it. Use `.magnitude_si` and say what unit you mean.
    fn __float__(&self) -> PyResult<f64> {
        Err(pyo3::exceptions::PyTypeError::new_err(format!(
            "refusing to convert a Qty to float: it carries the unit '{}'. \
             Use .magnitude_si to take the SI value explicitly.",
            self.unit
        )))
    }
}

/// A warning, transported.
#[pyclass(frozen, skip_from_py_object, module = "azoth._core", name = "Warning")]
#[derive(Debug, Clone, PartialEq)]
pub struct PyWarning {
    /// The machine-readable code, as its SCREAMING_SNAKE_CASE string.
    #[pyo3(get)]
    pub code: String,
    /// Human-readable explanation.
    #[pyo3(get)]
    pub message: String,
    /// The input or output it is about, if any.
    #[pyo3(get)]
    pub field: Option<String>,
}

#[pymethods]
impl PyWarning {
    fn __repr__(&self) -> String {
        match &self.field {
            Some(field) => format!("Warning({}, {field}: {})", self.code, self.message),
            None => format!("Warning({}, {})", self.code, self.message),
        }
    }
}

impl From<&Warning> for PyWarning {
    fn from(w: &Warning) -> Self {
        Self {
            code: w.code.as_str().to_string(),
            message: w.message.clone(),
            field: w.field.clone(),
        }
    }
}

pub(crate) fn transport(warnings: &[Warning]) -> Vec<PyWarning> {
    warnings.iter().map(PyWarning::from).collect()
}

/// One fitting's contribution, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "KComponent"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyKComponent {
    /// Fitting id as it appears in the registry.
    #[pyo3(get)]
    pub fitting_id: String,
    /// Equivalent length ratio (L_eq / D).
    #[pyo3(get)]
    pub n_ld: f64,
    /// This fitting's resistance coefficient.
    #[pyo3(get)]
    pub k: f64,
}

#[pymethods]
impl PyKComponent {
    fn __repr__(&self) -> String {
        format!(
            "KComponent({}, n_ld={}, k={})",
            self.fitting_id, self.n_ld, self.k
        )
    }
}

/// Every warning code the Rust side can emit.
///
/// Exists so the Python suite can assert the two code sets are identical
/// without parsing Rust source.
#[pyfunction]
#[must_use]
pub fn warning_codes() -> Vec<String> {
    WarningCode::all()
        .iter()
        .map(|c| c.as_str().to_string())
        .collect()
}

/// Every canonical unit string the spec schema permits.
///
/// The same job `warning_codes` does for the warning vocabulary, and it exists for
/// the same reason: the schema, `azoth.core.units.CANONICAL_UNITS` and this crate's
/// `UNIT_NAMES` all have to name one set, and a Python test can only assert that if
/// the third list is reachable from Python. Without this function the Rust half of
/// that contract would be unverifiable - which is exactly how `K` sat in the
/// vocabulary for the whole life of the project with no Rust conversion behind it.
#[pyfunction]
#[must_use]
pub fn unit_names() -> Vec<String> {
    UNIT_NAMES.iter().map(|name| (*name).to_string()).collect()
}

/// The dimension of one canonical unit, as exponents in `SLOTS` order.
///
/// `None` for a name not in the vocabulary, rather than an error: a caller asking
/// about a name that is not there is a question about the vocabulary, and the
/// answer is that it has no dimension because it is not in it.
///
/// This exists so the table's exponents can be compared from Python against the
/// ones the Rust side was generated with. The compile-time assertion in
/// `azoth_core::unit_vocab_gen` catches a wrong exponent only when it makes the
/// `uom` quantity differ, and two dimensions can share a quantity type - so this is
/// the check that sees all seven slots.
#[pyfunction]
#[must_use]
pub fn unit_dimensions(name: &str) -> Option<Vec<i8>> {
    azoth_core::unit_vocab_gen::dimension(name).map(|d| d.to_vec())
}

/// How many SI base units one of `name` is worth, by running the conversion a
/// calculation runs.
///
/// Not a stored number: this calls the generated conversion table, which is the
/// same function `to_si`'s Rust counterpart reaches. What it is compared against
/// is `pint`'s own answer for the same unit name, in
/// `python/tests/test_units_cross_library.py` - so the two units libraries are
/// compared and neither is restated.
#[pyfunction]
#[must_use]
pub fn unit_si_factor(name: &str) -> Option<f64> {
    azoth_core::unit_vocab_gen::si_factor(name)
}

/// One display unit's conversion to the SI base magnitude: `(value + offset) * factor`.
///
/// **The counterpart of [`unit_si_factor`] for the units a spec may not declare.** That function
/// reads a scale at one, which is exactly what an affine unit's answer at one is *not* - a degree
/// Celsius is 274.15 K at one - so the two are different calls and the check that holds each is
/// `python/tests/test_units_cross_library.py`, at one value for a scale and two for an affine
/// unit, because one point cannot tell a scale from a shifted one.
#[pyfunction]
#[must_use]
pub fn unit_affine_si(name: &str, value: f64) -> Option<f64> {
    azoth_core::unit_vocab_gen::affine_si(name, value)
}

/// The base dimensions, in the order the exponent tuples above are written in.
#[pyfunction]
#[must_use]
pub fn unit_slots() -> Vec<String> {
    azoth_core::unit_vocab_gen::SLOTS
        .iter()
        .map(|slot| (*slot).to_string())
        .collect()
}

/// The named unit sets, as `{set id: {dimension: unit}}`.
///
/// The same sets a front end switches a display between, and the same table
/// `azoth.core._units_gen.UNIT_SETS` is generated from: this function is what makes
/// the two halves of one table comparable from Python, so a set that reads one way
/// in Rust and another in Python is a failing test rather than a difference a
/// reader would have to notice.
#[pyfunction]
#[must_use]
pub fn unit_sets() -> std::collections::BTreeMap<String, std::collections::BTreeMap<String, String>>
{
    azoth_core::unit_vocab_gen::UNIT_SETS
        .iter()
        .map(|set| {
            let units = set
                .units
                .iter()
                .map(|(dimension, unit)| ((*dimension).to_string(), (*unit).to_string()))
                .collect();
            (set.id.to_string(), units)
        })
        .collect()
}

/// Every solver kind this crate implements, in the schema's spelling.
///
/// The third leg of the same contract `warning_codes` and `unit_names` each
/// provide one leg of: `specs/schema/calc.schema.json`'s `solver.kind` enum,
/// `azoth.core.solver.SolverKind` and this crate's `SolverKind::ALL` must name one
/// set, and a Python test can only assert that if the Rust list is reachable from
/// Python.
///
/// The schema's own description of `solver.kind` names this function as the
/// missing piece - "Nothing does that for solver kinds today, which is why this
/// enum is narrow rather than merely unchecked" - so its absence was the stated
/// reason no second solver kind could be added. It exists now, which is what
/// makes widening that enum a mechanical act rather than an unchecked one.
#[pyfunction]
#[must_use]
pub fn solver_kinds() -> Vec<String> {
    SolverKind::ALL
        .iter()
        .map(|kind| kind.as_str().to_string())
        .collect()
}

/// Version of the Rust core, so provenance can record which build answered.
#[pyfunction]
#[must_use]
pub fn version() -> String {
    azoth_core::VERSION.to_string()
}

#[cfg(test)]
mod transport_tests {
    //! The Rust side of the boundary, exercised without Python's help.
    //!
    //! `test_result_shapes_agree_across_languages` compares the *names* the extension
    //! exposes against the Rust `FIELDS` lists, reading them from Python. It cannot see
    //! a *mapping* fault: if the conversion wrote `gas_flow: r.liquid_flow`, both names
    //! would still be present, both would still be `f64`, and every shape check would
    //! pass while the two outlets reported each other's flow.
    //!
    //! So each result below is built with a distinct value in every field and read back
    //! through its Python attribute. Distinctness is the whole instrument: two fields
    //! that could be confused must not carry the same number.
    //!
    //! Building the struct as a literal is a second check for free. A field added to
    //! `PsFlashResult` and not to `PyPsFlashResult` stops these tests compiling, which
    //! is the only mechanism that turns a silently unexposed field into a visible one.

    use super::*;
    use crate::transport_gen::*;
    use azoth_core::CalcResult;
    use azoth_core::units::kelvins;
    use azoth_eos::Phase;
    use azoth_eos::results::PsFlashResult;

    /// Read a `f64` attribute from a converted result.
    fn number(result: &Bound<'_, PyAny>, name: &str) -> f64 {
        result
            .getattr(name)
            .unwrap_or_else(|e| panic!("no attribute {name}: {e}"))
            .extract()
            .unwrap_or_else(|e| panic!("{name} is not a number: {e}"))
    }

    /// Read a nested quantity's magnitude, e.g. `T.magnitude_si`.
    fn magnitude(result: &Bound<'_, PyAny>, name: &str) -> f64 {
        result
            .getattr(name)
            .unwrap_or_else(|e| panic!("no attribute {name}: {e}"))
            .getattr("magnitude_si")
            .unwrap_or_else(|e| panic!("{name}.magnitude_si: {e}"))
            .extract()
            .unwrap_or_else(|e| panic!("{name}.magnitude_si is not a number: {e}"))
    }

    fn numbers(result: &Bound<'_, PyAny>, name: &str) -> Vec<f64> {
        result
            .getattr(name)
            .unwrap_or_else(|e| panic!("no attribute {name}: {e}"))
            .extract()
            .unwrap_or_else(|e| panic!("{name} is not a list of numbers: {e}"))
    }

    /// Every name the Rust result declares must be readable on the converted object.
    fn assert_exposes_all_fields<T: CalcResult>(result: &Bound<'_, PyAny>) {
        for name in T::FIELDS {
            assert!(
                result.hasattr(*name).unwrap_or(false),
                "{name} is in {}::FIELDS but the converted object has no such \
                 attribute, so Python cannot see a field the Rust side declares",
                std::any::type_name::<T>()
            );
        }
    }

    fn ps_flash() -> PsFlashResult {
        PsFlashResult {
            temperature: kelvins(311.0),
            vapour_fraction: Some(0.5),
            x: vec![0.31, 0.32],
            y: vec![0.41, 0.42],
            k: vec![1.51, 1.52],
            phase: Phase::TwoPhase,
            z_liquid: 0.061,
            z_vapour: 0.062,
            iterations: 9,
            residual: 1.5e-11,
            warnings: Vec::new(),
        }
    }

    #[test]
    fn ps_flash_exposes_every_field_under_its_own_name() {
        Python::attach(|py| {
            let converted = Py::new(py, PyPsFlashResult::from(&ps_flash())).unwrap();
            let result = converted.bind(py);

            assert_exposes_all_fields::<PsFlashResult>(result);
            assert_eq!(magnitude(result, "T"), 311.0);
            assert_eq!(number(result, "vapour_fraction"), 0.5);
            assert_eq!(numbers(result, "x"), vec![0.31, 0.32]);
            assert_eq!(numbers(result, "y"), vec![0.41, 0.42]);
            assert_eq!(numbers(result, "k"), vec![1.51, 1.52]);
            assert_eq!(number(result, "z_liquid"), 0.061);
            assert_eq!(number(result, "z_vapour"), 0.062);
            assert_eq!(number(result, "iterations"), 9.0);
            // `ps_flash`'s residual is an entropy difference, so the spec declares it in
            // `J/(mol*K)` and it transports as a quantity rather than a bare number.
            assert_eq!(magnitude(result, "residual"), 1.5e-11);
        });
    }
}
