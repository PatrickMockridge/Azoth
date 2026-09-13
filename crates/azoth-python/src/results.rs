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

use azoth_core::CalcResult;
use azoth_core::units::UNIT_NAMES;
use azoth_core::warning::{Warning, WarningCode};
use azoth_thermal::results::ConductionPlaneWallResult;

use azoth_hydraulics::results::{
    ColebrookResult, DarcyWeisbachResult, HaalandResult, KComponent, KFactorsResult,
    OrificeFlowResult, ReynoldsNumberResult, SwameeJainResult,
};
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

fn transport(warnings: &[Warning]) -> Vec<PyWarning> {
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

/// Result of `hydraulics.reynolds_number`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ReynoldsNumberResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyReynoldsNumberResult {
    /// Reynolds number. Dimensionless.
    #[pyo3(get)]
    pub re: f64,
    /// Flow regime: `laminar`, `transitional` or `turbulent`.
    #[pyo3(get)]
    pub regime: String,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyReynoldsNumberResult {
    fn __repr__(&self) -> String {
        format!(
            "ReynoldsNumberResult(re={}, regime={}, {} warning(s))",
            self.re,
            self.regime,
            self.warnings.len()
        )
    }
}

impl From<&ReynoldsNumberResult> for PyReynoldsNumberResult {
    fn from(r: &ReynoldsNumberResult) -> Self {
        Self {
            re: r.re,
            regime: r.regime.as_str().to_string(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `hydraulics.friction_factor_colebrook`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ColebrookResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyColebrookResult {
    /// Darcy friction factor. Dimensionless.
    #[pyo3(get)]
    pub f: f64,
    /// Iterations performed.
    #[pyo3(get)]
    pub iterations: u32,
    /// Whether the iteration met its tolerance.
    #[pyo3(get)]
    pub converged: bool,
    /// Final change between iterates.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyColebrookResult {
    fn __repr__(&self) -> String {
        format!(
            "ColebrookResult(f={}, iterations={}, converged={})",
            self.f, self.iterations, self.converged
        )
    }
}

impl From<&ColebrookResult> for PyColebrookResult {
    fn from(r: &ColebrookResult) -> Self {
        Self {
            f: r.f,
            iterations: r.iterations,
            converged: r.converged,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `hydraulics.friction_factor_swamee_jain`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SwameeJainResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySwameeJainResult {
    /// Darcy friction factor. Dimensionless.
    #[pyo3(get)]
    pub f: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySwameeJainResult {
    fn __repr__(&self) -> String {
        format!("SwameeJainResult(f={})", self.f)
    }
}

impl From<&SwameeJainResult> for PySwameeJainResult {
    fn from(r: &SwameeJainResult) -> Self {
        Self {
            f: r.f,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `hydraulics.orifice_flow`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "OrificeFlowResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyOrificeFlowResult {
    /// Volumetric flow rate, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub q: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyOrificeFlowResult {
    fn __repr__(&self) -> String {
        format!(
            "OrificeFlowResult(q={} {})",
            self.q.magnitude_si, self.q.unit
        )
    }
}

impl From<&OrificeFlowResult> for PyOrificeFlowResult {
    fn from(r: &OrificeFlowResult) -> Self {
        Self {
            q: PyQty {
                magnitude_si: r.q.value,
                unit: "m**3/s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `thermal.conduction_plane_wall`, transported.
///
/// Carries a dimensioned output, so it transports a [`PyQty`] rather than a bare
/// float - the same shape `darcy_weisbach` uses, and the reason the transport layer
/// has a quantity type at all.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ConductionPlaneWallResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyConductionPlaneWallResult {
    /// Heat flow rate through the wall, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub q: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyConductionPlaneWallResult {
    fn __repr__(&self) -> String {
        format!(
            "ConductionPlaneWallResult(q={} {})",
            self.q.magnitude_si, self.q.unit
        )
    }
}

impl From<&ConductionPlaneWallResult> for PyConductionPlaneWallResult {
    fn from(r: &ConductionPlaneWallResult) -> Self {
        Self {
            q: PyQty {
                magnitude_si: r.q.value,
                unit: "W".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `hydraulics.friction_factor_haaland`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HaalandResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHaalandResult {
    /// Darcy friction factor. Dimensionless.
    #[pyo3(get)]
    pub f: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHaalandResult {
    fn __repr__(&self) -> String {
        format!("HaalandResult(f={})", self.f)
    }
}

impl From<&HaalandResult> for PyHaalandResult {
    fn from(r: &HaalandResult) -> Self {
        Self {
            f: r.f,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `hydraulics.crane_k_factors`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "KFactorsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyKFactorsResult {
    /// Total resistance coefficient.
    #[pyo3(get)]
    pub k_total: f64,
    /// The friction factor the coefficients were based on.
    #[pyo3(get)]
    pub f_t: f64,
    /// Per-fitting breakdown.
    #[pyo3(get)]
    pub components: Vec<PyKComponent>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyKFactorsResult {
    fn __repr__(&self) -> String {
        format!(
            "KFactorsResult(k_total={}, {} component(s))",
            self.k_total,
            self.components.len()
        )
    }
}

impl From<&KFactorsResult> for PyKFactorsResult {
    fn from(r: &KFactorsResult) -> Self {
        Self {
            k_total: r.k_total,
            f_t: r.f_t,
            components: r
                .components
                .iter()
                .map(|c: &KComponent| PyKComponent {
                    fitting_id: c.fitting_id.clone(),
                    n_ld: c.n_ld,
                    k: c.k,
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `hydraulics.darcy_weisbach`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "DarcyWeisbachResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyDarcyWeisbachResult {
    /// Pressure drop over the pipe length.
    #[pyo3(get)]
    pub dp: PyQty,
    /// The friction factor the drop was computed with.
    #[pyo3(get)]
    pub f: f64,
    /// Reynolds number, present only when viscosity was supplied.
    #[pyo3(get)]
    pub re: Option<f64>,
    /// Flow regime, present only when viscosity was supplied.
    #[pyo3(get)]
    pub regime: Option<String>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyDarcyWeisbachResult {
    fn __repr__(&self) -> String {
        format!(
            "DarcyWeisbachResult(dp={} Pa, f={}, re={:?})",
            self.dp.magnitude_si, self.f, self.re
        )
    }
}

impl From<&DarcyWeisbachResult> for PyDarcyWeisbachResult {
    fn from(r: &DarcyWeisbachResult) -> Self {
        Self {
            dp: PyQty {
                magnitude_si: r.dp.value,
                unit: "Pa".to_string(),
            },
            f: r.f,
            re: r.re,
            regime: r.regime.map(|regime| regime.as_str().to_string()),
            warnings: transport(&r.warnings),
        }
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

/// The public field names of a calc's result, in declaration order.
///
/// Returns an empty list for an unknown id rather than raising: this is an
/// introspection helper for tests and diagnostics, and a missing calc should
/// read as "no fields" rather than as an error to handle.
#[pyfunction]
#[must_use]
pub fn result_fields(calc_id: &str) -> Vec<String> {
    match calc_id {
        ReynoldsNumberResult::CALC_ID => ReynoldsNumberResult::FIELDS.to_vec(),
        ColebrookResult::CALC_ID => ColebrookResult::FIELDS.to_vec(),
        SwameeJainResult::CALC_ID => SwameeJainResult::FIELDS.to_vec(),
        HaalandResult::CALC_ID => HaalandResult::FIELDS.to_vec(),
        ConductionPlaneWallResult::CALC_ID => ConductionPlaneWallResult::FIELDS.to_vec(),
        OrificeFlowResult::CALC_ID => OrificeFlowResult::FIELDS.to_vec(),
        KFactorsResult::CALC_ID => KFactorsResult::FIELDS.to_vec(),
        DarcyWeisbachResult::CALC_ID => DarcyWeisbachResult::FIELDS.to_vec(),
        _ => Vec::new(),
    }
    .into_iter()
    .map(str::to_string)
    .collect()
}

/// The calc ids this extension implements.
#[pyfunction]
#[must_use]
pub fn calc_ids() -> Vec<String> {
    vec![
        ReynoldsNumberResult::CALC_ID.to_string(),
        ColebrookResult::CALC_ID.to_string(),
        SwameeJainResult::CALC_ID.to_string(),
        HaalandResult::CALC_ID.to_string(),
        ConductionPlaneWallResult::CALC_ID.to_string(),
        OrificeFlowResult::CALC_ID.to_string(),
        KFactorsResult::CALC_ID.to_string(),
        DarcyWeisbachResult::CALC_ID.to_string(),
    ]
}

/// Version of the Rust core, so provenance can record which build answered.
#[pyfunction]
#[must_use]
pub fn version() -> String {
    azoth_core::VERSION.to_string()
}
