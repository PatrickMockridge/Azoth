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
use azoth_core::solver::SolverKind;
use azoth_core::units::UNIT_NAMES;
use azoth_core::warning::{Warning, WarningCode};
use azoth_eos::results::{
    BubblePressureResult, DewPressureResult, IdealGasCpResult, PrAlphaAbResult, PrDepartureResult,
    PrKappaResult, PrMassDensityResult, PrMolarVolumeResult, PrZFactorResult, PrsvKappaResult,
    PtFlashResult, PureSaturationResult, RachfordRiceBinaryResult, Vdw1fMixBinaryResult,
};
use azoth_thermal::results::ConductionPlaneWallResult;

use azoth_hydraulics::results::{
    ChokedFlowAreaResult, ColebrookResult, ControlValveCvResult, DarcyWeisbachResult,
    HaalandResult, KComponent, KFactorsResult, OrificeFlowResult, PumpPowerResult,
    ReynoldsNumberResult, SwameeJainResult,
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

/// Result of `hydraulics.pump_power`, transported.
///
/// Carries a dimensioned output, so it transports a [`PyQty`] rather than a bare
/// float - the same shape `darcy_weisbach` and `conduction_plane_wall` use.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PumpPowerResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPumpPowerResult {
    /// Shaft power, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub power: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPumpPowerResult {
    fn __repr__(&self) -> String {
        format!(
            "PumpPowerResult(power={} {})",
            self.power.magnitude_si, self.power.unit
        )
    }
}

impl From<&PumpPowerResult> for PyPumpPowerResult {
    fn from(r: &PumpPowerResult) -> Self {
        Self {
            power: PyQty {
                magnitude_si: r.power.value,
                unit: "W".to_string(),
            },
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

/// Result of `hydraulics.control_valve_cv`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ControlValveCvResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyControlValveCvResult {
    /// Volumetric flow rate, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub q: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyControlValveCvResult {
    fn __repr__(&self) -> String {
        format!(
            "ControlValveCvResult(q={} {})",
            self.q.magnitude_si, self.q.unit
        )
    }
}

impl From<&ControlValveCvResult> for PyControlValveCvResult {
    fn from(r: &ControlValveCvResult) -> Self {
        Self {
            q: PyQty {
                magnitude_si: r.q.value,
                unit: "m**3/s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `hydraulics.choked_flow_area`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ChokedFlowAreaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyChokedFlowAreaResult {
    /// Throat area, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub a: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyChokedFlowAreaResult {
    fn __repr__(&self) -> String {
        format!(
            "ChokedFlowAreaResult(a={} {})",
            self.a.magnitude_si, self.a.unit
        )
    }
}

impl From<&ChokedFlowAreaResult> for PyChokedFlowAreaResult {
    fn from(r: &ChokedFlowAreaResult) -> Self {
        Self {
            a: PyQty {
                magnitude_si: r.a.value,
                unit: "m**2".to_string(),
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

/// Result of `eos.pr_kappa`, transported.
///
/// Both the input and the output are dimensionless, so this carries a bare `f64`
/// and no [`PyQty`] - the same shape `reynolds_number` uses, and the reason there
/// is no unit string here to keep in step with the spec.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrKappaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrKappaResult {
    /// The Peng-Robinson alpha-function coefficient. Dimensionless.
    #[pyo3(get)]
    pub kappa: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrKappaResult {
    fn __repr__(&self) -> String {
        format!(
            "PrKappaResult(kappa={}, {} warning(s))",
            self.kappa,
            self.warnings.len()
        )
    }
}

impl From<&PrKappaResult> for PyPrKappaResult {
    fn from(r: &PrKappaResult) -> Self {
        Self {
            kappa: r.kappa,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_alpha_ab`, transported.
///
/// Three dimensionless outputs, so three bare `f64`s and no [`PyQty`].
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrAlphaAbResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrAlphaAbResult {
    /// The alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// The cubic's `A`. Dimensionless.
    #[pyo3(get)]
    pub a_reduced: f64,
    /// The cubic's `B`. Dimensionless.
    #[pyo3(get)]
    pub b_reduced: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrAlphaAbResult {
    fn __repr__(&self) -> String {
        format!(
            "PrAlphaAbResult(alpha={}, a_reduced={}, b_reduced={}, {} warning(s))",
            self.alpha,
            self.a_reduced,
            self.b_reduced,
            self.warnings.len()
        )
    }
}

impl From<&PrAlphaAbResult> for PyPrAlphaAbResult {
    fn from(r: &PrAlphaAbResult) -> Self {
        Self {
            alpha: r.alpha,
            a_reduced: r.a_reduced,
            b_reduced: r.b_reduced,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_z_factor`, transported.
///
/// Carries `root_structure` as the spec's string rather than as an enum, so the
/// bridge rebuilds `azoth.core.result.RootStructure` - the same arrangement
/// `PyReynoldsNumberResult.regime` uses, and for the same reason: the enum is
/// Python's, and this layer must not grow a second definition of it.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrZFactorResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrZFactorResult {
    /// The smallest admissible root. Dimensionless.
    #[pyo3(get)]
    pub z_min: f64,
    /// The largest admissible root. Dimensionless.
    #[pyo3(get)]
    pub z_max: f64,
    /// `one_root` or `three_roots`.
    #[pyo3(get)]
    pub root_structure: String,
    /// Newton steps the polish took.
    #[pyo3(get)]
    pub iterations: u32,
    /// Whether the polish met its stopping rule.
    #[pyo3(get)]
    pub converged: bool,
    /// The largest `|x_k - x_{k-1}|` at the final polish step.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrZFactorResult {
    fn __repr__(&self) -> String {
        format!(
            "PrZFactorResult(z_min={}, z_max={}, root_structure={}, {} warning(s))",
            self.z_min,
            self.z_max,
            self.root_structure,
            self.warnings.len()
        )
    }
}

/// Result of `eos.prsv_kappa`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrsvKappaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrsvKappaResult {
    /// The PRSV alpha-function coefficient. Dimensionless.
    #[pyo3(get)]
    pub kappa: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrsvKappaResult {
    fn __repr__(&self) -> String {
        format!(
            "PrsvKappaResult(kappa={}, {} warning(s))",
            self.kappa,
            self.warnings.len()
        )
    }
}

impl From<&PrsvKappaResult> for PyPrsvKappaResult {
    fn from(r: &PrsvKappaResult) -> Self {
        Self {
            kappa: r.kappa,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_departure`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrDepartureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrDepartureResult {
    /// The logarithm of the fugacity coefficient. Dimensionless.
    #[pyo3(get)]
    pub ln_phi: f64,
    /// The departure enthalpy over `R*T`. Dimensionless.
    #[pyo3(get)]
    pub h_dep_rt: f64,
    /// The departure entropy over `R`. Dimensionless.
    #[pyo3(get)]
    pub s_dep_r: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrDepartureResult {
    fn __repr__(&self) -> String {
        format!(
            "PrDepartureResult(ln_phi={}, h_dep_rt={}, s_dep_r={}, {} warning(s))",
            self.ln_phi,
            self.h_dep_rt,
            self.s_dep_r,
            self.warnings.len()
        )
    }
}

impl From<&PrDepartureResult> for PyPrDepartureResult {
    fn from(r: &PrDepartureResult) -> Self {
        Self {
            ln_phi: r.ln_phi,
            h_dep_rt: r.h_dep_rt,
            s_dep_r: r.s_dep_r,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.vdw1f_mix_binary`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "Vdw1fMixBinaryResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyVdw1fMixBinaryResult {
    /// The mixture's attraction parameter. Dimensionless.
    #[pyo3(get)]
    pub a_mix: f64,
    /// The mixture's repulsion parameter. Dimensionless.
    #[pyo3(get)]
    pub b_mix: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyVdw1fMixBinaryResult {
    fn __repr__(&self) -> String {
        format!(
            "Vdw1fMixBinaryResult(a_mix={}, b_mix={}, {} warning(s))",
            self.a_mix,
            self.b_mix,
            self.warnings.len()
        )
    }
}

impl From<&Vdw1fMixBinaryResult> for PyVdw1fMixBinaryResult {
    fn from(r: &Vdw1fMixBinaryResult) -> Self {
        Self {
            a_mix: r.a_mix,
            b_mix: r.b_mix,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.rachford_rice_binary`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "RachfordRiceBinaryResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyRachfordRiceBinaryResult {
    /// The vapour fraction. Dimensionless, and outside `[0, 1]` when the feed is
    /// single phase - in which case the result carries a warning.
    #[pyo3(get)]
    pub beta: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyRachfordRiceBinaryResult {
    fn __repr__(&self) -> String {
        format!(
            "RachfordRiceBinaryResult(beta={}, {} warning(s))",
            self.beta,
            self.warnings.len()
        )
    }
}

impl From<&RachfordRiceBinaryResult> for PyRachfordRiceBinaryResult {
    fn from(r: &RachfordRiceBinaryResult) -> Self {
        Self {
            beta: r.beta,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_molar_volume`, transported.
///
/// The only transport class in this namespace carrying a `PyQty`, because it is the
/// only calc here that returns a dimensioned quantity.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrMolarVolumeResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrMolarVolumeResult {
    /// Molar volume, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub v: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrMolarVolumeResult {
    fn __repr__(&self) -> String {
        format!(
            "PrMolarVolumeResult(v={} {})",
            self.v.magnitude_si, self.v.unit
        )
    }
}

impl From<&PrMolarVolumeResult> for PyPrMolarVolumeResult {
    fn from(r: &PrMolarVolumeResult) -> Self {
        Self {
            v: PyQty {
                magnitude_si: r.v.value,
                unit: "m**3/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_mass_density`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrMassDensityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrMassDensityResult {
    /// Mass density, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub rho: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrMassDensityResult {
    fn __repr__(&self) -> String {
        format!(
            "PrMassDensityResult(rho={} {})",
            self.rho.magnitude_si, self.rho.unit
        )
    }
}

/// Result of `eos.pure_saturation`, transported.
///
/// A model's result, shaped like any other. The difference between a model and a
/// calculation is in how the answer was reached, not in what an answer is.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PureSaturationResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPureSaturationResult {
    /// The saturation pressure.
    #[pyo3(get)]
    pub p_sat: PyQty,
    /// The common `ln phi` at the converged pressure.
    #[pyo3(get)]
    pub ln_phi: f64,
    /// Bisection steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The final bracket's dimensionless half-width.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPureSaturationResult {
    fn __repr__(&self) -> String {
        format!(
            "PureSaturationResult(p_sat={} {}, {} iteration(s))",
            self.p_sat.magnitude_si, self.p_sat.unit, self.iterations
        )
    }
}

impl From<&PureSaturationResult> for PyPureSaturationResult {
    fn from(r: &PureSaturationResult) -> Self {
        Self {
            p_sat: PyQty {
                magnitude_si: r.p_sat.value,
                unit: "Pa".to_string(),
            },
            ln_phi: r.ln_phi,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pt_flash`, transported.
///
/// The first result here whose fields are vectors, and the first with an optional
/// scalar. `beta` crosses as `Option<f64>` and arrives as `None` - that is the whole
/// point of it being optional, so flattening it to a sentinel number at the boundary
/// would undo the design one layer below where it was made.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PtFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPtFlashResult {
    /// The vapour fraction, or `None` when there is none to report.
    #[pyo3(get)]
    pub beta: Option<f64>,
    /// Liquid-phase mole fractions.
    #[pyo3(get)]
    pub x: Vec<f64>,
    /// Vapour-phase mole fractions.
    #[pyo3(get)]
    pub y: Vec<f64>,
    /// `K_i = y_i / x_i`.
    #[pyo3(get)]
    pub k: Vec<f64>,
    /// `ln phi_i` in the liquid phase.
    #[pyo3(get)]
    pub ln_phi_liquid: Vec<f64>,
    /// `ln phi_i` in the vapour phase.
    #[pyo3(get)]
    pub ln_phi_vapour: Vec<f64>,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// The smallest `T / Tc_i` over the components.
    #[pyo3(get)]
    pub min_t_over_tc: f64,
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// Successive-substitution steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The rms change in `ln K` at the last completed step, or `NaN`.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPtFlashResult {
    fn __repr__(&self) -> String {
        match self.beta {
            Some(beta) => format!(
                "PtFlashResult(phase={}, beta={beta}, z_liquid={}, z_vapour={})",
                self.phase, self.z_liquid, self.z_vapour
            ),
            None => format!(
                "PtFlashResult(phase={}, beta=None, {} iteration(s))",
                self.phase, self.iterations
            ),
        }
    }
}

/// Result of `eos.ideal_gas_cp`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "IdealGasCpResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyIdealGasCpResult {
    /// The polynomial's value, `Cp/R`.
    #[pyo3(get)]
    pub cp_over_r: f64,
    /// The ideal-gas heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyIdealGasCpResult {
    fn __repr__(&self) -> String {
        format!(
            "IdealGasCpResult(cp={} {})",
            self.cp.magnitude_si, self.cp.unit
        )
    }
}

impl From<&IdealGasCpResult> for PyIdealGasCpResult {
    fn from(r: &IdealGasCpResult) -> Self {
        Self {
            cp_over_r: r.cp_over_r,
            cp: PyQty {
                magnitude_si: r.cp.value,
                unit: "J/(mol*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.bubble_pressure` or `eos.dew_pressure`, transported.
///
/// One transport type for two models, because the Rust results have the same shape
/// and the difference is only which phase is which - which the bridge resolves by
/// naming. Two `pyclass`es would be two copies of the same thirteen lines.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PhaseBoundaryResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPhaseBoundaryResult {
    /// The boundary pressure.
    #[pyo3(get)]
    pub pressure: PyQty,
    /// The incipient phase's composition.
    #[pyo3(get)]
    pub incipient: Vec<f64>,
    /// K-values at the converged pressure.
    #[pyo3(get)]
    pub k: Vec<f64>,
    /// The liquid root of the cubic at the converged state.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// The smallest `T / Tc_i` over the components.
    #[pyo3(get)]
    pub min_t_over_tc: f64,
    /// Pressure updates taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The final residual.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPhaseBoundaryResult {
    fn __repr__(&self) -> String {
        format!(
            "PhaseBoundaryResult(pressure={} {}, {} iteration(s))",
            self.pressure.magnitude_si, self.pressure.unit, self.iterations
        )
    }
}

impl From<&BubblePressureResult> for PyPhaseBoundaryResult {
    fn from(r: &BubblePressureResult) -> Self {
        Self {
            pressure: PyQty {
                magnitude_si: r.pressure.value,
                unit: "Pa".to_string(),
            },
            incipient: r.incipient.clone(),
            k: r.k.clone(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            min_t_over_tc: r.min_t_over_tc,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

impl From<&DewPressureResult> for PyPhaseBoundaryResult {
    fn from(r: &DewPressureResult) -> Self {
        Self {
            pressure: PyQty {
                magnitude_si: r.pressure.value,
                unit: "Pa".to_string(),
            },
            incipient: r.incipient.clone(),
            k: r.k.clone(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            min_t_over_tc: r.min_t_over_tc,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

impl From<&PtFlashResult> for PyPtFlashResult {
    fn from(r: &PtFlashResult) -> Self {
        Self {
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            ln_phi_liquid: r.ln_phi_liquid.clone(),
            ln_phi_vapour: r.ln_phi_vapour.clone(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            min_t_over_tc: r.min_t_over_tc,
            phase: r.phase.as_str().to_string(),
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

impl From<&PrMassDensityResult> for PyPrMassDensityResult {
    fn from(r: &PrMassDensityResult) -> Self {
        Self {
            rho: PyQty {
                magnitude_si: r.rho.value,
                unit: "kg/m**3".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

impl From<&PrZFactorResult> for PyPrZFactorResult {
    fn from(r: &PrZFactorResult) -> Self {
        Self {
            z_min: r.z_min,
            z_max: r.z_max,
            root_structure: r.root_structure.as_str().to_string(),
            iterations: r.iterations,
            converged: r.converged,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
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
        OrificeFlowResult::CALC_ID => OrificeFlowResult::FIELDS.to_vec(),
        ControlValveCvResult::CALC_ID => ControlValveCvResult::FIELDS.to_vec(),
        ChokedFlowAreaResult::CALC_ID => ChokedFlowAreaResult::FIELDS.to_vec(),
        ConductionPlaneWallResult::CALC_ID => ConductionPlaneWallResult::FIELDS.to_vec(),
        PrKappaResult::CALC_ID => PrKappaResult::FIELDS.to_vec(),
        PrAlphaAbResult::CALC_ID => PrAlphaAbResult::FIELDS.to_vec(),
        PrZFactorResult::CALC_ID => PrZFactorResult::FIELDS.to_vec(),
        PrsvKappaResult::CALC_ID => PrsvKappaResult::FIELDS.to_vec(),
        PrDepartureResult::CALC_ID => PrDepartureResult::FIELDS.to_vec(),
        Vdw1fMixBinaryResult::CALC_ID => Vdw1fMixBinaryResult::FIELDS.to_vec(),
        RachfordRiceBinaryResult::CALC_ID => RachfordRiceBinaryResult::FIELDS.to_vec(),
        PrMolarVolumeResult::CALC_ID => PrMolarVolumeResult::FIELDS.to_vec(),
        PrMassDensityResult::CALC_ID => PrMassDensityResult::FIELDS.to_vec(),
        // Models. Present here because a result's *shape* is a cross-language
        // contract whether or not its spec calls it a calculation, and before this
        // the model results were covered by no shape check at all.
        PureSaturationResult::CALC_ID => PureSaturationResult::FIELDS.to_vec(),
        PtFlashResult::CALC_ID => PtFlashResult::FIELDS.to_vec(),
        BubblePressureResult::CALC_ID => BubblePressureResult::FIELDS.to_vec(),
        DewPressureResult::CALC_ID => DewPressureResult::FIELDS.to_vec(),
        IdealGasCpResult::CALC_ID => IdealGasCpResult::FIELDS.to_vec(),
        PumpPowerResult::CALC_ID => PumpPowerResult::FIELDS.to_vec(),
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
        OrificeFlowResult::CALC_ID.to_string(),
        ControlValveCvResult::CALC_ID.to_string(),
        ChokedFlowAreaResult::CALC_ID.to_string(),
        ConductionPlaneWallResult::CALC_ID.to_string(),
        PrKappaResult::CALC_ID.to_string(),
        PrAlphaAbResult::CALC_ID.to_string(),
        PrZFactorResult::CALC_ID.to_string(),
        PrsvKappaResult::CALC_ID.to_string(),
        PrDepartureResult::CALC_ID.to_string(),
        Vdw1fMixBinaryResult::CALC_ID.to_string(),
        RachfordRiceBinaryResult::CALC_ID.to_string(),
        PrMolarVolumeResult::CALC_ID.to_string(),
        PrMassDensityResult::CALC_ID.to_string(),
        IdealGasCpResult::CALC_ID.to_string(),
        PumpPowerResult::CALC_ID.to_string(),
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
