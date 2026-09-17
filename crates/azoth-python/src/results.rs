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
    AmmoniaPhaseResult, AntoineVaporPressureResult, ArgonSolidPhaseResult, BubblePressureResult,
    BubbleTemperatureResult, BwrsPhaseResult, ChungConductivityResult, ChungViscosityResult,
    Co2PhaseResult, Co2WaterDiffusivityResult, CostaldMolarVolumeResult, CriticalPointResult,
    DewPressureResult, DewTemperatureResult, EosCgPhaseResult, Gerg2008PhaseResult,
    HaydukMinhasDiffusivityResult, HeatOfVaporizationResult, HeliumPhaseResult,
    HydrogenPhaseResult, IdealGasCpResult, LiquidHeatCapacityResult, MasonSaxenaConductivityResult,
    Matcop5PrumrAlphaResult, MatcopAlphaResult, MatcopPrAlphaResult, MatcopPrumrAlphaResult,
    MatcopPrumrNewAlphaResult, MolarEnthalpyEntropyResult, MollerupAlphaResult,
    NrtlActivityCoefficientsResult, ParachorSurfaceTensionResult, ParahydrogenSolidPhaseResult,
    PhFlashResult, Pr78KappaResult, PrAlphaAbResult, PrDaneshAlphaResult, PrDelft1998AlphaResult,
    PrDepartureResult, PrGassem2001AlphaResult, PrKappaResult, PrLeeKeslerAlphaResult,
    PrMassDensityResult, PrMolarVolumeResult, PrPenelouxShiftResult, PrZFactorResult,
    PrsvKappaResult, PsFlashResult, PtFlashResult, PtPhaseEnvelopeResult, PuFlashResult,
    PureSaturationResult, PvFlashResult, RachfordRiceBinaryResult, RackettMolarVolumeResult,
    RkAlphaAbResult, RkDepartureResult, SchwartzentruberAlphaResult, SiddiqiLucasDiffusivityResult,
    SoreideWhitsonAlphaResult, SrkAlphaAbResult, SrkDepartureResult, SrkKappaResult,
    SrkPenelouxShiftResult, SrkZFactorResult, StabilityTestResult, ThFlashResult,
    ThermalConductivityResult, TsFlashResult, TuFlashResult, TvFlashResult, TwuKappaResult,
    TwucoonAlphaResult, TwucoonParamAlphaResult, TwucoonStatoilAlphaResult,
    TynCalusDiffusivityResult, UmrprAlphaResult, UnifacActivityCoefficientsResult,
    UniquacActivityCoefficientsResult, VanLaarAcidActivityCoefficientsResult, Vdw1fMixBinaryResult,
    ViscosityResult, VuFlashResult, WaterPhaseResult, WilkeChangDiffusivityResult,
    WilkeViscosityResult, WilsonActivityCoefficientsResult,
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

/// Result of `eos.matcop_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MatcopAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMatcopAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMatcopAlphaResult {
    fn __repr__(&self) -> String {
        format!("MatcopAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&MatcopAlphaResult> for PyMatcopAlphaResult {
    fn from(r: &MatcopAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.matcop_pr_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MatcopPrAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMatcopPrAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMatcopPrAlphaResult {
    fn __repr__(&self) -> String {
        format!("MatcopPrAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&MatcopPrAlphaResult> for PyMatcopPrAlphaResult {
    fn from(r: &MatcopPrAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.matcop_prumr_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MatcopPrumrAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMatcopPrumrAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMatcopPrumrAlphaResult {
    fn __repr__(&self) -> String {
        format!("MatcopPrumrAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&MatcopPrumrAlphaResult> for PyMatcopPrumrAlphaResult {
    fn from(r: &MatcopPrumrAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.matcop_prumr_new_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MatcopPrumrNewAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMatcopPrumrNewAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMatcopPrumrNewAlphaResult {
    fn __repr__(&self) -> String {
        format!("MatcopPrumrNewAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&MatcopPrumrNewAlphaResult> for PyMatcopPrumrNewAlphaResult {
    fn from(r: &MatcopPrumrNewAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.matcop5_prumr_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "Matcop5PrumrAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMatcop5PrumrAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMatcop5PrumrAlphaResult {
    fn __repr__(&self) -> String {
        format!("Matcop5PrumrAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&Matcop5PrumrAlphaResult> for PyMatcop5PrumrAlphaResult {
    fn from(r: &Matcop5PrumrAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.mollerup_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MollerupAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMollerupAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMollerupAlphaResult {
    fn __repr__(&self) -> String {
        format!("MollerupAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&MollerupAlphaResult> for PyMollerupAlphaResult {
    fn from(r: &MollerupAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_danesh_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrDaneshAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrDaneshAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrDaneshAlphaResult {
    fn __repr__(&self) -> String {
        format!("PrDaneshAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&PrDaneshAlphaResult> for PyPrDaneshAlphaResult {
    fn from(r: &PrDaneshAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_delft1998_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrDelft1998AlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrDelft1998AlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrDelft1998AlphaResult {
    fn __repr__(&self) -> String {
        format!("PrDelft1998AlphaResult(alpha={})", self.alpha)
    }
}

impl From<&PrDelft1998AlphaResult> for PyPrDelft1998AlphaResult {
    fn from(r: &PrDelft1998AlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_gassem2001_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrGassem2001AlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrGassem2001AlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrGassem2001AlphaResult {
    fn __repr__(&self) -> String {
        format!("PrGassem2001AlphaResult(alpha={})", self.alpha)
    }
}

impl From<&PrGassem2001AlphaResult> for PyPrGassem2001AlphaResult {
    fn from(r: &PrGassem2001AlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_lee_kesler_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrLeeKeslerAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrLeeKeslerAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrLeeKeslerAlphaResult {
    fn __repr__(&self) -> String {
        format!("PrLeeKeslerAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&PrLeeKeslerAlphaResult> for PyPrLeeKeslerAlphaResult {
    fn from(r: &PrLeeKeslerAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
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

/// Result of `eos.pr78_kappa`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "Pr78KappaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPr78KappaResult {
    /// The 1978 Peng-Robinson alpha-function coefficient. Dimensionless.
    #[pyo3(get)]
    pub kappa: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPr78KappaResult {
    fn __repr__(&self) -> String {
        format!(
            "Pr78KappaResult(kappa={}, {} warning(s))",
            self.kappa,
            self.warnings.len()
        )
    }
}

impl From<&Pr78KappaResult> for PyPr78KappaResult {
    fn from(r: &Pr78KappaResult) -> Self {
        Self {
            kappa: r.kappa,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.twu_kappa`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TwuKappaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTwuKappaResult {
    /// Twu's alpha-function coefficient. Dimensionless.
    #[pyo3(get)]
    pub kappa: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTwuKappaResult {
    fn __repr__(&self) -> String {
        format!(
            "TwuKappaResult(kappa={}, {} warning(s))",
            self.kappa,
            self.warnings.len()
        )
    }
}

impl From<&TwuKappaResult> for PyTwuKappaResult {
    fn from(r: &TwuKappaResult) -> Self {
        Self {
            kappa: r.kappa,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.twucoon_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TwucoonAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTwucoonAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTwucoonAlphaResult {
    fn __repr__(&self) -> String {
        format!("TwucoonAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&TwucoonAlphaResult> for PyTwucoonAlphaResult {
    fn from(r: &TwucoonAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.twucoon_param_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TwucoonParamAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTwucoonParamAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTwucoonParamAlphaResult {
    fn __repr__(&self) -> String {
        format!("TwucoonParamAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&TwucoonParamAlphaResult> for PyTwucoonParamAlphaResult {
    fn from(r: &TwucoonParamAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.twucoon_statoil_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TwucoonStatoilAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTwucoonStatoilAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTwucoonStatoilAlphaResult {
    fn __repr__(&self) -> String {
        format!("TwucoonStatoilAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&TwucoonStatoilAlphaResult> for PyTwucoonStatoilAlphaResult {
    fn from(r: &TwucoonStatoilAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
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
    /// The departure heat capacity over `R`. Dimensionless.
    #[pyo3(get)]
    pub cp_dep_r: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrDepartureResult {
    fn __repr__(&self) -> String {
        format!(
            "PrDepartureResult(ln_phi={}, h_dep_rt={}, s_dep_r={}, cp_dep_r={}, {} warning(s))",
            self.ln_phi,
            self.h_dep_rt,
            self.s_dep_r,
            self.cp_dep_r,
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
            cp_dep_r: r.cp_dep_r,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.srk_kappa`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SrkKappaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySrkKappaResult {
    /// The Soave-Redlich-Kwong alpha-function coefficient. Dimensionless.
    #[pyo3(get)]
    pub kappa: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySrkKappaResult {
    fn __repr__(&self) -> String {
        format!(
            "SrkKappaResult(kappa={}, {} warning(s))",
            self.kappa,
            self.warnings.len()
        )
    }
}

impl From<&SrkKappaResult> for PySrkKappaResult {
    fn from(r: &SrkKappaResult) -> Self {
        Self {
            kappa: r.kappa,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_peneloux_shift`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrPenelouxShiftResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrPenelouxShiftResult {
    /// The Peneloux volume-translation parameter, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub c: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrPenelouxShiftResult {
    fn __repr__(&self) -> String {
        format!(
            "PrPenelouxShiftResult(c={} {})",
            self.c.magnitude_si, self.c.unit
        )
    }
}

impl From<&PrPenelouxShiftResult> for PyPrPenelouxShiftResult {
    fn from(r: &PrPenelouxShiftResult) -> Self {
        Self {
            c: PyQty {
                magnitude_si: r.c.value,
                unit: "m**3/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.srk_peneloux_shift`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SrkPenelouxShiftResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySrkPenelouxShiftResult {
    /// The Peneloux volume-translation parameter, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub c: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySrkPenelouxShiftResult {
    fn __repr__(&self) -> String {
        format!(
            "SrkPenelouxShiftResult(c={} {})",
            self.c.magnitude_si, self.c.unit
        )
    }
}

impl From<&SrkPenelouxShiftResult> for PySrkPenelouxShiftResult {
    fn from(r: &SrkPenelouxShiftResult) -> Self {
        Self {
            c: PyQty {
                magnitude_si: r.c.value,
                unit: "m**3/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.heat_of_vaporization`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HeatOfVaporizationResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHeatOfVaporizationResult {
    /// The pure-component heat of vaporisation, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub hov: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHeatOfVaporizationResult {
    fn __repr__(&self) -> String {
        format!(
            "HeatOfVaporizationResult(hov={} {})",
            self.hov.magnitude_si, self.hov.unit
        )
    }
}

impl From<&HeatOfVaporizationResult> for PyHeatOfVaporizationResult {
    fn from(r: &HeatOfVaporizationResult) -> Self {
        Self {
            hov: PyQty {
                magnitude_si: r.hov.value,
                unit: "J/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.liquid_heat_capacity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "LiquidHeatCapacityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyLiquidHeatCapacityResult {
    /// The pure-component liquid heat capacity, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub cp: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyLiquidHeatCapacityResult {
    fn __repr__(&self) -> String {
        format!(
            "LiquidHeatCapacityResult(cp={} {})",
            self.cp.magnitude_si, self.cp.unit
        )
    }
}

impl From<&LiquidHeatCapacityResult> for PyLiquidHeatCapacityResult {
    fn from(r: &LiquidHeatCapacityResult) -> Self {
        Self {
            cp: PyQty {
                magnitude_si: r.cp.value,
                unit: "J/(mol*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.antoine_vapor_pressure`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "AntoineVaporPressureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyAntoineVaporPressureResult {
    /// The pure-component vapour pressure, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub p_sat: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyAntoineVaporPressureResult {
    fn __repr__(&self) -> String {
        format!(
            "AntoineVaporPressureResult(p_sat={} {})",
            self.p_sat.magnitude_si, self.p_sat.unit
        )
    }
}

impl From<&AntoineVaporPressureResult> for PyAntoineVaporPressureResult {
    fn from(r: &AntoineVaporPressureResult) -> Self {
        Self {
            p_sat: PyQty {
                magnitude_si: r.p_sat.value,
                unit: "Pa".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.rackett_molar_volume`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "RackettMolarVolumeResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyRackettMolarVolumeResult {
    /// The saturated liquid molar volume, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub v: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyRackettMolarVolumeResult {
    fn __repr__(&self) -> String {
        format!(
            "RackettMolarVolumeResult(v={} {})",
            self.v.magnitude_si, self.v.unit
        )
    }
}

impl From<&RackettMolarVolumeResult> for PyRackettMolarVolumeResult {
    fn from(r: &RackettMolarVolumeResult) -> Self {
        Self {
            v: PyQty {
                magnitude_si: r.v.value,
                unit: "m**3/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.costald_molar_volume`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "CostaldMolarVolumeResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyCostaldMolarVolumeResult {
    /// The saturated liquid molar volume, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub v: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyCostaldMolarVolumeResult {
    fn __repr__(&self) -> String {
        format!(
            "CostaldMolarVolumeResult(v={} {})",
            self.v.magnitude_si, self.v.unit
        )
    }
}

impl From<&CostaldMolarVolumeResult> for PyCostaldMolarVolumeResult {
    fn from(r: &CostaldMolarVolumeResult) -> Self {
        Self {
            v: PyQty {
                magnitude_si: r.v.value,
                unit: "m**3/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.chung_viscosity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ChungViscosityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyChungViscosityResult {
    /// The gas dynamic viscosity, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub mu: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyChungViscosityResult {
    fn __repr__(&self) -> String {
        format!(
            "ChungViscosityResult(mu={} {})",
            self.mu.magnitude_si, self.mu.unit
        )
    }
}

impl From<&ChungViscosityResult> for PyChungViscosityResult {
    fn from(r: &ChungViscosityResult) -> Self {
        Self {
            mu: PyQty {
                magnitude_si: r.mu.value,
                unit: "Pa*s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.chung_conductivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ChungConductivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyChungConductivityResult {
    /// The gas thermal conductivity, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub k: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyChungConductivityResult {
    fn __repr__(&self) -> String {
        format!(
            "ChungConductivityResult(k={} {})",
            self.k.magnitude_si, self.k.unit
        )
    }
}

impl From<&ChungConductivityResult> for PyChungConductivityResult {
    fn from(r: &ChungConductivityResult) -> Self {
        Self {
            k: PyQty {
                magnitude_si: r.k.value,
                unit: "W/(m*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.mason_saxena_conductivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MasonSaxenaConductivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMasonSaxenaConductivityResult {
    /// The gas mixture thermal conductivity, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub k: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMasonSaxenaConductivityResult {
    fn __repr__(&self) -> String {
        format!(
            "MasonSaxenaConductivityResult(k={} {})",
            self.k.magnitude_si, self.k.unit
        )
    }
}

impl From<&MasonSaxenaConductivityResult> for PyMasonSaxenaConductivityResult {
    fn from(r: &MasonSaxenaConductivityResult) -> Self {
        Self {
            k: PyQty {
                magnitude_si: r.k.value,
                unit: "W/(m*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.tyn_calus_diffusivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TynCalusDiffusivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTynCalusDiffusivityResult {
    /// The binary diffusion coefficient, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub d: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTynCalusDiffusivityResult {
    fn __repr__(&self) -> String {
        format!(
            "TynCalusDiffusivityResult(d={} {})",
            self.d.magnitude_si, self.d.unit
        )
    }
}

impl From<&TynCalusDiffusivityResult> for PyTynCalusDiffusivityResult {
    fn from(r: &TynCalusDiffusivityResult) -> Self {
        Self {
            d: PyQty {
                magnitude_si: r.d.value,
                unit: "m**2/s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.wilke_chang_diffusivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "WilkeChangDiffusivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyWilkeChangDiffusivityResult {
    /// The binary diffusion coefficient, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub d: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyWilkeChangDiffusivityResult {
    fn __repr__(&self) -> String {
        format!(
            "WilkeChangDiffusivityResult(d={} {})",
            self.d.magnitude_si, self.d.unit
        )
    }
}

impl From<&WilkeChangDiffusivityResult> for PyWilkeChangDiffusivityResult {
    fn from(r: &WilkeChangDiffusivityResult) -> Self {
        Self {
            d: PyQty {
                magnitude_si: r.d.value,
                unit: "m**2/s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.hayduk_minhas_diffusivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HaydukMinhasDiffusivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHaydukMinhasDiffusivityResult {
    /// The binary diffusion coefficient, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub d: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHaydukMinhasDiffusivityResult {
    fn __repr__(&self) -> String {
        format!(
            "HaydukMinhasDiffusivityResult(d={} {})",
            self.d.magnitude_si, self.d.unit
        )
    }
}

impl From<&HaydukMinhasDiffusivityResult> for PyHaydukMinhasDiffusivityResult {
    fn from(r: &HaydukMinhasDiffusivityResult) -> Self {
        Self {
            d: PyQty {
                magnitude_si: r.d.value,
                unit: "m**2/s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.schwartzentruber_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SchwartzentruberAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySchwartzentruberAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySchwartzentruberAlphaResult {
    fn __repr__(&self) -> String {
        format!("SchwartzentruberAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&SchwartzentruberAlphaResult> for PySchwartzentruberAlphaResult {
    fn from(r: &SchwartzentruberAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.siddiqi_lucas_diffusivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SiddiqiLucasDiffusivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySiddiqiLucasDiffusivityResult {
    /// The binary diffusion coefficient, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub d: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySiddiqiLucasDiffusivityResult {
    fn __repr__(&self) -> String {
        format!(
            "SiddiqiLucasDiffusivityResult(d={} {})",
            self.d.magnitude_si, self.d.unit
        )
    }
}

impl From<&SiddiqiLucasDiffusivityResult> for PySiddiqiLucasDiffusivityResult {
    fn from(r: &SiddiqiLucasDiffusivityResult) -> Self {
        Self {
            d: PyQty {
                magnitude_si: r.d.value,
                unit: "m**2/s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.co2_water_diffusivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "Co2WaterDiffusivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyCo2WaterDiffusivityResult {
    /// The binary diffusion coefficient, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub d: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyCo2WaterDiffusivityResult {
    fn __repr__(&self) -> String {
        format!(
            "Co2WaterDiffusivityResult(d={} {})",
            self.d.magnitude_si, self.d.unit
        )
    }
}

impl From<&Co2WaterDiffusivityResult> for PyCo2WaterDiffusivityResult {
    fn from(r: &Co2WaterDiffusivityResult) -> Self {
        Self {
            d: PyQty {
                magnitude_si: r.d.value,
                unit: "m**2/s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.parachor_surface_tension`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ParachorSurfaceTensionResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyParachorSurfaceTensionResult {
    /// The surface tension, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub sigma: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyParachorSurfaceTensionResult {
    fn __repr__(&self) -> String {
        format!(
            "ParachorSurfaceTensionResult(sigma={} {})",
            self.sigma.magnitude_si, self.sigma.unit
        )
    }
}

impl From<&ParachorSurfaceTensionResult> for PyParachorSurfaceTensionResult {
    fn from(r: &ParachorSurfaceTensionResult) -> Self {
        Self {
            sigma: PyQty {
                magnitude_si: r.sigma.value,
                unit: "N/m".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.viscosity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ViscosityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyViscosityResult {
    /// The liquid viscosity, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub mu: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyViscosityResult {
    fn __repr__(&self) -> String {
        format!(
            "ViscosityResult(mu={} {})",
            self.mu.magnitude_si, self.mu.unit
        )
    }
}

impl From<&ViscosityResult> for PyViscosityResult {
    fn from(r: &ViscosityResult) -> Self {
        Self {
            mu: PyQty {
                magnitude_si: r.mu.value,
                unit: "Pa*s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.thermal_conductivity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ThermalConductivityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyThermalConductivityResult {
    /// The liquid thermal conductivity, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub k: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyThermalConductivityResult {
    fn __repr__(&self) -> String {
        format!(
            "ThermalConductivityResult(k={} {})",
            self.k.magnitude_si, self.k.unit
        )
    }
}

impl From<&ThermalConductivityResult> for PyThermalConductivityResult {
    fn from(r: &ThermalConductivityResult) -> Self {
        Self {
            k: PyQty {
                magnitude_si: r.k.value,
                unit: "W/(m*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.wilke_viscosity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "WilkeViscosityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyWilkeViscosityResult {
    /// The gas mixture dynamic viscosity, as an SI magnitude and display unit.
    #[pyo3(get)]
    pub mu: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyWilkeViscosityResult {
    fn __repr__(&self) -> String {
        format!(
            "WilkeViscosityResult(mu={} {})",
            self.mu.magnitude_si, self.mu.unit
        )
    }
}

impl From<&WilkeViscosityResult> for PyWilkeViscosityResult {
    fn from(r: &WilkeViscosityResult) -> Self {
        Self {
            mu: PyQty {
                magnitude_si: r.mu.value,
                unit: "Pa*s".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.nrtl_activity_coefficients`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "NrtlActivityCoefficientsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyNrtlActivityCoefficientsResult {
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyNrtlActivityCoefficientsResult {
    fn __repr__(&self) -> String {
        format!(
            "NrtlActivityCoefficientsResult(ln_gamma={:?}, gamma={:?})",
            self.ln_gamma, self.gamma
        )
    }
}

impl From<&NrtlActivityCoefficientsResult> for PyNrtlActivityCoefficientsResult {
    fn from(r: &NrtlActivityCoefficientsResult) -> Self {
        Self {
            ln_gamma: r.ln_gamma.clone(),
            gamma: r.gamma.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.umrpr_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "UmrprAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyUmrprAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyUmrprAlphaResult {
    fn __repr__(&self) -> String {
        format!("UmrprAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&UmrprAlphaResult> for PyUmrprAlphaResult {
    fn from(r: &UmrprAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.unifac_activity_coefficients`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "UnifacActivityCoefficientsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyUnifacActivityCoefficientsResult {
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyUnifacActivityCoefficientsResult {
    fn __repr__(&self) -> String {
        format!(
            "UnifacActivityCoefficientsResult(ln_gamma={:?}, gamma={:?})",
            self.ln_gamma, self.gamma
        )
    }
}

impl From<&UnifacActivityCoefficientsResult> for PyUnifacActivityCoefficientsResult {
    fn from(r: &UnifacActivityCoefficientsResult) -> Self {
        Self {
            ln_gamma: r.ln_gamma.clone(),
            gamma: r.gamma.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.van_laar_acid_activity_coefficients`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "VanLaarAcidActivityCoefficientsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyVanLaarAcidActivityCoefficientsResult {
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyVanLaarAcidActivityCoefficientsResult {
    fn __repr__(&self) -> String {
        format!(
            "VanLaarAcidActivityCoefficientsResult(ln_gamma={:?}, gamma={:?})",
            self.ln_gamma, self.gamma
        )
    }
}

impl From<&VanLaarAcidActivityCoefficientsResult> for PyVanLaarAcidActivityCoefficientsResult {
    fn from(r: &VanLaarAcidActivityCoefficientsResult) -> Self {
        Self {
            ln_gamma: r.ln_gamma.clone(),
            gamma: r.gamma.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.uniquac_activity_coefficients`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "UniquacActivityCoefficientsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyUniquacActivityCoefficientsResult {
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyUniquacActivityCoefficientsResult {
    fn __repr__(&self) -> String {
        format!(
            "UniquacActivityCoefficientsResult(ln_gamma={:?}, gamma={:?})",
            self.ln_gamma, self.gamma
        )
    }
}

impl From<&UniquacActivityCoefficientsResult> for PyUniquacActivityCoefficientsResult {
    fn from(r: &UniquacActivityCoefficientsResult) -> Self {
        Self {
            ln_gamma: r.ln_gamma.clone(),
            gamma: r.gamma.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.wilson_activity_coefficients`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "WilsonActivityCoefficientsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyWilsonActivityCoefficientsResult {
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyWilsonActivityCoefficientsResult {
    fn __repr__(&self) -> String {
        format!(
            "WilsonActivityCoefficientsResult(ln_gamma={:?}, gamma={:?})",
            self.ln_gamma, self.gamma
        )
    }
}

impl From<&WilsonActivityCoefficientsResult> for PyWilsonActivityCoefficientsResult {
    fn from(r: &WilsonActivityCoefficientsResult) -> Self {
        Self {
            ln_gamma: r.ln_gamma.clone(),
            gamma: r.gamma.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.soreide_whitson_alpha`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SoreideWhitsonAlphaResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySoreideWhitsonAlphaResult {
    /// The temperature-dependent alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySoreideWhitsonAlphaResult {
    fn __repr__(&self) -> String {
        format!("SoreideWhitsonAlphaResult(alpha={})", self.alpha)
    }
}

impl From<&SoreideWhitsonAlphaResult> for PySoreideWhitsonAlphaResult {
    fn from(r: &SoreideWhitsonAlphaResult) -> Self {
        Self {
            alpha: r.alpha,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.srk_alpha_ab`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SrkAlphaAbResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySrkAlphaAbResult {
    /// The Soave alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// `A = a*alpha*P/(R**2*T**2)`. Dimensionless.
    #[pyo3(get)]
    pub a_reduced: f64,
    /// `B = b*P/(R*T)`. Dimensionless.
    #[pyo3(get)]
    pub b_reduced: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySrkAlphaAbResult {
    fn __repr__(&self) -> String {
        format!(
            "SrkAlphaAbResult(alpha={}, a_reduced={}, b_reduced={}, {} warning(s))",
            self.alpha,
            self.a_reduced,
            self.b_reduced,
            self.warnings.len()
        )
    }
}

impl From<&SrkAlphaAbResult> for PySrkAlphaAbResult {
    fn from(r: &SrkAlphaAbResult) -> Self {
        Self {
            alpha: r.alpha,
            a_reduced: r.a_reduced,
            b_reduced: r.b_reduced,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.srk_z_factor`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SrkZFactorResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySrkZFactorResult {
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
impl PySrkZFactorResult {
    fn __repr__(&self) -> String {
        format!(
            "SrkZFactorResult(z_min={}, z_max={}, root_structure={}, {} warning(s))",
            self.z_min,
            self.z_max,
            self.root_structure,
            self.warnings.len()
        )
    }
}

impl From<&SrkZFactorResult> for PySrkZFactorResult {
    fn from(r: &SrkZFactorResult) -> Self {
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

/// Result of `eos.srk_departure`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SrkDepartureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySrkDepartureResult {
    /// The logarithm of the fugacity coefficient. Dimensionless.
    #[pyo3(get)]
    pub ln_phi: f64,
    /// The departure enthalpy over `R*T`. Dimensionless.
    #[pyo3(get)]
    pub h_dep_rt: f64,
    /// The departure entropy over `R`. Dimensionless.
    #[pyo3(get)]
    pub s_dep_r: f64,
    /// The departure heat capacity over `R`. Dimensionless.
    #[pyo3(get)]
    pub cp_dep_r: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySrkDepartureResult {
    fn __repr__(&self) -> String {
        format!(
            "SrkDepartureResult(ln_phi={}, h_dep_rt={}, s_dep_r={}, cp_dep_r={}, {} warning(s))",
            self.ln_phi,
            self.h_dep_rt,
            self.s_dep_r,
            self.cp_dep_r,
            self.warnings.len()
        )
    }
}

impl From<&SrkDepartureResult> for PySrkDepartureResult {
    fn from(r: &SrkDepartureResult) -> Self {
        Self {
            ln_phi: r.ln_phi,
            h_dep_rt: r.h_dep_rt,
            s_dep_r: r.s_dep_r,
            cp_dep_r: r.cp_dep_r,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.rk_alpha_ab`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "RkAlphaAbResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyRkAlphaAbResult {
    /// The Redlich-Kwong alpha function. Dimensionless.
    #[pyo3(get)]
    pub alpha: f64,
    /// `A = a*alpha*P/(R**2*T**2)`. Dimensionless.
    #[pyo3(get)]
    pub a_reduced: f64,
    /// `B = b*P/(R*T)`. Dimensionless.
    #[pyo3(get)]
    pub b_reduced: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyRkAlphaAbResult {
    fn __repr__(&self) -> String {
        format!(
            "RkAlphaAbResult(alpha={}, a_reduced={}, b_reduced={}, {} warning(s))",
            self.alpha,
            self.a_reduced,
            self.b_reduced,
            self.warnings.len()
        )
    }
}

impl From<&RkAlphaAbResult> for PyRkAlphaAbResult {
    fn from(r: &RkAlphaAbResult) -> Self {
        Self {
            alpha: r.alpha,
            a_reduced: r.a_reduced,
            b_reduced: r.b_reduced,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.rk_departure`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "RkDepartureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyRkDepartureResult {
    /// The logarithm of the fugacity coefficient. Dimensionless.
    #[pyo3(get)]
    pub ln_phi: f64,
    /// The departure enthalpy over `R*T`. Dimensionless.
    #[pyo3(get)]
    pub h_dep_rt: f64,
    /// The departure entropy over `R`. Dimensionless.
    #[pyo3(get)]
    pub s_dep_r: f64,
    /// The departure heat capacity over `R`. Dimensionless.
    #[pyo3(get)]
    pub cp_dep_r: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyRkDepartureResult {
    fn __repr__(&self) -> String {
        format!(
            "RkDepartureResult(ln_phi={}, h_dep_rt={}, s_dep_r={}, cp_dep_r={}, {} warning(s))",
            self.ln_phi,
            self.h_dep_rt,
            self.s_dep_r,
            self.cp_dep_r,
            self.warnings.len()
        )
    }
}

impl From<&RkDepartureResult> for PyRkDepartureResult {
    fn from(r: &RkDepartureResult) -> Self {
        Self {
            ln_phi: r.ln_phi,
            h_dep_rt: r.h_dep_rt,
            s_dep_r: r.s_dep_r,
            cp_dep_r: r.cp_dep_r,
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
            cp: PyQty {
                magnitude_si: r.cp.value,
                unit: "J/(mol*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.molar_enthalpy_entropy`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MolarEnthalpyEntropyResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMolarEnthalpyEntropyResult {
    /// The molar enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The molar entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The ideal-gas part of the enthalpy.
    #[pyo3(get)]
    pub h_ideal: PyQty,
    /// The ideal-gas part of the entropy.
    #[pyo3(get)]
    pub s_ideal: PyQty,
    /// The residual enthalpy.
    #[pyo3(get)]
    pub h_departure: PyQty,
    /// The residual entropy.
    #[pyo3(get)]
    pub s_departure: PyQty,
    /// The composition-weighted average of the components' `psi`.
    #[pyo3(get)]
    pub psi_bar: f64,
    /// The molar heat capacity at constant pressure.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The ideal-gas part of the heat capacity.
    #[pyo3(get)]
    pub cp_ideal: PyQty,
    /// The residual heat capacity.
    #[pyo3(get)]
    pub cp_departure: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMolarEnthalpyEntropyResult {
    fn __repr__(&self) -> String {
        format!(
            "MolarEnthalpyEntropyResult(h={} {}, s={} {})",
            self.h.magnitude_si, self.h.unit, self.s.magnitude_si, self.s.unit
        )
    }
}

impl From<&MolarEnthalpyEntropyResult> for PyMolarEnthalpyEntropyResult {
    fn from(r: &MolarEnthalpyEntropyResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            h_ideal: qty(r.h_ideal.value, "J/mol"),
            s_ideal: qty(r.s_ideal.value, "J/(mol*K)"),
            h_departure: qty(r.h_departure.value, "J/mol"),
            s_departure: qty(r.s_departure.value, "J/(mol*K)"),
            psi_bar: r.psi_bar,
            cp: qty(r.cp.value, "J/(mol*K)"),
            cp_ideal: qty(r.cp_ideal.value, "J/(mol*K)"),
            cp_departure: qty(r.cp_departure.value, "J/(mol*K)"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.bwrs_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "BwrsPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyBwrsPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The fugacity coefficients, as logarithms, one per component.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The residual enthalpy.
    #[pyo3(get)]
    pub h_res: PyQty,
    /// The residual entropy.
    #[pyo3(get)]
    pub s_res: PyQty,
    /// The residual isobaric heat capacity.
    #[pyo3(get)]
    pub cp_res: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyBwrsPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "BwrsPhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&BwrsPhaseResult> for PyBwrsPhaseResult {
    fn from(r: &BwrsPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            h_res: qty(r.h_res.value, "J/mol"),
            s_res: qty(r.s_res.value, "J/(mol*K)"),
            cp_res: qty(r.cp_res.value, "J/(mol*K)"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.ammonia_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "AmmoniaPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyAmmoniaPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyAmmoniaPhaseResult {
    fn __repr__(&self) -> String {
        format!("AmmoniaPhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&AmmoniaPhaseResult> for PyAmmoniaPhaseResult {
    fn from(r: &AmmoniaPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.co2_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "Co2PhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyCo2PhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyCo2PhaseResult {
    fn __repr__(&self) -> String {
        format!("Co2PhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&Co2PhaseResult> for PyCo2PhaseResult {
    fn from(r: &Co2PhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.water_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "WaterPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyWaterPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyWaterPhaseResult {
    fn __repr__(&self) -> String {
        format!("WaterPhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&WaterPhaseResult> for PyWaterPhaseResult {
    fn from(r: &WaterPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.argon_solid_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ArgonSolidPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyArgonSolidPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyArgonSolidPhaseResult {
    fn __repr__(&self) -> String {
        format!("ArgonSolidPhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&ArgonSolidPhaseResult> for PyArgonSolidPhaseResult {
    fn from(r: &ArgonSolidPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.parahydrogen_solid_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ParahydrogenSolidPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyParahydrogenSolidPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyParahydrogenSolidPhaseResult {
    fn __repr__(&self) -> String {
        format!("ParahydrogenSolidPhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&ParahydrogenSolidPhaseResult> for PyParahydrogenSolidPhaseResult {
    fn from(r: &ParahydrogenSolidPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.eos_cg_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "EosCgPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyEosCgPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyEosCgPhaseResult {
    fn __repr__(&self) -> String {
        format!("EosCgPhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&EosCgPhaseResult> for PyEosCgPhaseResult {
    fn from(r: &EosCgPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.gerg2008_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "Gerg2008PhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGerg2008PhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGerg2008PhaseResult {
    fn __repr__(&self) -> String {
        format!("Gerg2008PhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&Gerg2008PhaseResult> for PyGerg2008PhaseResult {
    fn from(r: &Gerg2008PhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.helium_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HeliumPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHeliumPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHeliumPhaseResult {
    fn __repr__(&self) -> String {
        format!("HeliumPhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&HeliumPhaseResult> for PyHeliumPhaseResult {
    fn from(r: &HeliumPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.hydrogen_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HydrogenPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrogenPhaseResult {
    /// The compressibility factor.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The internal energy.
    #[pyo3(get)]
    pub u: PyQty,
    /// The enthalpy.
    #[pyo3(get)]
    pub h: PyQty,
    /// The entropy.
    #[pyo3(get)]
    pub s: PyQty,
    /// The isochoric heat capacity.
    #[pyo3(get)]
    pub cv: PyQty,
    /// The isobaric heat capacity.
    #[pyo3(get)]
    pub cp: PyQty,
    /// The Gibbs energy.
    #[pyo3(get)]
    pub g: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHydrogenPhaseResult {
    fn __repr__(&self) -> String {
        format!("HydrogenPhaseResult(z_factor={})", self.z_factor)
    }
}

impl From<&HydrogenPhaseResult> for PyHydrogenPhaseResult {
    fn from(r: &HydrogenPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            u: qty(r.u.value, "J/mol"),
            h: qty(r.h.value, "J/mol"),
            s: qty(r.s.value, "J/(mol*K)"),
            cv: qty(r.cv.value, "J/(mol*K)"),
            cp: qty(r.cp.value, "J/(mol*K)"),
            g: qty(r.g.value, "J/mol"),
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

/// Result of `eos.bubble_temperature` or `eos.dew_temperature`, transported.
///
/// The temperature twin of [`PyPhaseBoundaryResult`]: same shape, but the answer is a
/// temperature rather than a pressure, so the scalar field is `temperature` in kelvin.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PhaseBoundaryTemperatureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPhaseBoundaryTemperatureResult {
    /// The boundary temperature.
    #[pyo3(get)]
    pub temperature: PyQty,
    /// The incipient phase's composition.
    #[pyo3(get)]
    pub incipient: Vec<f64>,
    /// K-values at the converged temperature.
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
    /// Temperature updates taken.
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
impl PyPhaseBoundaryTemperatureResult {
    fn __repr__(&self) -> String {
        format!(
            "PhaseBoundaryTemperatureResult(temperature={} {}, {} iteration(s))",
            self.temperature.magnitude_si, self.temperature.unit, self.iterations
        )
    }
}

impl From<&BubbleTemperatureResult> for PyPhaseBoundaryTemperatureResult {
    fn from(r: &BubbleTemperatureResult) -> Self {
        Self {
            temperature: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
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

impl From<&DewTemperatureResult> for PyPhaseBoundaryTemperatureResult {
    fn from(r: &DewTemperatureResult) -> Self {
        Self {
            temperature: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
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

/// Result of `eos.ph_flash`, transported.
///
/// The temperature crosses as a quantity like every other, and the phase as the spec's
/// spelling, so the adapter rebuilds the enum and a caller compares `Phase.TWO_PHASE`
/// regardless of which backend answered.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PhFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `T` is the symbol the spec and the Python result both use
pub struct PyPhFlashResult {
    /// The temperature that satisfies the enthalpy.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction, or `None` for a single-phase feed.
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
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Bisection steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `|H(T) - H_target| / max(|H_target|, 1)` at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPhFlashResult {
    fn __repr__(&self) -> String {
        match self.beta {
            Some(beta) => format!(
                "PhFlashResult(T={} K, phase={}, beta={beta})",
                self.T.magnitude_si, self.phase
            ),
            None => format!(
                "PhFlashResult(T={} K, phase={}, {} iteration(s))",
                self.T.magnitude_si, self.phase, self.iterations
            ),
        }
    }
}

impl From<&PhFlashResult> for PyPhFlashResult {
    fn from(r: &PhFlashResult) -> Self {
        Self {
            T: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            phase: r.phase.as_str().to_string(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.ps_flash`, transported.
///
/// The same fields as `PyPhFlashResult`, which is the point: the two models differ in
/// which property they invert and agree on everything else, so a caller who has read one
/// already knows how to read the other.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PsFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `T` is the symbol the spec and the Python result both use
pub struct PyPsFlashResult {
    /// The temperature that satisfies the entropy.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction, or `None` for a single-phase feed.
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
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Bisection steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `|S(T) - S_target| / max(|S_target|, 1)` at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPsFlashResult {
    fn __repr__(&self) -> String {
        match self.beta {
            Some(beta) => format!(
                "PsFlashResult(T={} K, phase={}, beta={beta})",
                self.T.magnitude_si, self.phase
            ),
            None => format!(
                "PsFlashResult(T={} K, phase={}, {} iteration(s))",
                self.T.magnitude_si, self.phase, self.iterations
            ),
        }
    }
}

impl From<&PsFlashResult> for PyPsFlashResult {
    fn from(r: &PsFlashResult) -> Self {
        Self {
            T: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            phase: r.phase.as_str().to_string(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.tv_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TvFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `P` is the symbol the spec and the Python result both use
pub struct PyTvFlashResult {
    /// The pressure that satisfies the volume.
    #[pyo3(get)]
    pub P: PyQty,
    /// The vapour fraction, or `None` for a single-phase feed.
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
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Newton steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `V(P) - V_target` at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTvFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "TvFlashResult(P={} Pa, phase={}, beta={:?})",
            self.P.magnitude_si, self.phase, self.beta
        )
    }
}

impl From<&TvFlashResult> for PyTvFlashResult {
    fn from(r: &TvFlashResult) -> Self {
        Self {
            P: PyQty {
                magnitude_si: r.pressure.value,
                unit: "Pa".to_string(),
            },
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            phase: r.phase.as_str().to_string(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pv_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PvFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `T` is the symbol the spec and the Python result both use
pub struct PyPvFlashResult {
    /// The temperature that satisfies the volume.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction, or `None` for a single-phase feed.
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
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Newton steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `V(T) - V_target` at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPvFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "PvFlashResult(T={} K, phase={}, beta={:?})",
            self.T.magnitude_si, self.phase, self.beta
        )
    }
}

impl From<&PvFlashResult> for PyPvFlashResult {
    fn from(r: &PvFlashResult) -> Self {
        Self {
            T: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            phase: r.phase.as_str().to_string(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// The pressure-specified flashes `eos.th_flash`, `eos.ts_flash` and `eos.tu_flash`
/// share this transport shape.
///
/// `#[rustfmt::skip]`: rustfmt cannot settle the indentation of the `#[pyclass]`
/// attribute's arguments here - the `name = $py_name` metavariable breaks its
/// indentation calculation, and each run re-indents by another level. Skipping the
/// whole macro keeps the three structs it expands to formatted once, by hand.
#[rustfmt::skip]
macro_rules! py_pressure_flash_result {
    ($name:ident, $py_name:literal) => {
        #[pyclass(
            frozen,
            skip_from_py_object,
            module = "azoth._core",
            name = $py_name
        )]
        #[derive(Debug, Clone, PartialEq)]
        #[allow(non_snake_case)] // `P` is the symbol the spec and the Python result both use
        pub struct $name {
            /// The pressure that satisfies the property.
            #[pyo3(get)]
            pub P: PyQty,
            /// The vapour fraction, or `None` for a single-phase feed.
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
            /// The spec's spelling of the phase.
            #[pyo3(get)]
            pub phase: String,
            /// The liquid root of the cubic.
            #[pyo3(get)]
            pub z_liquid: f64,
            /// The vapour root.
            #[pyo3(get)]
            pub z_vapour: f64,
            /// Newton steps taken.
            #[pyo3(get)]
            pub iterations: u32,
            /// The relative property residual at the answer.
            #[pyo3(get)]
            pub residual: f64,
            /// Caveats, deduplicated.
            #[pyo3(get)]
            pub warnings: Vec<PyWarning>,
        }

        #[pymethods]
        impl $name {
            fn __repr__(&self) -> String {
                format!(
                    concat!($py_name, "(P={} Pa, phase={}, beta={:?})"),
                    self.P.magnitude_si, self.phase, self.beta
                )
            }
        }
    };
}

py_pressure_flash_result!(PyThFlashResult, "ThFlashResult");
py_pressure_flash_result!(PyTsFlashResult, "TsFlashResult");
py_pressure_flash_result!(PyTuFlashResult, "TuFlashResult");

macro_rules! py_pressure_flash_from {
    ($name:ident, $src:ident) => {
        impl From<&$src> for $name {
            fn from(r: &$src) -> Self {
                Self {
                    P: PyQty {
                        magnitude_si: r.pressure.value,
                        unit: "Pa".to_string(),
                    },
                    beta: r.beta,
                    x: r.x.clone(),
                    y: r.y.clone(),
                    k: r.k.clone(),
                    phase: r.phase.as_str().to_string(),
                    z_liquid: r.z_liquid,
                    z_vapour: r.z_vapour,
                    iterations: r.iterations,
                    residual: r.residual,
                    warnings: transport(&r.warnings),
                }
            }
        }
    };
}

py_pressure_flash_from!(PyThFlashResult, ThFlashResult);
py_pressure_flash_from!(PyTsFlashResult, TsFlashResult);
py_pressure_flash_from!(PyTuFlashResult, TuFlashResult);

/// Result of `eos.pu_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PuFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `T` is the symbol the spec and the Python result both use
pub struct PyPuFlashResult {
    /// The temperature that satisfies the internal energy.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction, or `None` for a single-phase feed.
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
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Newton steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The relative internal-energy residual at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPuFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "PuFlashResult(T={} K, phase={}, beta={:?})",
            self.T.magnitude_si, self.phase, self.beta
        )
    }
}

impl From<&PuFlashResult> for PyPuFlashResult {
    fn from(r: &PuFlashResult) -> Self {
        Self {
            T: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            phase: r.phase.as_str().to_string(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.vu_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "VuFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `P` and `T` are the symbols the spec and the Python result both use
pub struct PyVuFlashResult {
    /// The pressure that satisfies the volume and internal energy.
    #[pyo3(get)]
    pub P: PyQty,
    /// The temperature that satisfies the volume and internal energy.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction, or `None` for a single-phase feed.
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
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Newton steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The larger of the relative volume and internal-energy residuals.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyVuFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "VuFlashResult(P={} Pa, T={} K, phase={}, beta={:?})",
            self.P.magnitude_si, self.T.magnitude_si, self.phase, self.beta
        )
    }
}

impl From<&VuFlashResult> for PyVuFlashResult {
    fn from(r: &VuFlashResult) -> Self {
        Self {
            P: PyQty {
                magnitude_si: r.pressure.value,
                unit: "Pa".to_string(),
            },
            T: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            phase: r.phase.as_str().to_string(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.stability_test`, transported.
///
/// The first result here whose fields are all vectors or lists, and the first whose
/// enum is *the* answer rather than a qualifier on one: `verdict` is what the model
/// was asked, and the two distances are the evidence for it. It crosses as the
/// spec's spelling, like `phase` does, and the adapter rebuilds the enum - so a
/// caller compares `StabilityVerdict.UNSTABLE` regardless of which backend answered.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "StabilityTestResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyStabilityTestResult {
    /// The spec's spelling of the verdict.
    #[pyo3(get)]
    pub verdict: String,
    /// The tangent-plane distance at each trial's stationary point, in trial order.
    #[pyo3(get)]
    pub tm: Vec<f64>,
    /// The stationary-point composition of each trial, in the same order.
    #[pyo3(get)]
    pub w: Vec<Vec<f64>>,
    /// Iterations each trial took, in the same order.
    #[pyo3(get)]
    pub iterations: Vec<u32>,
    /// The smallest `T / Tc_i` over the components.
    #[pyo3(get)]
    pub min_t_over_tc: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyStabilityTestResult {
    fn __repr__(&self) -> String {
        // The distances first: they are the evidence for the verdict, and a reader
        // seeing `stable` alone would not know whether the trials found nothing or
        // found separate stationary points.
        format!(
            "StabilityTestResult({}, tm={:?}, {} + {} trial(s))",
            self.verdict,
            self.tm,
            self.iterations.first().copied().unwrap_or_default(),
            self.iterations.get(1).copied().unwrap_or_default()
        )
    }
}

impl From<&StabilityTestResult> for PyStabilityTestResult {
    fn from(r: &StabilityTestResult) -> Self {
        Self {
            verdict: r.verdict.as_str().to_string(),
            tm: r.tm.clone(),
            w: r.w.clone(),
            iterations: r.iterations.clone(),
            min_t_over_tc: r.min_t_over_tc,
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

/// The base dimensions, in the order the exponent tuples above are written in.
#[pyfunction]
#[must_use]
pub fn unit_slots() -> Vec<String> {
    azoth_core::unit_vocab_gen::SLOTS
        .iter()
        .map(|slot| (*slot).to_string())
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

/// Result of `eos.pt_phase_envelope`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PhaseEnvelopeResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPhaseEnvelopeResult {
    /// The dew-point temperatures, in kelvin.
    #[pyo3(get)]
    pub dew_temperature: Vec<f64>,
    /// The dew-point pressures, in pascal.
    #[pyo3(get)]
    pub dew_pressure: Vec<f64>,
    /// The bubble-point temperatures, in kelvin.
    #[pyo3(get)]
    pub bubble_temperature: Vec<f64>,
    /// The bubble-point pressures, in pascal.
    #[pyo3(get)]
    pub bubble_pressure: Vec<f64>,
    /// The temperature at the cricondenbar.
    #[pyo3(get)]
    pub cricondenbar_temperature: PyQty,
    /// The cricondenbar pressure.
    #[pyo3(get)]
    pub cricondenbar_pressure: PyQty,
    /// The cricondentherm temperature.
    #[pyo3(get)]
    pub cricondentherm_temperature: PyQty,
    /// The pressure at the cricondentherm.
    #[pyo3(get)]
    pub cricondentherm_pressure: PyQty,
    /// The critical temperature.
    #[pyo3(get)]
    pub critical_temperature: PyQty,
    /// The critical pressure.
    #[pyo3(get)]
    pub critical_pressure: PyQty,
    /// The number of continuation points traced.
    #[pyo3(get)]
    pub iterations: u32,
    /// The distance between the two branches' endpoints.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

impl From<&PtPhaseEnvelopeResult> for PyPhaseEnvelopeResult {
    fn from(r: &PtPhaseEnvelopeResult) -> Self {
        Self {
            dew_temperature: r.dew_temperature.clone(),
            dew_pressure: r.dew_pressure.clone(),
            bubble_temperature: r.bubble_temperature.clone(),
            bubble_pressure: r.bubble_pressure.clone(),
            cricondenbar_temperature: PyQty {
                magnitude_si: r.cricondenbar_temperature.value,
                unit: "K".to_string(),
            },
            cricondenbar_pressure: PyQty {
                magnitude_si: r.cricondenbar_pressure.value,
                unit: "Pa".to_string(),
            },
            cricondentherm_temperature: PyQty {
                magnitude_si: r.cricondentherm_temperature.value,
                unit: "K".to_string(),
            },
            cricondentherm_pressure: PyQty {
                magnitude_si: r.cricondentherm_pressure.value,
                unit: "Pa".to_string(),
            },
            critical_temperature: PyQty {
                magnitude_si: r.critical_temperature.value,
                unit: "K".to_string(),
            },
            critical_pressure: PyQty {
                magnitude_si: r.critical_pressure.value,
                unit: "Pa".to_string(),
            },
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
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
        PrLeeKeslerAlphaResult::CALC_ID => PrLeeKeslerAlphaResult::FIELDS.to_vec(),
        Matcop5PrumrAlphaResult::CALC_ID => Matcop5PrumrAlphaResult::FIELDS.to_vec(),
        MatcopAlphaResult::CALC_ID => MatcopAlphaResult::FIELDS.to_vec(),
        MatcopPrAlphaResult::CALC_ID => MatcopPrAlphaResult::FIELDS.to_vec(),
        MatcopPrumrAlphaResult::CALC_ID => MatcopPrumrAlphaResult::FIELDS.to_vec(),
        MatcopPrumrNewAlphaResult::CALC_ID => MatcopPrumrNewAlphaResult::FIELDS.to_vec(),
        MollerupAlphaResult::CALC_ID => MollerupAlphaResult::FIELDS.to_vec(),
        PrAlphaAbResult::CALC_ID => PrAlphaAbResult::FIELDS.to_vec(),
        PrDaneshAlphaResult::CALC_ID => PrDaneshAlphaResult::FIELDS.to_vec(),
        PrDelft1998AlphaResult::CALC_ID => PrDelft1998AlphaResult::FIELDS.to_vec(),
        PrGassem2001AlphaResult::CALC_ID => PrGassem2001AlphaResult::FIELDS.to_vec(),
        PrZFactorResult::CALC_ID => PrZFactorResult::FIELDS.to_vec(),
        PrsvKappaResult::CALC_ID => PrsvKappaResult::FIELDS.to_vec(),
        PrDepartureResult::CALC_ID => PrDepartureResult::FIELDS.to_vec(),
        SrkKappaResult::CALC_ID => SrkKappaResult::FIELDS.to_vec(),
        SrkAlphaAbResult::CALC_ID => SrkAlphaAbResult::FIELDS.to_vec(),
        SrkZFactorResult::CALC_ID => SrkZFactorResult::FIELDS.to_vec(),
        SrkDepartureResult::CALC_ID => SrkDepartureResult::FIELDS.to_vec(),
        RkAlphaAbResult::CALC_ID => RkAlphaAbResult::FIELDS.to_vec(),
        RkDepartureResult::CALC_ID => RkDepartureResult::FIELDS.to_vec(),
        Pr78KappaResult::CALC_ID => Pr78KappaResult::FIELDS.to_vec(),
        TwuKappaResult::CALC_ID => TwuKappaResult::FIELDS.to_vec(),
        TwucoonAlphaResult::CALC_ID => TwucoonAlphaResult::FIELDS.to_vec(),
        TwucoonParamAlphaResult::CALC_ID => TwucoonParamAlphaResult::FIELDS.to_vec(),
        TwucoonStatoilAlphaResult::CALC_ID => TwucoonStatoilAlphaResult::FIELDS.to_vec(),
        Vdw1fMixBinaryResult::CALC_ID => Vdw1fMixBinaryResult::FIELDS.to_vec(),
        RachfordRiceBinaryResult::CALC_ID => RachfordRiceBinaryResult::FIELDS.to_vec(),
        PrMolarVolumeResult::CALC_ID => PrMolarVolumeResult::FIELDS.to_vec(),
        PrMassDensityResult::CALC_ID => PrMassDensityResult::FIELDS.to_vec(),
        PrPenelouxShiftResult::CALC_ID => PrPenelouxShiftResult::FIELDS.to_vec(),
        SrkPenelouxShiftResult::CALC_ID => SrkPenelouxShiftResult::FIELDS.to_vec(),
        HeatOfVaporizationResult::CALC_ID => HeatOfVaporizationResult::FIELDS.to_vec(),
        LiquidHeatCapacityResult::CALC_ID => LiquidHeatCapacityResult::FIELDS.to_vec(),
        AntoineVaporPressureResult::CALC_ID => AntoineVaporPressureResult::FIELDS.to_vec(),
        RackettMolarVolumeResult::CALC_ID => RackettMolarVolumeResult::FIELDS.to_vec(),
        CostaldMolarVolumeResult::CALC_ID => CostaldMolarVolumeResult::FIELDS.to_vec(),
        ChungViscosityResult::CALC_ID => ChungViscosityResult::FIELDS.to_vec(),
        ChungConductivityResult::CALC_ID => ChungConductivityResult::FIELDS.to_vec(),
        WilkeViscosityResult::CALC_ID => WilkeViscosityResult::FIELDS.to_vec(),
        MasonSaxenaConductivityResult::CALC_ID => MasonSaxenaConductivityResult::FIELDS.to_vec(),
        NrtlActivityCoefficientsResult::CALC_ID => NrtlActivityCoefficientsResult::FIELDS.to_vec(),
        UnifacActivityCoefficientsResult::CALC_ID => {
            UnifacActivityCoefficientsResult::FIELDS.to_vec()
        }
        UniquacActivityCoefficientsResult::CALC_ID => {
            UniquacActivityCoefficientsResult::FIELDS.to_vec()
        }
        VanLaarAcidActivityCoefficientsResult::CALC_ID => {
            VanLaarAcidActivityCoefficientsResult::FIELDS.to_vec()
        }
        WilsonActivityCoefficientsResult::CALC_ID => {
            WilsonActivityCoefficientsResult::FIELDS.to_vec()
        }
        TynCalusDiffusivityResult::CALC_ID => TynCalusDiffusivityResult::FIELDS.to_vec(),
        UmrprAlphaResult::CALC_ID => UmrprAlphaResult::FIELDS.to_vec(),
        WilkeChangDiffusivityResult::CALC_ID => WilkeChangDiffusivityResult::FIELDS.to_vec(),
        HaydukMinhasDiffusivityResult::CALC_ID => HaydukMinhasDiffusivityResult::FIELDS.to_vec(),
        SchwartzentruberAlphaResult::CALC_ID => SchwartzentruberAlphaResult::FIELDS.to_vec(),
        SoreideWhitsonAlphaResult::CALC_ID => SoreideWhitsonAlphaResult::FIELDS.to_vec(),
        SiddiqiLucasDiffusivityResult::CALC_ID => SiddiqiLucasDiffusivityResult::FIELDS.to_vec(),
        Co2WaterDiffusivityResult::CALC_ID => Co2WaterDiffusivityResult::FIELDS.to_vec(),
        ParachorSurfaceTensionResult::CALC_ID => ParachorSurfaceTensionResult::FIELDS.to_vec(),
        // Models. Present here because a result's *shape* is a cross-language
        // contract whether or not its spec calls it a calculation.
        PureSaturationResult::CALC_ID => PureSaturationResult::FIELDS.to_vec(),
        PtFlashResult::CALC_ID => PtFlashResult::FIELDS.to_vec(),
        PtPhaseEnvelopeResult::CALC_ID => PtPhaseEnvelopeResult::FIELDS.to_vec(),
        PhFlashResult::CALC_ID => PhFlashResult::FIELDS.to_vec(),
        PsFlashResult::CALC_ID => PsFlashResult::FIELDS.to_vec(),
        TvFlashResult::CALC_ID => TvFlashResult::FIELDS.to_vec(),
        PvFlashResult::CALC_ID => PvFlashResult::FIELDS.to_vec(),
        ThFlashResult::CALC_ID => ThFlashResult::FIELDS.to_vec(),
        TsFlashResult::CALC_ID => TsFlashResult::FIELDS.to_vec(),
        TuFlashResult::CALC_ID => TuFlashResult::FIELDS.to_vec(),
        PuFlashResult::CALC_ID => PuFlashResult::FIELDS.to_vec(),
        VuFlashResult::CALC_ID => VuFlashResult::FIELDS.to_vec(),
        StabilityTestResult::CALC_ID => StabilityTestResult::FIELDS.to_vec(),
        BubblePressureResult::CALC_ID => BubblePressureResult::FIELDS.to_vec(),
        BubbleTemperatureResult::CALC_ID => BubbleTemperatureResult::FIELDS.to_vec(),
        CriticalPointResult::CALC_ID => CriticalPointResult::FIELDS.to_vec(),
        BwrsPhaseResult::CALC_ID => BwrsPhaseResult::FIELDS.to_vec(),
        AmmoniaPhaseResult::CALC_ID => AmmoniaPhaseResult::FIELDS.to_vec(),
        Co2PhaseResult::CALC_ID => Co2PhaseResult::FIELDS.to_vec(),
        HeliumPhaseResult::CALC_ID => HeliumPhaseResult::FIELDS.to_vec(),
        HydrogenPhaseResult::CALC_ID => HydrogenPhaseResult::FIELDS.to_vec(),
        WaterPhaseResult::CALC_ID => WaterPhaseResult::FIELDS.to_vec(),
        ArgonSolidPhaseResult::CALC_ID => ArgonSolidPhaseResult::FIELDS.to_vec(),
        ParahydrogenSolidPhaseResult::CALC_ID => ParahydrogenSolidPhaseResult::FIELDS.to_vec(),
        EosCgPhaseResult::CALC_ID => EosCgPhaseResult::FIELDS.to_vec(),
        Gerg2008PhaseResult::CALC_ID => Gerg2008PhaseResult::FIELDS.to_vec(),
        DewPressureResult::CALC_ID => DewPressureResult::FIELDS.to_vec(),
        DewTemperatureResult::CALC_ID => DewTemperatureResult::FIELDS.to_vec(),
        IdealGasCpResult::CALC_ID => IdealGasCpResult::FIELDS.to_vec(),
        MolarEnthalpyEntropyResult::CALC_ID => MolarEnthalpyEntropyResult::FIELDS.to_vec(),
        ViscosityResult::CALC_ID => ViscosityResult::FIELDS.to_vec(),
        ThermalConductivityResult::CALC_ID => ThermalConductivityResult::FIELDS.to_vec(),
        // Unit operations. In the same table for the same reason the models are: a
        // result's shape is a cross-language contract whether or not its spec calls it
        // a calculation.
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
        PrLeeKeslerAlphaResult::CALC_ID.to_string(),
        Matcop5PrumrAlphaResult::CALC_ID.to_string(),
        MatcopAlphaResult::CALC_ID.to_string(),
        MatcopPrAlphaResult::CALC_ID.to_string(),
        MatcopPrumrAlphaResult::CALC_ID.to_string(),
        MatcopPrumrNewAlphaResult::CALC_ID.to_string(),
        MollerupAlphaResult::CALC_ID.to_string(),
        PrAlphaAbResult::CALC_ID.to_string(),
        PrDaneshAlphaResult::CALC_ID.to_string(),
        PrDelft1998AlphaResult::CALC_ID.to_string(),
        PrGassem2001AlphaResult::CALC_ID.to_string(),
        PrZFactorResult::CALC_ID.to_string(),
        PrsvKappaResult::CALC_ID.to_string(),
        PrDepartureResult::CALC_ID.to_string(),
        SrkKappaResult::CALC_ID.to_string(),
        SrkAlphaAbResult::CALC_ID.to_string(),
        SrkZFactorResult::CALC_ID.to_string(),
        SrkDepartureResult::CALC_ID.to_string(),
        RkAlphaAbResult::CALC_ID.to_string(),
        RkDepartureResult::CALC_ID.to_string(),
        Pr78KappaResult::CALC_ID.to_string(),
        TwuKappaResult::CALC_ID.to_string(),
        TwucoonAlphaResult::CALC_ID.to_string(),
        TwucoonParamAlphaResult::CALC_ID.to_string(),
        TwucoonStatoilAlphaResult::CALC_ID.to_string(),
        Vdw1fMixBinaryResult::CALC_ID.to_string(),
        RachfordRiceBinaryResult::CALC_ID.to_string(),
        PrMolarVolumeResult::CALC_ID.to_string(),
        PrMassDensityResult::CALC_ID.to_string(),
        PrPenelouxShiftResult::CALC_ID.to_string(),
        SrkPenelouxShiftResult::CALC_ID.to_string(),
        HeatOfVaporizationResult::CALC_ID.to_string(),
        LiquidHeatCapacityResult::CALC_ID.to_string(),
        AntoineVaporPressureResult::CALC_ID.to_string(),
        RackettMolarVolumeResult::CALC_ID.to_string(),
        CostaldMolarVolumeResult::CALC_ID.to_string(),
        ChungViscosityResult::CALC_ID.to_string(),
        ChungConductivityResult::CALC_ID.to_string(),
        TynCalusDiffusivityResult::CALC_ID.to_string(),
        UmrprAlphaResult::CALC_ID.to_string(),
        WilkeChangDiffusivityResult::CALC_ID.to_string(),
        HaydukMinhasDiffusivityResult::CALC_ID.to_string(),
        SchwartzentruberAlphaResult::CALC_ID.to_string(),
        SoreideWhitsonAlphaResult::CALC_ID.to_string(),
        SiddiqiLucasDiffusivityResult::CALC_ID.to_string(),
        Co2WaterDiffusivityResult::CALC_ID.to_string(),
        ParachorSurfaceTensionResult::CALC_ID.to_string(),
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

/// Result of `eos.critical_point`.
///
/// The four state variables of a mixture critical point. The units are strings here,
/// as everywhere on this boundary, and each is the unit the spec declares.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "CriticalPointResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyCriticalPointResult {
    /// The critical temperature.
    #[pyo3(get)]
    pub tc: PyQty,
    /// The critical pressure.
    #[pyo3(get)]
    pub pc: PyQty,
    /// The critical molar volume.
    #[pyo3(get)]
    pub vc: PyQty,
    /// `Pc Vc/(R Tc)`.
    #[pyo3(get)]
    pub z_c: f64,
    /// Outer iterations taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `max(|smallest eigenvalue|, |cubic form|)` at the returned state.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyCriticalPointResult {
    fn __repr__(&self) -> String {
        format!(
            "CriticalPointResult(tc={} K, pc={} Pa, vc={} m**3/mol, z_c={})",
            self.tc.magnitude_si, self.pc.magnitude_si, self.vc.magnitude_si, self.z_c
        )
    }
}

impl From<&CriticalPointResult> for PyCriticalPointResult {
    fn from(r: &CriticalPointResult) -> Self {
        Self {
            tc: PyQty {
                magnitude_si: r.tc.value,
                unit: "K".to_string(),
            },
            pc: PyQty {
                magnitude_si: r.pc.value,
                unit: "Pa".to_string(),
            },
            vc: PyQty {
                magnitude_si: r.vc.value,
                unit: "m**3/mol".to_string(),
            },
            z_c: r.z_c,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
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
    use azoth_core::units::kelvins;
    use azoth_eos::Phase;

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
            beta: Some(0.5),
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
            assert_eq!(number(result, "beta"), 0.5);
            assert_eq!(numbers(result, "x"), vec![0.31, 0.32]);
            assert_eq!(numbers(result, "y"), vec![0.41, 0.42]);
            assert_eq!(numbers(result, "k"), vec![1.51, 1.52]);
            assert_eq!(number(result, "z_liquid"), 0.061);
            assert_eq!(number(result, "z_vapour"), 0.062);
            assert_eq!(number(result, "iterations"), 9.0);
            assert_eq!(number(result, "residual"), 1.5e-11);
        });
    }
}
