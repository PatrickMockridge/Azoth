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
use azoth_eos::EffectiveDiffusionResult as KernelEffectiveDiffusionResult;
use azoth_eos::results::AqueousViscosityResult as KernelAqueousViscosityResult;
use azoth_eos::results::{
    AmmoniaPhaseResult, AntoineVaporPressureResult, ArgonSolidPhaseResult, BubblePressureResult,
    BubbleTemperatureResult, BwrsPhaseResult, CapillaryDewPointResult, ChungConductivityResult,
    ChungViscosityResult, Co2PhaseResult, Co2WaterDiffusivityResult, CostaldMolarVolumeResult,
    CriticalPointResult, DesmukhMatherPhaseResult, DewPressureResult, DewTemperatureResult,
    EosCgPhaseResult, FreezingPointResult, FurstElectrolyteMod2004PhaseResult,
    FurstElectrolytePhaseResult, GeFlashResult, GeNrtlFlashResult, GeNrtlPhaseResult,
    GeUnifacPhaseResult, GeUniquacPhaseResult, GeVanLaarAcidPhaseResult, GeWilsonPhaseResult,
    Gerg2008PhaseResult, HaydukMinhasDiffusivityResult, HeatOfVaporizationResult,
    HeliumPhaseResult, HybridEosGeFlashResult, HydrateEquilibriumLineResult,
    HydrateFormationPressureResult, HydrateFormationTemperatureResult, HydrateFractionResult,
    HydrateInhibitorConcentrationResult, HydrateInhibitorWtResult, HydrogenPhaseResult,
    IapwsHenryLawResult, IdealGasCpResult, KentEisenbergPhaseResult, LiquidHeatCapacityResult,
    MasonSaxenaConductivityResult, Matcop5PrumrAlphaResult, MatcopAlphaResult, MatcopPrAlphaResult,
    MatcopPrumrAlphaResult, MatcopPrumrNewAlphaResult, MolarEnthalpyEntropyResult,
    MollerupAlphaResult, NitricSulfuricAcidVaporPressureResult, NrtlActivityCoefficientsResult,
    ParachorSurfaceTensionResult, ParahydrogenSolidPhaseResult, PcsaftRahmatPhaseResult,
    PhFlashResult, PitzerPhaseResult, Pr78KappaResult, PrAlphaAbResult, PrCpaPhaseResult,
    PrDaneshAlphaResult, PrDelft1998AlphaResult, PrDepartureResult, PrGassem2001AlphaResult,
    PrKappaResult, PrLeeKeslerAlphaResult, PrMassDensityResult, PrMolarVolumeResult,
    PrPenelouxShiftResult, PrZFactorResult, PrsvKappaResult, PsFlashResult, PtFlashResult,
    PtPhaseEnvelopeResult, PuFlashResult, PureSaturationResult, PvFlashResult, PvRefluxFlashResult,
    PvfFlashResult, RachfordRiceBinaryResult, RachfordRiceResult, RackettMolarVolumeResult,
    RkAlphaAbResult, RkDepartureResult, SaftFlashResult, SaftVrMiePhaseResult,
    SaltPrecipitationResult, ScaleSaturationRatioResult, SchwartzentruberAlphaResult,
    SiddiqiLucasDiffusivityResult, SolidFugacityResult, SoreideWhitsonAlphaResult,
    SoreideWhitsonPhaseResult, SrkAlphaAbResult, SrkCpaPhaseResult, SrkDepartureResult,
    SrkKappaResult, SrkPenelouxShiftResult, SrkZFactorResult, StabilityTestResult,
    TbpFractionPropertiesResult, ThFlashResult, ThermalConductivityResult, TpMultiflashResult,
    TpMultiflashWaxResult, TpSolidFlashResult, TsFlashResult, TuFlashResult, TvFlashResult,
    TvFractionFlashResult, TwuKappaResult, TwucoonAlphaResult, TwucoonParamAlphaResult,
    TwucoonStatoilAlphaResult, TynCalusDiffusivityResult, UmrCpaPhaseResult, UmrprAlphaResult,
    UnifacActivityCoefficientsResult, UnifacPsrkActivityCoefficientsResult,
    UnifacUmrpruActivityCoefficientsResult, UniquacActivityCoefficientsResult,
    VanLaarAcidActivityCoefficientsResult, Vdw1fMixBinaryResult, VhFlashResult, ViscosityResult,
    VsFlashResult, VuFlashResult, VuFlashSingleCompResult, WaterPhaseResult,
    WaxSolidFugacityResult, WilkeChangDiffusivityResult, WilkeViscosityResult,
    WilsonActivityCoefficientsResult,
};
use azoth_process::{
    ComponentSplitterResult, CompressorResult, CoolerResult, DistillationColumnResult,
    EjectorResult, ExpanderResult, FilterResult, GasScrubberResult, HeatExchangerResult,
    HeaterResult, ManifoldResult, MixerResult, PumpResult, SeparatorResult,
    ShortcutDistillationColumnResult, SplitterResult, TankResult, ThreePhaseSeparatorResult,
    ThrottlingValveResult, pipe::PipeResult,
};
use azoth_reactions::chemical_equilibrium::ChemicalEquilibriumResult;
use azoth_reactions::equilibrium_constant::EquilibriumConstantResult;
use azoth_reactions::kinetic_rate_law::KineticRateLawResult as KernelKineticRateLawResult;
use azoth_reactions::kinetics::KineticsResult as KernelKineticsResult;
use azoth_reactions::reactive_hybrid_eos_ge_flash::ReactiveHybridEosGeFlashResult;
use azoth_reactions::reactive_ph_flash::ReactivePhFlashResult;
use azoth_reactions::reactive_phase_equilibrium::ReactivePhaseEquilibriumResult;
use azoth_reactions::reactive_tp_flash::ReactiveTpFlashResult;
use azoth_reactions::reference_potentials::ReferencePotentialsResult;
use azoth_standards::Iso6976Result;
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

/// Result of `process.heat_exchanger`, transported.
///
/// Ten fields under this unit operation's own port names, the shape `SeparatorResult` has.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HeatExchangerResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHeatExchangerResult {
    /// Hot outlet molar flow, mol/s.
    #[pyo3(get)]
    pub hot_out_n: PyQty,
    /// Hot outlet composition.
    #[pyo3(get)]
    pub hot_out_z: Vec<f64>,
    /// Hot outlet pressure.
    #[pyo3(get)]
    pub hot_out_p: PyQty,
    /// Hot outlet temperature.
    #[pyo3(get)]
    pub hot_out_t: PyQty,
    /// Hot outlet molar enthalpy.
    #[pyo3(get)]
    pub hot_out_h: PyQty,
    /// Cold outlet molar flow, mol/s.
    #[pyo3(get)]
    pub cold_out_n: PyQty,
    /// Cold outlet composition.
    #[pyo3(get)]
    pub cold_out_z: Vec<f64>,
    /// Cold outlet pressure.
    #[pyo3(get)]
    pub cold_out_p: PyQty,
    /// Cold outlet temperature.
    #[pyo3(get)]
    pub cold_out_t: PyQty,
    /// Cold outlet molar enthalpy.
    #[pyo3(get)]
    pub cold_out_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHeatExchangerResult {
    fn __repr__(&self) -> String {
        format!(
            "HeatExchangerResult(hot_out_t={} {}, cold_out_t={})",
            self.hot_out_t.magnitude_si, self.hot_out_t.unit, self.cold_out_t.magnitude_si
        )
    }
}

impl From<&HeatExchangerResult> for PyHeatExchangerResult {
    fn from(r: &HeatExchangerResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            hot_out_n: quantity(r.hot_out_n, "mol/s"),
            hot_out_z: r.hot_out_z.clone(),
            hot_out_p: quantity(r.hot_out_p.value, "Pa"),
            hot_out_t: quantity(r.hot_out_t.value, "K"),
            hot_out_h: quantity(r.hot_out_h.value, "J/mol"),
            cold_out_n: quantity(r.cold_out_n, "mol/s"),
            cold_out_z: r.cold_out_z.clone(),
            cold_out_p: quantity(r.cold_out_p.value, "Pa"),
            cold_out_t: quantity(r.cold_out_t.value, "K"),
            cold_out_h: quantity(r.cold_out_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.mixer`, transported.
///
/// The counterpart of `SplitterResult`: a `many` port on the way in and the record's five
/// fields on the way out, which is why this has no vectors in it.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "MixerResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyMixerResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub product_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub product_z: Vec<f64>,
    /// Outlet pressure.
    #[pyo3(get)]
    pub product_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub product_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub product_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyMixerResult {
    fn __repr__(&self) -> String {
        format!(
            "MixerResult(product_n={} {}, product_t={})",
            self.product_n.magnitude_si, self.product_n.unit, self.product_t.magnitude_si
        )
    }
}

impl From<&MixerResult> for PyMixerResult {
    fn from(r: &MixerResult) -> Self {
        Self {
            product_n: PyQty {
                magnitude_si: r.product_n,
                unit: "mol/s".to_string(),
            },
            product_z: r.product_z.clone(),
            product_p: PyQty {
                magnitude_si: r.product_p.value,
                unit: "Pa".to_string(),
            },
            product_t: PyQty {
                magnitude_si: r.product_t.value,
                unit: "K".to_string(),
            },
            product_h: PyQty {
                magnitude_si: r.product_h.value,
                unit: "J/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.throttling_valve`, transported.
///
/// Result of `process.pipe`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PipeResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPipeResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// The pressure drop the line implies.
    #[pyo3(get)]
    pub pressure_drop: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPipeResult {
    fn __repr__(&self) -> String {
        format!(
            "PipeResult(outlet_p={} {}, pressure_drop={} {})",
            self.outlet_p.magnitude_si,
            self.outlet_p.unit,
            self.pressure_drop.magnitude_si,
            self.pressure_drop.unit
        )
    }
}

impl From<&PipeResult> for PyPipeResult {
    fn from(r: &PipeResult) -> Self {
        Self {
            outlet_n: PyQty {
                magnitude_si: r.outlet_n,
                unit: "mol/s".to_string(),
            },
            outlet_z: r.outlet_z.clone(),
            outlet_p: PyQty {
                magnitude_si: r.outlet_p.value,
                unit: "Pa".to_string(),
            },
            outlet_t: PyQty {
                magnitude_si: r.outlet_t.value,
                unit: "K".to_string(),
            },
            outlet_h: PyQty {
                magnitude_si: r.outlet_h.value,
                unit: "J/mol".to_string(),
            },
            pressure_drop: PyQty {
                magnitude_si: r.pressure_drop.value,
                unit: "Pa".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.compressor`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "CompressorResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyCompressorResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyCompressorResult {
    fn __repr__(&self) -> String {
        format!(
            "CompressorResult(outlet_t={} {}, outlet_h={} {})",
            self.outlet_t.magnitude_si,
            self.outlet_t.unit,
            self.outlet_h.magnitude_si,
            self.outlet_h.unit
        )
    }
}

impl From<&CompressorResult> for PyCompressorResult {
    fn from(r: &CompressorResult) -> Self {
        Self {
            outlet_n: PyQty {
                magnitude_si: r.outlet_n,
                unit: "mol/s".to_string(),
            },
            outlet_z: r.outlet_z.clone(),
            outlet_p: PyQty {
                magnitude_si: r.outlet_p.value,
                unit: "Pa".to_string(),
            },
            outlet_t: PyQty {
                magnitude_si: r.outlet_t.value,
                unit: "K".to_string(),
            },
            outlet_h: PyQty {
                magnitude_si: r.outlet_h.value,
                unit: "J/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.expander`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ExpanderResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyExpanderResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyExpanderResult {
    fn __repr__(&self) -> String {
        format!(
            "ExpanderResult(outlet_t={} {}, outlet_h={} {})",
            self.outlet_t.magnitude_si,
            self.outlet_t.unit,
            self.outlet_h.magnitude_si,
            self.outlet_h.unit
        )
    }
}

impl From<&ExpanderResult> for PyExpanderResult {
    fn from(r: &ExpanderResult) -> Self {
        Self {
            outlet_n: PyQty {
                magnitude_si: r.outlet_n,
                unit: "mol/s".to_string(),
            },
            outlet_z: r.outlet_z.clone(),
            outlet_p: PyQty {
                magnitude_si: r.outlet_p.value,
                unit: "Pa".to_string(),
            },
            outlet_t: PyQty {
                magnitude_si: r.outlet_t.value,
                unit: "K".to_string(),
            },
            outlet_h: PyQty {
                magnitude_si: r.outlet_h.value,
                unit: "J/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.filter`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "FilterResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFilterResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// The drop actually applied, which differs from the request on a clamped row.
    #[pyo3(get)]
    pub applied_drop: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyFilterResult {
    fn __repr__(&self) -> String {
        format!(
            "FilterResult(outlet_p={} {}, applied_drop={} {})",
            self.outlet_p.magnitude_si,
            self.outlet_p.unit,
            self.applied_drop.magnitude_si,
            self.applied_drop.unit
        )
    }
}

impl From<&FilterResult> for PyFilterResult {
    fn from(r: &FilterResult) -> Self {
        Self {
            outlet_n: PyQty {
                magnitude_si: r.outlet_n,
                unit: "mol/s".to_string(),
            },
            outlet_z: r.outlet_z.clone(),
            outlet_p: PyQty {
                magnitude_si: r.outlet_p.value,
                unit: "Pa".to_string(),
            },
            outlet_t: PyQty {
                magnitude_si: r.outlet_t.value,
                unit: "K".to_string(),
            },
            outlet_h: PyQty {
                magnitude_si: r.outlet_h.value,
                unit: "J/mol".to_string(),
            },
            applied_drop: PyQty {
                magnitude_si: r.applied_drop.value,
                unit: "Pa".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.cooler`, transported.
///
/// Field for field `PyHeaterResult`, and its own class for the reason the model's result is
/// its own struct: `result_fields` answers by id, so one transport for two ids is one handle
/// for two things.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "CoolerResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyCoolerResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// The duty moved, negative when heat was removed.
    #[pyo3(get)]
    pub outlet_duty: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyCoolerResult {
    fn __repr__(&self) -> String {
        format!(
            "CoolerResult(outlet_t={} {}, outlet_duty={} {})",
            self.outlet_t.magnitude_si,
            self.outlet_t.unit,
            self.outlet_duty.magnitude_si,
            self.outlet_duty.unit
        )
    }
}

impl From<&CoolerResult> for PyCoolerResult {
    fn from(r: &CoolerResult) -> Self {
        Self {
            outlet_n: PyQty {
                magnitude_si: r.outlet_n,
                unit: "mol/s".to_string(),
            },
            outlet_z: r.outlet_z.clone(),
            outlet_p: PyQty {
                magnitude_si: r.outlet_p.value,
                unit: "Pa".to_string(),
            },
            outlet_t: PyQty {
                magnitude_si: r.outlet_t.value,
                unit: "K".to_string(),
            },
            outlet_h: PyQty {
                magnitude_si: r.outlet_h.value,
                unit: "J/mol".to_string(),
            },
            outlet_duty: PyQty {
                magnitude_si: r.outlet_duty.value,
                unit: "W".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.heater`, transported.
///
/// `PumpResult`'s five record fields plus the duty, which is the one number a machine that
/// moves heat adds: the outlet's own record cannot carry it, because a record describes a
/// *state* and a duty describes a step.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HeaterResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHeaterResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// The duty moved, as `outlet_n * (outlet_h - inlet_h)`.
    #[pyo3(get)]
    pub outlet_duty: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHeaterResult {
    fn __repr__(&self) -> String {
        format!(
            "HeaterResult(outlet_t={} {}, outlet_duty={} {})",
            self.outlet_t.magnitude_si,
            self.outlet_t.unit,
            self.outlet_duty.magnitude_si,
            self.outlet_duty.unit
        )
    }
}

impl From<&HeaterResult> for PyHeaterResult {
    fn from(r: &HeaterResult) -> Self {
        Self {
            outlet_n: PyQty {
                magnitude_si: r.outlet_n,
                unit: "mol/s".to_string(),
            },
            outlet_z: r.outlet_z.clone(),
            outlet_p: PyQty {
                magnitude_si: r.outlet_p.value,
                unit: "Pa".to_string(),
            },
            outlet_t: PyQty {
                magnitude_si: r.outlet_t.value,
                unit: "K".to_string(),
            },
            outlet_h: PyQty {
                magnitude_si: r.outlet_h.value,
                unit: "J/mol".to_string(),
            },
            outlet_duty: PyQty {
                magnitude_si: r.outlet_duty.value,
                unit: "W".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.manifold`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ManifoldResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyManifoldResult {
    /// Molar flow out, mol/s, one entry per outlet.
    #[pyo3(get)]
    pub products_n: Vec<PyQty>,
    /// Outlet compositions, one row per outlet.
    #[pyo3(get)]
    pub products_z: Vec<Vec<f64>>,
    /// Outlet pressures.
    #[pyo3(get)]
    pub products_p: Vec<PyQty>,
    /// Outlet temperatures.
    #[pyo3(get)]
    pub products_t: Vec<PyQty>,
    /// Outlet molar enthalpies.
    #[pyo3(get)]
    pub products_h: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyManifoldResult {
    fn __repr__(&self) -> String {
        format!(
            "ManifoldResult(products_n={} outlet(s))",
            self.products_n.len()
        )
    }
}

impl From<&ManifoldResult> for PyManifoldResult {
    fn from(r: &ManifoldResult) -> Self {
        Self {
            products_n: r
                .products_n
                .iter()
                .map(|n| PyQty {
                    magnitude_si: *n,
                    unit: "mol/s".to_string(),
                })
                .collect(),
            products_z: r.products_z.clone(),
            products_p: r
                .products_p
                .iter()
                .map(|p| PyQty {
                    magnitude_si: p.value,
                    unit: "Pa".to_string(),
                })
                .collect(),
            products_t: r
                .products_t
                .iter()
                .map(|t| PyQty {
                    magnitude_si: t.value,
                    unit: "K".to_string(),
                })
                .collect(),
            products_h: r
                .products_h
                .iter()
                .map(|h| PyQty {
                    magnitude_si: h.value,
                    unit: "J/mol".to_string(),
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}

/// The same shape as `PumpResult`, which is what a two-port unit operation's outlet is.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ThrottlingValveResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyThrottlingValveResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyThrottlingValveResult {
    fn __repr__(&self) -> String {
        format!(
            "ThrottlingValveResult(outlet_n={} {}, outlet_t={})",
            self.outlet_n.magnitude_si, self.outlet_n.unit, self.outlet_t.magnitude_si
        )
    }
}

impl From<&ThrottlingValveResult> for PyThrottlingValveResult {
    fn from(r: &ThrottlingValveResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            outlet_n: quantity(r.outlet_n, "mol/s"),
            outlet_z: r.outlet_z.clone(),
            outlet_p: quantity(r.outlet_p.value, "Pa"),
            outlet_t: quantity(r.outlet_t.value, "K"),
            outlet_h: quantity(r.outlet_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.separator`, transported.
///
/// **Two single-multiplicity ports, so ten fields.** The record is five fields per port
/// and a separator has two outlets, so this is the pair written out under the ports' own
/// names rather than a list - which is what the port rule gives a `one`-multiplicity port
/// even when a unit operation has two of them.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SeparatorResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySeparatorResult {
    /// Vapour outlet molar flow, mol/s.
    #[pyo3(get)]
    pub vapour_n: PyQty,
    /// Vapour outlet composition.
    #[pyo3(get)]
    pub vapour_z: Vec<f64>,
    /// Vapour outlet pressure.
    #[pyo3(get)]
    pub vapour_p: PyQty,
    /// Vapour outlet temperature.
    #[pyo3(get)]
    pub vapour_t: PyQty,
    /// Vapour outlet molar enthalpy.
    #[pyo3(get)]
    pub vapour_h: PyQty,
    /// Liquid outlet molar flow, mol/s.
    #[pyo3(get)]
    pub liquid_n: PyQty,
    /// Liquid outlet composition.
    #[pyo3(get)]
    pub liquid_z: Vec<f64>,
    /// Liquid outlet pressure.
    #[pyo3(get)]
    pub liquid_p: PyQty,
    /// Liquid outlet temperature.
    #[pyo3(get)]
    pub liquid_t: PyQty,
    /// Liquid outlet molar enthalpy.
    #[pyo3(get)]
    pub liquid_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySeparatorResult {
    fn __repr__(&self) -> String {
        format!(
            "SeparatorResult(vapour_n={} {}, liquid_n={} {})",
            self.vapour_n.magnitude_si,
            self.vapour_n.unit,
            self.liquid_n.magnitude_si,
            self.liquid_n.unit
        )
    }
}

impl From<&SeparatorResult> for PySeparatorResult {
    fn from(r: &SeparatorResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            vapour_n: quantity(r.vapour_n, "mol/s"),
            vapour_z: r.vapour_z.clone(),
            vapour_p: quantity(r.vapour_p.value, "Pa"),
            vapour_t: quantity(r.vapour_t.value, "K"),
            vapour_h: quantity(r.vapour_h.value, "J/mol"),
            liquid_n: quantity(r.liquid_n, "mol/s"),
            liquid_z: r.liquid_z.clone(),
            liquid_p: quantity(r.liquid_p.value, "Pa"),
            liquid_t: quantity(r.liquid_t.value, "K"),
            liquid_h: quantity(r.liquid_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}
/// Result of `process.component_splitter`, transported.
///
/// Two named outlets rather than a vector, because the class fixes the count at two.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ComponentSplitterResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyComponentSplitterResult {
    /// Overhead molar flow, mol/s.
    #[pyo3(get)]
    pub overhead_n: PyQty,
    /// Overhead composition.
    #[pyo3(get)]
    pub overhead_z: Vec<f64>,
    /// Overhead pressure.
    #[pyo3(get)]
    pub overhead_p: PyQty,
    /// Overhead temperature.
    #[pyo3(get)]
    pub overhead_t: PyQty,
    /// Overhead molar enthalpy.
    #[pyo3(get)]
    pub overhead_h: PyQty,
    /// Bottoms molar flow, mol/s.
    #[pyo3(get)]
    pub bottoms_n: PyQty,
    /// Bottoms composition.
    #[pyo3(get)]
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    #[pyo3(get)]
    pub bottoms_p: PyQty,
    /// Bottoms temperature.
    #[pyo3(get)]
    pub bottoms_t: PyQty,
    /// Bottoms molar enthalpy.
    #[pyo3(get)]
    pub bottoms_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyComponentSplitterResult {
    fn __repr__(&self) -> String {
        format!(
            "ComponentSplitterResult(overhead_n={} {}, bottoms_n={} {})",
            self.overhead_n.magnitude_si,
            self.overhead_n.unit,
            self.bottoms_n.magnitude_si,
            self.bottoms_n.unit
        )
    }
}

impl From<&ComponentSplitterResult> for PyComponentSplitterResult {
    fn from(r: &ComponentSplitterResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            overhead_n: quantity(r.overhead_n, "mol/s"),
            overhead_z: r.overhead_z.clone(),
            overhead_p: quantity(r.overhead_p.value, "Pa"),
            overhead_t: quantity(r.overhead_t.value, "K"),
            overhead_h: quantity(r.overhead_h.value, "J/mol"),
            bottoms_n: quantity(r.bottoms_n, "mol/s"),
            bottoms_z: r.bottoms_z.clone(),
            bottoms_p: quantity(r.bottoms_p.value, "Pa"),
            bottoms_t: quantity(r.bottoms_t.value, "K"),
            bottoms_h: quantity(r.bottoms_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.gas_scrubber`, transported.
///
/// **Two single-multiplicity ports, so ten fields.** The record is five fields per port
/// and a separator has two outlets, so this is the pair written out under the ports' own
/// names rather than a list - which is what the port rule gives a `one`-multiplicity port
/// even when a unit operation has two of them.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GasScrubberResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGasScrubberResult {
    /// Vapour outlet molar flow, mol/s.
    #[pyo3(get)]
    pub vapour_n: PyQty,
    /// Vapour outlet composition.
    #[pyo3(get)]
    pub vapour_z: Vec<f64>,
    /// Vapour outlet pressure.
    #[pyo3(get)]
    pub vapour_p: PyQty,
    /// Vapour outlet temperature.
    #[pyo3(get)]
    pub vapour_t: PyQty,
    /// Vapour outlet molar enthalpy.
    #[pyo3(get)]
    pub vapour_h: PyQty,
    /// Liquid outlet molar flow, mol/s.
    #[pyo3(get)]
    pub liquid_n: PyQty,
    /// Liquid outlet composition.
    #[pyo3(get)]
    pub liquid_z: Vec<f64>,
    /// Liquid outlet pressure.
    #[pyo3(get)]
    pub liquid_p: PyQty,
    /// Liquid outlet temperature.
    #[pyo3(get)]
    pub liquid_t: PyQty,
    /// Liquid outlet molar enthalpy.
    #[pyo3(get)]
    pub liquid_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGasScrubberResult {
    fn __repr__(&self) -> String {
        format!(
            "GasScrubberResult(vapour_n={} {}, liquid_n={} {})",
            self.vapour_n.magnitude_si,
            self.vapour_n.unit,
            self.liquid_n.magnitude_si,
            self.liquid_n.unit
        )
    }
}

impl From<&GasScrubberResult> for PyGasScrubberResult {
    fn from(r: &GasScrubberResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            vapour_n: quantity(r.vapour_n, "mol/s"),
            vapour_z: r.vapour_z.clone(),
            vapour_p: quantity(r.vapour_p.value, "Pa"),
            vapour_t: quantity(r.vapour_t.value, "K"),
            vapour_h: quantity(r.vapour_h.value, "J/mol"),
            liquid_n: quantity(r.liquid_n, "mol/s"),
            liquid_z: r.liquid_z.clone(),
            liquid_p: quantity(r.liquid_p.value, "Pa"),
            liquid_t: quantity(r.liquid_t.value, "K"),
            liquid_h: quantity(r.liquid_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.tank`, transported.
///
/// Two outlets named for what they carry - `gas` and `liquid` - rather than a separator's
/// `vapour` and `liquid`, so the result says which phase of the split it is reporting.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TankResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTankResult {
    /// Gas outlet molar flow, mol/s.
    #[pyo3(get)]
    pub gas_n: PyQty,
    /// Gas outlet composition.
    #[pyo3(get)]
    pub gas_z: Vec<f64>,
    /// Gas outlet pressure.
    #[pyo3(get)]
    pub gas_p: PyQty,
    /// Gas outlet temperature.
    #[pyo3(get)]
    pub gas_t: PyQty,
    /// Gas outlet molar enthalpy.
    #[pyo3(get)]
    pub gas_h: PyQty,
    /// Liquid outlet molar flow, mol/s.
    #[pyo3(get)]
    pub liquid_n: PyQty,
    /// Liquid outlet composition.
    #[pyo3(get)]
    pub liquid_z: Vec<f64>,
    /// Liquid outlet pressure.
    #[pyo3(get)]
    pub liquid_p: PyQty,
    /// Liquid outlet temperature.
    #[pyo3(get)]
    pub liquid_t: PyQty,
    /// Liquid outlet molar enthalpy.
    #[pyo3(get)]
    pub liquid_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTankResult {
    fn __repr__(&self) -> String {
        format!(
            "TankResult(gas_n={} {}, liquid_n={} {})",
            self.gas_n.magnitude_si,
            self.gas_n.unit,
            self.liquid_n.magnitude_si,
            self.liquid_n.unit
        )
    }
}

impl From<&TankResult> for PyTankResult {
    fn from(r: &TankResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            gas_n: quantity(r.gas_n, "mol/s"),
            gas_z: r.gas_z.clone(),
            gas_p: quantity(r.gas_p.value, "Pa"),
            gas_t: quantity(r.gas_t.value, "K"),
            gas_h: quantity(r.gas_h.value, "J/mol"),
            liquid_n: quantity(r.liquid_n, "mol/s"),
            liquid_z: r.liquid_z.clone(),
            liquid_p: quantity(r.liquid_p.value, "Pa"),
            liquid_t: quantity(r.liquid_t.value, "K"),
            liquid_h: quantity(r.liquid_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `standards.iso6976`, transported.
///
/// The standard's quantities in SI: molar quantities per mole and densities per cubic metre.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "Iso6976Result"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyIso6976Result {
    /// The mixture's molar mass, kg/mol.
    #[pyo3(get)]
    pub molar_mass: PyQty,
    /// The mixture's compression factor.
    #[pyo3(get)]
    pub compression_factor: f64,
    /// Density relative to dry air.
    #[pyo3(get)]
    pub relative_density: f64,
    /// The ideal-gas density, kg/m³.
    #[pyo3(get)]
    pub density_ideal: PyQty,
    /// The real-gas density, kg/m³.
    #[pyo3(get)]
    pub density_real: PyQty,
    /// Superior molar calorific value, J/mol.
    #[pyo3(get)]
    pub superior_calorific_value: PyQty,
    /// Inferior molar calorific value, J/mol.
    #[pyo3(get)]
    pub inferior_calorific_value: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyIso6976Result {
    fn __repr__(&self) -> String {
        format!(
            "Iso6976Result(molar_mass={} {}, inferior_calorific_value={} {})",
            self.molar_mass.magnitude_si,
            self.molar_mass.unit,
            self.inferior_calorific_value.magnitude_si,
            self.inferior_calorific_value.unit
        )
    }
}

impl From<&Iso6976Result> for PyIso6976Result {
    fn from(r: &Iso6976Result) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            molar_mass: quantity(r.molar_mass.value, "kg/mol"),
            compression_factor: r.compression_factor.value,
            relative_density: r.relative_density.value,
            density_ideal: quantity(r.density_ideal.value, "kg/m**3"),
            density_real: quantity(r.density_real.value, "kg/m**3"),
            superior_calorific_value: quantity(r.superior_calorific_value.value, "J/mol"),
            inferior_calorific_value: quantity(r.inferior_calorific_value.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.ejector`, transported.
///
/// One outlet, five fields: a two-inlet machine still discharges through one port.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "EjectorResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyEjectorResult {
    /// Outlet molar flow, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyEjectorResult {
    fn __repr__(&self) -> String {
        format!(
            "EjectorResult(outlet_n={} {}, outlet_t={} {})",
            self.outlet_n.magnitude_si,
            self.outlet_n.unit,
            self.outlet_t.magnitude_si,
            self.outlet_t.unit
        )
    }
}

impl From<&EjectorResult> for PyEjectorResult {
    fn from(r: &EjectorResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            outlet_n: quantity(r.outlet_n, "mol/s"),
            outlet_z: r.outlet_z.clone(),
            outlet_p: quantity(r.outlet_p.value, "Pa"),
            outlet_t: quantity(r.outlet_t.value, "K"),
            outlet_h: quantity(r.outlet_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.three_phase_separator`, transported.
///
/// Three outlets named for what they carry - `vapour`, `light_liquid` (the class's
/// `getOilOutStream`) and `heavy_liquid` (`getWaterOutStream`) - each the record's five
/// fields, so the result says which phase of the split each pair belongs to.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ThreePhaseSeparatorResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyThreePhaseSeparatorResult {
    /// Vapour outlet molar flow, mol/s.
    #[pyo3(get)]
    pub vapour_n: PyQty,
    /// Vapour outlet composition.
    #[pyo3(get)]
    pub vapour_z: Vec<f64>,
    /// Vapour outlet pressure.
    #[pyo3(get)]
    pub vapour_p: PyQty,
    /// Vapour outlet temperature.
    #[pyo3(get)]
    pub vapour_t: PyQty,
    /// Vapour outlet molar enthalpy.
    #[pyo3(get)]
    pub vapour_h: PyQty,
    /// Oil outlet molar flow, mol/s.
    #[pyo3(get)]
    pub light_liquid_n: PyQty,
    /// Oil outlet composition.
    #[pyo3(get)]
    pub light_liquid_z: Vec<f64>,
    /// Oil outlet pressure.
    #[pyo3(get)]
    pub light_liquid_p: PyQty,
    /// Oil outlet temperature.
    #[pyo3(get)]
    pub light_liquid_t: PyQty,
    /// Oil outlet molar enthalpy.
    #[pyo3(get)]
    pub light_liquid_h: PyQty,
    /// Aqueous outlet molar flow, mol/s.
    #[pyo3(get)]
    pub heavy_liquid_n: PyQty,
    /// Aqueous outlet composition.
    #[pyo3(get)]
    pub heavy_liquid_z: Vec<f64>,
    /// Aqueous outlet pressure.
    #[pyo3(get)]
    pub heavy_liquid_p: PyQty,
    /// Aqueous outlet temperature.
    #[pyo3(get)]
    pub heavy_liquid_t: PyQty,
    /// Aqueous outlet molar enthalpy.
    #[pyo3(get)]
    pub heavy_liquid_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyThreePhaseSeparatorResult {
    fn __repr__(&self) -> String {
        format!(
            "ThreePhaseSeparatorResult(vapour_n={} {}, light_liquid_n={} {}, \
             heavy_liquid_n={} {})",
            self.vapour_n.magnitude_si,
            self.vapour_n.unit,
            self.light_liquid_n.magnitude_si,
            self.light_liquid_n.unit,
            self.heavy_liquid_n.magnitude_si,
            self.heavy_liquid_n.unit
        )
    }
}

impl From<&ThreePhaseSeparatorResult> for PyThreePhaseSeparatorResult {
    fn from(r: &ThreePhaseSeparatorResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            vapour_n: quantity(r.vapour_n, "mol/s"),
            vapour_z: r.vapour_z.clone(),
            vapour_p: quantity(r.vapour_p.value, "Pa"),
            vapour_t: quantity(r.vapour_t.value, "K"),
            vapour_h: quantity(r.vapour_h.value, "J/mol"),
            light_liquid_n: quantity(r.light_liquid_n, "mol/s"),
            light_liquid_z: r.light_liquid_z.clone(),
            light_liquid_p: quantity(r.light_liquid_p.value, "Pa"),
            light_liquid_t: quantity(r.light_liquid_t.value, "K"),
            light_liquid_h: quantity(r.light_liquid_h.value, "J/mol"),
            heavy_liquid_n: quantity(r.heavy_liquid_n, "mol/s"),
            heavy_liquid_z: r.heavy_liquid_z.clone(),
            heavy_liquid_p: quantity(r.heavy_liquid_p.value, "Pa"),
            heavy_liquid_t: quantity(r.heavy_liquid_t.value, "K"),
            heavy_liquid_h: quantity(r.heavy_liquid_h.value, "J/mol"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.splitter`, transported.
///
/// **The `many` port's shape in full**, and the first result that carries one: a vector per
/// scalar record field with one entry per outlet, and `products_z` a matrix with one row per
/// outlet. The alternative - a nested list of outlet objects - is not available to a result,
/// which has to be a flat set of named fields on both sides of the boundary.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SplitterResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySplitterResult {
    /// Molar flow of each outlet, mol/s.
    #[pyo3(get)]
    pub products_n: Vec<PyQty>,
    /// Composition of each outlet, one row per outlet.
    #[pyo3(get)]
    pub products_z: Vec<Vec<f64>>,
    /// Pressure of each outlet.
    #[pyo3(get)]
    pub products_p: Vec<PyQty>,
    /// Temperature of each outlet.
    #[pyo3(get)]
    pub products_t: Vec<PyQty>,
    /// Molar enthalpy of each outlet.
    #[pyo3(get)]
    pub products_h: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySplitterResult {
    fn __repr__(&self) -> String {
        let flows: Vec<f64> = self.products_n.iter().map(|n| n.magnitude_si).collect();
        format!("SplitterResult(products_n={flows:?})")
    }
}

impl From<&SplitterResult> for PySplitterResult {
    fn from(r: &SplitterResult) -> Self {
        let quantity = |values: &[f64], unit: &str| -> Vec<PyQty> {
            values
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: unit.to_string(),
                })
                .collect()
        };
        Self {
            products_n: quantity(&r.products_n, "mol/s"),
            products_z: r.products_z.clone(),
            products_p: quantity(
                &r.products_p.iter().map(|p| p.value).collect::<Vec<_>>(),
                "Pa",
            ),
            products_t: quantity(
                &r.products_t.iter().map(|t| t.value).collect::<Vec<_>>(),
                "K",
            ),
            products_h: quantity(
                &r.products_h.iter().map(|h| h.value).collect::<Vec<_>>(),
                "J/mol",
            ),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.pump`, transported.
///
/// **The first result of a unit operation**, and its shape is the palette's own port
/// record spelled out: a `many` port crosses as a matrix and vectors, so a model with
/// several outlets has one set of fields per port rather than a nested object. There is no
/// `Stream` in a result because a result has to be a flat set of named fields on both sides
/// of the language boundary.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PumpResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPumpResult {
    /// Molar flow out, mol/s.
    #[pyo3(get)]
    pub outlet_n: PyQty,
    /// Outlet composition.
    #[pyo3(get)]
    pub outlet_z: Vec<f64>,
    /// Outlet pressure, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub outlet_p: PyQty,
    /// Outlet temperature.
    #[pyo3(get)]
    pub outlet_t: PyQty,
    /// Outlet molar enthalpy.
    #[pyo3(get)]
    pub outlet_h: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPumpResult {
    fn __repr__(&self) -> String {
        format!(
            "PumpResult(outlet_n={} {}, outlet_t={})",
            self.outlet_n.magnitude_si, self.outlet_n.unit, self.outlet_t.magnitude_si
        )
    }
}

impl From<&PumpResult> for PyPumpResult {
    fn from(r: &PumpResult) -> Self {
        Self {
            outlet_n: PyQty {
                magnitude_si: r.outlet_n,
                unit: "mol/s".to_string(),
            },
            outlet_z: r.outlet_z.clone(),
            outlet_p: PyQty {
                magnitude_si: r.outlet_p.value,
                unit: "Pa".to_string(),
            },
            outlet_t: PyQty {
                magnitude_si: r.outlet_t.value,
                unit: "K".to_string(),
            },
            outlet_h: PyQty {
                magnitude_si: r.outlet_h.value,
                unit: "J/mol".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `reactions.equilibrium_constant`, transported.
///
/// Two dimensionless outputs and two quantities, which is what a correlation whose
/// constants are fitted looks like at a boundary: the constant and its logarithm carry
/// no unit, and the derivative and the heat do.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "EquilibriumConstantResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyEquilibriumConstantResult {
    /// `ln K`. Dimensionless.
    #[pyo3(get)]
    pub ln_k: f64,
    /// `K`. Dimensionless.
    #[pyo3(get)]
    pub k: f64,
    /// `d(ln K)/dT`, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub ln_k_derivative: PyQty,
    /// The heat of reaction, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub reaction_heat: PyQty,
    /// The row's own citation.
    #[pyo3(get)]
    pub reference: String,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyEquilibriumConstantResult {
    fn __repr__(&self) -> String {
        format!(
            "EquilibriumConstantResult(ln_k={}, k={}, reference={})",
            self.ln_k, self.k, self.reference
        )
    }
}

impl From<&EquilibriumConstantResult> for PyEquilibriumConstantResult {
    fn from(r: &EquilibriumConstantResult) -> Self {
        Self {
            ln_k: r.ln_k,
            k: r.k,
            ln_k_derivative: PyQty {
                // **A bare `f64` and not a quantity**: uom has no reciprocal-temperature
                // type, so `1/K` converts by identity and the number is already SI.
                magnitude_si: r.ln_k_derivative,
                unit: "1/K".to_string(),
            },
            reaction_heat: PyQty {
                magnitude_si: r.reaction_heat.value,
                unit: "J/mol".to_string(),
            },
            reference: r.reference.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `reactions.chemical_equilibrium`, transported.
///
/// Four fields, one of them a flag: **an unconverged solve is a result and not an
/// exception**, because NeqSim returns its last iterate and its callers read it.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ChemicalEquilibriumResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyChemicalEquilibriumResult {
    /// The moles of each species, at the answer or where the solve gave up, each as an SI
    /// magnitude and a display unit.
    #[pyo3(get)]
    pub moles: Vec<PyQty>,
    /// Passes taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The final error.
    #[pyo3(get)]
    pub error: f64,
    /// Whether the error came in under the tolerance.
    #[pyo3(get)]
    pub converged: bool,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyChemicalEquilibriumResult {
    fn __repr__(&self) -> String {
        format!(
            "ChemicalEquilibriumResult(converged={}, iterations={}, error={})",
            self.converged, self.iterations, self.error
        )
    }
}

impl From<&ChemicalEquilibriumResult> for PyChemicalEquilibriumResult {
    fn from(r: &ChemicalEquilibriumResult) -> Self {
        Self {
            moles: r
                .moles
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "mol".to_string(),
                })
                .collect(),
            iterations: r.iterations,
            error: r.error,
            converged: r.converged,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `reactions.reactive_phase_equilibrium`, transported.
///
/// The first result carrying a **matrix**, and the first carrying a flag that says whether
/// the solve ran at all. `skipped` crosses as a boolean and not as a sentinel: the
/// composition on a skip is the caller's own, so a flattened flag would make an untouched
/// answer indistinguishable from a converged one.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ReactivePhaseEquilibriumResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyReactivePhaseEquilibriumResult {
    /// Whether the phase was one the solve runs in.
    #[pyo3(get)]
    pub skipped: bool,
    /// The element matrix, one row per element and the electroneutrality row last.
    #[pyo3(get)]
    pub a_matrix: Vec<Vec<f64>>,
    /// The element amounts the solve conserves, the charge correction last.
    #[pyo3(get)]
    pub b: Vec<PyQty>,
    /// Each component's standard-state reference potential.
    #[pyo3(get)]
    pub chem_ref: Vec<PyQty>,
    /// The composition the phase is left holding.
    #[pyo3(get)]
    pub moles: Vec<PyQty>,
    /// Passes the solve took, zero on a skip.
    #[pyo3(get)]
    pub iterations: u32,
    /// The solve's final error, zero on a skip.
    #[pyo3(get)]
    pub error: f64,
    /// The solve's own convergence flag, false on a skip.
    #[pyo3(get)]
    pub converged: bool,
    /// NeqSim's refinement loop: how many ran.
    #[pyo3(get)]
    pub refinements: u32,
    /// Whether all three residuals came in under their tolerances.
    #[pyo3(get)]
    pub certified: bool,
    /// `max |ln Q - ln K|` over the reactions the fluid can run.
    #[pyo3(get)]
    pub max_reaction_log_residual: f64,
    /// The phase's net charge, `nan` on a skip.
    #[pyo3(get)]
    pub net_charge_moles: PyQty,
    /// `max |A n - b|` over the element rows.
    #[pyo3(get)]
    pub max_element_residual: PyQty,
    /// Whether the linear program's estimate became the starting composition.
    #[pyo3(get)]
    pub seed_applied: bool,
    /// The composition the solve started from, one `mol` quantity per component.
    #[pyo3(get)]
    pub seed_moles: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyReactivePhaseEquilibriumResult {
    fn __repr__(&self) -> String {
        format!(
            "ReactivePhaseEquilibriumResult(skipped={}, converged={}, iterations={})",
            self.skipped, self.converged, self.iterations
        )
    }
}

/// Result of `reactions.reactive_ph_flash`, transported.
///
/// The temperature a reactive fluid's enthalpy asks for, and the passes it cost. **The two
/// counts are path quantities** and cross as they are.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ReactivePhFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyReactivePhFlashResult {
    /// The temperature the loop stopped at.
    #[pyo3(get)]
    pub temperature: PyQty,
    /// `isConverged`, also true where the bracket closed rather than the residual.
    #[pyo3(get)]
    pub converged: bool,
    /// The temperature steps taken.
    #[pyo3(get)]
    pub outer_iterations: u32,
    /// Every inner flash's passes, summed.
    #[pyo3(get)]
    pub total_inner_iterations: u32,
    /// Caveats, as `(code, field, message)` triples.
    #[pyo3(get)]
    pub warnings: Vec<(String, String, String)>,
}

#[pymethods]
impl PyReactivePhFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "ReactivePhFlashResult(temperature={}, converged={}, outer_iterations={})",
            self.temperature.magnitude_si, self.converged, self.outer_iterations
        )
    }
}

impl From<&ReactivePhFlashResult> for PyReactivePhFlashResult {
    fn from(r: &ReactivePhFlashResult) -> Self {
        Self {
            temperature: PyQty {
                magnitude_si: r.temperature,
                unit: "K".to_string(),
            },
            converged: r.converged,
            outer_iterations: r.outer_iterations,
            total_inner_iterations: r.total_inner_iterations,
            warnings: r
                .warnings
                .iter()
                .map(|w| {
                    (
                        format!("{:?}", w.code),
                        w.field.clone().unwrap_or_default(),
                        w.message.clone(),
                    )
                })
                .collect(),
        }
    }
}

/// Result of `reactions.reactive_hybrid_eos_ge_flash`, transported.
///
/// The first result whose two mole vectors are in **different orders**: `coupled_moles` is one
/// entry per component, in the fluid's, and `aqueous_moles` is one per *reactive* component, in
/// the order `reactive_components` reports - which is the same fluid order with the spectators
/// dropped, so the two are not aligned entry by entry.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ReactiveHybridEosGeFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyReactiveHybridEosGeFlashResult {
    /// Each role's mole fraction of the feed, in `[gas, oil, aqueous]` order.
    #[pyo3(get)]
    pub beta: Vec<f64>,
    /// Each role's composition, one row per role and one column per component.
    #[pyo3(get)]
    pub x: Vec<Vec<f64>>,
    /// The reaction-adjusted overall inventory the last pass solved at.
    #[pyo3(get)]
    pub coupled_moles: Vec<PyQty>,
    /// The brine's species amounts at the answer, in the reactive set's order.
    #[pyo3(get)]
    pub aqueous_moles: Vec<PyQty>,
    /// How many coupled passes the loop made.
    #[pyo3(get)]
    pub passes: u32,
    /// The last pass's composition deviation.
    #[pyo3(get)]
    pub chemical_deviation: f64,
    /// The last fraction solve's own residual.
    #[pyo3(get)]
    pub residual: f64,
    /// The worst material-balance deviation over the coupled inventory.
    #[pyo3(get)]
    pub max_material_balance_residual: f64,
    /// The worst cross-role `ln(x_i phi_i P)` spread.
    #[pyo3(get)]
    pub max_log_fugacity_residual: f64,
    /// The worst element residual of the split against the coupled inventory.
    #[pyo3(get)]
    pub element_residual: f64,
    /// The net charge the phases fail to account for.
    #[pyo3(get)]
    pub charge_residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyReactiveHybridEosGeFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "ReactiveHybridEosGeFlashResult(passes={}, deviation={:e}, residual={:e})",
            self.passes, self.chemical_deviation, self.residual
        )
    }
}

fn moles(values: &[f64]) -> Vec<PyQty> {
    values
        .iter()
        .map(|value| PyQty {
            magnitude_si: *value,
            unit: "mol".to_string(),
        })
        .collect()
}

impl From<&ReactiveHybridEosGeFlashResult> for PyReactiveHybridEosGeFlashResult {
    fn from(r: &ReactiveHybridEosGeFlashResult) -> Self {
        Self {
            beta: r.beta.clone(),
            x: r.x.clone(),
            coupled_moles: moles(&r.coupled_moles),
            aqueous_moles: moles(&r.aqueous_moles),
            passes: r.passes,
            chemical_deviation: r.chemical_deviation,
            residual: r.residual,
            max_material_balance_residual: r.max_material_balance_residual,
            max_log_fugacity_residual: r.max_log_fugacity_residual,
            element_residual: r.element_residual,
            charge_residual: r.charge_residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `reactions.reactive_tp_flash`, transported.
///
/// The first result carrying a **matrix of moles**, one row per phase, and the first whose
/// rows are the driver's own phase order - with **no phase type**, because NeqSim's types
/// are its system's bookkeeping and not a state the model computes.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ReactiveTpFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyReactiveTpFlashResult {
    /// How many phases the driver stopped on.
    #[pyo3(get)]
    pub phase_count: u32,
    /// Each phase's mole numbers, one row per phase and one column per component.
    #[pyo3(get)]
    pub phase_moles: Vec<Vec<PyQty>>,
    /// Each phase's share of the fluid.
    #[pyo3(get)]
    pub phase_fraction: Vec<f64>,
    /// `isConverged`, true on the branches that accept an answer without a converged solve.
    #[pyo3(get)]
    pub converged: bool,
    /// Every solve's passes, summed. Zero on the `NR = 0` fallback.
    #[pyo3(get)]
    pub total_iterations: u32,
    /// The last solve's total moles.
    #[pyo3(get)]
    pub equilibrium_total_moles: PyQty,
    /// `computeGibbsEnergy`, dimensionless.
    #[pyo3(get)]
    pub gibbs_energy: f64,
    /// What the solve stopped on.
    #[pyo3(get)]
    pub residual: f64,
    /// The element residual on its own.
    #[pyo3(get)]
    pub element_residual: f64,
    /// Caveats, as `(code, field, message)` triples.
    #[pyo3(get)]
    pub warnings: Vec<(String, String, String)>,
}

#[pymethods]
impl PyReactiveTpFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "ReactiveTpFlashResult(phase_count={}, converged={}, iterations={})",
            self.phase_count, self.converged, self.total_iterations
        )
    }
}

impl From<&ReactiveTpFlashResult> for PyReactiveTpFlashResult {
    fn from(r: &ReactiveTpFlashResult) -> Self {
        Self {
            phase_count: r.phase_count as u32,
            phase_moles: r
                .phase_moles
                .iter()
                .map(|phase| {
                    phase
                        .iter()
                        .map(|value| PyQty {
                            magnitude_si: *value,
                            unit: "mol".to_string(),
                        })
                        .collect()
                })
                .collect(),
            phase_fraction: r.phase_fraction.clone(),
            converged: r.converged,
            total_iterations: r.total_iterations,
            equilibrium_total_moles: PyQty {
                magnitude_si: r.equilibrium_total_moles,
                unit: "mol".to_string(),
            },
            gibbs_energy: r.gibbs_energy,
            residual: r.residual,
            element_residual: r.element_residual,
            warnings: r
                .warnings
                .iter()
                .map(|w| {
                    (
                        format!("{:?}", w.code),
                        w.field.clone().unwrap_or_default(),
                        w.message.clone(),
                    )
                })
                .collect(),
        }
    }
}

impl From<&ReactivePhaseEquilibriumResult> for PyReactivePhaseEquilibriumResult {
    fn from(r: &ReactivePhaseEquilibriumResult) -> Self {
        Self {
            skipped: r.skipped,
            a_matrix: r.a_matrix.clone(),
            b: r.b
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "mol".to_string(),
                })
                .collect(),
            chem_ref: r
                .chem_ref
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "J/mol".to_string(),
                })
                .collect(),
            moles: r
                .moles
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "mol".to_string(),
                })
                .collect(),
            iterations: r.iterations,
            error: r.error,
            converged: r.converged,
            refinements: r.refinements,
            certified: r.certified,
            max_reaction_log_residual: r.max_reaction_log_residual,
            net_charge_moles: PyQty {
                magnitude_si: r.net_charge_moles,
                unit: "mol".to_string(),
            },
            max_element_residual: PyQty {
                magnitude_si: r.max_element_residual,
                unit: "mol".to_string(),
            },
            seed_applied: r.seed_applied,
            seed_moles: r
                .seed_moles
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "mol".to_string(),
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `reactions.reference_potentials`, transported.
///
/// Two masks travel beside the potentials: they are what says which components the basis
/// solved for and which reactions survived, and a port that chose a different basis would
/// land on plausible numbers without them.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ReferencePotentialsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyReferencePotentialsResult {
    /// The standard-state reference potentials, in the caller's order, each as an SI
    /// magnitude and a display unit.
    #[pyo3(get)]
    pub potentials: Vec<PyQty>,
    /// A mask over the components: 1.0 where the basis solved directly.
    #[pyo3(get)]
    pub independent: Vec<f64>,
    /// A mask over the source's loaded reactions, in table order.
    #[pyo3(get)]
    pub survivors: Vec<f64>,
    /// The rank the reaction basis reached.
    #[pyo3(get)]
    pub rank: usize,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyReferencePotentialsResult {
    fn __repr__(&self) -> String {
        format!(
            "ReferencePotentialsResult(rank={}, {} component(s))",
            self.rank,
            self.potentials.len()
        )
    }
}

impl From<&ReferencePotentialsResult> for PyReferencePotentialsResult {
    fn from(r: &ReferencePotentialsResult) -> Self {
        Self {
            potentials: r
                .potentials
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "J/mol".to_string(),
                })
                .collect(),
            independent: r.independent.clone(),
            survivors: r.survivors.clone(),
            rank: r.rank,
            warnings: transport(&r.warnings),
        }
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

/// Result of `eos.aqueous_viscosity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "AqueousViscosityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyAqueousViscosityResult {
    /// The phase's dynamic viscosity, as an SI magnitude and a display unit.
    #[pyo3(get)]
    pub viscosity: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyAqueousViscosityResult {
    fn __repr__(&self) -> String {
        format!(
            "AqueousViscosityResult(viscosity={} {})",
            self.viscosity.magnitude_si, self.viscosity.unit
        )
    }
}

impl From<&KernelAqueousViscosityResult> for PyAqueousViscosityResult {
    fn from(r: &KernelAqueousViscosityResult) -> Self {
        Self {
            viscosity: PyQty {
                magnitude_si: r.viscosity.value,
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

/// Result of `eos.nitric_sulfuric_acid_vapor_pressure`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "NitricSulfuricAcidVaporPressureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyNitricSulfuricAcidVaporPressureResult {
    /// The pure-component vapour pressure of water.
    #[pyo3(get)]
    pub p_water: PyQty,
    /// The pure-component vapour pressure of nitric acid.
    #[pyo3(get)]
    pub p_nitric_acid: PyQty,
    /// The pure-component vapour pressure of sulfuric acid.
    #[pyo3(get)]
    pub p_sulfuric_acid: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyNitricSulfuricAcidVaporPressureResult {
    fn __repr__(&self) -> String {
        format!(
            "NitricSulfuricAcidVaporPressureResult(p_water={} Pa, p_nitric_acid={} Pa)",
            self.p_water.magnitude_si, self.p_nitric_acid.magnitude_si
        )
    }
}

impl From<&NitricSulfuricAcidVaporPressureResult> for PyNitricSulfuricAcidVaporPressureResult {
    fn from(r: &NitricSulfuricAcidVaporPressureResult) -> Self {
        Self {
            p_water: PyQty {
                magnitude_si: r.p_water.value,
                unit: "Pa".to_string(),
            },
            p_nitric_acid: PyQty {
                magnitude_si: r.p_nitric_acid.value,
                unit: "Pa".to_string(),
            },
            p_sulfuric_acid: PyQty {
                magnitude_si: r.p_sulfuric_acid.value,
                unit: "Pa".to_string(),
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

/// Result of `eos.ge_nrtl_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GeNrtlPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGeNrtlPhaseResult {
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The pure-component saturation pressure of each component, as an SI magnitude
    /// and display unit.
    #[pyo3(get)]
    pub p_sat: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGeNrtlPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "GeNrtlPhaseResult(gamma={:?}, ln_phi={:?})",
            self.gamma, self.ln_phi
        )
    }
}

impl From<&GeNrtlPhaseResult> for PyGeNrtlPhaseResult {
    fn from(r: &GeNrtlPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            p_sat: r
                .p_sat
                .iter()
                .map(|p| PyQty {
                    magnitude_si: p.value,
                    unit: "Pa".to_string(),
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.ge_nrtl_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GeNrtlFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGeNrtlFlashResult {
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
    /// `ln phi_i` in the liquid.
    #[pyo3(get)]
    pub ln_phi_liquid: Vec<f64>,
    /// `ln phi_i` in the vapour.
    #[pyo3(get)]
    pub ln_phi_vapour: Vec<f64>,
    /// The vapour root of the cubic.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// The smallest `T / Tc_i`.
    #[pyo3(get)]
    pub min_t_over_tc: f64,
    /// The converged state, as its spec spelling.
    #[pyo3(get)]
    pub phase: String,
    /// Successive-substitution steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The convergence residual.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGeNrtlFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "GeNrtlFlashResult(phase={}, beta={:?}, x={:?}, y={:?})",
            self.phase, self.beta, self.x, self.y
        )
    }
}

impl From<&GeNrtlFlashResult> for PyGeNrtlFlashResult {
    fn from(r: &GeNrtlFlashResult) -> Self {
        Self {
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            ln_phi_liquid: r.ln_phi_liquid.clone(),
            ln_phi_vapour: r.ln_phi_vapour.clone(),
            z_vapour: r.z_vapour,
            min_t_over_tc: r.min_t_over_tc,
            phase: r.phase.as_str().to_string(),
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.ge_nrtl_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GeUnifacPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGeUnifacPhaseResult {
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The pure-component saturation pressure of each component, as an SI magnitude
    /// and display unit.
    #[pyo3(get)]
    pub p_sat: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGeUnifacPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "GeUnifacPhaseResult(gamma={:?}, ln_phi={:?})",
            self.gamma, self.ln_phi
        )
    }
}

impl From<&GeUnifacPhaseResult> for PyGeUnifacPhaseResult {
    fn from(r: &GeUnifacPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            p_sat: r
                .p_sat
                .iter()
                .map(|p| PyQty {
                    magnitude_si: p.value,
                    unit: "Pa".to_string(),
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.ge_wilson_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GeWilsonPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGeWilsonPhaseResult {
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The pure-component saturation pressure of each component, as an SI magnitude
    /// and display unit.
    #[pyo3(get)]
    pub p_sat: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGeWilsonPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "GeWilsonPhaseResult(gamma={:?}, ln_phi={:?})",
            self.gamma, self.ln_phi
        )
    }
}

impl From<&GeWilsonPhaseResult> for PyGeWilsonPhaseResult {
    fn from(r: &GeWilsonPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            p_sat: r
                .p_sat
                .iter()
                .map(|p| PyQty {
                    magnitude_si: p.value,
                    unit: "Pa".to_string(),
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.desmukh_mather_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "DesmukhMatherPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyDesmukhMatherPhaseResult {
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// Each component's molality, in mol/kg of solvent.
    #[pyo3(get)]
    pub molality: Vec<f64>,
    /// The ionic strength, in mol/kg.
    #[pyo3(get)]
    pub ionic_strength: f64,
    /// The mean molar mass of the `solvent`-reference components, in kg/mol.
    #[pyo3(get)]
    pub solvent_molar_mass: f64,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyDesmukhMatherPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "DesmukhMatherPhaseResult(I={}, ln_phi={:?})",
            self.ionic_strength, self.ln_phi
        )
    }
}

impl From<&DesmukhMatherPhaseResult> for PyDesmukhMatherPhaseResult {
    fn from(r: &DesmukhMatherPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            molality: r.molality.clone(),
            ionic_strength: r.ionic_strength,
            solvent_molar_mass: r.solvent_molar_mass,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.kent_eisenberg_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "KentEisenbergPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyKentEisenbergPhaseResult {
    /// The activity coefficient of each component, identically one.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient, identically zero.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyKentEisenbergPhaseResult {
    fn __repr__(&self) -> String {
        format!("KentEisenbergPhaseResult(ln_phi={:?})", self.ln_phi)
    }
}

impl From<&KentEisenbergPhaseResult> for PyKentEisenbergPhaseResult {
    fn from(r: &KentEisenbergPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pitzer_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PitzerPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPitzerPhaseResult {
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The Henry coefficient each component's arm read.
    #[pyo3(get)]
    pub henry: Vec<PyQty>,
    /// The infinite-dilution activity coefficient each arm divided by.
    #[pyo3(get)]
    pub gamma_inf: Vec<f64>,
    /// Each component's molality, in mol/kg of solvent.
    #[pyo3(get)]
    pub molality: Vec<f64>,
    /// The ionic strength, in mol/kg.
    #[pyo3(get)]
    pub ionic_strength: f64,
    /// The Pitzer osmotic coefficient of the water.
    #[pyo3(get)]
    pub osmotic_coefficient: f64,
    /// The water activity.
    #[pyo3(get)]
    pub water_activity: f64,
    /// Which parameter dataset answered: `phreeqc` or `legacy`.
    #[pyo3(get)]
    pub dataset: String,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPitzerPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "PitzerPhaseResult(dataset={}, I={}, phi={})",
            self.dataset, self.ionic_strength, self.osmotic_coefficient
        )
    }
}

impl From<&PitzerPhaseResult> for PyPitzerPhaseResult {
    fn from(r: &PitzerPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            henry: r
                .henry
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "Pa".to_string(),
                })
                .collect(),
            gamma_inf: r.gamma_inf.clone(),
            molality: r.molality.clone(),
            ionic_strength: r.ionic_strength,
            osmotic_coefficient: r.osmotic_coefficient,
            water_activity: r.water_activity,
            dataset: r.dataset.name().to_string(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.ge_uniquac_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GeUniquacPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGeUniquacPhaseResult {
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The pure-component saturation pressure of each component, as an SI magnitude
    /// and display unit.
    #[pyo3(get)]
    pub p_sat: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGeUniquacPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "GeUniquacPhaseResult(gamma={:?}, ln_phi={:?})",
            self.gamma, self.ln_phi
        )
    }
}

impl From<&GeUniquacPhaseResult> for PyGeUniquacPhaseResult {
    fn from(r: &GeUniquacPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            p_sat: r
                .p_sat
                .iter()
                .map(|p| PyQty {
                    magnitude_si: p.value,
                    unit: "Pa".to_string(),
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}
/// Result of `eos.ge_van_laar_acid_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GeVanLaarAcidPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGeVanLaarAcidPhaseResult {
    /// The activity coefficient of each component.
    #[pyo3(get)]
    pub gamma: Vec<f64>,
    /// The natural logarithm of each activity coefficient.
    #[pyo3(get)]
    pub ln_gamma: Vec<f64>,
    /// The natural logarithm of each fugacity coefficient.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The pure-component saturation pressure of each component, as an SI magnitude
    /// and display unit.
    #[pyo3(get)]
    pub p_sat: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGeVanLaarAcidPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "GeVanLaarAcidPhaseResult(gamma={:?}, ln_phi={:?})",
            self.gamma, self.ln_phi
        )
    }
}

impl From<&GeVanLaarAcidPhaseResult> for PyGeVanLaarAcidPhaseResult {
    fn from(r: &GeVanLaarAcidPhaseResult) -> Self {
        Self {
            gamma: r.gamma.clone(),
            ln_gamma: r.ln_gamma.clone(),
            ln_phi: r.ln_phi.clone(),
            p_sat: r
                .p_sat
                .iter()
                .map(|p| PyQty {
                    magnitude_si: p.value,
                    unit: "Pa".to_string(),
                })
                .collect(),
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

/// Result of `eos.unifac_umrpru_activity_coefficients`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "UnifacUmrpruActivityCoefficientsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyUnifacUmrpruActivityCoefficientsResult {
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
impl PyUnifacUmrpruActivityCoefficientsResult {
    fn __repr__(&self) -> String {
        format!(
            "UnifacUmrpruActivityCoefficientsResult(ln_gamma={:?}, gamma={:?})",
            self.ln_gamma, self.gamma
        )
    }
}

impl From<&UnifacUmrpruActivityCoefficientsResult> for PyUnifacUmrpruActivityCoefficientsResult {
    fn from(r: &UnifacUmrpruActivityCoefficientsResult) -> Self {
        Self {
            ln_gamma: r.ln_gamma.clone(),
            gamma: r.gamma.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.unifac_psrk_activity_coefficients`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "UnifacPsrkActivityCoefficientsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyUnifacPsrkActivityCoefficientsResult {
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
impl PyUnifacPsrkActivityCoefficientsResult {
    fn __repr__(&self) -> String {
        format!(
            "UnifacPsrkActivityCoefficientsResult(ln_gamma={:?}, gamma={:?})",
            self.ln_gamma, self.gamma
        )
    }
}

impl From<&UnifacPsrkActivityCoefficientsResult> for PyUnifacPsrkActivityCoefficientsResult {
    fn from(r: &UnifacPsrkActivityCoefficientsResult) -> Self {
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

/// Result of `eos.rachford_rice`, transported.
///
/// The model rather than the binary closed form, so `beta` is a root the solver
/// found: outside `[0, 1]` it is the negative flash, which crosses as a value and
/// not as an error.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "RachfordRiceResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyRachfordRiceResult {
    /// The vapour fraction. Dimensionless.
    #[pyo3(get)]
    pub beta: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyRachfordRiceResult {
    fn __repr__(&self) -> String {
        format!(
            "RachfordRiceResult(beta={}, {} warning(s))",
            self.beta,
            self.warnings.len()
        )
    }
}

impl From<&RachfordRiceResult> for PyRachfordRiceResult {
    fn from(r: &RachfordRiceResult) -> Self {
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
    name = "HydrateFormationTemperatureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrateFormationTemperatureResult {
    /// The temperature at which the hydrate's water fugacity meets the fluid's.
    #[pyo3(get)]
    pub temperature: PyQty,
    /// The stable structure, as its spec spelling.
    #[pyo3(get)]
    pub structure: String,
    /// Flash evaluations taken, the scan's and the bisection's together.
    #[pyo3(get)]
    pub iterations: u32,
    /// `f_w^hydrate / f_w^fluid - 1` at the reported temperature.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHydrateFormationTemperatureResult {
    fn __repr__(&self) -> String {
        format!(
            "HydrateFormationTemperatureResult(temperature={} {}, structure={}, {} iteration(s))",
            self.temperature.magnitude_si, self.temperature.unit, self.structure, self.iterations
        )
    }
}

impl From<&HydrateFormationTemperatureResult> for PyHydrateFormationTemperatureResult {
    fn from(r: &HydrateFormationTemperatureResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            temperature: qty(r.temperature.value, "K"),
            structure: r.structure.as_str().to_string(),
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.tp_multiflash_wax`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TpMultiflashWaxResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTpMultiflashWaxResult {
    /// The fraction of the feed's moles in the wax phase.
    #[pyo3(get)]
    pub wax_fraction: f64,
    /// How many phases the feed splits into.
    #[pyo3(get)]
    pub phase_count: u32,
    /// The mole fraction of the feed in each phase, the wax last.
    #[pyo3(get)]
    pub beta: Vec<f64>,
    /// The composition of each phase, one vector per phase.
    #[pyo3(get)]
    pub x: Vec<Vec<f64>>,
    /// Fraction-solve steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The norm of the last fraction correction.
    #[pyo3(get)]
    pub residual: f64,
    /// Whether the residual met the tolerance.
    #[pyo3(get)]
    pub converged: bool,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTpMultiflashWaxResult {
    fn __repr__(&self) -> String {
        format!(
            "TpMultiflashWaxResult(phase_count={}, wax_fraction={}, {} iteration(s), converged={})",
            self.phase_count, self.wax_fraction, self.iterations, self.converged
        )
    }
}

impl From<&TpMultiflashWaxResult> for PyTpMultiflashWaxResult {
    fn from(r: &TpMultiflashWaxResult) -> Self {
        Self {
            wax_fraction: r.wax_fraction,
            phase_count: r.phase_count,
            beta: r.beta.clone(),
            x: r.x.clone(),
            iterations: r.iterations,
            residual: r.residual,
            converged: r.converged,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.tp_solid_flash`, transported.
#[pyclass(
    frozen,
    module = "azoth._core",
    name = "TpSolidFlashResult",
    eq,
    skip_from_py_object
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTpSolidFlashResult {
    /// The fraction of the feed's moles in the pure solid.
    #[pyo3(get)]
    pub solid_fraction: f64,
    /// How many phases the feed splits into, the solid counted when it is there.
    #[pyo3(get)]
    pub phase_count: u32,
    /// The mole fraction of the feed in each phase, the solid last.
    #[pyo3(get)]
    pub beta: Vec<f64>,
    /// The composition of each phase, one vector per phase.
    #[pyo3(get)]
    pub x: Vec<Vec<f64>>,
    /// The pure solid's fugacity coefficient.
    #[pyo3(get)]
    pub solid_fugacity_coefficient: f64,
    /// Fraction-solve steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The norm of the last fraction correction.
    #[pyo3(get)]
    pub residual: f64,
    /// Whether the residual met the tolerance.
    #[pyo3(get)]
    pub converged: bool,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTpSolidFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "TpSolidFlashResult(phase_count={}, solid_fraction={}, solid_fugacity_coefficient={}, {} iteration(s), converged={})",
            self.phase_count,
            self.solid_fraction,
            self.solid_fugacity_coefficient,
            self.iterations,
            self.converged
        )
    }
}

impl From<&TpSolidFlashResult> for PyTpSolidFlashResult {
    fn from(r: &TpSolidFlashResult) -> Self {
        Self {
            solid_fraction: r.solid_fraction,
            phase_count: r.phase_count,
            beta: r.beta.clone(),
            x: r.x.clone(),
            solid_fugacity_coefficient: r.solid_fugacity_coefficient,
            iterations: r.iterations,
            residual: r.residual,
            converged: r.converged,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.salt_precipitation`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SaltPrecipitationResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySaltPrecipitationResult {
    /// The solid taken, in moles per mole of feed.
    #[pyo3(get)]
    pub precipitated_moles: f64,
    /// `IAP/Ksp` before anything was taken.
    #[pyo3(get)]
    pub initial_saturation_ratio: f64,
    /// `IAP/Ksp` at the answer.
    #[pyo3(get)]
    pub final_saturation_ratio: f64,
    /// Bisection steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The extent reached against the bracket's upper end.
    #[pyo3(get)]
    pub extent_of_maximum: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySaltPrecipitationResult {
    fn __repr__(&self) -> String {
        format!(
            "SaltPrecipitationResult(precipitated_moles={}, initial_saturation_ratio={}, {} iteration(s))",
            self.precipitated_moles, self.initial_saturation_ratio, self.iterations
        )
    }
}

impl From<&SaltPrecipitationResult> for PySaltPrecipitationResult {
    fn from(r: &SaltPrecipitationResult) -> Self {
        Self {
            precipitated_moles: r.precipitated_moles,
            initial_saturation_ratio: r.initial_saturation_ratio,
            final_saturation_ratio: r.final_saturation_ratio,
            iterations: r.iterations,
            extent_of_maximum: r.extent_of_maximum,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.solid_fugacity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SolidFugacityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySolidFugacityResult {
    /// The solid's fugacity coefficient for one component.
    #[pyo3(get)]
    pub fugacity_coefficient: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySolidFugacityResult {
    fn __repr__(&self) -> String {
        format!(
            "SolidFugacityResult(fugacity_coefficient={})",
            self.fugacity_coefficient
        )
    }
}

impl From<&SolidFugacityResult> for PySolidFugacityResult {
    fn from(r: &SolidFugacityResult) -> Self {
        Self {
            fugacity_coefficient: r.fugacity_coefficient,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.scale_saturation_ratio`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "ScaleSaturationRatioResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyScaleSaturationRatioResult {
    /// `IAP/Ksp`, one at saturation.
    #[pyo3(get)]
    pub saturation_ratio: f64,
    /// The ion activity product.
    #[pyo3(get)]
    pub ion_activity_product: f64,
    /// `Ksp(T, P)` after the override and the pressure term.
    #[pyo3(get)]
    pub solubility_product: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyScaleSaturationRatioResult {
    fn __repr__(&self) -> String {
        format!(
            "ScaleSaturationRatioResult(saturation_ratio={}, solubility_product={})",
            self.saturation_ratio, self.solubility_product
        )
    }
}

impl From<&ScaleSaturationRatioResult> for PyScaleSaturationRatioResult {
    fn from(r: &ScaleSaturationRatioResult) -> Self {
        Self {
            saturation_ratio: r.saturation_ratio,
            ion_activity_product: r.ion_activity_product,
            solubility_product: r.solubility_product,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.wax_solid_fugacity`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "WaxSolidFugacityResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyWaxSolidFugacityResult {
    /// The wax phase's fugacity coefficient for one component.
    #[pyo3(get)]
    pub fugacity_coefficient: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyWaxSolidFugacityResult {
    fn __repr__(&self) -> String {
        format!(
            "WaxSolidFugacityResult(fugacity_coefficient={})",
            self.fugacity_coefficient
        )
    }
}

impl From<&WaxSolidFugacityResult> for PyWaxSolidFugacityResult {
    fn from(r: &WaxSolidFugacityResult) -> Self {
        Self {
            fugacity_coefficient: r.fugacity_coefficient,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.tbp_fraction_properties`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TbpFractionPropertiesResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTbpFractionPropertiesResult {
    /// Critical temperature.
    #[pyo3(get)]
    pub tc: PyQty,
    /// Critical pressure.
    #[pyo3(get)]
    pub pc: PyQty,
    /// Normal boiling point.
    #[pyo3(get)]
    pub boiling_temperature: PyQty,
    /// Pitzer's acentric factor.
    #[pyo3(get)]
    pub acentric_factor: f64,
    /// The `m` of the cubic's alpha function.
    #[pyo3(get)]
    pub attraction_exponent: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTbpFractionPropertiesResult {
    fn __repr__(&self) -> String {
        format!(
            "TbpFractionPropertiesResult(tc={} {}, pc={} {}, acentric_factor={}, m={})",
            self.tc.magnitude_si,
            self.tc.unit,
            self.pc.magnitude_si,
            self.pc.unit,
            self.acentric_factor,
            self.attraction_exponent
        )
    }
}

impl From<&TbpFractionPropertiesResult> for PyTbpFractionPropertiesResult {
    fn from(r: &TbpFractionPropertiesResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            tc: qty(r.tc.value, "K"),
            pc: qty(r.pc.value, "Pa"),
            boiling_temperature: qty(r.boiling_temperature.value, "K"),
            acentric_factor: r.acentric_factor,
            attraction_exponent: r.attraction_exponent,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.hydrate_formation_pressure`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HydrateFormationPressureResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrateFormationPressureResult {
    /// The pressure at which the hydrate's water fugacity meets the fluid's.
    #[pyo3(get)]
    pub pressure: PyQty,
    /// The stable structure, as its spec spelling.
    #[pyo3(get)]
    pub structure: String,
    /// Flash evaluations taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `f_w^hydrate / f_w^fluid - 1` at the reported pressure.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHydrateFormationPressureResult {
    fn __repr__(&self) -> String {
        format!(
            "HydrateFormationPressureResult(pressure={} {}, structure={}, {} iteration(s))",
            self.pressure.magnitude_si, self.pressure.unit, self.structure, self.iterations
        )
    }
}

impl From<&HydrateFormationPressureResult> for PyHydrateFormationPressureResult {
    fn from(r: &HydrateFormationPressureResult) -> Self {
        Self {
            pressure: PyQty {
                magnitude_si: r.pressure.value,
                unit: "Pa".to_string(),
            },
            structure: r.structure.as_str().to_string(),
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.hydrate_inhibitor_wt`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HydrateInhibitorWtResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrateInhibitorWtResult {
    /// The inhibitor's moles at the answer, the feed's own included.
    #[pyo3(get)]
    pub inhibitor_moles: f64,
    /// The aqueous phase's inhibitor mass fraction.
    #[pyo3(get)]
    pub weight_fraction: f64,
    /// How many phases the answer's state has.
    #[pyo3(get)]
    pub phases: u32,
    /// Secant steps taken, with a floor of three.
    #[pyo3(get)]
    pub iterations: u32,
    /// `-(wtp - wt_target)` at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHydrateInhibitorWtResult {
    fn __repr__(&self) -> String {
        format!(
            "HydrateInhibitorWtResult({} mol, {:.1} wt% aqueous, {} phase(s))",
            self.inhibitor_moles,
            self.weight_fraction * 100.0,
            self.phases
        )
    }
}

impl From<&HydrateInhibitorWtResult> for PyHydrateInhibitorWtResult {
    fn from(r: &HydrateInhibitorWtResult) -> Self {
        Self {
            inhibitor_moles: r.inhibitor_moles,
            weight_fraction: r.weight_fraction,
            phases: r.phases,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.hydrate_inhibitor_concentration`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HydrateInhibitorConcentrationResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrateInhibitorConcentrationResult {
    /// The inhibitor's moles at the answer, the feed's own included.
    #[pyo3(get)]
    pub inhibitor_moles: f64,
    /// The inhibitor's mass fraction of the inhibitor-and-water pair.
    #[pyo3(get)]
    pub weight_fraction: f64,
    /// The hydrate temperature the answer's composition gives.
    #[pyo3(get)]
    pub hydrate_temperature: PyQty,
    /// Secant steps taken, with a floor of three.
    #[pyo3(get)]
    pub iterations: u32,
    /// `T_hydrate - T_target` at the answer, in kelvin.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHydrateInhibitorConcentrationResult {
    fn __repr__(&self) -> String {
        format!(
            "HydrateInhibitorConcentrationResult({} mol, {:.1} wt%, {} iteration(s))",
            self.inhibitor_moles,
            self.weight_fraction * 100.0,
            self.iterations
        )
    }
}

impl From<&HydrateInhibitorConcentrationResult> for PyHydrateInhibitorConcentrationResult {
    fn from(r: &HydrateInhibitorConcentrationResult) -> Self {
        Self {
            inhibitor_moles: r.inhibitor_moles,
            weight_fraction: r.weight_fraction,
            hydrate_temperature: PyQty {
                magnitude_si: r.hydrate_temperature.value,
                unit: "K".to_string(),
            },
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.hydrate_equilibrium_line`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HydrateEquilibriumLineResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrateEquilibriumLineResult {
    /// The formation temperatures along the grid, in kelvin.
    #[pyo3(get)]
    pub temperature: Vec<f64>,
    /// The pressures the points were solved at, in pascals, parallel to `temperature`.
    #[pyo3(get)]
    pub pressure: Vec<f64>,
    /// Caveats, from every point of the grid.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHydrateEquilibriumLineResult {
    fn __repr__(&self) -> String {
        format!(
            "HydrateEquilibriumLineResult({} point(s), {} to {} K)",
            self.temperature.len(),
            self.temperature.first().copied().unwrap_or(f64::NAN),
            self.temperature.last().copied().unwrap_or(f64::NAN)
        )
    }
}

impl From<&HydrateEquilibriumLineResult> for PyHydrateEquilibriumLineResult {
    fn from(r: &HydrateEquilibriumLineResult) -> Self {
        Self {
            temperature: r.temperature.clone(),
            pressure: r.pressure.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.hydrate_fraction`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HydrateFractionResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHydrateFractionResult {
    /// The fraction of the feed's moles that is hydrate.
    #[pyo3(get)]
    pub beta: f64,
    /// The stable structure, as its spec spelling.
    #[pyo3(get)]
    pub structure: String,
    /// The largest `|sum_p beta_p x_ip - z_i|` over the components at the answer.
    #[pyo3(get)]
    pub balance_error: f64,
    /// Flash evaluations taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `ln(f_w^hydrate/f_w^fluid)` at the feed.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHydrateFractionResult {
    fn __repr__(&self) -> String {
        format!(
            "HydrateFractionResult(beta={}, structure={}, balance_error={:.3e}, {} iteration(s))",
            self.beta, self.structure, self.balance_error, self.iterations
        )
    }
}

impl From<&HydrateFractionResult> for PyHydrateFractionResult {
    fn from(r: &HydrateFractionResult) -> Self {
        Self {
            beta: r.beta,
            structure: r.structure.as_str().to_string(),
            balance_error: r.balance_error,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "FreezingPointResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFreezingPointResult {
    /// The temperature at which the calibrated solid's Gibbs energy meets the fluid's.
    #[pyo3(get)]
    pub temperature: PyQty,
    /// Which substance's freezing point this is, of the fluid's candidates.
    #[pyo3(get)]
    pub component: String,
    /// Bracket expansions and bisection steps together.
    #[pyo3(get)]
    pub iterations: u32,
    /// The dimensionless Gibbs difference at the reported temperature.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyFreezingPointResult {
    fn __repr__(&self) -> String {
        format!(
            "FreezingPointResult(temperature={} {}, {} iteration(s))",
            self.temperature.magnitude_si, self.temperature.unit, self.iterations
        )
    }
}

impl From<&FreezingPointResult> for PyFreezingPointResult {
    fn from(r: &FreezingPointResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            component: r.component.clone(),
            temperature: qty(r.temperature.value, "K"),
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

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

/// Result of `eos.hybrid_eos_ge_flash`, transported.
///
/// The three roles are fixed before the fractions are solved, so the vector order is the
/// model's own - `[gas, oil, aqueous]` - and no role field crosses beside it.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "HybridEosGeFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyHybridEosGeFlashResult {
    /// The mole fraction of the feed in each role.
    #[pyo3(get)]
    pub beta: Vec<f64>,
    /// The composition of each role, one row per role.
    #[pyo3(get)]
    pub x: Vec<Vec<f64>>,
    /// `ln phi_i` in each role.
    #[pyo3(get)]
    pub ln_phi: Vec<Vec<f64>>,
    /// Newton steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The solver's own convergence measure.
    #[pyo3(get)]
    pub residual: f64,
    /// The worst material-balance residual.
    #[pyo3(get)]
    pub max_material_balance_residual: f64,
    /// The worst cross-role log-fugacity spread.
    #[pyo3(get)]
    pub max_log_fugacity_residual: f64,
    /// The smallest `T / Tc_i`.
    #[pyo3(get)]
    pub min_t_over_tc: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyHybridEosGeFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "HybridEosGeFlashResult(beta={:?}, iterations={})",
            self.beta, self.iterations
        )
    }
}

impl From<&HybridEosGeFlashResult> for PyHybridEosGeFlashResult {
    fn from(r: &HybridEosGeFlashResult) -> Self {
        Self {
            beta: r.beta.clone(),
            x: r.x.clone(),
            ln_phi: r.ln_phi.clone(),
            iterations: r.iterations,
            residual: r.residual,
            max_material_balance_residual: r.max_material_balance_residual,
            max_log_fugacity_residual: r.max_log_fugacity_residual,
            min_t_over_tc: r.min_t_over_tc,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.iapws_henry_law`, transported.
///
/// The status crosses as the spec's spelling, so the adapter rebuilds the enum and a caller
/// compares `HenryStatus.WITHIN_FITTED_RANGE` regardless of which backend answered.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "IapwsHenryLawResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyIapwsHenryLawResult {
    /// `kH`, the Henry constant.
    #[pyo3(get)]
    pub henry: PyQty,
    /// `ln kH`.
    #[pyo3(get)]
    pub ln_henry: f64,
    /// `d(ln kH)/dT`.
    #[pyo3(get)]
    pub d_ln_henry_d_t: f64,
    /// The spec's spelling of the fitted-range status.
    #[pyo3(get)]
    pub status: String,
    /// The row's reported root-mean-square residual in `ln kH`.
    #[pyo3(get)]
    pub rms_log_residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyIapwsHenryLawResult {
    fn __repr__(&self) -> String {
        format!(
            "IapwsHenryLawResult(henry={}, status={})",
            self.henry.magnitude_si, self.status
        )
    }
}

impl From<&IapwsHenryLawResult> for PyIapwsHenryLawResult {
    fn from(r: &IapwsHenryLawResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            henry: qty(r.henry.value, "Pa"),
            ln_henry: r.ln_henry,
            d_ln_henry_d_t: r.d_ln_henry_d_t,
            status: r.status.as_str().to_string(),
            rms_log_residual: r.rms_log_residual,
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

/// Result of `eos.capillary_dew_point`, transported.
///
/// Its own type rather than a reuse of `PyPhaseBoundaryTemperatureResult`, because it carries
/// one field that one does not - the Young-Laplace pressure - and a shared struct with an
/// `Option` in it would make every caller of the flat model handle a curvature it cannot have.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "CapillaryDewPointResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyCapillaryDewPointResult {
    /// The dew-point temperature.
    #[pyo3(get)]
    pub temperature: PyQty,
    /// The composition of the liquid that first appears.
    #[pyo3(get)]
    pub incipient: Vec<f64>,
    /// K-values at the converged temperature, with the Kelvin shift applied.
    #[pyo3(get)]
    pub k: Vec<f64>,
    /// The liquid root of the cubic at the converged state.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// The Young-Laplace pressure across the interface.
    #[pyo3(get)]
    pub capillary_pressure: PyQty,
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
impl PyCapillaryDewPointResult {
    fn __repr__(&self) -> String {
        format!(
            "CapillaryDewPointResult(T={} K, dP_cap={} Pa, {} iteration(s))",
            self.temperature.magnitude_si, self.capillary_pressure.magnitude_si, self.iterations
        )
    }
}

impl From<&CapillaryDewPointResult> for PyCapillaryDewPointResult {
    fn from(r: &CapillaryDewPointResult) -> Self {
        Self {
            temperature: PyQty {
                magnitude_si: r.temperature.value,
                unit: "K".to_string(),
            },
            incipient: r.incipient.clone(),
            k: r.k.clone(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            capillary_pressure: PyQty {
                magnitude_si: r.capillary_pressure.value,
                unit: "Pa".to_string(),
            },
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

/// Result of `eos.vu_flash_single_comp`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "VuFlashSingleCompResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `T` and `V` are the symbols the spec and the Python result both use
pub struct PyVuFlashSingleCompResult {
    /// The saturation temperature at the pressure asked for.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction, from the lever rule.
    #[pyo3(get)]
    pub beta: f64,
    /// The molar volume the split implies.
    #[pyo3(get)]
    pub V: PyQty,
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyVuFlashSingleCompResult {
    fn __repr__(&self) -> String {
        format!(
            "VuFlashSingleCompResult(T={} K, beta={}, phase={})",
            self.T.magnitude_si, self.beta, self.phase
        )
    }
}

impl From<&VuFlashSingleCompResult> for PyVuFlashSingleCompResult {
    fn from(r: &VuFlashSingleCompResult) -> Self {
        Self {
            T: PyQty {
                magnitude_si: r.t.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            V: PyQty {
                magnitude_si: r.v.value,
                unit: "m**3/mol".to_string(),
            },
            phase: r.phase.as_str().to_string(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pv_reflux_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PvRefluxFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `T` is the symbol the spec and the Python result both use
pub struct PyPvRefluxFlashResult {
    /// The temperature at which the phase ratio is the one asked for.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction at the answer, or `None` for a single-phase feed.
    #[pyo3(get)]
    pub beta: Option<f64>,
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// Liquid-phase mole fractions.
    #[pyo3(get)]
    pub x: Vec<f64>,
    /// Vapour-phase mole fractions.
    #[pyo3(get)]
    pub y: Vec<f64>,
    /// `K_i = y_i / x_i`.
    #[pyo3(get)]
    pub k: Vec<f64>,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Secant steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `|reflux_ratio(T) - reflux_spec|` at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPvRefluxFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "PvRefluxFlashResult(T={} K, beta={:?}, phase={})",
            self.T.magnitude_si, self.beta, self.phase
        )
    }
}

impl From<&PvRefluxFlashResult> for PyPvRefluxFlashResult {
    fn from(r: &PvRefluxFlashResult) -> Self {
        Self {
            T: PyQty {
                magnitude_si: r.t.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            phase: r.phase.as_str().to_string(),
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pvf_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PvfFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `T` is the symbol the spec and the Python result both use
pub struct PyPvfFlashResult {
    /// The temperature at which the feed's vapour fraction is the one asked for.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction at the answer.
    #[pyo3(get)]
    pub beta: f64,
    /// The spec's spelling of the phase.
    #[pyo3(get)]
    pub phase: String,
    /// Liquid-phase mole fractions.
    #[pyo3(get)]
    pub x: Vec<f64>,
    /// Vapour-phase mole fractions.
    #[pyo3(get)]
    pub y: Vec<f64>,
    /// `K_i = y_i / x_i`.
    #[pyo3(get)]
    pub k: Vec<f64>,
    /// The liquid root of the cubic.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The vapour root.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// Illinois steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// `|beta(T) - beta_spec|` at the answer.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPvfFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "PvfFlashResult(T={} K, beta={}, phase={})",
            self.T.magnitude_si, self.beta, self.phase
        )
    }
}

impl From<&PvfFlashResult> for PyPvfFlashResult {
    fn from(r: &PvfFlashResult) -> Self {
        Self {
            T: PyQty {
                magnitude_si: r.t.value,
                unit: "K".to_string(),
            },
            beta: r.beta,
            phase: r.phase.as_str().to_string(),
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.vs_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "VsFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `P` and `T` are the symbols the spec and the Python result both use
pub struct PyVsFlashResult {
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
    /// The larger of the relative volume and entropy residuals.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyVsFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "VsFlashResult(P={} Pa, T={} K, phase={}, beta={:?})",
            self.P.magnitude_si, self.T.magnitude_si, self.phase, self.beta
        )
    }
}

impl From<&VsFlashResult> for PyVsFlashResult {
    fn from(r: &VsFlashResult) -> Self {
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

/// Result of `eos.vh_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "VhFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `P` and `T` are the symbols the spec and the Python result both use
pub struct PyVhFlashResult {
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
    /// The larger of the relative volume and enthalpy residuals.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats, deduplicated.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyVhFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "VhFlashResult(P={} Pa, T={} K, phase={}, beta={:?})",
            self.P.magnitude_si, self.T.magnitude_si, self.phase, self.beta
        )
    }
}

impl From<&VhFlashResult> for PyVhFlashResult {
    fn from(r: &VhFlashResult) -> Self {
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

/// Result of `eos.tv_fraction_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TvFractionFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
#[allow(non_snake_case)] // `P` and `T` are the symbols the spec and the Python result both use
pub struct PyTvFractionFlashResult {
    /// The pressure that satisfies the temperature and volume fraction.
    #[pyo3(get)]
    pub P: PyQty,
    /// The temperature that satisfies the temperature and volume fraction.
    #[pyo3(get)]
    pub T: PyQty,
    /// The vapour fraction, or `None` for a single-phase feed.
    #[pyo3(get)]
    pub beta: Option<f64>,
    /// The gas phase's volume share at the answer.
    #[pyo3(get)]
    pub volume_fraction: f64,
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
impl PyTvFractionFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "TvFractionFlashResult(P={} Pa, T={} K, phase={}, beta={:?})",
            self.P.magnitude_si, self.T.magnitude_si, self.phase, self.beta
        )
    }
}

impl From<&TvFractionFlashResult> for PyTvFractionFlashResult {
    fn from(r: &TvFractionFlashResult) -> Self {
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
            volume_fraction: r.volume_fraction,
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

/// Result of `eos.tp_multiflash`, transported.
///
/// `x` and `ln_phi` cross as rows of tuples, like every other per-phase matrix here, and
/// `seeded` crosses as the spec's spelling and is rebuilt as the enum - so a caller compares
/// `TpMultiflashSeed.STABILITY_SEEDED` whichever backend answered.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TpMultiflashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTpMultiflashResult {
    /// How many phases the feed splits into.
    #[pyo3(get)]
    pub phase_count: u32,
    /// The mole fraction of the feed in each phase.
    #[pyo3(get)]
    pub beta: Vec<f64>,
    /// The composition of each phase, one row per phase.
    #[pyo3(get)]
    pub x: Vec<Vec<f64>>,
    /// The root of the cubic each phase sits on.
    #[pyo3(get)]
    pub z_factor: Vec<f64>,
    /// `ln phi_i` in each phase, one row per phase.
    #[pyo3(get)]
    pub ln_phi: Vec<Vec<f64>>,
    /// The spec's spelling of which path produced the answer.
    #[pyo3(get)]
    pub seeded: String,
    /// The tangent-plane distance at each trial's stationary point.
    #[pyo3(get)]
    pub tm: Vec<f64>,
    /// Fraction-solve steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The larger of the last step's norm and the gradient norm.
    #[pyo3(get)]
    pub residual: f64,
    /// The smallest `T / Tc_i` over the components.
    #[pyo3(get)]
    pub min_t_over_tc: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTpMultiflashResult {
    fn __repr__(&self) -> String {
        format!(
            "TpMultiflashResult({} phase(s), {}, beta={:?})",
            self.phase_count, self.seeded, self.beta
        )
    }
}

impl From<&TpMultiflashResult> for PyTpMultiflashResult {
    fn from(r: &TpMultiflashResult) -> Self {
        Self {
            phase_count: r.phase_count,
            beta: r.beta.clone(),
            x: r.x.clone(),
            z_factor: r.z_factor.clone(),
            ln_phi: r.ln_phi.clone(),
            seeded: r.seeded.as_str().to_string(),
            tm: r.tm.clone(),
            iterations: r.iterations,
            residual: r.residual,
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
        PumpResult::CALC_ID => PumpResult::FIELDS.to_vec(),
        HeaterResult::CALC_ID => HeaterResult::FIELDS.to_vec(),
        CoolerResult::CALC_ID => CoolerResult::FIELDS.to_vec(),
        FilterResult::CALC_ID => FilterResult::FIELDS.to_vec(),
        CompressorResult::CALC_ID => CompressorResult::FIELDS.to_vec(),
        ExpanderResult::CALC_ID => ExpanderResult::FIELDS.to_vec(),
        PipeResult::CALC_ID => PipeResult::FIELDS.to_vec(),
        HeatExchangerResult::CALC_ID => HeatExchangerResult::FIELDS.to_vec(),
        MixerResult::CALC_ID => MixerResult::FIELDS.to_vec(),
        ManifoldResult::CALC_ID => ManifoldResult::FIELDS.to_vec(),
        SeparatorResult::CALC_ID => SeparatorResult::FIELDS.to_vec(),
        GasScrubberResult::CALC_ID => GasScrubberResult::FIELDS.to_vec(),
        ComponentSplitterResult::CALC_ID => ComponentSplitterResult::FIELDS.to_vec(),
        ShortcutDistillationColumnResult::CALC_ID => {
            ShortcutDistillationColumnResult::FIELDS.to_vec()
        }
        DistillationColumnResult::CALC_ID => DistillationColumnResult::FIELDS.to_vec(),
        ThrottlingValveResult::CALC_ID => ThrottlingValveResult::FIELDS.to_vec(),
        SplitterResult::CALC_ID => SplitterResult::FIELDS.to_vec(),
        TankResult::CALC_ID => TankResult::FIELDS.to_vec(),
        ThreePhaseSeparatorResult::CALC_ID => ThreePhaseSeparatorResult::FIELDS.to_vec(),
        EjectorResult::CALC_ID => EjectorResult::FIELDS.to_vec(),
        Iso6976Result::CALC_ID => Iso6976Result::FIELDS.to_vec(),
        EquilibriumConstantResult::CALC_ID => EquilibriumConstantResult::FIELDS.to_vec(),
        ChemicalEquilibriumResult::CALC_ID => ChemicalEquilibriumResult::FIELDS.to_vec(),
        ReactivePhaseEquilibriumResult::CALC_ID => ReactivePhaseEquilibriumResult::FIELDS.to_vec(),
        ReferencePotentialsResult::CALC_ID => ReferencePotentialsResult::FIELDS.to_vec(),
        ReactiveHybridEosGeFlashResult::CALC_ID => ReactiveHybridEosGeFlashResult::FIELDS.to_vec(),
        ReactiveTpFlashResult::CALC_ID => ReactiveTpFlashResult::FIELDS.to_vec(),
        ReactivePhFlashResult::CALC_ID => ReactivePhFlashResult::FIELDS.to_vec(),
        KernelKineticRateLawResult::CALC_ID => KernelKineticRateLawResult::FIELDS.to_vec(),
        KernelKineticsResult::CALC_ID => KernelKineticsResult::FIELDS.to_vec(),
        KernelEffectiveDiffusionResult::CALC_ID => KernelEffectiveDiffusionResult::FIELDS.to_vec(),
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
        RachfordRiceResult::CALC_ID => RachfordRiceResult::FIELDS.to_vec(),
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
        NitricSulfuricAcidVaporPressureResult::CALC_ID => {
            NitricSulfuricAcidVaporPressureResult::FIELDS.to_vec()
        }
        UnifacActivityCoefficientsResult::CALC_ID => {
            UnifacActivityCoefficientsResult::FIELDS.to_vec()
        }
        UniquacActivityCoefficientsResult::CALC_ID => {
            UniquacActivityCoefficientsResult::FIELDS.to_vec()
        }
        UnifacUmrpruActivityCoefficientsResult::CALC_ID => {
            UnifacUmrpruActivityCoefficientsResult::FIELDS.to_vec()
        }
        UnifacPsrkActivityCoefficientsResult::CALC_ID => {
            UnifacPsrkActivityCoefficientsResult::FIELDS.to_vec()
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
        TvFractionFlashResult::CALC_ID => TvFractionFlashResult::FIELDS.to_vec(),
        PuFlashResult::CALC_ID => PuFlashResult::FIELDS.to_vec(),
        PvRefluxFlashResult::CALC_ID => PvRefluxFlashResult::FIELDS.to_vec(),
        PvfFlashResult::CALC_ID => PvfFlashResult::FIELDS.to_vec(),
        VhFlashResult::CALC_ID => VhFlashResult::FIELDS.to_vec(),
        VsFlashResult::CALC_ID => VsFlashResult::FIELDS.to_vec(),
        VuFlashResult::CALC_ID => VuFlashResult::FIELDS.to_vec(),
        VuFlashSingleCompResult::CALC_ID => VuFlashSingleCompResult::FIELDS.to_vec(),
        StabilityTestResult::CALC_ID => StabilityTestResult::FIELDS.to_vec(),
        TpMultiflashResult::CALC_ID => TpMultiflashResult::FIELDS.to_vec(),
        BubblePressureResult::CALC_ID => BubblePressureResult::FIELDS.to_vec(),
        BubbleTemperatureResult::CALC_ID => BubbleTemperatureResult::FIELDS.to_vec(),
        CriticalPointResult::CALC_ID => CriticalPointResult::FIELDS.to_vec(),
        BwrsPhaseResult::CALC_ID => BwrsPhaseResult::FIELDS.to_vec(),
        SrkCpaPhaseResult::CALC_ID => SrkCpaPhaseResult::FIELDS.to_vec(),
        PcsaftRahmatPhaseResult::CALC_ID => PcsaftRahmatPhaseResult::FIELDS.to_vec(),
        SaftVrMiePhaseResult::CALC_ID => SaftVrMiePhaseResult::FIELDS.to_vec(),
        SaftFlashResult::CALC_ID => SaftFlashResult::FIELDS.to_vec(),
        PrCpaPhaseResult::CALC_ID => PrCpaPhaseResult::FIELDS.to_vec(),
        UmrCpaPhaseResult::CALC_ID => UmrCpaPhaseResult::FIELDS.to_vec(),
        SoreideWhitsonPhaseResult::CALC_ID => SoreideWhitsonPhaseResult::FIELDS.to_vec(),
        FurstElectrolytePhaseResult::CALC_ID => FurstElectrolytePhaseResult::FIELDS.to_vec(),
        FurstElectrolyteMod2004PhaseResult::CALC_ID => {
            FurstElectrolyteMod2004PhaseResult::FIELDS.to_vec()
        }
        AmmoniaPhaseResult::CALC_ID => AmmoniaPhaseResult::FIELDS.to_vec(),
        Co2PhaseResult::CALC_ID => Co2PhaseResult::FIELDS.to_vec(),
        HeliumPhaseResult::CALC_ID => HeliumPhaseResult::FIELDS.to_vec(),
        HydrogenPhaseResult::CALC_ID => HydrogenPhaseResult::FIELDS.to_vec(),
        WaterPhaseResult::CALC_ID => WaterPhaseResult::FIELDS.to_vec(),
        ArgonSolidPhaseResult::CALC_ID => ArgonSolidPhaseResult::FIELDS.to_vec(),
        FreezingPointResult::CALC_ID => FreezingPointResult::FIELDS.to_vec(),
        HydrateFormationTemperatureResult::CALC_ID => {
            HydrateFormationTemperatureResult::FIELDS.to_vec()
        }
        HydrateEquilibriumLineResult::CALC_ID => HydrateEquilibriumLineResult::FIELDS.to_vec(),
        HydrateInhibitorConcentrationResult::CALC_ID => {
            HydrateInhibitorConcentrationResult::FIELDS.to_vec()
        }
        HydrateInhibitorWtResult::CALC_ID => HydrateInhibitorWtResult::FIELDS.to_vec(),
        HydrateFractionResult::CALC_ID => HydrateFractionResult::FIELDS.to_vec(),
        HydrateFormationPressureResult::CALC_ID => HydrateFormationPressureResult::FIELDS.to_vec(),
        TbpFractionPropertiesResult::CALC_ID => TbpFractionPropertiesResult::FIELDS.to_vec(),
        ScaleSaturationRatioResult::CALC_ID => ScaleSaturationRatioResult::FIELDS.to_vec(),
        SolidFugacityResult::CALC_ID => SolidFugacityResult::FIELDS.to_vec(),
        SaltPrecipitationResult::CALC_ID => SaltPrecipitationResult::FIELDS.to_vec(),
        WaxSolidFugacityResult::CALC_ID => WaxSolidFugacityResult::FIELDS.to_vec(),
        TpMultiflashWaxResult::CALC_ID => TpMultiflashWaxResult::FIELDS.to_vec(),
        TpSolidFlashResult::CALC_ID => TpSolidFlashResult::FIELDS.to_vec(),
        ParahydrogenSolidPhaseResult::CALC_ID => ParahydrogenSolidPhaseResult::FIELDS.to_vec(),
        EosCgPhaseResult::CALC_ID => EosCgPhaseResult::FIELDS.to_vec(),
        Gerg2008PhaseResult::CALC_ID => Gerg2008PhaseResult::FIELDS.to_vec(),
        GeNrtlPhaseResult::CALC_ID => GeNrtlPhaseResult::FIELDS.to_vec(),
        GeNrtlFlashResult::CALC_ID => GeNrtlFlashResult::FIELDS.to_vec(),
        GeFlashResult::CALC_ID => GeFlashResult::FIELDS.to_vec(),
        GeUnifacPhaseResult::CALC_ID => GeUnifacPhaseResult::FIELDS.to_vec(),
        GeUniquacPhaseResult::CALC_ID => GeUniquacPhaseResult::FIELDS.to_vec(),
        GeVanLaarAcidPhaseResult::CALC_ID => GeVanLaarAcidPhaseResult::FIELDS.to_vec(),
        GeWilsonPhaseResult::CALC_ID => GeWilsonPhaseResult::FIELDS.to_vec(),
        PitzerPhaseResult::CALC_ID => PitzerPhaseResult::FIELDS.to_vec(),
        KentEisenbergPhaseResult::CALC_ID => KentEisenbergPhaseResult::FIELDS.to_vec(),
        DesmukhMatherPhaseResult::CALC_ID => DesmukhMatherPhaseResult::FIELDS.to_vec(),
        DewPressureResult::CALC_ID => DewPressureResult::FIELDS.to_vec(),
        CapillaryDewPointResult::CALC_ID => CapillaryDewPointResult::FIELDS.to_vec(),
        DewTemperatureResult::CALC_ID => DewTemperatureResult::FIELDS.to_vec(),
        HybridEosGeFlashResult::CALC_ID => HybridEosGeFlashResult::FIELDS.to_vec(),
        IapwsHenryLawResult::CALC_ID => IapwsHenryLawResult::FIELDS.to_vec(),
        IdealGasCpResult::CALC_ID => IdealGasCpResult::FIELDS.to_vec(),
        MolarEnthalpyEntropyResult::CALC_ID => MolarEnthalpyEntropyResult::FIELDS.to_vec(),
        ViscosityResult::CALC_ID => ViscosityResult::FIELDS.to_vec(),
        KernelAqueousViscosityResult::CALC_ID => KernelAqueousViscosityResult::FIELDS.to_vec(),
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
        EquilibriumConstantResult::CALC_ID.to_string(),
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
        NitricSulfuricAcidVaporPressureResult::CALC_ID.to_string(),
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
        IapwsHenryLawResult::CALC_ID.to_string(),
        IdealGasCpResult::CALC_ID.to_string(),
        PumpPowerResult::CALC_ID.to_string(),
        KFactorsResult::CALC_ID.to_string(),
        DarcyWeisbachResult::CALC_ID.to_string(),
        TbpFractionPropertiesResult::CALC_ID.to_string(),
        ScaleSaturationRatioResult::CALC_ID.to_string(),
        SolidFugacityResult::CALC_ID.to_string(),
        WaxSolidFugacityResult::CALC_ID.to_string(),
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

/// Result of `eos.srk_cpa_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SrkCpaPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySrkCpaPhaseResult {
    /// The compressibility factor at the chosen root.
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
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySrkCpaPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "SrkCpaPhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&azoth_eos::results::SrkCpaPhaseResult> for PySrkCpaPhaseResult {
    fn from(r: &azoth_eos::results::SrkCpaPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            h_res: qty(r.h_res.value, "J/mol"),
            s_res: qty(r.s_res.value, "J/(mol*K)"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.tp_flash_saft`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "TpFlashSaftResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyTpFlashSaftResult {
    /// The vapour fraction, absent when the flash reports one phase.
    #[pyo3(get)]
    pub beta: Option<f64>,
    /// Liquid-phase mole fractions.
    #[pyo3(get)]
    pub x: Vec<f64>,
    /// Vapour-phase mole fractions.
    #[pyo3(get)]
    pub y: Vec<f64>,
    /// The K-values the loop converged on.
    #[pyo3(get)]
    pub k: Vec<f64>,
    /// `ln phi_i` in the liquid phase.
    #[pyo3(get)]
    pub ln_phi_liquid: Vec<f64>,
    /// `ln phi_i` in the vapour phase.
    #[pyo3(get)]
    pub ln_phi_vapour: Vec<f64>,
    /// The compressibility factor of the liquid solve.
    #[pyo3(get)]
    pub z_liquid: f64,
    /// The compressibility factor of the vapour solve.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// What the converged state is.
    #[pyo3(get)]
    pub phase: String,
    /// Successive-substitution steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The largest relative change in a K-value at the last step.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyTpFlashSaftResult {
    fn __repr__(&self) -> String {
        format!(
            "TpFlashSaftResult(phase={}, beta={:?})",
            self.phase, self.beta
        )
    }
}

impl From<&azoth_eos::results::SaftFlashResult> for PyTpFlashSaftResult {
    fn from(r: &azoth_eos::results::SaftFlashResult) -> Self {
        Self {
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            ln_phi_liquid: r.ln_phi_liquid.clone(),
            ln_phi_vapour: r.ln_phi_vapour.clone(),
            z_liquid: r.z_liquid,
            z_vapour: r.z_vapour,
            phase: r.phase.as_str().to_string(),
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.saft_vr_mie_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SaftVrMiePhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySaftVrMiePhaseResult {
    /// The compressibility factor at the chosen root.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The fugacity coefficients, as logarithms, one per component.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The molar volume at the chosen root.
    #[pyo3(get)]
    pub v: PyQty,
    /// The residual enthalpy.
    #[pyo3(get)]
    pub h_res: PyQty,
    /// The residual entropy.
    #[pyo3(get)]
    pub s_res: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PySaftVrMiePhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "SaftVrMiePhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&azoth_eos::results::SaftVrMiePhaseResult> for PySaftVrMiePhaseResult {
    fn from(r: &azoth_eos::results::SaftVrMiePhaseResult) -> Self {
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            v: PyQty {
                magnitude_si: r.v.value,
                unit: "m**3/mol".to_string(),
            },
            h_res: PyQty {
                magnitude_si: r.h_res.value,
                unit: "J/mol".to_string(),
            },
            s_res: PyQty {
                magnitude_si: r.s_res.value,
                unit: "J/(mol*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pcsaft_rahmat_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PcsaftRahmatPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPcsaftRahmatPhaseResult {
    /// The compressibility factor at the chosen root.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The fugacity coefficients, as logarithms, one per component.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// The molar volume at the chosen root.
    #[pyo3(get)]
    pub v: PyQty,
    /// The residual enthalpy.
    #[pyo3(get)]
    pub h_res: PyQty,
    /// The residual entropy.
    #[pyo3(get)]
    pub s_res: PyQty,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPcsaftRahmatPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "PcsaftRahmatPhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&azoth_eos::results::PcsaftRahmatPhaseResult> for PyPcsaftRahmatPhaseResult {
    fn from(r: &azoth_eos::results::PcsaftRahmatPhaseResult) -> Self {
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            v: PyQty {
                magnitude_si: r.v.value,
                unit: "m**3/mol".to_string(),
            },
            h_res: PyQty {
                magnitude_si: r.h_res.value,
                unit: "J/mol".to_string(),
            },
            s_res: PyQty {
                magnitude_si: r.s_res.value,
                unit: "J/(mol*K)".to_string(),
            },
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.pr_cpa_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "PrCpaPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyPrCpaPhaseResult {
    /// The compressibility factor at the chosen root.
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
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyPrCpaPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "PrCpaPhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&azoth_eos::results::PrCpaPhaseResult> for PyPrCpaPhaseResult {
    fn from(r: &azoth_eos::results::PrCpaPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            h_res: qty(r.h_res.value, "J/mol"),
            s_res: qty(r.s_res.value, "J/(mol*K)"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.umr_cpa_phase`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "UmrCpaPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyUmrCpaPhaseResult {
    /// The compressibility factor at the chosen root.
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
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "SoreideWhitsonPhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PySoreideWhitsonPhaseResult {
    /// The compressibility factor at the chosen root.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The fugacity coefficients, as logarithms, one per component.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "FurstElectrolytePhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFurstElectrolytePhaseResult {
    /// The compressibility factor at the chosen root.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The fugacity coefficients, as logarithms, one per component.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "FurstElectrolyteMod2004PhaseResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyFurstElectrolyteMod2004PhaseResult {
    /// The compressibility factor at the chosen root.
    #[pyo3(get)]
    pub z_factor: f64,
    /// The fugacity coefficients, as logarithms, one per component.
    #[pyo3(get)]
    pub ln_phi: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyFurstElectrolyteMod2004PhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "FurstElectrolyteMod2004PhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&FurstElectrolyteMod2004PhaseResult> for PyFurstElectrolyteMod2004PhaseResult {
    fn from(r: &FurstElectrolyteMod2004PhaseResult) -> Self {
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

#[pymethods]
impl PyFurstElectrolytePhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "FurstElectrolytePhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&FurstElectrolytePhaseResult> for PyFurstElectrolytePhaseResult {
    fn from(r: &FurstElectrolytePhaseResult) -> Self {
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

#[pymethods]
impl PySoreideWhitsonPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "SoreideWhitsonPhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&azoth_eos::results::SoreideWhitsonPhaseResult> for PySoreideWhitsonPhaseResult {
    fn from(r: &azoth_eos::results::SoreideWhitsonPhaseResult) -> Self {
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

#[pymethods]
impl PyUmrCpaPhaseResult {
    fn __repr__(&self) -> String {
        format!(
            "UmrCpaPhaseResult(z_factor={}, ln_phi={:?})",
            self.z_factor, self.ln_phi
        )
    }
}

impl From<&azoth_eos::results::UmrCpaPhaseResult> for PyUmrCpaPhaseResult {
    fn from(r: &azoth_eos::results::UmrCpaPhaseResult) -> Self {
        let qty = |v: f64, unit: &str| PyQty {
            magnitude_si: v,
            unit: unit.to_string(),
        };
        Self {
            z_factor: r.z_factor,
            ln_phi: r.ln_phi.clone(),
            h_res: qty(r.h_res.value, "J/mol"),
            s_res: qty(r.s_res.value, "J/(mol*K)"),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `reactions.kinetic_rate_law`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "KineticRateLawResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyKineticRateLawResult {
    /// The reaction's rate factor at `T`, by the selected law.
    #[pyo3(get)]
    pub rate_factor: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyKineticRateLawResult {
    fn __repr__(&self) -> String {
        format!("KineticRateLawResult(rate_factor={})", self.rate_factor)
    }
}

impl From<&azoth_reactions::KineticRateLawResult> for PyKineticRateLawResult {
    fn from(r: &azoth_reactions::KineticRateLawResult) -> Self {
        Self {
            rate_factor: r.rate_factor,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `reactions.kinetics`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "KineticsResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyKineticsResult {
    /// `reacCoef` per component.
    #[pyo3(get)]
    pub coefficient: Vec<f64>,
    /// `getPhiInfinite` per component, zero where no reaction produced one.
    #[pyo3(get)]
    pub phi_infinite: Vec<f64>,
    /// Whether the irreversibility test fired while each component's row was built.
    #[pyo3(get)]
    pub irreversible: Vec<f64>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyKineticsResult {
    fn __repr__(&self) -> String {
        format!(
            "KineticsResult(coefficient={:?}, irreversible={:?})",
            self.coefficient, self.irreversible
        )
    }
}

impl From<&azoth_reactions::KineticsResult> for PyKineticsResult {
    fn from(r: &azoth_reactions::KineticsResult) -> Self {
        Self {
            coefficient: r.coefficient.clone(),
            phi_infinite: r.phi_infinite.clone(),
            irreversible: r.irreversible.clone(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.effective_diffusion`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "EffectiveDiffusionResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyEffectiveDiffusionResult {
    /// One effective coefficient per component, in m²/s.
    #[pyo3(get)]
    pub effective_diffusion: Vec<PyQty>,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyEffectiveDiffusionResult {
    fn __repr__(&self) -> String {
        format!(
            "EffectiveDiffusionResult(effective_diffusion={:?})",
            self.effective_diffusion
                .iter()
                .map(|value| value.magnitude_si)
                .collect::<Vec<f64>>()
        )
    }
}

impl From<&azoth_eos::EffectiveDiffusionResult> for PyEffectiveDiffusionResult {
    fn from(r: &azoth_eos::EffectiveDiffusionResult) -> Self {
        Self {
            effective_diffusion: r
                .effective_diffusion
                .iter()
                .map(|value| PyQty {
                    magnitude_si: *value,
                    unit: "m**2/s".to_string(),
                })
                .collect(),
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `eos.ge_flash`, transported.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "GeFlashResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyGeFlashResult {
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
    /// `ln phi_i` in the liquid.
    #[pyo3(get)]
    pub ln_phi_liquid: Vec<f64>,
    /// `ln phi_i` in the vapour.
    #[pyo3(get)]
    pub ln_phi_vapour: Vec<f64>,
    /// The vapour root of the cubic.
    #[pyo3(get)]
    pub z_vapour: f64,
    /// The smallest `T / Tc_i`.
    #[pyo3(get)]
    pub min_t_over_tc: f64,
    /// The converged state, as its spec spelling.
    #[pyo3(get)]
    pub phase: String,
    /// Successive-substitution steps taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The convergence residual.
    #[pyo3(get)]
    pub residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyGeFlashResult {
    fn __repr__(&self) -> String {
        format!(
            "GeFlashResult(phase={}, beta={:?}, x={:?}, y={:?})",
            self.phase, self.beta, self.x, self.y
        )
    }
}

impl From<&azoth_eos::GeFlashResult> for PyGeFlashResult {
    fn from(r: &azoth_eos::GeFlashResult) -> Self {
        Self {
            beta: r.beta,
            x: r.x.clone(),
            y: r.y.clone(),
            k: r.k.clone(),
            ln_phi_liquid: r.ln_phi_liquid.clone(),
            ln_phi_vapour: r.ln_phi_vapour.clone(),
            z_vapour: r.z_vapour,
            min_t_over_tc: r.min_t_over_tc,
            phase: r.phase.as_str().to_string(),
            iterations: r.iterations,
            residual: r.residual,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.shortcut_distillation_column`, transported.
///
/// **Ten stream fields and eight scalars.** The two products cross as `process.separator`'s
/// pair does; the eight scalars are the class's whole answer, because a shortcut column's
/// output is not a profile but a stage count, a reflux ratio, a feed tray and two duties.
/// `relative_volatility` is among them because it is the one intermediate the NeqSim class
/// exposes, and it is what makes the capture's K-values checkable against the class's flash.
#[pyclass(module = "azoth._core")]
pub struct PyShortcutDistillationColumnResult {
    /// Distillate molar flow, mol/s.
    #[pyo3(get)]
    pub distillate_n: PyQty,
    /// Distillate composition.
    #[pyo3(get)]
    pub distillate_z: Vec<f64>,
    /// Distillate pressure.
    #[pyo3(get)]
    pub distillate_p: PyQty,
    /// Distillate temperature.
    #[pyo3(get)]
    pub distillate_t: PyQty,
    /// Distillate molar enthalpy.
    #[pyo3(get)]
    pub distillate_h: PyQty,
    /// Bottoms molar flow, mol/s.
    #[pyo3(get)]
    pub bottoms_n: PyQty,
    /// Bottoms composition.
    #[pyo3(get)]
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    #[pyo3(get)]
    pub bottoms_p: PyQty,
    /// Bottoms temperature.
    #[pyo3(get)]
    pub bottoms_t: PyQty,
    /// Bottoms molar enthalpy.
    #[pyo3(get)]
    pub bottoms_h: PyQty,
    /// Fenske's minimum stages.
    #[pyo3(get)]
    pub minimum_stages: f64,
    /// Underwood's minimum reflux ratio.
    #[pyo3(get)]
    pub minimum_reflux_ratio: f64,
    /// Molokanov's actual stage count.
    #[pyo3(get)]
    pub actual_stages: f64,
    /// The reflux ratio the correlation was evaluated at.
    #[pyo3(get)]
    pub actual_reflux_ratio: f64,
    /// The feed stage counted from the top.
    #[pyo3(get)]
    pub feed_tray_number: i64,
    /// The condenser duty the class reports.
    #[pyo3(get)]
    pub condenser_duty: PyQty,
    /// The reboiler duty the class reports.
    #[pyo3(get)]
    pub reboiler_duty: PyQty,
    /// `alpha_LK/HK`.
    #[pyo3(get)]
    pub relative_volatility: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyShortcutDistillationColumnResult {
    fn __repr__(&self) -> String {
        format!(
            "ShortcutDistillationColumnResult(actual_stages={}, feed_tray_number={})",
            self.actual_stages, self.feed_tray_number
        )
    }
}

impl From<&ShortcutDistillationColumnResult> for PyShortcutDistillationColumnResult {
    fn from(r: &ShortcutDistillationColumnResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            distillate_n: quantity(r.distillate_n, "mol/s"),
            distillate_z: r.distillate_z.clone(),
            distillate_p: quantity(r.distillate_p.value, "Pa"),
            distillate_t: quantity(r.distillate_t.value, "K"),
            distillate_h: quantity(r.distillate_h.value, "J/mol"),
            bottoms_n: quantity(r.bottoms_n, "mol/s"),
            bottoms_z: r.bottoms_z.clone(),
            bottoms_p: quantity(r.bottoms_p.value, "Pa"),
            bottoms_t: quantity(r.bottoms_t.value, "K"),
            bottoms_h: quantity(r.bottoms_h.value, "J/mol"),
            minimum_stages: r.minimum_stages,
            minimum_reflux_ratio: r.minimum_reflux_ratio,
            actual_stages: r.actual_stages,
            actual_reflux_ratio: r.actual_reflux_ratio,
            feed_tray_number: r.feed_tray_number,
            condenser_duty: quantity(r.condenser_duty.value, "W"),
            reboiler_duty: quantity(r.reboiler_duty.value, "W"),
            relative_volatility: r.relative_volatility,
            warnings: transport(&r.warnings),
        }
    }
}

/// Result of `process.distillation_column`, transported.
///
/// **Four vectors whose length is the tray count**, which the inputs decide: the profile is
/// what a column's answer is, and a fixed-length shape could not carry it. The two pressure
/// and temperature vectors cross as quantities, one `PyQty` per tray.
#[pyclass(module = "azoth._core")]
pub struct PyDistillationColumnResult {
    /// Each tray's temperature.
    #[pyo3(get)]
    pub tray_temperature: Vec<PyQty>,
    /// Each tray's pressure.
    #[pyo3(get)]
    pub tray_pressure: Vec<PyQty>,
    /// Each tray's vapour traffic.
    #[pyo3(get)]
    pub tray_gas_n: Vec<PyQty>,
    /// Each tray's liquid traffic.
    #[pyo3(get)]
    pub tray_liquid_n: Vec<PyQty>,
    /// Distillate molar flow, mol/s.
    #[pyo3(get)]
    pub distillate_n: PyQty,
    /// Distillate composition.
    #[pyo3(get)]
    pub distillate_z: Vec<f64>,
    /// Distillate pressure.
    #[pyo3(get)]
    pub distillate_p: PyQty,
    /// Distillate temperature.
    #[pyo3(get)]
    pub distillate_t: PyQty,
    /// Distillate molar enthalpy.
    #[pyo3(get)]
    pub distillate_h: PyQty,
    /// Bottoms molar flow, mol/s.
    #[pyo3(get)]
    pub bottoms_n: PyQty,
    /// Bottoms composition.
    #[pyo3(get)]
    pub bottoms_z: Vec<f64>,
    /// Bottoms pressure.
    #[pyo3(get)]
    pub bottoms_p: PyQty,
    /// Bottoms temperature.
    #[pyo3(get)]
    pub bottoms_t: PyQty,
    /// Bottoms molar enthalpy.
    #[pyo3(get)]
    pub bottoms_h: PyQty,
    /// The condenser's duty, W.
    #[pyo3(get)]
    pub condenser_duty: PyQty,
    /// The reboiler's duty, W.
    #[pyo3(get)]
    pub reboiler_duty: PyQty,
    /// Iterations taken.
    #[pyo3(get)]
    pub iterations: u32,
    /// The mean tray-temperature change at the last iteration, K.
    #[pyo3(get)]
    pub temperature_residual: f64,
    /// The mass closure.
    #[pyo3(get)]
    pub mass_residual: f64,
    /// The enthalpy closure.
    #[pyo3(get)]
    pub energy_residual: f64,
    /// Caveats.
    #[pyo3(get)]
    pub warnings: Vec<PyWarning>,
}

#[pymethods]
impl PyDistillationColumnResult {
    fn __repr__(&self) -> String {
        format!(
            "DistillationColumnResult(trays={}, condenserduty={} W)",
            self.tray_temperature.len(),
            self.condenser_duty.magnitude_si
        )
    }
}

impl From<&DistillationColumnResult> for PyDistillationColumnResult {
    fn from(r: &DistillationColumnResult) -> Self {
        let quantity = |magnitude_si: f64, unit: &str| PyQty {
            magnitude_si,
            unit: unit.to_string(),
        };
        Self {
            tray_temperature: r
                .tray_temperature
                .iter()
                .map(|t| quantity(t.value, "K"))
                .collect(),
            tray_pressure: r
                .tray_pressure
                .iter()
                .map(|p| quantity(p.value, "Pa"))
                .collect(),
            tray_gas_n: r.tray_gas_n.iter().map(|n| quantity(*n, "mol/s")).collect(),
            tray_liquid_n: r
                .tray_liquid_n
                .iter()
                .map(|n| quantity(*n, "mol/s"))
                .collect(),
            distillate_n: quantity(r.distillate_n, "mol/s"),
            distillate_z: r.distillate_z.clone(),
            distillate_p: quantity(r.distillate_p.value, "Pa"),
            distillate_t: quantity(r.distillate_t.value, "K"),
            distillate_h: quantity(r.distillate_h.value, "J/mol"),
            bottoms_n: quantity(r.bottoms_n, "mol/s"),
            bottoms_z: r.bottoms_z.clone(),
            bottoms_p: quantity(r.bottoms_p.value, "Pa"),
            bottoms_t: quantity(r.bottoms_t.value, "K"),
            bottoms_h: quantity(r.bottoms_h.value, "J/mol"),
            condenser_duty: quantity(r.condenser_duty.value, "W"),
            reboiler_duty: quantity(r.reboiler_duty.value, "W"),
            iterations: r.iterations,
            temperature_residual: r.temperature_residual,
            mass_residual: r.mass_residual,
            energy_residual: r.energy_residual,
            warnings: transport(&r.warnings),
        }
    }
}
