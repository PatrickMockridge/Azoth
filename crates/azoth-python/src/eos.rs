//! The equations-of-state calculations, exposed to Python.
//!
//! Same rule as `hydraulics.rs` and `thermal.rs`: every parameter below is an SI
//! magnitude, not a `uom` quantity. Here that rule costs nothing, because the
//! quantities in this namespace are dimensionless by construction - the reduced
//! variables and constitutive coefficients an equation of state is written in -
//! so there is no conversion to do at the boundary and no unit string to keep in
//! step with the spec.

use azoth_core::units::{
    cubic_meters_per_mole, joules_per_mole, joules_per_mole_kelvin, kelvins,
    kilograms_per_cubic_meter, kilograms_per_mole, pascal_seconds, pascals,
};
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::{
    PyAmmoniaPhaseResult, PyAntoineVaporPressureResult, PyAqueousViscosityResult,
    PyArgonSolidPhaseResult, PyBwrsPhaseResult, PyCapillaryDewPointResult,
    PyChungConductivityResult, PyChungViscosityResult, PyCo2PhaseResult,
    PyCo2WaterDiffusivityResult, PyCostaldMolarVolumeResult, PyCriticalPointResult,
    PyDesmukhMatherPhaseResult, PyEosCgPhaseResult, PyGeNrtlFlashResult, PyGeNrtlPhaseResult,
    PyGeUnifacPhaseResult, PyGeUniquacPhaseResult, PyGeVanLaarAcidPhaseResult,
    PyGeWilsonPhaseResult, PyGerg2008PhaseResult, PyHaydukMinhasDiffusivityResult,
    PyHeatOfVaporizationResult, PyHeliumPhaseResult, PyHybridEosGeFlashResult,
    PyHydrogenPhaseResult, PyIapwsHenryLawResult, PyIdealGasCpResult, PyKentEisenbergPhaseResult,
    PyLiquidHeatCapacityResult, PyMasonSaxenaConductivityResult, PyMatcop5PrumrAlphaResult,
    PyMatcopAlphaResult, PyMatcopPrAlphaResult, PyMatcopPrumrAlphaResult,
    PyMatcopPrumrNewAlphaResult, PyMolarEnthalpyEntropyResult, PyMollerupAlphaResult,
    PyNitricSulfuricAcidVaporPressureResult, PyNrtlActivityCoefficientsResult,
    PyParachorSurfaceTensionResult, PyParahydrogenSolidPhaseResult, PyPhFlashResult,
    PyPhaseBoundaryResult, PyPhaseBoundaryTemperatureResult, PyPhaseEnvelopeResult,
    PyPitzerPhaseResult, PyPr78KappaResult, PyPrAlphaAbResult, PyPrDaneshAlphaResult,
    PyPrDelft1998AlphaResult, PyPrDepartureResult, PyPrGassem2001AlphaResult, PyPrKappaResult,
    PyPrLeeKeslerAlphaResult, PyPrMassDensityResult, PyPrMolarVolumeResult,
    PyPrPenelouxShiftResult, PyPrZFactorResult, PyPrsvKappaResult, PyPsFlashResult,
    PyPtFlashResult, PyPuFlashResult, PyPureSaturationResult, PyPvFlashResult,
    PyPvRefluxFlashResult, PyPvfFlashResult, PyRachfordRiceBinaryResult, PyRachfordRiceResult,
    PyRackettMolarVolumeResult, PyRkAlphaAbResult, PyRkDepartureResult,
    PySchwartzentruberAlphaResult, PySiddiqiLucasDiffusivityResult, PySoreideWhitsonAlphaResult,
    PySrkAlphaAbResult, PySrkDepartureResult, PySrkKappaResult, PySrkPenelouxShiftResult,
    PySrkZFactorResult, PyStabilityTestResult, PyThFlashResult, PyThermalConductivityResult,
    PyTpMultiflashResult, PyTsFlashResult, PyTuFlashResult, PyTvFlashResult,
    PyTvFractionFlashResult, PyTwuKappaResult, PyTwucoonAlphaResult, PyTwucoonParamAlphaResult,
    PyTwucoonStatoilAlphaResult, PyTynCalusDiffusivityResult, PyUmrprAlphaResult,
    PyUnifacActivityCoefficientsResult, PyUnifacPsrkActivityCoefficientsResult,
    PyUnifacUmrpruActivityCoefficientsResult, PyUniquacActivityCoefficientsResult,
    PyVanLaarAcidActivityCoefficientsResult, PyVdw1fMixBinaryResult, PyVhFlashResult,
    PyViscosityResult, PyVsFlashResult, PyVuFlashResult, PyVuFlashSingleCompResult,
    PyWaterPhaseResult, PyWilkeChangDiffusivityResult, PyWilkeViscosityResult,
    PyWilsonActivityCoefficientsResult,
};

/// The Peng-Robinson alpha-function coefficient.
///
/// `omega` is the acentric factor. A `kappa` below zero is returned carrying an
/// `OUT_OF_VALID_RANGE` warning rather than refused - see the calc's own
/// documentation for why that is the honest answer rather than a refusal.
#[pyfunction]
#[pyo3(signature = (omega))]
#[pyo3(text_signature = "(omega)")]
pub fn pr_kappa(py: Python<'_>, omega: f64) -> PyResult<PyPrKappaResult> {
    azoth_eos::pr_kappa(omega)
        .map(|r| PyPrKappaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Mathias-Copeman alpha function.
#[pyfunction]
#[pyo3(signature = (mc1, mc2, mc3, Tr))]
#[pyo3(text_signature = "(mc1, mc2, mc3, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn matcop_alpha(
    py: Python<'_>,
    mc1: f64,
    mc2: f64,
    mc3: f64,
    Tr: f64,
) -> PyResult<PyMatcopAlphaResult> {
    azoth_eos::matcop_alpha(mc1, mc2, mc3, Tr)
        .map(|r| PyMatcopAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Mathias-Copeman alpha function with a Peng-Robinson supercritical fallback.
#[pyfunction]
#[pyo3(signature = (omega, mc1, mc2, mc3, Tr))]
#[pyo3(text_signature = "(omega, mc1, mc2, mc3, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn matcop_pr_alpha(
    py: Python<'_>,
    omega: f64,
    mc1: f64,
    mc2: f64,
    mc3: f64,
    Tr: f64,
) -> PyResult<PyMatcopPrAlphaResult> {
    azoth_eos::matcop_pr_alpha(omega, mc1, mc2, mc3, Tr)
        .map(|r| PyMatcopPrAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Mathias-Copeman alpha function with the UMR-PR fallback.
#[pyfunction]
#[pyo3(signature = (omega, mc1, mc2, mc3, Tr))]
#[pyo3(text_signature = "(omega, mc1, mc2, mc3, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn matcop_prumr_alpha(
    py: Python<'_>,
    omega: f64,
    mc1: f64,
    mc2: f64,
    mc3: f64,
    Tr: f64,
) -> PyResult<PyMatcopPrumrAlphaResult> {
    azoth_eos::matcop_prumr_alpha(omega, mc1, mc2, mc3, Tr)
        .map(|r| PyMatcopPrumrAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The five-parameter Mathias-Copeman alpha function, UMR-PR new variant.
#[pyfunction]
#[pyo3(signature = (omega, mc1, mc2, mc3, mc4, mc5, Tr))]
#[pyo3(text_signature = "(omega, mc1, mc2, mc3, mc4, mc5, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn matcop_prumr_new_alpha(
    py: Python<'_>,
    omega: f64,
    mc1: f64,
    mc2: f64,
    mc3: f64,
    mc4: f64,
    mc5: f64,
    Tr: f64,
) -> PyResult<PyMatcopPrumrNewAlphaResult> {
    azoth_eos::matcop_prumr_new_alpha(omega, mc1, mc2, mc3, mc4, mc5, Tr)
        .map(|r| PyMatcopPrumrNewAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The five-parameter Mathias-Copeman alpha function.
#[pyfunction]
#[pyo3(signature = (omega, mc1, mc2, mc3, mc4, mc5, Tr))]
#[pyo3(text_signature = "(omega, mc1, mc2, mc3, mc4, mc5, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn matcop5_prumr_alpha(
    py: Python<'_>,
    omega: f64,
    mc1: f64,
    mc2: f64,
    mc3: f64,
    mc4: f64,
    mc5: f64,
    Tr: f64,
) -> PyResult<PyMatcop5PrumrAlphaResult> {
    azoth_eos::matcop5_prumr_alpha(omega, mc1, mc2, mc3, mc4, mc5, Tr)
        .map(|r| PyMatcop5PrumrAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Mollerup alpha function.
#[pyfunction]
#[pyo3(signature = (p1, p2, p3, Tr))]
#[pyo3(text_signature = "(p1, p2, p3, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn mollerup_alpha(
    py: Python<'_>,
    p1: f64,
    p2: f64,
    p3: f64,
    Tr: f64,
) -> PyResult<PyMollerupAlphaResult> {
    azoth_eos::mollerup_alpha(p1, p2, p3, Tr)
        .map(|r| PyMollerupAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Danesh alpha function.
#[pyfunction]
#[pyo3(signature = (omega, Tr))]
#[pyo3(text_signature = "(omega, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_danesh_alpha(py: Python<'_>, omega: f64, Tr: f64) -> PyResult<PyPrDaneshAlphaResult> {
    azoth_eos::pr_danesh_alpha(omega, Tr)
        .map(|r| PyPrDaneshAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Peng-Robinson alpha function, Delft (1998).
#[pyfunction]
#[pyo3(signature = (omega, Tr))]
#[pyo3(text_signature = "(omega, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_delft1998_alpha(
    py: Python<'_>,
    omega: f64,
    Tr: f64,
) -> PyResult<PyPrDelft1998AlphaResult> {
    azoth_eos::pr_delft1998_alpha(omega, Tr)
        .map(|r| PyPrDelft1998AlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Gassem (2001) alpha function.
#[pyfunction]
#[pyo3(signature = (omega, Tr))]
#[pyo3(text_signature = "(omega, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_gassem2001_alpha(
    py: Python<'_>,
    omega: f64,
    Tr: f64,
) -> PyResult<PyPrGassem2001AlphaResult> {
    azoth_eos::pr_gassem2001_alpha(omega, Tr)
        .map(|r| PyPrGassem2001AlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Peng-Robinson alpha function with a Soave-form m-factor.
#[pyfunction]
#[pyo3(signature = (omega, Tr))]
#[pyo3(text_signature = "(omega, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_lee_kesler_alpha(
    py: Python<'_>,
    omega: f64,
    Tr: f64,
) -> PyResult<PyPrLeeKeslerAlphaResult> {
    azoth_eos::pr_lee_kesler_alpha(omega, Tr)
        .map(|r| PyPrLeeKeslerAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Peng-Robinson alpha function and the reduced attraction parameters.
///
/// `kappa` comes from `pr_kappa`; `Tr` and `Pr` are reduced against the caller's
/// critical point. All six quantities are dimensionless, so nothing here touches
/// units in either direction.
#[pyfunction]
#[pyo3(signature = (kappa, Tr, Pr))]
#[pyo3(text_signature = "(kappa, Tr, Pr)")]
#[allow(non_snake_case)] // `Tr` and `Pr` are the symbols in the published equation
pub fn pr_alpha_ab(py: Python<'_>, kappa: f64, Tr: f64, Pr: f64) -> PyResult<PyPrAlphaAbResult> {
    azoth_eos::pr_alpha_ab(kappa, Tr, Pr)
        .map(|r| PyPrAlphaAbResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Peng-Robinson compressibility factor.
///
/// Both arguments are the cubic's dimensionless parameters, so nothing here
/// touches units. `root_structure` crosses as the spec's spelling and the bridge
/// rebuilds the enum - the same arrangement `reynolds_number`'s `regime` uses.
#[pyfunction]
#[pyo3(signature = (a_reduced, b_reduced))]
#[pyo3(text_signature = "(a_reduced, b_reduced)")]
pub fn pr_z_factor(py: Python<'_>, a_reduced: f64, b_reduced: f64) -> PyResult<PyPrZFactorResult> {
    azoth_eos::pr_z_factor(a_reduced, b_reduced)
        .map(|r| PyPrZFactorResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The PRSV alpha-function coefficient.
///
/// All three arguments are dimensionless, and `kappa1` is the caller's - this
/// library ships no fitted values for it, which is the one thing about PRSV a
/// caller is most likely to expect and not get.
#[pyfunction]
#[pyo3(signature = (omega, Tr, kappa1))]
#[pyo3(text_signature = "(omega, Tr, kappa1)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn prsv_kappa(py: Python<'_>, omega: f64, Tr: f64, kappa1: f64) -> PyResult<PyPrsvKappaResult> {
    azoth_eos::prsv_kappa(omega, Tr, kappa1)
        .map(|r| PyPrsvKappaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The 1978 Peng-Robinson alpha-function coefficient.
#[pyfunction]
#[pyo3(signature = (omega))]
#[pyo3(text_signature = "(omega)")]
pub fn pr78_kappa(py: Python<'_>, omega: f64) -> PyResult<PyPr78KappaResult> {
    azoth_eos::pr78_kappa(omega)
        .map(|r| PyPr78KappaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Twu's alpha-function coefficient.
#[pyfunction]
#[pyo3(signature = (omega))]
#[pyo3(text_signature = "(omega)")]
pub fn twu_kappa(py: Python<'_>, omega: f64) -> PyResult<PyTwuKappaResult> {
    azoth_eos::twu_kappa(omega)
        .map(|r| PyTwuKappaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Twu-Coon alpha function.
#[pyfunction]
#[pyo3(signature = (omega, Tr))]
#[pyo3(text_signature = "(omega, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn twucoon_alpha(py: Python<'_>, omega: f64, Tr: f64) -> PyResult<PyTwucoonAlphaResult> {
    azoth_eos::twucoon_alpha(omega, Tr)
        .map(|r| PyTwucoonAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Twu-Coon parameter alpha function.
#[pyfunction]
#[pyo3(signature = (a, b, c, Tr))]
#[pyo3(text_signature = "(a, b, c, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn twucoon_param_alpha(
    py: Python<'_>,
    a: f64,
    b: f64,
    c: f64,
    Tr: f64,
) -> PyResult<PyTwucoonParamAlphaResult> {
    azoth_eos::twucoon_param_alpha(a, b, c, Tr)
        .map(|r| PyTwucoonParamAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Twu-Coon Statoil alpha function.
#[pyfunction]
#[pyo3(signature = (a, b, c, Tr))]
#[pyo3(text_signature = "(a, b, c, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn twucoon_statoil_alpha(
    py: Python<'_>,
    a: f64,
    b: f64,
    c: f64,
    Tr: f64,
) -> PyResult<PyTwucoonStatoilAlphaResult> {
    azoth_eos::twucoon_statoil_alpha(a, b, c, Tr)
        .map(|r| PyTwucoonStatoilAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Peng-Robinson fugacity coefficient and departure functions.
///
/// Five dimensionless arguments and three dimensionless outputs. `kappa` may come
/// from either `pr_kappa` or `prsv_kappa` - this calc uses it only through the
/// logarithmic derivative of the alpha function, which both correlations feed.
#[pyfunction]
#[pyo3(signature = (a_reduced, b_reduced, z, kappa, Tr))]
#[pyo3(text_signature = "(a_reduced, b_reduced, z, kappa, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn pr_departure(
    py: Python<'_>,
    a_reduced: f64,
    b_reduced: f64,
    z: f64,
    kappa: f64,
    Tr: f64,
) -> PyResult<PyPrDepartureResult> {
    azoth_eos::pr_departure(a_reduced, b_reduced, z, kappa, Tr)
        .map(|r| PyPrDepartureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Soave-Redlich-Kwong alpha-function coefficient.
#[pyfunction]
#[pyo3(signature = (omega))]
#[pyo3(text_signature = "(omega)")]
pub fn srk_kappa(py: Python<'_>, omega: f64) -> PyResult<PySrkKappaResult> {
    azoth_eos::srk_kappa(omega)
        .map(|r| PySrkKappaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Soave-Redlich-Kwong alpha function and the reduced attraction parameters.
#[pyfunction]
#[pyo3(signature = (kappa, Tr, Pr))]
#[pyo3(text_signature = "(kappa, Tr, Pr)")]
#[allow(non_snake_case)] // `Tr` and `Pr` are the symbols in the published equation
pub fn srk_alpha_ab(py: Python<'_>, kappa: f64, Tr: f64, Pr: f64) -> PyResult<PySrkAlphaAbResult> {
    azoth_eos::srk_alpha_ab(kappa, Tr, Pr)
        .map(|r| PySrkAlphaAbResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Soave-Redlich-Kwong compressibility factor.
#[pyfunction]
#[pyo3(signature = (a_reduced, b_reduced))]
#[pyo3(text_signature = "(a_reduced, b_reduced)")]
pub fn srk_z_factor(
    py: Python<'_>,
    a_reduced: f64,
    b_reduced: f64,
) -> PyResult<PySrkZFactorResult> {
    azoth_eos::srk_z_factor(a_reduced, b_reduced)
        .map(|r| PySrkZFactorResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Soave-Redlich-Kwong fugacity coefficient and departure functions.
#[pyfunction]
#[pyo3(signature = (a_reduced, b_reduced, z, kappa, Tr))]
#[pyo3(text_signature = "(a_reduced, b_reduced, z, kappa, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn srk_departure(
    py: Python<'_>,
    a_reduced: f64,
    b_reduced: f64,
    z: f64,
    kappa: f64,
    Tr: f64,
) -> PyResult<PySrkDepartureResult> {
    azoth_eos::srk_departure(a_reduced, b_reduced, z, kappa, Tr)
        .map(|r| PySrkDepartureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Redlich-Kwong alpha function and the reduced attraction parameters.
#[pyfunction]
#[pyo3(signature = (Tr, Pr))]
#[pyo3(text_signature = "(Tr, Pr)")]
#[allow(non_snake_case)] // `Tr` and `Pr` are the symbols in the published equation
pub fn rk_alpha_ab(py: Python<'_>, Tr: f64, Pr: f64) -> PyResult<PyRkAlphaAbResult> {
    azoth_eos::rk_alpha_ab(Tr, Pr)
        .map(|r| PyRkAlphaAbResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Redlich-Kwong fugacity coefficient and departure functions.
#[pyfunction]
#[pyo3(signature = (a_reduced, b_reduced, z))]
#[pyo3(text_signature = "(a_reduced, b_reduced, z)")]
pub fn rk_departure(
    py: Python<'_>,
    a_reduced: f64,
    b_reduced: f64,
    z: f64,
) -> PyResult<PyRkDepartureResult> {
    azoth_eos::rk_departure(a_reduced, b_reduced, z)
        .map(|r| PyRkDepartureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The van der Waals one-fluid mixture parameters for a binary.
///
/// Six dimensionless arguments and two dimensionless outputs. `k12` is the
/// caller's - this library ships no fitted binary parameters.
#[pyfunction]
#[pyo3(signature = (z1, a1, a2, b1, b2, k12))]
#[pyo3(text_signature = "(z1, a1, a2, b1, b2, k12)")]
pub fn vdw1f_mix_binary(
    py: Python<'_>,
    z1: f64,
    a1: f64,
    a2: f64,
    b1: f64,
    b2: f64,
    k12: f64,
) -> PyResult<PyVdw1fMixBinaryResult> {
    azoth_eos::vdw1f_mix_binary(z1, a1, a2, b1, b2, k12)
        .map(|r| PyVdw1fMixBinaryResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The binary Rachford-Rice vapour fraction.
///
/// A `beta` outside `[0, 1]` crosses as a value carrying an `OUT_OF_VALID_RANGE`
/// warning rather than as an error - the feed is single phase, which is a real
/// answer and not a failure.
#[pyfunction]
#[pyo3(signature = (z1, K1, K2))]
#[pyo3(text_signature = "(z1, K1, K2)")]
#[allow(non_snake_case)] // `K1` and `K2` are the symbols in the published equation
pub fn rachford_rice_binary(
    py: Python<'_>,
    z1: f64,
    K1: f64,
    K2: f64,
) -> PyResult<PyRachfordRiceBinaryResult> {
    azoth_eos::rachford_rice_binary(z1, K1, K2)
        .map(|r| PyRachfordRiceBinaryResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Rachford-Rice vapour fraction, by NeqSim's own solver.
///
/// A *model* rather than a calculation, so the procedure comes from the generated
/// model table and this function only carries it across the boundary. The root comes
/// back as the equation gives it: outside `[0, 1]` it is the negative flash, which is
/// a reading rather than a failure.
#[pyfunction]
#[pyo3(signature = (z, K))]
#[pyo3(text_signature = "(z, K)")]
#[allow(non_snake_case)] // `K` is the symbol in the published equation
pub fn rachford_rice(py: Python<'_>, z: Vec<f64>, K: Vec<f64>) -> PyResult<PyRachfordRiceResult> {
    azoth_eos::rachford_rice::rachford_rice(&z, &K)
        .map(|r| PyRachfordRiceResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Molar volume from a compressibility factor.
///
/// The one dimensional calc in this namespace, so unlike its neighbours this takes
/// and returns SI magnitudes with units attached - `T` in kelvin and `P` in pascals
/// in, a molar volume out.
#[pyfunction]
#[pyo3(signature = (z, T, P))]
#[pyo3(text_signature = "(z, T, P)")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the published equation
pub fn pr_molar_volume(py: Python<'_>, z: f64, T: f64, P: f64) -> PyResult<PyPrMolarVolumeResult> {
    azoth_eos::pr_molar_volume(z, kelvins(T), pascals(P))
        .map(|r| PyPrMolarVolumeResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Mass density from a molar mass and a molar volume.
///
/// `M` is in kg/mol - not the g/mol a table would quote - and `v` in m**3/mol.
/// Both are SI magnitudes crossing this boundary, converted once in Python.
#[pyfunction]
#[pyo3(signature = (M, v))]
#[pyo3(text_signature = "(M, v)")]
#[allow(non_snake_case)] // `M` is the symbol in the equation
pub fn pr_mass_density(py: Python<'_>, M: f64, v: f64) -> PyResult<PyPrMassDensityResult> {
    azoth_eos::pr_mass_density(kilograms_per_mole(M), cubic_meters_per_mole(v))
        .map(|r| PyPrMassDensityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Peng-Robinson Peneloux volume-translation parameter.
///
/// `Tc` in kelvin and `Pc` in pascals in, a volume shift in `m**3/mol` out - the
/// quantity subtracted from `eos.pr_molar_volume` to correct a reported volume.
#[pyfunction]
#[pyo3(signature = (omega, Tc, Pc, z_ra = None))]
#[pyo3(text_signature = "(omega, Tc, Pc, z_ra=None)")]
#[allow(non_snake_case)] // `Tc` and `Pc` are the symbols in the published equation
pub fn pr_peneloux_shift(
    py: Python<'_>,
    omega: f64,
    Tc: f64,
    Pc: f64,
    z_ra: Option<f64>,
) -> PyResult<PyPrPenelouxShiftResult> {
    azoth_eos::pr_peneloux_shift(omega, kelvins(Tc), pascals(Pc), z_ra)
        .map(|r| PyPrPenelouxShiftResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Soave-Redlich-Kwong Peneloux volume-translation parameter.
///
/// `Tc` in kelvin and `Pc` in pascals in, a volume shift in `m**3/mol` out.
#[pyfunction]
#[pyo3(signature = (omega, Tc, Pc))]
#[pyo3(text_signature = "(omega, Tc, Pc)")]
#[allow(non_snake_case)] // `Tc` and `Pc` are the symbols in the published equation
pub fn srk_peneloux_shift(
    py: Python<'_>,
    omega: f64,
    Tc: f64,
    Pc: f64,
) -> PyResult<PySrkPenelouxShiftResult> {
    azoth_eos::srk_peneloux_shift(omega, kelvins(Tc), pascals(Pc))
        .map(|r| PySrkPenelouxShiftResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pure-component heat of vaporisation, from NeqSim's correlation.
#[pyfunction]
#[pyo3(signature = (c0, c1, c2, c3, Tc, T))]
#[pyo3(text_signature = "(c0, c1, c2, c3, Tc, T)")]
#[allow(non_snake_case)] // `Tc` and `T` are the symbols in the published equation
pub fn heat_of_vaporization(
    py: Python<'_>,
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    Tc: f64,
    T: f64,
) -> PyResult<PyHeatOfVaporizationResult> {
    azoth_eos::heat_of_vaporization(c0, c1, c2, c3, kelvins(Tc), kelvins(T))
        .map(|r| PyHeatOfVaporizationResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pure-component liquid heat capacity, from NeqSim's polynomial.
#[pyfunction]
#[pyo3(signature = (c0, c1, c2, c3, c4, T))]
#[pyo3(text_signature = "(c0, c1, c2, c3, c4, T)")]
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn liquid_heat_capacity(
    py: Python<'_>,
    c0: f64,
    c1: f64,
    c2: f64,
    c3: f64,
    c4: f64,
    T: f64,
) -> PyResult<PyLiquidHeatCapacityResult> {
    azoth_eos::liquid_heat_capacity(c0, c1, c2, c3, c4, kelvins(T))
        .map(|r| PyLiquidHeatCapacityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pure-component vapour pressure, from NeqSim's Antoine correlation.
/// The three pure-component vapour pressures of the water-nitric-sulfuric acid system.
#[pyfunction]
#[pyo3(signature = (T))]
#[pyo3(text_signature = "(T)")]
#[allow(non_snake_case)] // `T` is the symbol in the chemistry
pub fn nitric_sulfuric_acid_vapor_pressure(
    py: Python<'_>,
    T: f64,
) -> PyResult<PyNitricSulfuricAcidVaporPressureResult> {
    azoth_eos::nitric_sulfuric_acid_vapor_pressure(kelvins(T))
        .map(|r| PyNitricSulfuricAcidVaporPressureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

#[pyfunction]
#[pyo3(signature = (A, B, C, D, E, form, Tc, Pc, T))]
#[pyo3(text_signature = "(A, B, C, D, E, form, Tc, Pc, T)")]
#[allow(non_snake_case)] // `A`-`E`, `Tc`, `Pc` and `T` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn antoine_vapor_pressure(
    py: Python<'_>,
    A: f64,
    B: f64,
    C: f64,
    D: f64,
    E: f64,
    form: &str,
    Tc: f64,
    Pc: f64,
    T: f64,
) -> PyResult<PyAntoineVaporPressureResult> {
    let form: azoth_eos::AntoineForm = form
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    azoth_eos::antoine_vapor_pressure(A, B, C, D, E, form, kelvins(Tc), pascals(Pc), kelvins(T))
        .map(|r| PyAntoineVaporPressureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The saturated liquid molar volume, from the Spencer-Danner Rackett equation.
#[pyfunction]
#[pyo3(signature = (omega, Tc, Pc, T))]
#[pyo3(text_signature = "(omega, Tc, Pc, T)")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the published equation
pub fn rackett_molar_volume(
    py: Python<'_>,
    omega: f64,
    Tc: f64,
    Pc: f64,
    T: f64,
) -> PyResult<PyRackettMolarVolumeResult> {
    azoth_eos::rackett_molar_volume(omega, kelvins(Tc), pascals(Pc), kelvins(T))
        .map(|r| PyRackettMolarVolumeResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The saturated liquid molar volume, from the Hankinson-Thomson COSTALD equation.
#[pyfunction]
#[pyo3(signature = (omega, Tc, Vc, M, rho_normal, T))]
#[pyo3(text_signature = "(omega, Tc, Vc, M, rho_normal, T)")]
#[allow(non_snake_case)] // `Tc`, `Vc`, `M` and `T` are the symbols in the published equation
pub fn costald_molar_volume(
    py: Python<'_>,
    omega: f64,
    Tc: f64,
    Vc: f64,
    M: f64,
    rho_normal: f64,
    T: f64,
) -> PyResult<PyCostaldMolarVolumeResult> {
    azoth_eos::costald_molar_volume(
        omega,
        kelvins(Tc),
        cubic_meters_per_mole(Vc),
        kilograms_per_mole(M),
        kilograms_per_cubic_meter(rho_normal),
        kelvins(T),
    )
    .map(|r| PyCostaldMolarVolumeResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The gas dynamic viscosity, from the Chung correlation.
#[pyfunction]
#[pyo3(signature = (omega, Tc, Vc, M, dipole, kappa, T, V))]
#[pyo3(text_signature = "(omega, Tc, Vc, M, dipole, kappa, T, V)")]
#[allow(non_snake_case)] // `Tc`, `Vc`, `M`, `T` and `V` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn chung_viscosity(
    py: Python<'_>,
    omega: f64,
    Tc: f64,
    Vc: f64,
    M: f64,
    dipole: f64,
    kappa: f64,
    T: f64,
    V: f64,
) -> PyResult<PyChungViscosityResult> {
    azoth_eos::chung_viscosity(
        omega,
        kelvins(Tc),
        cubic_meters_per_mole(Vc),
        kilograms_per_mole(M),
        dipole,
        kappa,
        kelvins(T),
        cubic_meters_per_mole(V),
    )
    .map(|r| PyChungViscosityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The gas thermal conductivity, from the Chung correlation.
#[pyfunction]
#[pyo3(signature = (Cv0, M, omega, Tc, Vc, dipole, kappa, T))]
#[pyo3(text_signature = "(Cv0, M, omega, Tc, Vc, dipole, kappa, T)")]
#[allow(non_snake_case)] // `Cv0`, `Tc`, `Vc`, `M` and `T` are the symbols in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn chung_conductivity(
    py: Python<'_>,
    Cv0: f64,
    M: f64,
    omega: f64,
    Tc: f64,
    Vc: f64,
    dipole: f64,
    kappa: f64,
    T: f64,
) -> PyResult<PyChungConductivityResult> {
    azoth_eos::chung_conductivity(
        joules_per_mole_kelvin(Cv0),
        kilograms_per_mole(M),
        omega,
        kelvins(Tc),
        cubic_meters_per_mole(Vc),
        dipole,
        kappa,
        kelvins(T),
    )
    .map(|r| PyChungConductivityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The liquid binary diffusivity, from the Tyn-Calus correlation.
#[pyfunction]
#[pyo3(signature = (VA, VB, T, eta))]
#[pyo3(text_signature = "(VA, VB, T, eta)")]
#[allow(non_snake_case)] // `VA`, `VB` and `T` are the symbols in the published equation
pub fn tyn_calus_diffusivity(
    py: Python<'_>,
    VA: f64,
    VB: f64,
    T: f64,
    eta: f64,
) -> PyResult<PyTynCalusDiffusivityResult> {
    azoth_eos::tyn_calus_diffusivity(
        cubic_meters_per_mole(VA),
        cubic_meters_per_mole(VB),
        kelvins(T),
        pascal_seconds(eta),
    )
    .map(|r| PyTynCalusDiffusivityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The liquid binary diffusivity, from the Wilke-Chang correlation.
#[pyfunction]
#[pyo3(signature = (phi, M, T, eta, VA))]
#[pyo3(text_signature = "(phi, M, T, eta, VA)")]
#[allow(non_snake_case)] // `M`, `T` and `VA` are the symbols in the published equation
pub fn wilke_chang_diffusivity(
    py: Python<'_>,
    phi: f64,
    M: f64,
    T: f64,
    eta: f64,
    VA: f64,
) -> PyResult<PyWilkeChangDiffusivityResult> {
    azoth_eos::wilke_chang_diffusivity(
        phi,
        kilograms_per_mole(M),
        kelvins(T),
        pascal_seconds(eta),
        cubic_meters_per_mole(VA),
    )
    .map(|r| PyWilkeChangDiffusivityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The liquid binary diffusivity, from the Hayduk-Minhas correlation.
#[pyfunction]
#[pyo3(signature = (form, VA, T, eta))]
#[pyo3(text_signature = "(form, VA, T, eta)")]
#[allow(non_snake_case)] // `VA` and `T` are the symbols in the published equation
pub fn hayduk_minhas_diffusivity(
    py: Python<'_>,
    form: &str,
    VA: f64,
    T: f64,
    eta: f64,
) -> PyResult<PyHaydukMinhasDiffusivityResult> {
    let form: azoth_eos::HaydukMinhasForm = form
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    azoth_eos::hayduk_minhas_diffusivity(
        form,
        cubic_meters_per_mole(VA),
        kelvins(T),
        pascal_seconds(eta),
    )
    .map(|r| PyHaydukMinhasDiffusivityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The Schwartzentruber-Renon alpha function.
#[pyfunction]
#[pyo3(signature = (omega, p1, p2, p3, Tr))]
#[pyo3(text_signature = "(omega, p1, p2, p3, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn schwartzentruber_alpha(
    py: Python<'_>,
    omega: f64,
    p1: f64,
    p2: f64,
    p3: f64,
    Tr: f64,
) -> PyResult<PySchwartzentruberAlphaResult> {
    azoth_eos::schwartzentruber_alpha(omega, p1, p2, p3, Tr)
        .map(|r| PySchwartzentruberAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Soreide-Whitson alpha function for water.
#[pyfunction]
#[pyo3(signature = (salinity, Tr))]
#[pyo3(text_signature = "(salinity, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn soreide_whitson_alpha(
    py: Python<'_>,
    salinity: f64,
    Tr: f64,
) -> PyResult<PySoreideWhitsonAlphaResult> {
    azoth_eos::soreide_whitson_alpha(salinity, Tr)
        .map(|r| PySoreideWhitsonAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The liquid binary diffusivity, from the Siddiqi-Lucas correlation.
#[pyfunction]
#[pyo3(signature = (form, VA, VB, T, eta))]
#[pyo3(text_signature = "(form, VA, VB, T, eta)")]
#[allow(non_snake_case)] // `VA`, `VB` and `T` are the symbols in the published equation
pub fn siddiqi_lucas_diffusivity(
    py: Python<'_>,
    form: &str,
    VA: f64,
    VB: f64,
    T: f64,
    eta: f64,
) -> PyResult<PySiddiqiLucasDiffusivityResult> {
    let form: azoth_eos::SiddiqiLucasForm = form
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    azoth_eos::siddiqi_lucas_diffusivity(
        form,
        cubic_meters_per_mole(VA),
        cubic_meters_per_mole(VB),
        kelvins(T),
        pascal_seconds(eta),
    )
    .map(|r| PySiddiqiLucasDiffusivityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The CO2-in-water binary diffusivity, from NeqSim's `CO2water` correlation.
#[pyfunction]
#[pyo3(signature = (T))]
#[pyo3(text_signature = "(T)")]
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn co2_water_diffusivity(py: Python<'_>, T: f64) -> PyResult<PyCo2WaterDiffusivityResult> {
    azoth_eos::co2_water_diffusivity(kelvins(T))
        .map(|r| PyCo2WaterDiffusivityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The interface surface tension between a gas and a liquid, from the parachor.
///
/// The mixture form of the pure-component correlation, and a *model* rather than a
/// calculation: both compositions are vectors.
#[pyfunction]
#[pyo3(signature = (parachors, rho_gas, M_gas, x_gas, rho_liquid, M_liquid, x_liquid))]
#[pyo3(text_signature = "(parachors, rho_gas, M_gas, x_gas, rho_liquid, M_liquid, x_liquid)")]
#[allow(non_snake_case)] // `M_gas` and `M_liquid` are the symbols in the equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn parachor_mixture_surface_tension(
    py: Python<'_>,
    parachors: Vec<f64>,
    rho_gas: f64,
    M_gas: f64,
    x_gas: Vec<f64>,
    rho_liquid: f64,
    M_liquid: f64,
    x_liquid: Vec<f64>,
) -> PyResult<crate::results::PyParachorMixtureSurfaceTensionResult> {
    azoth_eos::parachor_mixture_surface_tension(
        &parachors,
        kilograms_per_cubic_meter(rho_gas),
        kilograms_per_mole(M_gas),
        &x_gas,
        kilograms_per_cubic_meter(rho_liquid),
        kilograms_per_mole(M_liquid),
        &x_liquid,
    )
    .map(|r| crate::results::PyParachorMixtureSurfaceTensionResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The surface tension from the parachor (Macleod-Sugden) correlation.
#[pyfunction]
#[pyo3(signature = (parachor, rho_l, rho_v, M))]
#[pyo3(text_signature = "(parachor, rho_l, rho_v, M)")]
#[allow(non_snake_case)] // `M` is the symbol in the published equation
pub fn parachor_surface_tension(
    py: Python<'_>,
    parachor: f64,
    rho_l: f64,
    rho_v: f64,
    M: f64,
) -> PyResult<PyParachorSurfaceTensionResult> {
    azoth_eos::parachor_surface_tension(
        parachor,
        kilograms_per_cubic_meter(rho_l),
        kilograms_per_cubic_meter(rho_v),
        kilograms_per_mole(M),
    )
    .map(|r| PyParachorSurfaceTensionResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The gas mixture dynamic viscosity, from Wilke's rule over the pure Chung
/// viscosities. A *model* rather than a calculation: its arguments are vectors.
#[pyfunction]
#[pyo3(signature = (Tc, Vc, M, omega, dipole, kappa, T, V, z))]
#[pyo3(text_signature = "(Tc, Vc, M, omega, dipole, kappa, T, V, z)")]
#[allow(non_snake_case)] // `Tc`, `Vc`, `M`, `T` and `V` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn wilke_viscosity(
    py: Python<'_>,
    Tc: Vec<f64>,
    Vc: Vec<f64>,
    M: Vec<f64>,
    omega: Vec<f64>,
    dipole: Vec<f64>,
    kappa: Vec<f64>,
    T: f64,
    V: f64,
    z: Vec<f64>,
) -> PyResult<PyWilkeViscosityResult> {
    azoth_eos::wilke_viscosity(&Tc, &Vc, &M, &omega, &dipole, &kappa, T, V, &z)
        .map(|r| PyWilkeViscosityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The gas mixture thermal conductivity, from Mason-Saxena mixing over the pure
/// Chung conductivities. A *model* rather than a calculation: its arguments are
/// vectors.
#[pyfunction]
#[pyo3(signature = (Cv0, M, omega, Tc, Vc, dipole, kappa, T, z))]
#[pyo3(text_signature = "(Cv0, M, omega, Tc, Vc, dipole, kappa, T, z)")]
#[allow(non_snake_case)] // `Cv0`, `Tc`, `Vc`, `M` and `T` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn mason_saxena_conductivity(
    py: Python<'_>,
    Cv0: Vec<f64>,
    M: Vec<f64>,
    omega: Vec<f64>,
    Tc: Vec<f64>,
    Vc: Vec<f64>,
    dipole: Vec<f64>,
    kappa: Vec<f64>,
    T: f64,
    z: Vec<f64>,
) -> PyResult<PyMasonSaxenaConductivityResult> {
    azoth_eos::mason_saxena_conductivity(&Cv0, &M, &omega, &Tc, &Vc, &dipole, &kappa, T, &z)
        .map(|r| PyMasonSaxenaConductivityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from NRTL. A *model* rather than a
/// calculation: its resolved parameters cross the boundary flattened row-major.
#[pyfunction]
#[pyo3(signature = (alpha, dij, T, x))]
#[pyo3(text_signature = "(alpha, dij, T, x)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn nrtl_activity_coefficients(
    py: Python<'_>,
    alpha: Vec<f64>,
    dij: Vec<f64>,
    T: f64,
    x: Vec<f64>,
) -> PyResult<PyNrtlActivityCoefficientsResult> {
    let params = azoth_eos::databank::NrtlParameters { alpha, dij };
    azoth_eos::nrtl_activity_coefficients(&params, T, &x)
        .map(|r| PyNrtlActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The fugacity coefficients of an NRTL activity-coefficient liquid. A *model* rather
/// than a calculation: the resolved phase parameters cross flattened, one list per
/// field of `GeNrtlPhaseParameters`, and the Antoine record's own fields with them.
#[pyfunction]
#[pyo3(signature = (alpha, dij, antoine_type, antoine_coefficients, antoine_tc, antoine_pc, T, P, x))]
#[pyo3(
    text_signature = "(alpha, dij, antoine_type, antoine_coefficients, antoine_tc, antoine_pc, T, P, x)"
)]
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the resolved record's own fields.
pub fn ge_nrtl_phase(
    py: Python<'_>,
    alpha: Vec<f64>,
    dij: Vec<f64>,
    antoine_type: Vec<String>,
    antoine_coefficients: Vec<f64>,
    antoine_tc: Vec<f64>,
    antoine_pc: Vec<f64>,
    T: f64,
    P: f64,
    x: Vec<f64>,
) -> PyResult<PyGeNrtlPhaseResult> {
    let params = azoth_eos::databank::GeNrtlPhaseParameters {
        alpha,
        dij,
        antoine: antoine_records(
            antoine_type,
            &antoine_coefficients,
            &antoine_tc,
            &antoine_pc,
        ),
    };
    azoth_eos::ge_nrtl_phase(&params, T, P, &x)
        .map(|r| PyGeNrtlPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// One [`azoth_eos::databank::AntoineRecord`] per component, from the transport's four
/// parallel lists.
///
/// The coefficients are component-major, five per component, which is the order the
/// parameter record carries them in on the Python side.
fn antoine_records(
    antoine_type: Vec<String>,
    antoine_coefficients: &[f64],
    antoine_tc: &[f64],
    antoine_pc: &[f64],
) -> Vec<azoth_eos::databank::AntoineRecord> {
    antoine_type
        .into_iter()
        .enumerate()
        .map(|(i, antoine_type)| azoth_eos::databank::AntoineRecord {
            antoine_type,
            coefficients: antoine_coefficients[i * 5..i * 5 + 5]
                .try_into()
                .expect("5 coefficients per component"),
            tc: antoine_tc[i],
            pc: antoine_pc[i],
        })
        .collect()
}

/// The isothermal flash of an SRK vapour over an NRTL liquid. A *model* rather than a
/// calculation: the mixture crosses flattened, the parameter record's own fields with
/// it, and the NRTL matrix keeps its record name while the cubic's alpha correlation
/// is spelled `cubic_alpha`.
#[pyfunction]
#[pyo3(signature = (
    Tc,
    Pc,
    omega,
    kij,
    association,
    alpha,
    dij,
    antoine_type,
    antoine_coefficients,
    antoine_tc,
    antoine_pc,
    T,
    P,
    z,
    eos = "pr",
    cubic_alpha = "pr",
    alpha_params = None
))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, alpha, dij, antoine_type, antoine_coefficients, \
                         antoine_tc, antoine_pc, T, P, z, eos = \"pr\", cubic_alpha = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T`, `P` and `z` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the resolved record's own fields.
pub fn ge_nrtl_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    alpha: Vec<f64>,
    dij: Vec<f64>,
    antoine_type: Vec<String>,
    antoine_coefficients: Vec<f64>,
    antoine_tc: Vec<f64>,
    antoine_pc: Vec<f64>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    cubic_alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyGeNrtlFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        cubic_alpha,
        alpha_params.as_deref(),
    )?;
    let params = azoth_eos::databank::GeNrtlPhaseParameters {
        alpha,
        dij,
        antoine: antoine_records(
            antoine_type,
            &antoine_coefficients,
            &antoine_tc,
            &antoine_pc,
        ),
    };
    azoth_eos::ge_nrtl_flash::ge_nrtl_flash(&params, &mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyGeNrtlFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The fugacity coefficients of a UNIFAC activity-coefficient liquid. A *model* rather
/// than a calculation: the resolved group tables and the per-component vapour-pressure
/// columns cross flattened, one list per field of `GeUnifacPhaseParameters`.
#[pyfunction]
#[pyo3(signature = (
    groups,
    group_r,
    group_q,
    aij,
    antoine_type,
    antoine_coefficients,
    antoine_tc,
    antoine_pc,
    T,
    P,
    x
))]
#[pyo3(
    text_signature = "(groups, group_r, group_q, aij, antoine_type, antoine_coefficients, \
                         antoine_tc, antoine_pc, T, P, x)"
)]
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the resolved record's own fields.
pub fn ge_unifac_phase(
    py: Python<'_>,
    groups: Vec<f64>,
    group_r: Vec<f64>,
    group_q: Vec<f64>,
    aij: Vec<f64>,
    antoine_type: Vec<String>,
    antoine_coefficients: Vec<f64>,
    antoine_tc: Vec<f64>,
    antoine_pc: Vec<f64>,
    T: f64,
    P: f64,
    x: Vec<f64>,
) -> PyResult<PyGeUnifacPhaseResult> {
    let params = azoth_eos::databank::GeUnifacPhaseParameters {
        groups,
        group_r,
        group_q,
        aij,
        antoine: antoine_records(
            antoine_type,
            &antoine_coefficients,
            &antoine_tc,
            &antoine_pc,
        ),
    };
    azoth_eos::ge_unifac_phase::ge_unifac_phase(&params, T, P, &x)
        .map(|r| PyGeUnifacPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The fugacity coefficients of a Wilson activity-coefficient liquid. A *model* rather
/// than a calculation: the mixture crosses flattened with each component's molar mass,
/// which is what the Wilson correlation reads, and the vapour-pressure columns with it.
#[pyfunction]
#[pyo3(signature = (
    Tc,
    Pc,
    omega,
    kij,
    association,
    molar_mass,
    antoine_type,
    antoine_coefficients,
    antoine_tc,
    antoine_pc,
    T,
    P,
    x,
    eos = "pr",
    alpha = "pr",
    alpha_params = None
))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, molar_mass, antoine_type, antoine_coefficients, \
                         antoine_tc, antoine_pc, T, P, x)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T`, `P` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the resolved record's own fields.
pub fn ge_wilson_phase(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    molar_mass: Vec<f64>,
    antoine_type: Vec<String>,
    antoine_coefficients: Vec<f64>,
    antoine_tc: Vec<f64>,
    antoine_pc: Vec<f64>,
    T: f64,
    P: f64,
    x: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyGeWilsonPhaseResult> {
    let mixture = build_mixture_with_mass(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        &molar_mass,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let params = azoth_eos::databank::GeWilsonPhaseParameters {
        antoine: antoine_records(
            antoine_type,
            &antoine_coefficients,
            &antoine_tc,
            &antoine_pc,
        ),
    };
    azoth_eos::ge_wilson_phase::ge_wilson_phase(&params, &mixture, T, P, &x)
        .map(|r| PyGeWilsonPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a Desmukh-Mather electrolyte phase.
///
/// The extended Debye-Huckel plus pair-sum expression, then NeqSim's mole-fraction
/// conversion, and the `(gamma / gamma^infinity) H / P` fugacity branch. The solvent is
/// whatever `REFERENCESTATETYPE` says it is, not whatever is named `water`.
#[pyfunction]
#[pyo3(signature = (components, T, P, x))]
#[pyo3(text_signature = "(components, T, P, x)")]
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn desmukh_mather_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    x: Vec<f64>,
) -> PyResult<PyDesmukhMatherPhaseResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    azoth_eos::desmukh_mather_phase::desmukh_mather_phase(&names, T, P, &x)
        .map(|r| PyDesmukhMatherPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The fugacity coefficients of a Kent-Eisenberg phase, whose activity coefficients are one.
///
/// Takes the component names, because the branch each component takes is its databank
/// `REFERENCESTATETYPE` and its charge - a `solvent` gets `P0_i(T)/P`, a neutral solute
/// `H_i(T)/P`, and an ion the constant `1e8`.
#[pyfunction]
#[pyo3(signature = (components, T, P, x))]
#[pyo3(text_signature = "(components, T, P, x)")]
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn kent_eisenberg_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    x: Vec<f64>,
) -> PyResult<PyKentEisenbergPhaseResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    azoth_eos::kent_eisenberg_phase::kent_eisenberg_phase(&names, T, P, &x)
        .map(|r| PyKentEisenbergPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity and fugacity coefficients of an electrolyte phase whose non-ideality is
/// Pitzer's.
///
/// Takes the component **names** rather than resolved arrays, because the names are what
/// the parameter datasets are keyed by: the model resolves each one's charge, molar mass
/// and reference state against the databank itself. `P` moves no activity coefficient -
/// `getGamma`'s own pressure argument is unused - and the fugacity coefficients divide by
/// it, in bar.
#[pyfunction]
#[pyo3(signature = (components, T, P, x))]
#[pyo3(text_signature = "(components, T, P, x)")]
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn pitzer_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    x: Vec<f64>,
) -> PyResult<PyPitzerPhaseResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    azoth_eos::pitzer_phase(&names, T, P, &x)
        .map(|r| PyPitzerPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// A mixture whose components carry their molar mass, for the models that read it.
///
/// One mixture's association, as the Python side holds it.
///
/// **A record, and a required argument on every model whose Python side takes a `Mixture`.**
/// The two fitted cubic sets are what make an associating fluid a different fluid — water's
/// fitted covolume is `1.4515e-5 m³/mol` against `0.08664 R Tc/Pc`'s `2.11e-5` — so a mixture
/// whose association does not cross converges to a plausible answer for something else. That
/// is not hypothetical: `eos.pt_flash` through this boundary ran a classical SRK flash on a
/// fluid carrying the CPA interaction column and reported `all_liquid` where the associating
/// model splits at `beta = 0.208383589`.
///
/// The scheme crosses as a **name** (`1A`, `2A`, `2B`, `4C`) because the site count cannot
/// reconstruct it: upstream states `1A` at zero sites, and `2A` and `2B` are different
/// two-site schemes.
///
/// The numbers are in NeqSim's internal scale, which is the scale `AssociationRecord` holds
/// them in because it is the scale the table states them in. The conversion belongs with the
/// model that reads them, exactly as it does on the Python side.
#[pyclass(name = "AssociationSpec")]
#[derive(Debug)]
pub struct PyAssociationSpec {
    /// Whether the phase model runs the association at all — the model's decision, not the
    /// components'.
    #[pyo3(get)]
    pub associating: bool,
    /// One scheme name per component; empty for a component that carries none.
    #[pyo3(get)]
    pub schemes: Vec<String>,
    /// One row per component, in `AssociationRecord`'s field order after the scheme:
    /// `sites`, `energy`, `volume_srk`, `a_srk`, `b_srk`, `m_srk`, `volume_pr`, `a_pr`,
    /// `b_pr`, `m_pr`, `racket_z`, `volume_correction`.
    #[pyo3(get)]
    pub values: Vec<Vec<f64>>,
}

#[pymethods]
impl PyAssociationSpec {
    #[new]
    fn new(associating: bool, schemes: Vec<String>, values: Vec<Vec<f64>>) -> Self {
        Self {
            associating,
            schemes,
            values,
        }
    }
}

impl PyAssociationSpec {
    /// Component `i`'s record, or `None` when its scheme is empty.
    fn record(&self, i: usize) -> PyResult<Option<azoth_eos::association::AssociationRecord>> {
        let Some(name) = self.schemes.get(i).filter(|name| !name.is_empty()) else {
            return Ok(None);
        };
        let scheme =
            azoth_eos::association::SiteScheme::from_databank_name(name).ok_or_else(|| {
                pyo3::exceptions::PyValueError::new_err(format!(
                    "unknown association scheme `{name}`; expected `1A`, `2A`, `2B` or `4C`"
                ))
            })?;
        let row = self.values.get(i).ok_or_else(|| {
            pyo3::exceptions::PyValueError::new_err(format!(
                "component {i} names the scheme `{name}` but carries no parameters"
            ))
        })?;
        if row.len() < 12 {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "component {i}: an association row carries 12 numbers, got {}",
                row.len()
            )));
        }
        Ok(Some(azoth_eos::association::AssociationRecord {
            scheme,
            sites: row[0] as u32,
            energy: row[1],
            volume_srk: row[2],
            a_srk: row[3],
            b_srk: row[4],
            m_srk: row[5],
            volume_pr: row[6],
            a_pr: row[7],
            b_pr: row[8],
            m_pr: row[9],
            racket_z: row[10],
            volume_correction: row[11],
            // The boundary carries the cubic families' fitted sets only: a caller
            // building an association by hand is stating a CPA fluid, and the
            // `UMRCPA_*` set is read from the databank by the model that pairs it with
            // the UMR mixing rule.
            umr_cpa: None,
        }))
    }
}

/// [`build_mixture`] is for the models that read only the cubic's constants; this is the
/// same construction with `Component::with_molar_mass` applied, which `eos.viscosity`,
/// `eos.thermal_conductivity`, `eos.wilson_activity_coefficients` and
/// `eos.ge_wilson_phase` need.
#[allow(non_snake_case)] // `Tc` and `Pc` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The transport's own shape, not a design.
fn build_mixture_with_mass(
    py: Python<'_>,
    Tc: &[f64],
    Pc: &[f64],
    omega: &[f64],
    kij: Vec<f64>,
    association: &PyAssociationSpec,
    molar_mass: &[f64],
    eos: &str,
    alpha: &str,
    alpha_params: Option<&[Vec<f64>]>,
) -> PyResult<azoth_eos::Mixture> {
    let n = Tc.len();
    if Pc.len() != n || omega.len() != n || molar_mass.len() != n {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Tc, Pc, omega and molar_mass must be the same length; got {}, {}, {} and {}",
            Tc.len(),
            Pc.len(),
            omega.len(),
            molar_mass.len()
        )));
    }
    let associations: Vec<Option<azoth_eos::association::AssociationRecord>> = (0..n)
        .map(|i| association.record(i))
        .collect::<PyResult<Vec<_>>>()?;
    let components = (0..n)
        .map(|i| {
            let mut component =
                azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i])?
                    .with_molar_mass(Some(molar_mass[i]));
            if let Some(params) = alpha_params.and_then(|all| all.get(i)) {
                component = component.with_alpha_params(params.clone());
            }
            component = component.with_association(associations[i].clone());
            Ok(component)
        })
        .collect::<azoth_core::Result<Vec<_>>>()
        .map_err(|e| to_pyerr(py, e))?;
    let cubic: azoth_eos::Cubic = eos
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    let alpha: azoth_eos::Alpha = alpha
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    azoth_eos::Mixture::new(components, kij)
        .and_then(|m| {
            if association.associating {
                m.with_association()
            } else {
                Ok(m)
            }
        })
        .map(|m| m.with_cubic(cubic).with_alpha(alpha))
        .map_err(|e| to_pyerr(py, e))
}

/// The fugacity coefficients of a UNIQUAC activity-coefficient liquid. A *model* rather
/// than a calculation: the resolved record's fields cross flattened, and `aij` with them,
/// because no upstream table carries a UNIQUAC interaction matrix.
#[pyfunction]
#[pyo3(signature = (
    r,
    q,
    antoine_type,
    antoine_coefficients,
    antoine_tc,
    antoine_pc,
    T,
    P,
    x,
    aij
))]
#[pyo3(
    text_signature = "(r, q, antoine_type, antoine_coefficients, antoine_tc, antoine_pc, \
                         T, P, x, aij)"
)]
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the resolved record's own fields.
pub fn ge_uniquac_phase(
    py: Python<'_>,
    r: Vec<f64>,
    q: Vec<f64>,
    antoine_type: Vec<String>,
    antoine_coefficients: Vec<f64>,
    antoine_tc: Vec<f64>,
    antoine_pc: Vec<f64>,
    T: f64,
    P: f64,
    x: Vec<f64>,
    aij: Vec<Vec<f64>>,
) -> PyResult<PyGeUniquacPhaseResult> {
    let params = azoth_eos::databank::GeUniquacPhaseParameters {
        r,
        q,
        antoine: antoine_records(
            antoine_type,
            &antoine_coefficients,
            &antoine_tc,
            &antoine_pc,
        ),
    };
    azoth_eos::ge_uniquac_phase::ge_uniquac_phase(&params, T, P, &x, &aij)
        .map(|res| PyGeUniquacPhaseResult::from(&res))
        .map_err(|e| to_pyerr(py, e))
}

/// The fugacity coefficients of the water-nitric-sulfuric acid liquid. A *model* rather
/// than a calculation: the resolved record's fields cross flattened, the acid identities
/// first and the vapour-pressure columns beside them.
#[pyfunction]
#[pyo3(signature = (acid_index, antoine_type, antoine_coefficients, antoine_tc, antoine_pc, T, P, x))]
#[pyo3(
    text_signature = "(acid_index, antoine_type, antoine_coefficients, antoine_tc, antoine_pc, T, P, x)"
)]
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the resolved record's own fields.
pub fn ge_van_laar_acid_phase(
    py: Python<'_>,
    acid_index: Vec<u8>,
    antoine_type: Vec<String>,
    antoine_coefficients: Vec<f64>,
    antoine_tc: Vec<f64>,
    antoine_pc: Vec<f64>,
    T: f64,
    P: f64,
    x: Vec<f64>,
) -> PyResult<PyGeVanLaarAcidPhaseResult> {
    let params = azoth_eos::databank::GeVanLaarAcidPhaseParameters {
        acid_index,
        antoine: antoine_records(
            antoine_type,
            &antoine_coefficients,
            &antoine_tc,
            &antoine_pc,
        ),
    };
    azoth_eos::ge_van_laar_acid_phase::ge_van_laar_acid_phase(&params, T, P, &x)
        .map(|r| PyGeVanLaarAcidPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture containing water, nitric acid and sulfuric
/// acid, from the Taleb-Ponche-Mirabel Van Laar model. A *model* rather than a
/// calculation: the resolved acid identities cross flattened.
#[pyfunction]
#[pyo3(signature = (acid_index, T, x))]
#[pyo3(text_signature = "(acid_index, T, x)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn van_laar_acid_activity_coefficients(
    py: Python<'_>,
    acid_index: Vec<u8>,
    T: f64,
    x: Vec<f64>,
) -> PyResult<PyVanLaarAcidActivityCoefficientsResult> {
    let params = azoth_eos::databank::VanLaarAcidParameters { acid_index };
    azoth_eos::van_laar_acid_activity_coefficients(&params, T, &x)
        .map(|r| PyVanLaarAcidActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from UNIFAC with the PSRK interaction
/// parameters. A *model* rather than a calculation: the resolved group basis and the
/// three interaction matrices cross flattened.
#[pyfunction]
#[pyo3(signature = (groups, group_r, group_q, aij, bij, cij, T, x))]
#[pyo3(text_signature = "(groups, group_r, group_q, aij, bij, cij, T, x)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn unifac_psrk_activity_coefficients(
    py: Python<'_>,
    groups: Vec<f64>,
    group_r: Vec<f64>,
    group_q: Vec<f64>,
    aij: Vec<f64>,
    bij: Vec<f64>,
    cij: Vec<f64>,
    T: f64,
    x: Vec<f64>,
) -> PyResult<PyUnifacPsrkActivityCoefficientsResult> {
    let params = azoth_eos::databank::UnifacPsrkParameters {
        groups,
        group_r,
        group_q,
        aij,
        bij,
        cij,
    };
    azoth_eos::unifac_psrk_activity_coefficients(&params, T, &x)
        .map(|r| PyUnifacPsrkActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from UNIFAC with the UMR-PRU parameters.
/// A *model* rather than a calculation: the resolved basis and the three interaction
/// matrices of the chosen set cross flattened.
#[pyfunction]
#[pyo3(signature = (groups, group_r, group_q, aij, bij, cij, T, x))]
#[pyo3(text_signature = "(groups, group_r, group_q, aij, bij, cij, T, x)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn unifac_umrpru_activity_coefficients(
    py: Python<'_>,
    groups: Vec<f64>,
    group_r: Vec<f64>,
    group_q: Vec<f64>,
    aij: Vec<f64>,
    bij: Vec<f64>,
    cij: Vec<f64>,
    T: f64,
    x: Vec<f64>,
) -> PyResult<PyUnifacUmrpruActivityCoefficientsResult> {
    let params = azoth_eos::databank::UnifacUmrpruParameters {
        groups,
        group_r,
        group_q,
        aij,
        bij,
        cij,
    };
    azoth_eos::unifac_umrpru_activity_coefficients(&params, T, &x)
        .map(|r| PyUnifacUmrpruActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The UMR-PR alpha function.
#[pyfunction]
#[pyo3(signature = (omega, Tr))]
#[pyo3(text_signature = "(omega, Tr)")]
#[allow(non_snake_case)] // `Tr` is the symbol in the published equation
pub fn umrpr_alpha(py: Python<'_>, omega: f64, Tr: f64) -> PyResult<PyUmrprAlphaResult> {
    azoth_eos::umrpr_alpha(omega, Tr)
        .map(|r| PyUmrprAlphaResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from UNIFAC. A *model* rather than a
/// calculation: its resolved parameters cross the boundary flattened row-major.
#[pyfunction]
#[pyo3(signature = (groups, group_r, group_q, aij, T, x))]
#[pyo3(text_signature = "(groups, group_r, group_q, aij, T, x)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn unifac_activity_coefficients(
    py: Python<'_>,
    groups: Vec<f64>,
    group_r: Vec<f64>,
    group_q: Vec<f64>,
    aij: Vec<f64>,
    T: f64,
    x: Vec<f64>,
) -> PyResult<PyUnifacActivityCoefficientsResult> {
    let params = azoth_eos::databank::UnifacParameters {
        groups,
        group_r,
        group_q,
        aij,
    };
    azoth_eos::unifac_activity_coefficients(&params, T, &x)
        .map(|r| PyUnifacActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from UNIQUAC. A *model* rather than a
/// calculation: its resolved `r`/`q` cross flattened, and `aij` stays the caller's.
#[pyfunction]
#[pyo3(signature = (r, q, T, x, aij))]
#[pyo3(text_signature = "(r, q, T, x, aij)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn uniquac_activity_coefficients(
    py: Python<'_>,
    r: Vec<f64>,
    q: Vec<f64>,
    T: f64,
    x: Vec<f64>,
    aij: Vec<Vec<f64>>,
) -> PyResult<PyUniquacActivityCoefficientsResult> {
    let params = azoth_eos::databank::UniquacParameters { r, q };
    azoth_eos::uniquac_activity_coefficients(&params, T, &x, &aij)
        .map(|result| PyUniquacActivityCoefficientsResult::from(&result))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from the paraffin-wax Wilson model. A
/// *model* rather than a calculation: it reads the resolved critical constants and
/// molar mass, the same prefix the cubic models cross with.
///
/// `eos`, `alpha` and `alpha_params` are accepted for the boundary's uniformity with the
/// other mixture models and ignored: the Coutinho correlation is an activity model, not
/// a cubic, so the mixture's own cubic and alpha are not consulted.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, molar_mass, T, x, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, molar_mass, T, x, eos = \"pr\", alpha = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `x` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
#[allow(unused_variables)] // `eos`, `alpha` and `alpha_params` are boundary-only, see above.
pub fn wilson_activity_coefficients(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    molar_mass: Vec<f64>,
    T: f64,
    x: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyWilsonActivityCoefficientsResult> {
    let n = Tc.len();
    let components = (0..n)
        .map(|i| {
            azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i])
                .map(|c| c.with_molar_mass(Some(molar_mass[i])))
        })
        .collect::<Result<Vec<_>, _>>()
        .map_err(|e| to_pyerr(py, e))?;
    let mixture = azoth_eos::Mixture::new(components, kij).map_err(|e| to_pyerr(py, e))?;
    azoth_eos::wilson_activity_coefficients(&mixture, T, &x)
        .map(|r| PyWilsonActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The saturation pressure of a pure component.
///
/// A *model* rather than a calculation - its spec fixes a procedure and it composes
/// the kernels above rather than adding arithmetic of its own. It crosses the
/// boundary like any other function: SI magnitudes in, a quantity out.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, T))]
#[pyo3(text_signature = "(Tc, Pc, omega, T)")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the chemistry
pub fn pure_saturation(
    py: Python<'_>,
    Tc: f64,
    Pc: f64,
    omega: f64,
    T: f64,
) -> PyResult<PyPureSaturationResult> {
    azoth_eos::pure_saturation(kelvins(Tc), pascals(Pc), omega, kelvins(T))
        .map(|r| PyPureSaturationResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The two-phase flash of a mixture at a temperature and pressure.
///
/// The first function here whose arguments are vectors, and the first whose result
/// has any. That is the whole reason the flash is a *model* rather than a
/// calculation: the registry is scalar, so a composition vector has nowhere to live
/// in a spec's `inputs` - and even if it did, a registered calc's guarantee is a
/// worked example a human can retrace, which "40 components, converged in 21
/// iterations" is not.
///
/// `kij` crosses flattened row-major, with `N = Tc.len()`. The Python side builds
/// the nested shape the caller wrote and flattens it here, so the boundary carries
/// one list rather than a list of lists, and the component order is the one thing
/// the two sides have to agree about.
///
/// `beta` is `Option<f64>` on purpose: `None` crosses as `None`. A sentinel number
/// at this boundary would undo, one layer below where it was decided, the design
/// that stops a caller mistaking a trivial solution for a phase split.
/// A mixture from the three per-component vectors and a flattened `kij` matrix.
///
/// Shared by the flash and the two phase-boundary models: all three take the same
/// component arguments in the same order, and building the object in one place is
/// what keeps the component ordering from being decided three times.
///
/// `pub(crate)` because `process.rs` builds the same object from the same arguments.
/// A unit operation takes a mixture exactly as a flash does, so a second construction
/// here would be a second place for the component order to be decided.
#[allow(non_snake_case)] // `Tc` and `Pc` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub(crate) fn build_mixture(
    py: Python<'_>,
    Tc: &[f64],
    Pc: &[f64],
    omega: &[f64],
    kij: Vec<f64>,
    association: &PyAssociationSpec,
    eos: &str,
    alpha: &str,
    alpha_params: Option<&[Vec<f64>]>,
) -> PyResult<azoth_eos::Mixture> {
    let n = Tc.len();
    if Pc.len() != n || omega.len() != n {
        return Err(pyo3::exceptions::PyValueError::new_err(format!(
            "Tc, Pc and omega must be the same length; got {}, {} and {}",
            Tc.len(),
            Pc.len(),
            omega.len()
        )));
    }
    let associations: Vec<Option<azoth_eos::association::AssociationRecord>> = (0..n)
        .map(|i| association.record(i))
        .collect::<PyResult<Vec<_>>>()?;
    let components = (0..n)
        .map(|i| {
            let mut component =
                azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i])?;
            if let Some(params) = alpha_params.and_then(|all| all.get(i)) {
                component = component.with_alpha_params(params.clone());
            }
            component = component.with_association(associations[i].clone());
            Ok(component)
        })
        .collect::<azoth_core::Result<Vec<_>>>()
        .map_err(|e| to_pyerr(py, e))?;
    let cubic: azoth_eos::Cubic = eos
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    let alpha: azoth_eos::Alpha = alpha
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    azoth_eos::Mixture::new(components, kij)
        .and_then(|m| {
            if association.associating {
                m.with_association()
            } else {
                Ok(m)
            }
        })
        .map(|m| m.with_cubic(cubic).with_alpha(alpha))
        .map_err(|e| to_pyerr(py, e))
}

#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, T, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pt_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPtFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::pt_flash(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyPtFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The PT phase envelope of a mixture of composition `z`, traced from a low pressure.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pt_phase_envelope(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPhaseEnvelopeResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::pt_phase_envelope(&mixture, pascals(P), &z)
        .map(|r| PyPhaseEnvelopeResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The temperature at which a mixture has a given molar enthalpy at a pressure.
///
/// The ideal-gas vectors are the same six `eos.molar_enthalpy_entropy` takes, and they
/// are required for the same reason: the requested enthalpy is a *difference* from the
/// datum they carry, and a coefficient set without a reference state is not a
/// thermodynamic model. `s_ref` and `P_ref` do not enter this model's answer - an
/// isenthalpic flash does not need an entropy - and are taken because they belong to
/// the same model.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, H, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, H, z, eos = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T_ref` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn ph_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    P: f64,
    H: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPhFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::ph_flash(&mixture, &ideal_gas, pascals(P), joules_per_mole(H), &z)
        .map(|r| PyPhFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The temperature at which a mixture has a given molar entropy at a pressure.
///
/// The isentropic companion to `ph_flash`, and the same fifteen arguments: a compressor
/// or an expander assumed ideal knows the pressure it leaves at and the entropy it
/// arrived with, and not the temperature that results.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, S, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, S, z, eos = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T_ref` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn ps_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    P: f64,
    S: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPsFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::ps_flash(
        &mixture,
        &ideal_gas,
        pascals(P),
        joules_per_mole_kelvin(S),
        &z,
    )
    .map(|r| PyPsFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The pressure at which a mixture has a given molar volume at a temperature.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, V, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, V, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn tv_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    T: f64,
    V: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyTvFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::tv_flash(
        &mixture,
        &ideal_gas,
        kelvins(T),
        cubic_meters_per_mole(V),
        &z,
    )
    .map(|r| PyTvFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The temperature at which a mixture has a given molar volume at a pressure.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, V, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, V, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `P` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pv_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    P: f64,
    V: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPvFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::pv_flash(
        &mixture,
        &ideal_gas,
        pascals(P),
        cubic_meters_per_mole(V),
        &z,
    )
    .map(|r| PyPvFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The th-flash flash (T,H -> P).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, H, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, H, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn th_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    T: f64,
    H: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyThFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::th_flash(&mixture, &ideal_gas, kelvins(T), joules_per_mole(H), &z)
        .map(|r| PyThFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The ts-flash flash (T,S -> P).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, S, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, S, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn ts_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    T: f64,
    S: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyTsFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::ts_flash(
        &mixture,
        &ideal_gas,
        kelvins(T),
        joules_per_mole_kelvin(S),
        &z,
    )
    .map(|r| PyTsFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The tu-flash flash (T,U -> P).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, U, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, U, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn tu_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    T: f64,
    U: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyTuFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::tu_flash(&mixture, &ideal_gas, kelvins(T), joules_per_mole(U), &z)
        .map(|r| PyTuFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pu-flash flash (P,U -> T).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, U, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, U, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn pu_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    P: f64,
    U: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPuFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::pu_flash(&mixture, &ideal_gas, pascals(P), joules_per_mole(U), &z)
        .map(|r| PyPuFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pressure/reflux-ratio flash of a mixture (P,ratio -> T).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, P, reflux, phase, temperature, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, P, reflux, phase, temperature, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)] // the signature is the flash's inputs
pub fn pv_reflux_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    P: f64,
    reflux: f64,
    phase: &str,
    temperature: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPvRefluxFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let which = match phase {
        "liquid" => azoth_eos::pv_reflux_flash::RefluxPhase::Liquid,
        _ => azoth_eos::pv_reflux_flash::RefluxPhase::Vapour,
    };
    azoth_eos::pv_reflux_flash::pv_reflux_flash(
        &mixture,
        pascals(P),
        reflux,
        which,
        kelvins(temperature),
        &z,
    )
    .map(|r| PyPvRefluxFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The pressure/vapour-fraction flash of a mixture (P,beta -> T).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, P, beta, temperature, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, P, beta, temperature, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)] // the signature is the flash's inputs
pub fn pvf_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    P: f64,
    beta: f64,
    temperature: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPvfFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::pvf_flash::pvf_flash(&mixture, pascals(P), beta, kelvins(temperature), &z)
        .map(|r| PyPvfFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The temperature and vapour-volume-fraction flash of a mixture (T,fraction -> P).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, T, fraction, P, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, T, fraction, P, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)] // the signature is the flash's inputs
pub fn tv_fraction_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    T: f64,
    fraction: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyTvFractionFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::tv_fraction_flash::tv_fraction_flash(&mixture, kelvins(T), fraction, pascals(P), &z)
        .map(|r| PyTvFractionFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The volume-internal-energy flash of a mixture (V,U -> P,T).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, V, U, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, V, U, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn vu_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    V: f64,
    U: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyVuFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::vu_flash(
        &mixture,
        &ideal_gas,
        cubic_meters_per_mole(V),
        joules_per_mole(U),
        &z,
    )
    .map(|r| PyVuFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The volume-enthalpy flash of a mixture (V,H -> P,T).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, V, H, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, V, H, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn vh_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    V: f64,
    H: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyVhFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::vh_flash::vh_flash(
        &mixture,
        &ideal_gas,
        cubic_meters_per_mole(V),
        joules_per_mole(H),
        &z,
    )
    .map(|r| PyVhFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The volume-entropy flash of a mixture (V,S -> P,T).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, V, S, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, V, S, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn vs_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    V: f64,
    S: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyVsFlashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::vs_flash::vs_flash(
        &mixture,
        &ideal_gas,
        cubic_meters_per_mole(V),
        joules_per_mole_kelvin(S),
        &z,
    )
    .map(|r| PyVsFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The volume-internal-energy state of a pure component (P,V,U -> T,beta).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, V, U,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, P, V, U, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn vu_flash_single_comp(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    P: f64,
    V: f64,
    U: f64,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyVuFlashSingleCompResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::vu_flash_single_comp::vu_flash_single_comp(
        &mixture,
        &ideal_gas,
        pascals(P),
        cubic_meters_per_mole(V),
        joules_per_mole(U),
    )
    .map(|r| PyVuFlashSingleCompResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// Whether a feed at a temperature and pressure is stable as a single phase.
///
/// The same seven arguments as the flash, because it is the same state asked a
/// different question, and the two results are meant to be read together: a feed the
/// flash splits is `unstable` here, and one it reports as a single phase is `stable`.
///
/// `verdict` crosses as the spec's spelling and the bridge rebuilds the enum, the
/// same arrangement `phase` uses. `tm` and `w` are one entry per trial, always two,
/// in a fixed order - vapour-like first - so a caller reads `tm[0]` as the vapour-like
/// trial rather than having to look it up.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, T, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn stability_test(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyStabilityTestResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::stability_test(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyStabilityTestResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// How many phases a feed splits into at a temperature and pressure, and how much of each.
///
/// The two-phase flash with the stability seeding: a tangent-plane trial whose stationary
/// point is at a composition that is not already a phase is a phase the flash has not found,
/// and it is added before the fractions are solved. The same seven arguments as
/// `stability_test`, and `min_t_over_tc` is the smaller of the two models' - whichever of
/// them is nearer a critical point is the one the answer is least trustworthy at.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, T, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn tp_multiflash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyTpMultiflashResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::tp_multiflash(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyTpMultiflashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pressure at which a liquid of composition `x` first gives off vapour.
///
/// The first of the two phase-boundary models. They take the same six arguments as
/// the flash - three per-component vectors, a flattened interaction matrix,
/// temperature and a composition - and differ only in which composition it is,
/// which is why they are two functions rather than one with a switch.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, T, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, T, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn bubble_pressure(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    T: f64,
    held: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPhaseBoundaryResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::bubble_pressure(&mixture, kelvins(T), &held)
        .map(|r| PyPhaseBoundaryResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The critical point of a mixture of composition `z`.
///
/// The one model here that solves for a *state* rather than for a phase split, and the
/// only one whose answer a caller cannot check against a phase they can see. It takes
/// the same six arguments as the boundary models and returns four state variables
/// instead of one pressure.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `z` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn critical_point(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyCriticalPointResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::critical_point(&mixture, &z)
        .map(|r| PyCriticalPointResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pressure at which a vapour of composition `held` first condenses.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, T, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, T, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn dew_pressure(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    T: f64,
    held: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPhaseBoundaryResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::dew_pressure(&mixture, kelvins(T), &held)
        .map(|r| PyPhaseBoundaryResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The temperature at which a liquid of composition `held` first gives off vapour.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, P, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, P, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn bubble_temperature(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    P: f64,
    held: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPhaseBoundaryTemperatureResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::bubble_temperature(&mixture, pascals(P), &held)
        .map(|r| PyPhaseBoundaryTemperatureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The dew point of a vapour held in a pore, with the Kelvin shift that curvature puts on it.
///
/// The same boundary as `dew_temperature` on a curved interface: Young-Laplace puts the liquid
/// at `2 sigma cos(theta) / r` above the vapour, and the K-values shift by the Kelvin equation.
/// `pore_radius` is in metres, `contact_angle` in radians and `surface_tension` in N/m - the
/// tension is an argument here and an interphase property upstream, which is in the spec.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, P, held, pore_radius, contact_angle, surface_tension, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, P, held, pore_radius, contact_angle, surface_tension, eos = \"pr\", alpha = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn capillary_dew_point(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    P: f64,
    held: Vec<f64>,
    pore_radius: f64,
    contact_angle: f64,
    surface_tension: f64,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyCapillaryDewPointResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::capillary_dew_point(
        &mixture,
        pascals(P),
        &held,
        pore_radius,
        contact_angle,
        surface_tension,
    )
    .map(|r| PyCapillaryDewPointResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The temperature at which a vapour of composition `held` first gives off liquid.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, P, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, association, P, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn dew_temperature(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    P: f64,
    held: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyPhaseBoundaryTemperatureResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::dew_temperature(&mixture, pascals(P), &held)
        .map(|r| PyPhaseBoundaryTemperatureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The ideal-gas heat capacity from a four-term polynomial.
///
/// The first calc in this namespace whose constants are fitted data rather than
/// coefficients of a published equation, and the first whose arguments include a
/// caller-supplied coefficient set at all. All four coefficients cross as plain
/// floats: they are dimensionless by construction, because the polynomial is written
/// against `T/(1000 K)` and divided through by `R`.
#[pyfunction]
#[pyo3(signature = (cp_a, cp_b, cp_c, cp_d, cp_e, T))]
#[pyo3(text_signature = "(cp_a, cp_b, cp_c, cp_d, cp_e, T)")]
#[allow(non_snake_case)] // `T` is the symbol in the chemistry
pub fn ideal_gas_cp(
    py: Python<'_>,
    cp_a: f64,
    cp_b: f64,
    cp_c: f64,
    cp_d: f64,
    cp_e: f64,
    T: f64,
) -> PyResult<PyIdealGasCpResult> {
    azoth_eos::ideal_gas_cp(cp_a, cp_b, cp_c, cp_d, cp_e, kelvins(T))
        .map(|r| PyIdealGasCpResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The absolute molar enthalpy and entropy of a mixture at a state.
///
/// A `direct` model - vectors in, values out, no iteration - which is why it has no
/// algorithm block in its spec and why its arguments are the one place in this
/// namespace where sixteen of them arrive at once. They are flattened here rather
/// than passed as an object because the boundary carries numbers, and the Python side
/// is where the caller's `IdealGasModel` is unpacked.
#[pyfunction]
#[pyo3(signature = (
    Tc, Pc, omega, kij, association, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z,
    compressibility, eos = "pr", alpha = "pr", alpha_params = None
))]
#[pyo3(text_signature = "(
    Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z,
    compressibility, eos = \"pr\", alpha = \"pr\"
)")]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn molar_enthalpy_entropy(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    compressibility: f64,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyMolarEnthalpyEntropyResult> {
    let mixture = build_mixture(
        py,
        &Tc,
        &Pc,
        &omega,
        kij,
        &association,
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::molar_enthalpy_entropy(
        &mixture,
        &ideal_gas,
        kelvins(T),
        pascals(P),
        &z,
        compressibility,
    )
    .map(|r| PyMolarEnthalpyEntropyResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The BWRS (MBWR-32) phase state, computed in Rust.
///
/// `a` is the 32 coefficients per component, flattened component-major; `rhoc` the
/// per-component critical density. The bridge unpacks the named coefficients into
/// these two vectors the same way it does the cubic's `Tc`, `Pc` and `omega`.
#[pyfunction]
#[pyo3(signature = (a, rhoc, T, P, z))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn bwrs_phase(
    py: Python<'_>,
    a: Vec<f64>,
    rhoc: Vec<f64>,
    T: f64,
    P: f64,
    z: Vec<f64>,
) -> PyResult<PyBwrsPhaseResult> {
    let n = rhoc.len();
    let mut coeffs = Vec::with_capacity(n);
    for (i, &r) in rhoc.iter().enumerate() {
        let mut arr = [0.0; 32];
        arr.copy_from_slice(&a[i * 32..(i + 1) * 32]);
        coeffs.push(azoth_eos::bwrs::BwrsCoefficients { a: arr, rhoc: r });
    }
    azoth_eos::bwrs_phase(&coeffs, kelvins(T), pascals(P), &z)
        .map(|r| PyBwrsPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The ammonia reference phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (T, P))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn ammonia_phase(py: Python<'_>, T: f64, P: f64) -> PyResult<PyAmmoniaPhaseResult> {
    azoth_eos::ammonia_phase(kelvins(T), pascals(P))
        .map(|r| PyAmmoniaPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Span-Wagner CO2 phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (T, P))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn co2_phase(py: Python<'_>, T: f64, P: f64) -> PyResult<PyCo2PhaseResult> {
    azoth_eos::co2_phase(kelvins(T), pascals(P))
        .map(|r| PyCo2PhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Vega helium phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (T, P))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn helium_phase(py: Python<'_>, T: f64, P: f64) -> PyResult<PyHeliumPhaseResult> {
    azoth_eos::helium_phase(kelvins(T), pascals(P))
        .map(|r| PyHeliumPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The fixed-role gas-oil-brine flash, computed in Rust.
///
/// The component **names** cross unresolved: the two EoS roles' constants, the seeding's
/// classes and the brine's ion mask all resolve here from the same databank, so the two
/// languages cannot disagree about which substance is which.
#[pyfunction]
#[pyo3(signature = (components, cubic, T, P, moles))]
#[pyo3(text_signature = "(components, cubic, T, P, moles)")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn hybrid_eos_ge_flash(
    py: Python<'_>,
    components: Vec<String>,
    cubic: &str,
    T: f64,
    P: f64,
    moles: Vec<f64>,
) -> PyResult<PyHybridEosGeFlashResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    azoth_eos::hybrid_eos_ge_flash(
        &names,
        cubic
            .parse()
            .map_err(pyo3::exceptions::PyValueError::new_err)?,
        T,
        P,
        &moles,
    )
    .map(|r| PyHybridEosGeFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The Henry constant of a gas in water, computed in Rust.
///
/// `gas` crosses and is resolved by name, so the two backends accept the same vocabulary:
/// the reference twin resolves it the same way, and a caller who has `methane` rather than
/// `ch4` gets one answer and not two. The refusal is built here rather than left as a
/// `ValueError` from the enum's own parser, because an unknown gas is the reference's
/// `InvalidInputError` and the two backends must raise the same class.
#[pyfunction]
#[pyo3(signature = (gas, T))]
#[allow(non_snake_case)] // `T` is the symbol in the guideline
pub fn iapws_henry_law(py: Python<'_>, gas: &str, T: f64) -> PyResult<PyIapwsHenryLawResult> {
    let row = azoth_eos::gas_by_name(gas).map_err(|e| to_pyerr(py, e))?;
    azoth_eos::iapws_henry_law(row, kelvins(T))
        .map(|r| PyIapwsHenryLawResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Leachman hydrogen phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (T, P, hydrogen_type = "normal", compressed_phase = "vapour"))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn hydrogen_phase(
    py: Python<'_>,
    T: f64,
    P: f64,
    hydrogen_type: &str,
    compressed_phase: &str,
) -> PyResult<PyHydrogenPhaseResult> {
    azoth_eos::hydrogen_phase(kelvins(T), pascals(P), hydrogen_type, compressed_phase)
        .map(|r| PyHydrogenPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The hydrate formation temperature of a fluid, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the hydrate's guest
/// tables are keyed by name, so a mixture built from constants alone could not carry them.
#[pyfunction]
#[pyo3(signature = (components, P, z, eos = "srk", hydrate_model = "pvtsim"))]
#[allow(non_snake_case)] // `P` is the symbol in the chemistry
pub fn hydrate_formation_temperature(
    py: Python<'_>,
    components: Vec<String>,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    hydrate_model: &str,
) -> PyResult<crate::results::PyHydrateFormationTemperatureResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = azoth_eos::hydrate::hydrate_mixture_of(
        &names,
        eos.parse().unwrap_or(azoth_eos::Cubic::Srk),
        None,
        hydrate_model
            .parse()
            .unwrap_or(azoth_eos::hydrate::HydrateModel::Pvtsim),
    )
    .map_err(|e| to_pyerr(py, e))?;
    azoth_eos::hydrate_formation_temperature(&mixture, pascals(P), &z)
        .map(|r| crate::results::PyHydrateFormationTemperatureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The hydrate formation pressure of a fluid at a temperature, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the hydrate's guest
/// tables are keyed by name, so a mixture built from constants alone could not carry them.
#[pyfunction]
#[pyo3(signature = (components, T, z, eos = "srk", hydrate_model = "pvtsim"))]
#[allow(non_snake_case)] // `T` is the symbol in the chemistry
pub fn hydrate_formation_pressure(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    z: Vec<f64>,
    eos: &str,
    hydrate_model: &str,
) -> PyResult<crate::results::PyHydrateFormationPressureResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = azoth_eos::hydrate::hydrate_mixture_of(
        &names,
        eos.parse().unwrap_or(azoth_eos::Cubic::Srk),
        None,
        hydrate_model
            .parse()
            .unwrap_or(azoth_eos::hydrate::HydrateModel::Pvtsim),
    )
    .map_err(|e| to_pyerr(py, e))?;
    azoth_eos::hydrate_formation_pressure(&mixture, kelvins(T), &z)
        .map(|r| crate::results::PyHydrateFormationPressureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// A hydrate curve over a pressure grid, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the hydrate's guest
/// tables are keyed by name, so a mixture built from constants alone could not carry them.
#[pyfunction]
#[pyo3(signature = (components, P_min, P_max, z, eos = "srk", hydrate_model = "pvtsim"))]
#[allow(non_snake_case)] // `P_min` and `P_max` are the bounds in the chemistry
pub fn hydrate_equilibrium_line(
    py: Python<'_>,
    components: Vec<String>,
    P_min: f64,
    P_max: f64,
    z: Vec<f64>,
    eos: &str,
    hydrate_model: &str,
) -> PyResult<crate::results::PyHydrateEquilibriumLineResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = azoth_eos::hydrate::hydrate_mixture_of(
        &names,
        eos.parse().unwrap_or(azoth_eos::Cubic::Srk),
        None,
        hydrate_model
            .parse()
            .unwrap_or(azoth_eos::hydrate::HydrateModel::Pvtsim),
    )
    .map_err(|e| to_pyerr(py, e))?;
    azoth_eos::hydrate_equilibrium_line(&mixture, pascals(P_min), pascals(P_max), &z)
        .map(|r| crate::results::PyHydrateEquilibriumLineResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The inhibitor moles that hold a hydrate temperature down to a target, computed in Rust.
///
/// **The feed crosses in moles**, which is the one hydrate model here that does: the secant
/// adds an absolute amount to the inhibitor's entry, so a normalised feed would reproduce the
/// equation and not the path.
#[pyfunction]
#[pyo3(signature = (components, moles, inhibitor, T_target, P, eos = "srk", hydrate_model = "pvtsim"))]
#[allow(non_snake_case)] // `T_target`, `P` and `eos` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // the names, the feed, the inhibitor and the state
pub fn hydrate_inhibitor_concentration(
    py: Python<'_>,
    components: Vec<String>,
    moles: Vec<f64>,
    inhibitor: &str,
    T_target: f64,
    P: f64,
    eos: &str,
    hydrate_model: &str,
) -> PyResult<crate::results::PyHydrateInhibitorConcentrationResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = azoth_eos::hydrate::hydrate_mixture_of(
        &names,
        eos.parse().unwrap_or(azoth_eos::Cubic::Srk),
        None,
        hydrate_model
            .parse()
            .unwrap_or(azoth_eos::hydrate::HydrateModel::Pvtsim),
    )
    .map_err(|e| to_pyerr(py, e))?;
    azoth_eos::hydrate_inhibitor_concentration(
        &mixture,
        inhibitor,
        &moles,
        kelvins(T_target),
        pascals(P),
    )
    .map(|r| crate::results::PyHydrateInhibitorConcentrationResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The inhibitor dose that reaches a target aqueous mass fraction, computed in Rust.
///
/// **The feed crosses in moles**, as for its sibling: the secant adds an absolute amount to
/// the inhibitor's entry.
#[pyfunction]
#[pyo3(signature = (components, moles, inhibitor, wt_target, T, P, eos = "srk"))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // the names, the feed, the inhibitor and the state
pub fn hydrate_inhibitor_wt(
    py: Python<'_>,
    components: Vec<String>,
    moles: Vec<f64>,
    inhibitor: &str,
    wt_target: f64,
    T: f64,
    P: f64,
    eos: &str,
) -> PyResult<crate::results::PyHydrateInhibitorWtResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) =
        azoth_eos::databank::mixture_of(&names, eos.parse().unwrap_or(azoth_eos::Cubic::Srk), None)
            .map_err(|e| to_pyerr(py, e))?;
    azoth_eos::hydrate_inhibitor_wt(
        &mixture,
        inhibitor,
        &moles,
        wt_target,
        kelvins(T),
        pascals(P),
    )
    .map(|r| crate::results::PyHydrateInhibitorWtResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The fraction of a feed that is hydrate at a state, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the hydrate's guest
/// tables are keyed by name, so a mixture built from constants alone could not carry them.
#[pyfunction]
#[pyo3(signature = (components, T, P, z, eos = "srk", hydrate_model = "pvtsim"))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn hydrate_fraction(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    hydrate_model: &str,
) -> PyResult<crate::results::PyHydrateFractionResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) = azoth_eos::hydrate::hydrate_mixture_of(
        &names,
        eos.parse().unwrap_or(azoth_eos::Cubic::Srk),
        None,
        hydrate_model
            .parse()
            .unwrap_or(azoth_eos::hydrate::HydrateModel::Pvtsim),
    )
    .map_err(|e| to_pyerr(py, e))?;
    azoth_eos::hydrate_fraction(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| crate::results::PyHydrateFractionResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// A TBP pseudo-component's properties, computed in Rust.
#[pyfunction]
#[pyo3(signature = (molar_mass, density))]
pub fn tbp_fraction_properties(
    py: Python<'_>,
    molar_mass: f64,
    density: f64,
) -> PyResult<crate::results::PyTbpFractionPropertiesResult> {
    azoth_eos::tbp_fraction_properties(molar_mass, density)
        .map(|r| crate::results::PyTbpFractionPropertiesResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// A wax cut's solid fugacity coefficient, computed in Rust.
#[pyfunction]
#[pyo3(signature = (molar_mass, tc, pc, omega, heat_of_fusion, triple_point_temperature, T, P, eos = "srk"))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)]
pub fn wax_solid_fugacity(
    py: Python<'_>,
    molar_mass: f64,
    tc: f64,
    pc: f64,
    omega: f64,
    heat_of_fusion: f64,
    triple_point_temperature: f64,
    T: f64,
    P: f64,
    eos: &str,
) -> PyResult<crate::results::PyWaxSolidFugacityResult> {
    azoth_eos::wax_solid_fugacity(
        molar_mass,
        kelvins(tc),
        pascals(pc),
        omega,
        heat_of_fusion,
        kelvins(triple_point_temperature),
        kelvins(T),
        pascals(P),
        eos,
    )
    .map(|r| crate::results::PyWaxSolidFugacityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The wax fraction of a feed at a state, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the wax flag and
/// the melt data are the databank's own columns.
#[pyfunction]
#[pyo3(signature = (components, T, P, z, eos = "srk"))]
#[allow(non_snake_case)] // `T`, `P` and `z` are the symbols in the chemistry
pub fn tp_multiflash_wax(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
) -> PyResult<crate::results::PyTpMultiflashWaxResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let (mixture, _) =
        azoth_eos::databank::mixture_of(&names, eos.parse().unwrap_or(azoth_eos::Cubic::Srk), None)
            .map_err(|e| to_pyerr(py, e))?;
    azoth_eos::tp_multiflash_wax(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| crate::results::PyTpMultiflashWaxResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The fraction of a feed that has frozen out as one pure solid, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the melt data and
/// the density correlations are the databank's own columns.
#[pyfunction]
#[pyo3(signature = (components, solid, T, P, z, eos = "srk"))]
#[allow(non_snake_case)] // `T`, `P` and `z` are the symbols in the chemistry
pub fn tp_solid_flash(
    py: Python<'_>,
    components: Vec<String>,
    solid: &str,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
) -> PyResult<crate::results::PyTpSolidFlashResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    azoth_eos::tp_solid_flash(&names, solid, kelvins(T), pascals(P), &z, eos)
        .map(|r| crate::results::PyTpSolidFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// One salt's saturation ratio in a brine, computed in Rust.
#[pyfunction]
#[pyo3(signature = (salt, x1, x2, x_water, gamma1, gamma2, water_activity, T, P, h3o_molality))]
#[allow(clippy::too_many_arguments)]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn scale_saturation_ratio(
    py: Python<'_>,
    salt: &str,
    x1: f64,
    x2: f64,
    x_water: f64,
    gamma1: f64,
    gamma2: f64,
    water_activity: f64,
    T: f64,
    P: f64,
    h3o_molality: Option<f64>,
) -> PyResult<crate::results::PyScaleSaturationRatioResult> {
    azoth_eos::scale_saturation_ratio(
        salt,
        x1,
        x2,
        x_water,
        gamma1,
        gamma2,
        water_activity,
        h3o_molality,
        kelvins(T),
        pascals(P),
    )
    .map(|r| crate::results::PyScaleSaturationRatioResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The solid one mineral takes from a brine, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the brine is an
/// electrolyte and its coefficients come from the Pitzer phase over the same names.
#[pyfunction]
#[pyo3(signature = (components, z, salt, T, P))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn salt_precipitation(
    py: Python<'_>,
    components: Vec<String>,
    z: Vec<f64>,
    salt: &str,
    T: f64,
    P: f64,
) -> PyResult<crate::results::PySaltPrecipitationResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    azoth_eos::salt_precipitation(&names, salt, kelvins(T), pascals(P), &z)
        .map(|r| crate::results::PySaltPrecipitationResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// A pure solid's fugacity coefficient, computed in Rust.
#[pyfunction]
#[pyo3(signature = (heat_of_fusion, triple_point_temperature, delta_cp_sl, delta_solid_volume, tc, pc, omega, T, P, eos = "srk"))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)]
pub fn solid_fugacity(
    py: Python<'_>,
    heat_of_fusion: f64,
    triple_point_temperature: f64,
    delta_cp_sl: f64,
    delta_solid_volume: f64,
    tc: f64,
    pc: f64,
    omega: f64,
    T: f64,
    P: f64,
    eos: &str,
) -> PyResult<crate::results::PySolidFugacityResult> {
    azoth_eos::solid_fugacity(
        heat_of_fusion,
        triple_point_temperature,
        delta_cp_sl,
        delta_solid_volume,
        kelvins(tc),
        pascals(pc),
        omega,
        kelvins(T),
        pascals(P),
        eos,
    )
    .map(|r| crate::results::PySolidFugacityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The gas binary diffusivity, from the Chapman-Enskog theory.
#[pyfunction]
#[pyo3(signature = (MA, MB, sigma, eps, T, P))]
#[pyo3(text_signature = "(MA, MB, sigma, eps, T, P)")]
#[allow(non_snake_case)] // `MA`, `MB` and `T` are the symbols in the published equation
pub fn chapman_enskog_diffusivity(
    py: Python<'_>,
    MA: f64,
    MB: f64,
    sigma: f64,
    eps: f64,
    T: f64,
    P: f64,
) -> PyResult<crate::results::PyChapmanEnskogDiffusivityResult> {
    azoth_eos::chapman_enskog_diffusivity(
        kilograms_per_mole(MA),
        kilograms_per_mole(MB),
        sigma,
        kelvins(eps),
        kelvins(T),
        pascals(P),
    )
    .map(|r| crate::results::PyChapmanEnskogDiffusivityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The gas binary diffusivity, from the Fuller-Schettler-Giddings correlation.
#[pyfunction]
#[pyo3(signature = (MA, MB, VA, VB, T, P))]
#[pyo3(text_signature = "(MA, MB, VA, VB, T, P)")]
#[allow(non_snake_case)] // `MA`, `MB`, `VA`, `VB`, `T` and `P` are the equation's symbols
pub fn fuller_schettler_giddings_diffusivity(
    py: Python<'_>,
    MA: f64,
    MB: f64,
    VA: f64,
    VB: f64,
    T: f64,
    P: f64,
) -> PyResult<crate::results::PyFullerSchettlerGiddingsDiffusivityResult> {
    azoth_eos::fuller_schettler_giddings_diffusivity(
        kilograms_per_mole(MA),
        kilograms_per_mole(MB),
        cubic_meters_per_mole(VA),
        cubic_meters_per_mole(VB),
        kelvins(T),
        pascals(P),
    )
    .map(|r| crate::results::PyFullerSchettlerGiddingsDiffusivityResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// A fluid's freezing-point temperature, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves them: the tabulated solid
/// reads the databank's melting point, heat of fusion and density correlations, and the
/// para-hydrogen route reads the reference equation this crate carries.
#[pyfunction]
#[pyo3(signature = (components, z, solid, P))]
#[allow(non_snake_case)] // `P` is the symbol in the chemistry
pub fn freezing_point(
    py: Python<'_>,
    components: Vec<String>,
    z: Vec<f64>,
    solid: String,
    P: f64,
) -> PyResult<crate::results::PyFreezingPointResult> {
    azoth_eos::freezing_point(&components, &z, &solid, pascals(P))
        .map(|r| crate::results::PyFreezingPointResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The IAPWS-IF97 water phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (T, P))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn water_phase(py: Python<'_>, T: f64, P: f64) -> PyResult<PyWaterPhaseResult> {
    azoth_eos::water_phase(kelvins(T), pascals(P))
        .map(|r| PyWaterPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The solid argon phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (T, P))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn argon_solid_phase(py: Python<'_>, T: f64, P: f64) -> PyResult<PyArgonSolidPhaseResult> {
    azoth_eos::argon_solid_phase(kelvins(T), pascals(P))
        .map(|r| PyArgonSolidPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The solid para-hydrogen phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (T, P))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn parahydrogen_solid_phase(
    py: Python<'_>,
    T: f64,
    P: f64,
) -> PyResult<PyParahydrogenSolidPhaseResult> {
    azoth_eos::parahydrogen_solid_phase(kelvins(T), pascals(P))
        .map(|r| PyParahydrogenSolidPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The EOS-CG phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (components, T, P, z))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn eos_cg_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
) -> PyResult<PyEosCgPhaseResult> {
    azoth_eos::eos_cg_phase(&components, kelvins(T), pascals(P), &z)
        .map(|r| PyEosCgPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The GERG-2008 phase state, computed in Rust.
#[pyfunction]
#[pyo3(signature = (components, T, P, z))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn gerg2008_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
) -> PyResult<PyGerg2008PhaseResult> {
    azoth_eos::gerg2008_phase(&components, kelvins(T), pascals(P), &z)
        .map(|r| PyGerg2008PhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The liquid viscosity from the Pedersen (PFCT) heavy-oil correlation.
///
/// `eos`, `alpha` and `alpha_params` are accepted for the boundary's uniformity with the
/// other mixture models and ignored: the reference flash is pure methane SRK, so the
/// mixture's own cubic and alpha are not consulted.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, molar_mass, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, molar_mass, T, P, z, eos = \"pr\", alpha = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
#[allow(unused_variables)] // `eos`, `alpha` and `alpha_params` are boundary-only, see above.
pub fn viscosity(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    molar_mass: Vec<f64>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyViscosityResult> {
    let n = Tc.len();
    let components = (0..n)
        .map(|i| {
            azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i])
                .map(|c| c.with_molar_mass(Some(molar_mass[i])))
        })
        .collect::<azoth_core::Result<Vec<_>>>()
        .map_err(|e| to_pyerr(py, e))?;
    let mixture = azoth_eos::Mixture::new(components, kij).map_err(|e| to_pyerr(py, e))?;
    azoth_eos::viscosity(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyViscosityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The liquid viscosity NeqSim gives an **aqueous** phase.
///
/// `eos`, `alpha` and `alpha_params` are accepted for the boundary's uniformity and
/// ignored: this correlation reads the components' own `LIQVISC` sets and nothing about
/// the mixture's cubic.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, molar_mass, liqvisc, liqvisc_model, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, molar_mass, liqvisc, liqvisc_model, T, P, z, eos = \"pr\", alpha = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
#[allow(unused_variables)] // `eos`, `alpha` and `alpha_params` are boundary-only, see above.
pub fn aqueous_viscosity(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    molar_mass: Vec<f64>,
    liqvisc: Vec<f64>,
    liqvisc_model: Vec<u32>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyAqueousViscosityResult> {
    let n = Tc.len();
    let components = (0..n)
        .map(|i| {
            azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i]).map(|c| {
                c.with_molar_mass(Some(molar_mass[i]))
                    .with_liquid_viscosity(
                        [
                            liqvisc[4 * i],
                            liqvisc[4 * i + 1],
                            liqvisc[4 * i + 2],
                            liqvisc[4 * i + 3],
                        ],
                        liqvisc_model[i],
                    )
            })
        })
        .collect::<azoth_core::Result<Vec<_>>>()
        .map_err(|e| to_pyerr(py, e))?;
    let mixture = azoth_eos::Mixture::new(components, kij).map_err(|e| to_pyerr(py, e))?;
    azoth_eos::aqueous_viscosity(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyAqueousViscosityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The liquid thermal conductivity from the Pedersen (PFCT) correlation.
///
/// `eos`, `alpha` and `alpha_params` are accepted for the boundary's uniformity and
/// ignored: the reference flash is pure methane SRK.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, association, molar_mass, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, association, molar_mass, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z, eos = \"pr\", alpha = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
#[allow(unused_variables)] // `eos`, `alpha` and `alpha_params` are boundary-only, see above.
pub fn thermal_conductivity(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    association: PyRef<'_, PyAssociationSpec>,
    molar_mass: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    cp_e: Vec<f64>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    eos: &str,
    alpha: &str,
    alpha_params: Option<Vec<Vec<f64>>>,
) -> PyResult<PyThermalConductivityResult> {
    let n = Tc.len();
    let components = (0..n)
        .map(|i| {
            azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i])
                .map(|c| c.with_molar_mass(Some(molar_mass[i])))
        })
        .collect::<azoth_core::Result<Vec<_>>>()
        .map_err(|e| to_pyerr(py, e))?;
    let mixture = azoth_eos::Mixture::new(components, kij).map_err(|e| to_pyerr(py, e))?;
    let ideal_gas = azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        cp_e,
    };
    azoth_eos::thermal_conductivity(&mixture, &ideal_gas, kelvins(T), pascals(P), &z)
        .map(|r| PyThermalConductivityResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// Every model in the workspace, across every namespace that has one.
///
/// The table is generated into the crate that owns its namespace, so a second
/// namespace means a second table chained in here. There are three: `eos` and
/// `reactions`, and `process`, which is the unit-operation tier and is what P11 builds.
///
/// **This is the function that has to change when a namespace is added**, and a missing
/// chain is not self-announcing: `process.pump` landed without one, and the contract test
/// that would have said so is `@pytest.mark.requires_rust`, so a Python-backend run
/// skipped it. A namespace is added here in the same commit as its first model.
fn all_models() -> impl Iterator<Item = &'static azoth_core::ModelSpec> {
    azoth_eos::model_gen::models()
        .iter()
        .chain(azoth_reactions::model_gen::models().iter())
        .chain(azoth_standards::model_gen::models().iter())
        .chain(azoth_process::model_gen::models().iter())
        .copied()
}

/// Every model id this extension implements.
///
/// Separate from `calc_ids()` on purpose: the calc registry's id list is asserted to
/// be *exactly* the specs under `specs/calcs/`, so a model appearing there would
/// break that contract rather than extend it. Models have their own list, generated
/// from their own tree.
#[pyfunction]
#[must_use]
pub fn model_ids() -> Vec<String> {
    all_models().map(|m| m.id.to_string()).collect()
}

/// What a model's spec fixes, by model id: `procedure` or `direct`.
///
/// Held to the specs by the same contract test that holds `model_schemes`: a `kind`
/// field nothing reads would be a spec field nothing checks, which is the defect this
/// project is organised against. A `direct` model has no scheme, so this is the one
/// thing that distinguishes it before its algorithm is looked for.
#[pyfunction]
#[must_use]
pub fn model_kind(model_id: &str) -> String {
    match all_models().find(|spec| spec.id == model_id) {
        Some(spec) => spec.kind.to_string(),
        None => String::new(),
    }
}

/// The algorithm scheme each model runs, by model id.
///
/// The third leg of the same contract `solver_kinds()` provides for solvers: the
/// schema's `algorithm.scheme`, the implementations' own names, and this list must
/// agree, and `test_model_contract.py` asserts it rather than trusting three
/// hand-edited lists to stay in step.
#[pyfunction]
#[must_use]
pub fn model_schemes(model_id: &str) -> Vec<String> {
    match all_models().find(|spec| spec.id == model_id) {
        // A `direct` model has no scheme, and an empty list is the honest answer:
        // there is no procedure whose name could be compared, and the contract test
        // asserts the emptiness rather than papering over it with a placeholder.
        Some(spec) => spec
            .algorithm
            .map(|a| vec![a.scheme.to_string()])
            .unwrap_or_default(),
        None => Vec::new(),
    }
}

/// The SAFT-VR-Mie flash, computed in Rust.
///
/// **The component names cross unresolved**, and this side looks them up in the Rust
/// databank - the Mie set and the cubic constants the Wilson seed is built from.
#[pyfunction]
#[pyo3(signature = (components, T, P, z))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn tp_flash_saft(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
) -> PyResult<crate::results::PyTpFlashSaftResult> {
    azoth_eos::tp_flash_saft::tp_flash_saft(&components, kelvins(T), pascals(P), &z)
        .map(|r| crate::results::PyTpFlashSaftResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The SAFT-VR-Mie phase state, computed in Rust.
///
/// **The component names cross unresolved**, and this side looks them up in the Rust
/// databank - the five Mie columns, whose absence the table spells as zeros in three of
/// them and not in the exponents.
#[pyfunction]
#[pyo3(signature = (components, T, P, z, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn saft_vr_mie_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    compressed_phase: &str,
) -> PyResult<crate::results::PySaftVrMiePhaseResult> {
    azoth_eos::saft_vr_mie_phase::saft_vr_mie_phase(
        &components,
        kelvins(T),
        pascals(P),
        &z,
        compressed_phase,
    )
    .map(|r| crate::results::PySaftVrMiePhaseResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The PC-SAFT phase state, computed in Rust.
///
/// **The component names cross unresolved**, and this side looks them up in the Rust
/// databank - which for this model means the `mSAFT`/`sigmaSAFT`/`epsikSAFT` set, whose
/// absence the table spells as zeros, and the `KIJPCSAFT` column.
#[pyfunction]
#[pyo3(signature = (components, T, P, z, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn pcsaft_rahmat_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    compressed_phase: &str,
) -> PyResult<crate::results::PyPcsaftRahmatPhaseResult> {
    azoth_eos::pcsaft_rahmat_phase(&components, kelvins(T), pascals(P), &z, compressed_phase)
        .map(|r| crate::results::PyPcsaftRahmatPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The PR-CPA phase state, computed in Rust.
///
/// **The component names cross unresolved**, as they do for the SRK twin, and this side
/// looks them up in the Rust databank. What that buys here is the family: `Cubic::Pr`
/// selects the fitted `aCPA_PR`/`bCPA_PR`/`mCPA_PR` set and the `cpakij_PR` column, both of
/// which this is the first shipped model to read.
#[pyfunction]
#[pyo3(signature = (components, T, P, z, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn pr_cpa_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    compressed_phase: &str,
) -> PyResult<crate::results::PyPrCpaPhaseResult> {
    azoth_eos::pr_cpa_phase(&components, kelvins(T), pascals(P), &z, compressed_phase)
        .map(|r| crate::results::PyPrCpaPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The UMR-CPA phase state, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves the UMR-CPA parameter
/// set, the `UMRCPA_MC1..5` coefficients and the `UNIFACcompUMRPRU` group decomposition
/// itself. What that buys here is everything the model is: a third fitted set, a
/// five-parameter alpha and a universal mixing rule, none of which a caller could state as
/// numbers without reimplementing the model.
#[pyfunction]
#[pyo3(signature = (components, T, P, z, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn umr_cpa_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    compressed_phase: &str,
) -> PyResult<crate::results::PyUmrCpaPhaseResult> {
    azoth_eos::umr_cpa_phase(&components, kelvins(T), pascals(P), &z, compressed_phase)
        .map(|r| crate::results::PyUmrCpaPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The Soreide-Whitson phase state, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves each name's role in the
/// aqueous correlation and reads its row of `KIJWhitsonSoriede`. Both are `name` lookups on
/// NeqSim's side too, so a caller could not state them as numbers without reimplementing the
/// model - and the roles are what decide which component gets the salinity-dependent alpha.
#[pyfunction]
#[pyo3(signature = (components, T, P, x, salinity, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn soreide_whitson_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    x: Vec<f64>,
    salinity: f64,
    compressed_phase: &str,
) -> PyResult<crate::results::PySoreideWhitsonPhaseResult> {
    azoth_eos::soreide_whitson_phase(
        &components,
        kelvins(T),
        pascals(P),
        &x,
        salinity,
        compressed_phase,
    )
    .map(|r| crate::results::PySoreideWhitsonPhaseResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The Fürst electrolyte phase state, computed in Rust.
///
/// **The component names cross unresolved**, and this side resolves each one's dielectric
/// coefficients, fitted covolume, Schwartzentruber parameters and short-range pair table.
/// The salt is a component and not a scalar, so an ion crosses as a name like any other.
#[pyfunction]
#[pyo3(signature = (components, T, P, x, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn furst_electrolyte_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    x: Vec<f64>,
    compressed_phase: &str,
) -> PyResult<crate::results::PyFurstElectrolytePhaseResult> {
    azoth_eos::furst_electrolyte_phase(&components, kelvins(T), pascals(P), &x, compressed_phase)
        .map(|r| crate::results::PyFurstElectrolytePhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The 2004 revision of the Fürst electrolyte phase state, computed in Rust.
///
/// The same kernels as `furst_electrolyte_phase` with five quantities zeroed and one term
/// added; the resolver carries the difference, so this is the same call.
#[pyfunction]
#[pyo3(signature = (components, T, P, x, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn furst_electrolyte_mod2004_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    x: Vec<f64>,
    compressed_phase: &str,
) -> PyResult<crate::results::PyFurstElectrolyteMod2004PhaseResult> {
    azoth_eos::furst_electrolyte_mod2004_phase(
        &components,
        kelvins(T),
        pascals(P),
        &x,
        compressed_phase,
    )
    .map(|r| crate::results::PyFurstElectrolyteMod2004PhaseResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The SRK-CPA phase state, computed in Rust.
///
/// **The component names cross unresolved**, and this side looks them up in the Rust
/// databank. That is the `eos.eos_cg_phase` precedent, and for this model it is what makes
/// the two-kernel comparison cover the *resolution* as well as the arithmetic: an
/// associating mixture mixes with `cpakij_SRK` and a classical one with `KIJPR`, and on
/// water/methanol those differ by a factor of two.
#[pyfunction]
#[pyo3(signature = (components, T, P, z, compressed_phase))]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the chemistry
pub fn srk_cpa_phase(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    z: Vec<f64>,
    compressed_phase: &str,
) -> PyResult<crate::results::PySrkCpaPhaseResult> {
    azoth_eos::srk_cpa_phase(&components, kelvins(T), pascals(P), &z, compressed_phase)
        .map(|r| crate::results::PySrkCpaPhaseResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The effective diffusion coefficients of a phase, from its binary matrix.
///
/// The matrix crosses row-major and unresolved: which pair coefficient is which is the
/// caller's matrix, and this is the eight-line assembly and nothing else.
#[pyfunction]
#[pyo3(signature = (binary_diffusion, x))]
#[pyo3(text_signature = "(binary_diffusion, x)")]
pub fn effective_diffusion(
    py: Python<'_>,
    binary_diffusion: Vec<Vec<f64>>,
    x: Vec<f64>,
) -> PyResult<crate::results::PyEffectiveDiffusionResult> {
    azoth_eos::effective_diffusion::effective_diffusion(&binary_diffusion, &x)
        .map(|r| crate::results::PyEffectiveDiffusionResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The generalised gamma-phi flash: a cubic vapour over a named activity-coefficient liquid.
///
/// The names cross **unresolved** and both halves resolve on the Rust side - the liquid's
/// parameters through `liquid_model` and the vapour's constants through `cubic` - so the two
/// languages cannot disagree about which row answered, and a caller cannot pair a liquid
/// with a vapour built from a different component list.
#[pyfunction]
#[pyo3(signature = (components, liquid_model, T, P, z, cubic))]
#[pyo3(text_signature = "(components, liquid_model, T, P, z, cubic)")]
#[allow(non_snake_case)] // `T`, `P` and `z` are the symbols in the chemistry
pub fn ge_flash(
    py: Python<'_>,
    components: Vec<String>,
    liquid_model: &str,
    T: f64,
    P: f64,
    z: Vec<f64>,
    cubic: &str,
) -> PyResult<crate::results::PyGeFlashResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    azoth_eos::ge_flash::ge_flash(&names, cubic, liquid_model, kelvins(T), pascals(P), &z)
        .map(|r| crate::results::PyGeFlashResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}
