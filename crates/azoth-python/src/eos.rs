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
    PyAmmoniaPhaseResult, PyAntoineVaporPressureResult, PyBwrsPhaseResult,
    PyChungConductivityResult, PyChungViscosityResult, PyCo2WaterDiffusivityResult,
    PyCostaldMolarVolumeResult, PyCriticalPointResult, PyHaydukMinhasDiffusivityResult,
    PyHeatOfVaporizationResult, PyIdealGasCpResult, PyLiquidHeatCapacityResult,
    PyMasonSaxenaConductivityResult, PyMatcop5PrumrAlphaResult, PyMatcopAlphaResult,
    PyMatcopPrAlphaResult, PyMatcopPrumrAlphaResult, PyMatcopPrumrNewAlphaResult,
    PyMolarEnthalpyEntropyResult, PyMollerupAlphaResult, PyNrtlActivityCoefficientsResult,
    PyParachorSurfaceTensionResult, PyPhFlashResult, PyPhaseBoundaryResult,
    PyPhaseBoundaryTemperatureResult, PyPhaseEnvelopeResult, PyPr78KappaResult, PyPrAlphaAbResult,
    PyPrDaneshAlphaResult, PyPrDelft1998AlphaResult, PyPrDepartureResult,
    PyPrGassem2001AlphaResult, PyPrKappaResult, PyPrLeeKeslerAlphaResult, PyPrMassDensityResult,
    PyPrMolarVolumeResult, PyPrPenelouxShiftResult, PyPrZFactorResult, PyPrsvKappaResult,
    PyPsFlashResult, PyPtFlashResult, PyPuFlashResult, PyPureSaturationResult, PyPvFlashResult,
    PyRachfordRiceBinaryResult, PyRackettMolarVolumeResult, PyRkAlphaAbResult, PyRkDepartureResult,
    PySchwartzentruberAlphaResult, PySiddiqiLucasDiffusivityResult, PySoreideWhitsonAlphaResult,
    PySrkAlphaAbResult, PySrkDepartureResult, PySrkKappaResult, PySrkPenelouxShiftResult,
    PySrkZFactorResult, PyStabilityTestResult, PyThFlashResult, PyThermalConductivityResult,
    PyTsFlashResult, PyTuFlashResult, PyTvFlashResult, PyTwuKappaResult, PyTwucoonAlphaResult,
    PyTwucoonParamAlphaResult, PyTwucoonStatoilAlphaResult, PyTynCalusDiffusivityResult,
    PyUmrprAlphaResult, PyUnifacActivityCoefficientsResult, PyUniquacActivityCoefficientsResult,
    PyVdw1fMixBinaryResult, PyViscosityResult, PyVuFlashResult, PyWilkeChangDiffusivityResult,
    PyWilkeViscosityResult, PyWilsonActivityCoefficientsResult,
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
#[pyo3(signature = (omega, Tc, Pc))]
#[pyo3(text_signature = "(omega, Tc, Pc)")]
#[allow(non_snake_case)] // `Tc` and `Pc` are the symbols in the published equation
pub fn pr_peneloux_shift(
    py: Python<'_>,
    omega: f64,
    Tc: f64,
    Pc: f64,
) -> PyResult<PyPrPenelouxShiftResult> {
    azoth_eos::pr_peneloux_shift(omega, kelvins(Tc), pascals(Pc))
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
/// calculation: its arguments are vectors and matrices.
#[pyfunction]
#[pyo3(signature = (T, x, Dij, alpha))]
#[pyo3(text_signature = "(T, x, Dij, alpha)")]
#[allow(non_snake_case)] // `Dij` and `T` are the symbols in the chemistry
pub fn nrtl_activity_coefficients(
    py: Python<'_>,
    T: f64,
    x: Vec<f64>,
    Dij: Vec<Vec<f64>>,
    alpha: Vec<Vec<f64>>,
) -> PyResult<PyNrtlActivityCoefficientsResult> {
    azoth_eos::nrtl_activity_coefficients(T, &x, &Dij, &alpha)
        .map(|r| PyNrtlActivityCoefficientsResult::from(&r))
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
/// calculation: its arguments are vectors and matrices.
#[pyfunction]
#[pyo3(signature = (T, x, groups, group_r, group_q, aij))]
#[pyo3(text_signature = "(T, x, groups, group_r, group_q, aij)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn unifac_activity_coefficients(
    py: Python<'_>,
    T: f64,
    x: Vec<f64>,
    groups: Vec<Vec<f64>>,
    group_r: Vec<f64>,
    group_q: Vec<f64>,
    aij: Vec<Vec<f64>>,
) -> PyResult<PyUnifacActivityCoefficientsResult> {
    azoth_eos::unifac_activity_coefficients(T, &x, &groups, &group_r, &group_q, &aij)
        .map(|r| PyUnifacActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from UNIQUAC. A *model* rather than a
/// calculation: its arguments are vectors and a matrix.
#[pyfunction]
#[pyo3(signature = (T, x, r, q, aij))]
#[pyo3(text_signature = "(T, x, r, q, aij)")]
#[allow(non_snake_case)] // `T` and `x` are the symbols in the chemistry
pub fn uniquac_activity_coefficients(
    py: Python<'_>,
    T: f64,
    x: Vec<f64>,
    r: Vec<f64>,
    q: Vec<f64>,
    aij: Vec<Vec<f64>>,
) -> PyResult<PyUniquacActivityCoefficientsResult> {
    azoth_eos::uniquac_activity_coefficients(T, &x, &r, &q, &aij)
        .map(|r| PyUniquacActivityCoefficientsResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The activity coefficients of a mixture, from the paraffin-wax Wilson model. A
/// *model* rather than a calculation: its arguments are vectors.
#[pyfunction]
#[pyo3(signature = (T, x, M, Tc))]
#[pyo3(text_signature = "(T, x, M, Tc)")]
#[allow(non_snake_case)] // `M`, `Tc`, `T` and `x` are the symbols in the chemistry
pub fn wilson_activity_coefficients(
    py: Python<'_>,
    T: f64,
    x: Vec<f64>,
    M: Vec<f64>,
    Tc: Vec<f64>,
) -> PyResult<PyWilsonActivityCoefficientsResult> {
    azoth_eos::wilson_activity_coefficients(T, &x, &M, &Tc)
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
    let components = (0..n)
        .map(|i| {
            let mut component =
                azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i])?;
            if let Some(params) = alpha_params.and_then(|all| all.get(i)) {
                component = component.with_alpha_params(params.clone());
            }
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
        .map(|m| m.with_cubic(cubic).with_alpha(alpha))
        .map_err(|e| to_pyerr(py, e))
}

#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pt_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pt_phase_envelope(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, H, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, H, z, eos = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T_ref` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn ph_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, S, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, S, z, eos = \"pr\")"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T_ref` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn ps_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, V, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, V, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn tv_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, V, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, V, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc`, `P` and the rest are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pv_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, H, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, H, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn th_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, S, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, S, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn ts_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, U, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, U, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn tu_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, U, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, U, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn pu_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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

/// The volume-internal-energy flash of a mixture (V,U -> P,T).
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, V, U, z,
    eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, V, U, z, eos = \"pr\", alpha = \"pr\", alpha_params = None)"
)]
#[allow(non_snake_case)]
#[allow(clippy::too_many_arguments)]
pub fn vu_flash(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn stability_test(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::stability_test(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyStabilityTestResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pressure at which a liquid of composition `x` first gives off vapour.
///
/// The first of the two phase-boundary models. They take the same six arguments as
/// the flash - three per-component vectors, a flattened interaction matrix,
/// temperature and a composition - and differ only in which composition it is,
/// which is why they are two functions rather than one with a switch.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, T, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn bubble_pressure(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `z` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn critical_point(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, T, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn dew_pressure(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
#[pyo3(signature = (Tc, Pc, omega, kij, P, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, P, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn bubble_temperature(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
        eos,
        alpha,
        alpha_params.as_deref(),
    )?;
    azoth_eos::bubble_temperature(&mixture, pascals(P), &held)
        .map(|r| PyPhaseBoundaryTemperatureResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The temperature at which a vapour of composition `held` first gives off liquid.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, P, held, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, P, held, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn dew_temperature(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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
    Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z,
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

/// The liquid viscosity from the Pedersen (PFCT) heavy-oil correlation.
///
/// `eos`, `alpha` and `alpha_params` are accepted for the boundary's uniformity with the
/// other mixture models and ignored: the reference flash is pure methane SRK, so the
/// mixture's own cubic and alpha are not consulted.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, molar_mass, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, molar_mass, T, P, z, eos = \"pr\", alpha = \"pr\")")]
#[allow(non_snake_case)] // `Tc`, `Pc`, `T` and `P` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
#[allow(unused_variables)] // `eos`, `alpha` and `alpha_params` are boundary-only, see above.
pub fn viscosity(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
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

/// The liquid thermal conductivity from the Pedersen (PFCT) correlation.
///
/// `eos`, `alpha` and `alpha_params` are accepted for the boundary's uniformity and
/// ignored: the reference flash is pure methane SRK.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, molar_mass, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z, eos = "pr", alpha = "pr", alpha_params = None))]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, molar_mass, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z, eos = \"pr\", alpha = \"pr\")"
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
/// namespace means a second table chained in here. `process` was one - a unit-operation
/// layer that asserted a process simulator this library does not have - and it was
/// deleted rather than repaired; `eos` is the whole of what remains.
///
/// **This is the function that has to change when a namespace is added.** It read
/// `azoth_eos::model_gen` alone until `process` arrived, and the three `model_*` functions
/// below stayed green for exactly as long as nobody asked them about a unit operation.
fn all_models() -> impl Iterator<Item = &'static azoth_core::ModelSpec> {
    azoth_eos::model_gen::models().iter().copied()
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
