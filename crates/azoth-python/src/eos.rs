//! The equations-of-state calculations, exposed to Python.
//!
//! Same rule as `hydraulics.rs` and `thermal.rs`: every parameter below is an SI
//! magnitude, not a `uom` quantity. Here that rule costs nothing, because the
//! quantities in this namespace are dimensionless by construction - the reduced
//! variables and constitutive coefficients an equation of state is written in -
//! so there is no conversion to do at the boundary and no unit string to keep in
//! step with the spec.

use azoth_core::units::{
    cubic_meters_per_mole, joules_per_mole, joules_per_mole_kelvin, kelvins, kilograms_per_mole,
    pascals,
};
use pyo3::prelude::*;

use crate::errors::to_pyerr;
use crate::results::{
    PyCriticalPointResult, PyIdealGasCpResult, PyMolarEnthalpyEntropyResult, PyPhFlashResult,
    PyPhaseBoundaryResult, PyPrAlphaAbResult, PyPrDepartureResult, PyPrKappaResult,
    PyPrMassDensityResult, PyPrMolarVolumeResult, PyPrZFactorResult, PyPrsvKappaResult,
    PyPsFlashResult, PyPtFlashResult, PyPureSaturationResult, PyRachfordRiceBinaryResult,
    PySrkAlphaAbResult, PySrkDepartureResult, PySrkKappaResult, PySrkZFactorResult,
    PyStabilityTestResult, PyVdw1fMixBinaryResult,
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
pub(crate) fn build_mixture(
    py: Python<'_>,
    Tc: &[f64],
    Pc: &[f64],
    omega: &[f64],
    kij: Vec<f64>,
    eos: &str,
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
        .map(|i| azoth_eos::Component::new(kelvins(Tc[i]), pascals(Pc[i]), omega[i]))
        .collect::<azoth_core::Result<Vec<_>>>()
        .map_err(|e| to_pyerr(py, e))?;
    let cubic: azoth_eos::Cubic = eos
        .parse()
        .map_err(pyo3::exceptions::PyValueError::new_err)?;
    azoth_eos::Mixture::new(components, kij)
        .map(|m| m.with_cubic(cubic))
        .map_err(|e| to_pyerr(py, e))
}

#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, T, P, z, eos = "pr"))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, P, z, eos = \"pr\")")]
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
) -> PyResult<PyPtFlashResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
    azoth_eos::pt_flash(&mixture, kelvins(T), pascals(P), &z)
        .map(|r| PyPtFlashResult::from(&r))
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, H, z, eos = "pr"))]
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
) -> PyResult<PyPhFlashResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
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
#[pyo3(signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, P, S, z, eos = "pr"))]
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
) -> PyResult<PyPsFlashResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
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
#[pyo3(signature = (Tc, Pc, omega, kij, T, P, z, eos = "pr"))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, P, z, eos = \"pr\")")]
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
) -> PyResult<PyStabilityTestResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
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
#[pyo3(signature = (Tc, Pc, omega, kij, T, held, eos = "pr"))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, held, eos = \"pr\")")]
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
) -> PyResult<PyPhaseBoundaryResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
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
#[pyo3(signature = (Tc, Pc, omega, kij, z, eos = "pr"))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, z, eos = \"pr\")")]
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
) -> PyResult<PyCriticalPointResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
    azoth_eos::critical_point(&mixture, &z)
        .map(|r| PyCriticalPointResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// The pressure at which a vapour of composition `held` first condenses.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, T, held, eos = "pr"))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, held, eos = \"pr\")")]
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
) -> PyResult<PyPhaseBoundaryResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
    azoth_eos::dew_pressure(&mixture, kelvins(T), &held)
        .map(|r| PyPhaseBoundaryResult::from(&r))
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
    compressibility, eos = "pr"
))]
#[pyo3(text_signature = "(
    Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, cp_e, T, P, z,
    compressibility, eos = \"pr\"
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
) -> PyResult<PyMolarEnthalpyEntropyResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij, eos)?;
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
