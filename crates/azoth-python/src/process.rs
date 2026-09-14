//! The unit operations, exposed to Python.
//!
//! Same rule as `eos.rs`, `hydraulics.rs` and `thermal.rs`: every parameter below is an
//! SI magnitude, not a `uom` quantity. A unit operation's arguments are the same
//! fifteen an equation-of-state model takes - the mixture, the ideal-gas block, a
//! temperature, a pressure and a composition - plus the two that make it a stream: a
//! molar flow and whatever the unit does to it.
//!
//! `n` is the one dimensional argument here with no `uom` type behind it, because
//! `uom` has no molar-flow quantity. It crosses as mol/s, which is already SI base, so
//! nothing is lost - see `crates/azoth-core/src/units.rs`, where the vocabulary says
//! the same thing.
//!
//! # One shape, eight times
//!
//! Every function below is the same five lines: build the mixture, build the ideal-gas
//! model, call the crate, adapt the result, map the error. The repetition is deliberate
//! rather than a missing abstraction - each `#[pyfunction]` is its own signature and its
//! own docstring, and the two helpers that *are* shared are the two things that would
//! otherwise be decided eight times: the component ordering and the ideal-gas block.
//! The same decision is recorded on the Python side, in `_rust_bridge.py`.

use azoth_core::units::{kelvins, pascals, watts};
use azoth_process as process;
use pyo3::prelude::*;

use crate::eos::build_mixture;
use crate::errors::to_pyerr;
use crate::results::{
    PyCompressorResult, PyExpanderResult, PyHeaterResult, PyMixerResult, PyPumpResult,
    PySeparatorResult, PySplitterResult, PyThrottlingValveResult,
};

/// The ideal-gas model the enthalpy and the entropy are measured from, built from the
/// same eight arguments `eos.ph_flash` and `eos.ps_flash` take.
///
/// A free function rather than a closure so every unit operation that needs one - and
/// most of them do, because most of them move energy - builds it the same way.
#[allow(non_snake_case)] // `T_ref` and `P_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub(crate) fn eos_ideal_gas(
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
) -> azoth_eos::IdealGasModel {
    azoth_eos::IdealGasModel {
        cp_a,
        cp_b,
        cp_c,
        cp_d,
        h_ref,
        s_ref,
        t_ref: kelvins(T_ref),
        p_ref: pascals(P_ref),
    }
}

/// One feed split into a gas and a liquid at a single temperature and pressure.
///
/// The first unit operation ported from NeqSim. Its physics is one flash and the
/// arithmetic that follows from it: `TPflash` at the reduced pressure when no duty is
/// set, `PHflash` at the feed's enthalpy plus the duty when one is, and the vapour
/// fraction splitting the molar flow.
#[pyfunction]
#[pyo3(
    signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, pressure_drop, heat_duty)
)]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, pressure_drop, heat_duty)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn separator(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
    T: f64,
    P: f64,
    n: f64,
    z: Vec<f64>,
    pressure_drop: f64,
    heat_duty: f64,
) -> PyResult<PySeparatorResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    let ideal_gas = eos_ideal_gas(cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref);
    process::separator(
        &mixture,
        &ideal_gas,
        kelvins(T),
        pascals(P),
        n,
        &z,
        pascals(pressure_drop),
        watts(heat_duty),
    )
    .map(|r| PySeparatorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// Several feeds blended into one, at the lowest inlet pressure.
///
/// The first unit whose state is a *set* of streams: `T`, `P` and `n` cross as vectors
/// with one entry per inlet, and `z` crosses **flattened row-major** - inlet `s` occupies
/// `z[s * N .. (s + 1) * N]`, with `N` derived from the mixture. That is the same
/// flattening `kij` uses, and for the same reason: the stream count is then derivable
/// from the arguments rather than being a fourth vector that could disagree with them.
#[pyfunction]
#[pyo3(
    signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z)
)]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn mixer(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
    T: Vec<f64>,
    P: Vec<f64>,
    n: Vec<f64>,
    z: Vec<f64>,
) -> PyResult<PyMixerResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    let ideal_gas = eos_ideal_gas(cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref);
    let temperatures: Vec<_> = T.into_iter().map(kelvins).collect();
    let pressures: Vec<_> = P.into_iter().map(pascals).collect();
    process::mixer(&mixture, &ideal_gas, &temperatures, &pressures, &n, &z)
        .map(|r| PyMixerResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// One feed divided into branches at the feed's own temperature and pressure.
///
/// The only unit operation that takes **no ideal-gas block**, because the only one that
/// does no energy balance. Its single thermodynamic call is a `TPflash`, which is a
/// function of the temperature, the pressure and the composition.
#[pyfunction]
#[pyo3(signature = (Tc, Pc, omega, kij, T, P, n, z, fractions))]
#[pyo3(text_signature = "(Tc, Pc, omega, kij, T, P, n, z, fractions)")]
#[allow(non_snake_case)] // `Tc` and `Pc` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn splitter(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    T: f64,
    P: f64,
    n: f64,
    z: Vec<f64>,
    fractions: Vec<f64>,
) -> PyResult<PySplitterResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    process::splitter(&mixture, kelvins(T), pascals(P), n, &z, &fractions)
        .map(|r| PySplitterResult::from(&r))
        .map_err(|e| to_pyerr(py, e))
}

/// A stream's pressure dropped at constant enthalpy.
///
/// **No molar flow crosses**, because the model takes none: an isenthalpic flash is a
/// molar property, so the outlet state does not depend on the flow at all.
#[pyfunction]
#[pyo3(
    signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, z, pressure_drop)
)]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, z, pressure_drop)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn throttling_valve(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
    T: f64,
    P: f64,
    z: Vec<f64>,
    pressure_drop: f64,
) -> PyResult<PyThrottlingValveResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    let ideal_gas = eos_ideal_gas(cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref);
    process::throttling_valve(
        &mixture,
        &ideal_gas,
        kelvins(T),
        pascals(P),
        &z,
        pascals(pressure_drop),
    )
    .map(|r| PyThrottlingValveResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// A duty applied to a stream at a fixed pressure. A negative duty is a cooler.
#[pyfunction]
#[pyo3(
    signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, pressure_drop, heat_duty)
)]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, pressure_drop, heat_duty)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn heater(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
    T: f64,
    P: f64,
    n: f64,
    z: Vec<f64>,
    pressure_drop: f64,
    heat_duty: f64,
) -> PyResult<PyHeaterResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    let ideal_gas = eos_ideal_gas(cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref);
    process::heater(
        &mixture,
        &ideal_gas,
        kelvins(T),
        pascals(P),
        n,
        &z,
        pascals(pressure_drop),
        watts(heat_duty),
    )
    .map(|r| PyHeaterResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// A pressure rise at a stated isentropic efficiency.
///
/// The first unit operation whose answer depends on an *entropy*, and so the first to
/// need `eos.ps_flash`. `efficiency` crosses as a plain fraction, not a percentage:
/// `0.75`, not `75`.
#[pyfunction]
#[pyo3(
    signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, outlet_pressure, efficiency)
)]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, outlet_pressure, efficiency)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn compressor(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
    T: f64,
    P: f64,
    n: f64,
    z: Vec<f64>,
    outlet_pressure: f64,
    efficiency: f64,
) -> PyResult<PyCompressorResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    let ideal_gas = eos_ideal_gas(cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref);
    process::compressor(
        &mixture,
        &ideal_gas,
        kelvins(T),
        pascals(P),
        n,
        &z,
        pascals(outlet_pressure),
        efficiency,
    )
    .map(|r| PyCompressorResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// A pressure rise in a liquid, at a stated isentropic efficiency.
#[pyfunction]
#[pyo3(
    signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, outlet_pressure, efficiency)
)]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, outlet_pressure, efficiency)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pump(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
    T: f64,
    P: f64,
    n: f64,
    z: Vec<f64>,
    outlet_pressure: f64,
    efficiency: f64,
) -> PyResult<PyPumpResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    let ideal_gas = eos_ideal_gas(cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref);
    process::pump(
        &mixture,
        &ideal_gas,
        kelvins(T),
        pascals(P),
        n,
        &z,
        pascals(outlet_pressure),
        efficiency,
    )
    .map(|r| PyPumpResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// A pressure drop that produces work.
///
/// The same procedure as a compressor with the efficiency multiplying rather than
/// dividing, so `power` crosses **negative** - the fluid is doing the work.
#[pyfunction]
#[pyo3(
    signature = (Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, outlet_pressure, efficiency)
)]
#[pyo3(
    text_signature = "(Tc, Pc, omega, kij, cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref, T, P, n, z, outlet_pressure, efficiency)"
)]
#[allow(non_snake_case)] // `Tc`, `Pc` and `T_ref` are the symbols in the chemistry
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn expander(
    py: Python<'_>,
    Tc: Vec<f64>,
    Pc: Vec<f64>,
    omega: Vec<f64>,
    kij: Vec<f64>,
    cp_a: Vec<f64>,
    cp_b: Vec<f64>,
    cp_c: Vec<f64>,
    cp_d: Vec<f64>,
    h_ref: Vec<f64>,
    s_ref: Vec<f64>,
    T_ref: f64,
    P_ref: f64,
    T: f64,
    P: f64,
    n: f64,
    z: Vec<f64>,
    outlet_pressure: f64,
    efficiency: f64,
) -> PyResult<PyExpanderResult> {
    let mixture = build_mixture(py, &Tc, &Pc, &omega, kij)?;
    let ideal_gas = eos_ideal_gas(cp_a, cp_b, cp_c, cp_d, h_ref, s_ref, T_ref, P_ref);
    process::expander(
        &mixture,
        &ideal_gas,
        kelvins(T),
        pascals(P),
        n,
        &z,
        pascals(outlet_pressure),
        efficiency,
    )
    .map(|r| PyExpanderResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
