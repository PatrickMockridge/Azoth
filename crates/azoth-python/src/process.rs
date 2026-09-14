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

use azoth_core::units::{kelvins, pascals, watts};
use azoth_process as process;
use pyo3::prelude::*;

use crate::eos::build_mixture;
use crate::errors::to_pyerr;
use crate::results::PySeparatorResult;

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

/// The ideal-gas model the enthalpy is measured from, built from the same six vectors
/// `eos.ph_flash` takes.
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
