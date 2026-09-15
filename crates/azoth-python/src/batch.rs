//! The batch entry point: one call, N evaluations, arrays in and arrays out.
//!
//! # The decision this module embodies
//!
//! **This is a loop over the existing scalar kernels, not a vectorised kernel.**
//! `README.md` and `docs/src/architecture/specification.md` both stake the project on there being exactly two
//! implementations that check each other. A vectorised Rust kernel would be a third -
//! a different shape and a different arithmetic order, cross-checked by nothing - so
//! the claim that two independent implementations agree would quietly stop being true
//! while every test still passed.
//!
//! The optimisation does not need it either. The README's own argument is that the PyO3
//! call overhead exceeds the cost of the arithmetic, so removing N-1 boundary crossings
//! is where the win is, and looping here captures exactly that.
//!
//! # Why one entry point rather than one per calc
//!
//! `batch_run` takes a calc id and a name-to-array map, because the per-calc argument
//! names and units belong on the Python side where the spec already describes them. The
//! alternative - a `#[pyfunction]` per calc - would restate every signature here in a
//! second place, and the two would drift.
//!
//! The `match` below is the dispatch, and it is exhaustive by construction: a calc added
//! to the registry without a batch arm fails a test rather than silently doing nothing.

use azoth_core::Warning;
use azoth_core::units::{
    cubic_meters_per_mole, cubic_meters_per_second, kelvin_intervals, kelvins,
    kilograms_per_cubic_meter, kilograms_per_mole, kilograms_per_second, meters, meters_per_second,
    pascal_seconds, pascals, square_meters, watts_per_meter_kelvin,
};
use azoth_eos as eos;
use azoth_hydraulics as hyd;
use azoth_thermal as therm;
use pyo3::prelude::*;

use crate::results::PyWarning;

/// The inputs of one batch call, as SI magnitudes keyed by the spec's input names.
type Inputs = std::collections::HashMap<String, Vec<f64>>;

/// One output column: its name, its unit for display, and its values.
///
/// Exactly one of `values` and `labels` is `Some`. Two variants rather than one
/// generic column because they are genuinely different things: a numeric column is an
/// array of SI magnitudes, and an enum column is a sequence of labels - `regime` is
/// `laminar`, `transitional` or `turbulent`, and there is no number it should be.
///
/// Encoding the enum as an index into a table would make the caller look the mapping up
/// to read a value, and would need a sentinel for "not available" that is not one of the
/// three real answers. `None` in a label column says that instead, and says it
/// unambiguously.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "BatchColumn"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyBatchColumn {
    /// Output field name, as the spec declares it.
    #[pyo3(get)]
    pub name: String,
    /// Canonical unit string, e.g. `Pa`, or `dimensionless`.
    #[pyo3(get)]
    pub unit: String,
    /// The values, in SI base magnitudes. `None` for an enum column.
    #[pyo3(get)]
    pub values: Option<Vec<f64>>,
    /// The labels, `None` where the value is absent. `None` for a numeric column.
    #[pyo3(get)]
    pub labels: Option<Vec<Option<String>>>,
}

/// The result of one batch call.
///
/// Struct-of-arrays, with warnings positional: `warnings[i]` is element `i`'s warnings,
/// empty for a clean element. Doing it the other way - a table of rows - is what the
/// batch shape exists to avoid, and an array of result objects would allocate N of them
/// plus N quantities, giving back most of what crossing the boundary once saved.
#[pyclass(
    frozen,
    skip_from_py_object,
    module = "azoth._core",
    name = "BatchResult"
)]
#[derive(Debug, Clone, PartialEq)]
pub struct PyBatchResult {
    /// One entry per output field, in the order the spec declares them.
    #[pyo3(get)]
    pub columns: Vec<PyBatchColumn>,
    /// Per element, in order. Empty tuples for clean elements.
    #[pyo3(get)]
    pub warnings: Vec<Vec<PyWarning>>,
}

/// Collect one numeric output column.
fn push_values(columns: &mut Vec<PyBatchColumn>, name: &str, unit: &str, values: Vec<f64>) {
    columns.push(PyBatchColumn {
        name: name.to_string(),
        unit: unit.to_string(),
        values: Some(values),
        labels: None,
    });
}

/// Collect one enum output column.
fn push_labels(columns: &mut Vec<PyBatchColumn>, name: &str, labels: Vec<Option<String>>) {
    columns.push(PyBatchColumn {
        name: name.to_string(),
        unit: "dimensionless".to_string(),
        values: None,
        labels: Some(labels),
    });
}

/// Take a named input array, or fail.
///
/// # Errors
/// Returns `InvalidInput` naming the missing key. A missing array is a programming
/// error on the Python side rather than a user condition, but it is reported as an
/// error rather than panicking, and the message names the key so the fault is obvious.
fn take(inputs: &Inputs, name: &str) -> PyResult<Vec<f64>> {
    inputs
        .get(name)
        .cloned()
        .ok_or_else(|| pyo3::exceptions::PyKeyError::new_err(format!("missing input `{name}`")))
}

fn take_optional(inputs: &Inputs, name: &str) -> Option<Vec<f64>> {
    inputs.get(name).cloned()
}

/// Warned values from a scalar result, transported.
fn transport(warnings: &[Warning]) -> Vec<PyWarning> {
    warnings.iter().map(PyWarning::from).collect()
}

/// One element's outcome: its warnings transported, or its error mapped.
///
/// Every arm below is `n` iterations of one scalar call, and this is the part of
/// each iteration that is identical. Writing it once means the arms differ only in
/// which kernel they call and which fields they read, so an arm cannot quietly
/// forget to transport a warning - the one omission no numerical test would catch.
fn element<T: azoth_core::CalcResult>(
    py: Python<'_>,
    outcome: Result<T, azoth_core::AzothError>,
    warnings: &mut Vec<Vec<PyWarning>>,
) -> PyResult<T> {
    match outcome {
        Ok(value) => {
            warnings.push(transport(value.warnings()));
            Ok(value)
        }
        Err(error) => Err(crate::errors::to_pyerr(py, error)),
    }
}

/// A boolean output as a column value.
///
/// A batch column is an array of numbers, so `converged` arrives as 1.0 or 0.0
/// rather than as `True` or `False`. The information is exact and the conversion
/// is lossless - unlike an enum, which is why an enum gets a label column instead.
/// Stated in the Python-side result docs, where a caller will meet it.
fn flag(value: bool) -> f64 {
    if value { 1.0 } else { 0.0 }
}

/// Run one calculation over arrays of inputs.
///
/// # Errors
/// Propagates whatever the scalar kernel raises, **fail-fast**: if element 7 of 1000
/// violates an error-severity bound the whole call raises rather than returning 999
/// values with a silent gap. A partial result is the failure this library is organised
/// against - a caller who does not notice which element is missing has a wrong answer
/// that looks complete.
///
/// Also raises for a calc id with no batch arm, which is a programming error: the
/// registry and this dispatch are supposed to cover the same calcs, and a test asserts
/// they do.
#[pyfunction]
pub fn batch_run(py: Python<'_>, calc_id: &str, inputs: Inputs) -> PyResult<PyBatchResult> {
    let mut columns: Vec<PyBatchColumn> = Vec::new();
    let mut warnings: Vec<Vec<PyWarning>> = Vec::new();

    // The length every input must share. Checked here as well as on the Python side,
    // because "the arrays were all the same length" is the assumption the whole loop
    // rests on and an index panic is a poor way to discover it is false.
    let n = inputs.values().map(Vec::len).max().unwrap_or(0);
    for (name, values) in &inputs {
        if values.len() != n {
            return Err(pyo3::exceptions::PyValueError::new_err(format!(
                "input `{name}` has {} value(s) but the batch is {n} long",
                values.len()
            )));
        }
    }

    match calc_id {
        "hydraulics.reynolds_number" => {
            let (rho, v, d, mu) = (
                take(&inputs, "rho")?,
                take(&inputs, "v")?,
                take(&inputs, "D")?,
                take(&inputs, "mu")?,
            );
            let mut re = Vec::with_capacity(n);
            let mut regime = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    hyd::reynolds_number(
                        kilograms_per_cubic_meter(rho[i]),
                        meters_per_second(v[i]),
                        meters(d[i]),
                        pascal_seconds(mu[i]),
                    ),
                    &mut warnings,
                )?;
                re.push(r.re);
                regime.push(Some(r.regime.as_str().to_string()));
            }
            push_values(&mut columns, "re", "dimensionless", re);
            push_labels(&mut columns, "regime", regime);
        }

        "hydraulics.darcy_weisbach" => {
            let (f, l, d, rho, v) = (
                take(&inputs, "f")?,
                take(&inputs, "L")?,
                take(&inputs, "D")?,
                take(&inputs, "rho")?,
                take(&inputs, "v")?,
            );
            // Optional per *batch*, not per element. `mu` is present for every element
            // or for none: a partly-supplied optional input has no sensible meaning
            // here, and inventing one would need a per-element absent-value encoding
            // for a case no caller has asked for. Absent for the whole batch means
            // every element carries RANGE_CHECK_SKIPPED, which is the honest answer.
            let mu = take_optional(&inputs, "mu");
            let mut dp = Vec::with_capacity(n);
            let mut re = Vec::with_capacity(n);
            let mut regime = Vec::with_capacity(n);
            for i in 0..n {
                let outcome = element(
                    py,
                    hyd::darcy_weisbach(
                        f[i],
                        meters(l[i]),
                        meters(d[i]),
                        kilograms_per_cubic_meter(rho[i]),
                        meters_per_second(v[i]),
                        mu.as_ref().map(|m| pascal_seconds(m[i])),
                    ),
                    &mut warnings,
                )?;
                dp.push(outcome.dp.value);
                // Absent rather than zero when viscosity was omitted: a Reynolds number
                // of zero is a physical claim, and the caller omitted the input that
                // would have let this calc make one. The RANGE_CHECK_SKIPPED warning on
                // the same element carries the reason.
                re.push(outcome.re.unwrap_or(f64::NAN));
                regime.push(outcome.regime.map(|r| r.as_str().to_string()));
            }
            push_values(&mut columns, "dp", "Pa", dp);
            push_values(&mut columns, "re", "dimensionless", re);
            push_labels(&mut columns, "regime", regime);
        }

        "hydraulics.friction_factor_colebrook" => {
            let (re, rr) = (take(&inputs, "re")?, take(&inputs, "relative_roughness")?);
            let (mut f, mut iterations, mut converged, mut residual) = (
                Vec::with_capacity(n),
                Vec::with_capacity(n),
                Vec::with_capacity(n),
                Vec::with_capacity(n),
            );
            // The solver report travels with the answer because the answer is the last
            // iterate: without it a caller cannot tell a solution from a failure to
            // converge, which is exactly what the scalar result says.
            for i in 0..n {
                let r = element(
                    py,
                    hyd::friction_factor_colebrook(re[i], rr[i]),
                    &mut warnings,
                )?;
                f.push(r.f);
                iterations.push(f64::from(r.iterations));
                converged.push(flag(r.converged));
                residual.push(r.residual);
            }
            push_values(&mut columns, "f", "dimensionless", f);
            push_values(&mut columns, "iterations", "dimensionless", iterations);
            push_values(&mut columns, "converged", "dimensionless", converged);
            push_values(&mut columns, "residual", "dimensionless", residual);
        }

        "hydraulics.friction_factor_swamee_jain" => {
            let (re, rr) = (take(&inputs, "re")?, take(&inputs, "relative_roughness")?);
            let mut f = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    hyd::friction_factor_swamee_jain(re[i], rr[i]),
                    &mut warnings,
                )?;
                f.push(r.f);
            }
            push_values(&mut columns, "f", "dimensionless", f);
        }

        "hydraulics.friction_factor_haaland" => {
            let (re, rr) = (take(&inputs, "re")?, take(&inputs, "relative_roughness")?);
            let mut f = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    hyd::friction_factor_haaland(re[i], rr[i]),
                    &mut warnings,
                )?;
                f.push(r.f);
            }
            push_values(&mut columns, "f", "dimensionless", f);
        }

        "hydraulics.pump_power" => {
            let (rho, q, h, eta) = (
                take(&inputs, "rho")?,
                take(&inputs, "q")?,
                take(&inputs, "H")?,
                take(&inputs, "eta")?,
            );
            let mut power = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    hyd::pump_power(
                        kilograms_per_cubic_meter(rho[i]),
                        cubic_meters_per_second(q[i]),
                        meters(h[i]),
                        eta[i],
                    ),
                    &mut warnings,
                )?;
                power.push(r.power.value);
            }
            push_values(&mut columns, "power", "W", power);
        }

        "hydraulics.orifice_flow" => {
            let (d, dp, rho, cd) = (
                take(&inputs, "d")?,
                take(&inputs, "dP")?,
                take(&inputs, "rho")?,
                take(&inputs, "Cd")?,
            );
            let mut q = Vec::with_capacity(n);
            for i in 0..n {
                // `d` arrives in metres: the spec declares millimetres and the Python
                // side applies the factor once for the whole array. Converting here as
                // well would apply it twice.
                let r = element(
                    py,
                    hyd::orifice_flow(
                        meters(d[i]),
                        pascals(dp[i]),
                        kilograms_per_cubic_meter(rho[i]),
                        cd[i],
                    ),
                    &mut warnings,
                )?;
                q.push(r.q.value);
            }
            push_values(&mut columns, "q", "m**3/s", q);
        }

        "hydraulics.control_valve_cv" => {
            let (cv, dp, sg) = (
                take(&inputs, "Cv")?,
                take(&inputs, "dP")?,
                take(&inputs, "SG")?,
            );
            let mut q = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    hyd::control_valve_cv(cv[i], pascals(dp[i]), sg[i]),
                    &mut warnings,
                )?;
                q.push(r.q.value);
            }
            push_values(&mut columns, "q", "m**3/s", q);
        }

        "hydraulics.choked_flow_area" => {
            let (m_dot, p0, rho0, k) = (
                take(&inputs, "m_dot")?,
                take(&inputs, "P0")?,
                take(&inputs, "rho0")?,
                take(&inputs, "k")?,
            );
            let mut a = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    hyd::choked_flow_area(
                        kilograms_per_second(m_dot[i]),
                        pascals(p0[i]),
                        kilograms_per_cubic_meter(rho0[i]),
                        k[i],
                    ),
                    &mut warnings,
                )?;
                a.push(r.a.value);
            }
            push_values(&mut columns, "a", "m**2", a);
        }

        "thermal.conduction_plane_wall" => {
            let (k, a, dt, l) = (
                take(&inputs, "k")?,
                take(&inputs, "A")?,
                take(&inputs, "dT")?,
                take(&inputs, "L")?,
            );
            let mut q = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    therm::conduction_plane_wall(
                        watts_per_meter_kelvin(k[i]),
                        square_meters(a[i]),
                        // An interval, not an absolute temperature: `dT` is a difference
                        // and adding 273.15 to it would be silently wrong.
                        kelvin_intervals(dt[i]),
                        meters(l[i]),
                    ),
                    &mut warnings,
                )?;
                q.push(r.q.value);
            }
            push_values(&mut columns, "q", "W", q);
        }

        "eos.pr_kappa" => {
            let omega = take(&inputs, "omega")?;
            let mut kappa = Vec::with_capacity(n);
            for value in &omega {
                // No unit wrapping in either direction: `omega` is a genuine
                // dimensionless quantity, so it crosses as the number it is. That
                // is the same rule that makes `f` and `re` plain floats in the
                // hydraulics arms.
                let r = element(py, eos::pr_kappa(*value), &mut warnings)?;
                kappa.push(r.kappa);
            }
            push_values(&mut columns, "kappa", "dimensionless", kappa);
        }

        "eos.pr_alpha_ab" => {
            let (kappa, tr, pr) = (
                take(&inputs, "kappa")?,
                take(&inputs, "Tr")?,
                take(&inputs, "Pr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            let mut a_reduced = Vec::with_capacity(n);
            let mut b_reduced = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::pr_alpha_ab(kappa[i], tr[i], pr[i]), &mut warnings)?;
                alpha.push(r.alpha);
                a_reduced.push(r.a_reduced);
                b_reduced.push(r.b_reduced);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
            push_values(&mut columns, "a_reduced", "dimensionless", a_reduced);
            push_values(&mut columns, "b_reduced", "dimensionless", b_reduced);
        }

        "eos.pr_z_factor" => {
            let (a, b) = (take(&inputs, "a_reduced")?, take(&inputs, "b_reduced")?);
            let mut z_min = Vec::with_capacity(n);
            let mut z_max = Vec::with_capacity(n);
            let mut root_structure = Vec::with_capacity(n);
            let mut iterations = Vec::with_capacity(n);
            let mut converged = Vec::with_capacity(n);
            let mut residual = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::pr_z_factor(a[i], b[i]), &mut warnings)?;
                z_min.push(r.z_min);
                z_max.push(r.z_max);
                // A label column, like `reynolds_number`'s regime: an enum has no
                // numeric form, and inventing one would make a caller look the
                // mapping up to read a value.
                root_structure.push(Some(r.root_structure.as_str().to_string()));
                iterations.push(f64::from(r.iterations));
                // A batch column is an array of numbers, so the flag is one;
                // `bool(...)` recovers it and the conversion is exact.
                converged.push(flag(r.converged));
                residual.push(r.residual);
            }
            push_values(&mut columns, "z_min", "dimensionless", z_min);
            push_values(&mut columns, "z_max", "dimensionless", z_max);
            push_labels(&mut columns, "root_structure", root_structure);
            push_values(&mut columns, "iterations", "dimensionless", iterations);
            push_values(&mut columns, "converged", "dimensionless", converged);
            push_values(&mut columns, "residual", "dimensionless", residual);
        }

        "eos.srk_kappa" => {
            let omega = take(&inputs, "omega")?;
            let mut kappa = Vec::with_capacity(n);
            for value in &omega {
                let r = element(py, eos::srk_kappa(*value), &mut warnings)?;
                kappa.push(r.kappa);
            }
            push_values(&mut columns, "kappa", "dimensionless", kappa);
        }

        "eos.pr78_kappa" => {
            let omega = take(&inputs, "omega")?;
            let mut kappa = Vec::with_capacity(n);
            for value in &omega {
                let r = element(py, eos::pr78_kappa(*value), &mut warnings)?;
                kappa.push(r.kappa);
            }
            push_values(&mut columns, "kappa", "dimensionless", kappa);
        }

        "eos.twu_kappa" => {
            let omega = take(&inputs, "omega")?;
            let mut kappa = Vec::with_capacity(n);
            for value in &omega {
                let r = element(py, eos::twu_kappa(*value), &mut warnings)?;
                kappa.push(r.kappa);
            }
            push_values(&mut columns, "kappa", "dimensionless", kappa);
        }

        "eos.srk_alpha_ab" => {
            let (kappa, tr, pr) = (
                take(&inputs, "kappa")?,
                take(&inputs, "Tr")?,
                take(&inputs, "Pr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            let mut a_reduced = Vec::with_capacity(n);
            let mut b_reduced = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::srk_alpha_ab(kappa[i], tr[i], pr[i]), &mut warnings)?;
                alpha.push(r.alpha);
                a_reduced.push(r.a_reduced);
                b_reduced.push(r.b_reduced);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
            push_values(&mut columns, "a_reduced", "dimensionless", a_reduced);
            push_values(&mut columns, "b_reduced", "dimensionless", b_reduced);
        }

        "eos.srk_z_factor" => {
            let (a, b) = (take(&inputs, "a_reduced")?, take(&inputs, "b_reduced")?);
            let mut z_min = Vec::with_capacity(n);
            let mut z_max = Vec::with_capacity(n);
            let mut root_structure = Vec::with_capacity(n);
            let mut iterations = Vec::with_capacity(n);
            let mut converged = Vec::with_capacity(n);
            let mut residual = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::srk_z_factor(a[i], b[i]), &mut warnings)?;
                z_min.push(r.z_min);
                z_max.push(r.z_max);
                root_structure.push(Some(r.root_structure.as_str().to_string()));
                iterations.push(f64::from(r.iterations));
                converged.push(flag(r.converged));
                residual.push(r.residual);
            }
            push_values(&mut columns, "z_min", "dimensionless", z_min);
            push_values(&mut columns, "z_max", "dimensionless", z_max);
            push_labels(&mut columns, "root_structure", root_structure);
            push_values(&mut columns, "iterations", "dimensionless", iterations);
            push_values(&mut columns, "converged", "dimensionless", converged);
            push_values(&mut columns, "residual", "dimensionless", residual);
        }

        "eos.srk_departure" => {
            let (a, b, z, kappa, tr) = (
                take(&inputs, "a_reduced")?,
                take(&inputs, "b_reduced")?,
                take(&inputs, "z")?,
                take(&inputs, "kappa")?,
                take(&inputs, "Tr")?,
            );
            let mut ln_phi = Vec::with_capacity(n);
            let mut h_dep_rt = Vec::with_capacity(n);
            let mut s_dep_r = Vec::with_capacity(n);
            let mut cp_dep_r = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::srk_departure(a[i], b[i], z[i], kappa[i], tr[i]),
                    &mut warnings,
                )?;
                ln_phi.push(r.ln_phi);
                h_dep_rt.push(r.h_dep_rt);
                s_dep_r.push(r.s_dep_r);
                cp_dep_r.push(r.cp_dep_r);
            }
            push_values(&mut columns, "ln_phi", "dimensionless", ln_phi);
            push_values(&mut columns, "h_dep_rt", "dimensionless", h_dep_rt);
            push_values(&mut columns, "s_dep_r", "dimensionless", s_dep_r);
            push_values(&mut columns, "cp_dep_r", "dimensionless", cp_dep_r);
        }

        "eos.rk_alpha_ab" => {
            let (tr, pr) = (take(&inputs, "Tr")?, take(&inputs, "Pr")?);
            let mut alpha = Vec::with_capacity(n);
            let mut a_reduced = Vec::with_capacity(n);
            let mut b_reduced = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::rk_alpha_ab(tr[i], pr[i]), &mut warnings)?;
                alpha.push(r.alpha);
                a_reduced.push(r.a_reduced);
                b_reduced.push(r.b_reduced);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
            push_values(&mut columns, "a_reduced", "dimensionless", a_reduced);
            push_values(&mut columns, "b_reduced", "dimensionless", b_reduced);
        }

        "eos.rk_departure" => {
            let (a, b, z) = (
                take(&inputs, "a_reduced")?,
                take(&inputs, "b_reduced")?,
                take(&inputs, "z")?,
            );
            let mut ln_phi = Vec::with_capacity(n);
            let mut h_dep_rt = Vec::with_capacity(n);
            let mut s_dep_r = Vec::with_capacity(n);
            let mut cp_dep_r = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::rk_departure(a[i], b[i], z[i]), &mut warnings)?;
                ln_phi.push(r.ln_phi);
                h_dep_rt.push(r.h_dep_rt);
                s_dep_r.push(r.s_dep_r);
                cp_dep_r.push(r.cp_dep_r);
            }
            push_values(&mut columns, "ln_phi", "dimensionless", ln_phi);
            push_values(&mut columns, "h_dep_rt", "dimensionless", h_dep_rt);
            push_values(&mut columns, "s_dep_r", "dimensionless", s_dep_r);
            push_values(&mut columns, "cp_dep_r", "dimensionless", cp_dep_r);
        }

        "eos.prsv_kappa" => {
            let (omega, tr, kappa1) = (
                take(&inputs, "omega")?,
                take(&inputs, "Tr")?,
                take(&inputs, "kappa1")?,
            );
            let mut kappa = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::prsv_kappa(omega[i], tr[i], kappa1[i]),
                    &mut warnings,
                )?;
                kappa.push(r.kappa);
            }
            push_values(&mut columns, "kappa", "dimensionless", kappa);
        }

        "eos.pr_departure" => {
            let (a, b, z, kappa, tr) = (
                take(&inputs, "a_reduced")?,
                take(&inputs, "b_reduced")?,
                take(&inputs, "z")?,
                take(&inputs, "kappa")?,
                take(&inputs, "Tr")?,
            );
            let mut ln_phi = Vec::with_capacity(n);
            let mut h_dep_rt = Vec::with_capacity(n);
            let mut s_dep_r = Vec::with_capacity(n);
            let mut cp_dep_r = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::pr_departure(a[i], b[i], z[i], kappa[i], tr[i]),
                    &mut warnings,
                )?;
                ln_phi.push(r.ln_phi);
                h_dep_rt.push(r.h_dep_rt);
                s_dep_r.push(r.s_dep_r);
                cp_dep_r.push(r.cp_dep_r);
            }
            push_values(&mut columns, "ln_phi", "dimensionless", ln_phi);
            push_values(&mut columns, "h_dep_rt", "dimensionless", h_dep_rt);
            push_values(&mut columns, "s_dep_r", "dimensionless", s_dep_r);
            push_values(&mut columns, "cp_dep_r", "dimensionless", cp_dep_r);
        }

        "eos.vdw1f_mix_binary" => {
            let (z1, a1, a2, b1, b2, k12) = (
                take(&inputs, "z1")?,
                take(&inputs, "a1")?,
                take(&inputs, "a2")?,
                take(&inputs, "b1")?,
                take(&inputs, "b2")?,
                take(&inputs, "k12")?,
            );
            let mut a_mix = Vec::with_capacity(n);
            let mut b_mix = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::vdw1f_mix_binary(z1[i], a1[i], a2[i], b1[i], b2[i], k12[i]),
                    &mut warnings,
                )?;
                a_mix.push(r.a_mix);
                b_mix.push(r.b_mix);
            }
            push_values(&mut columns, "a_mix", "dimensionless", a_mix);
            push_values(&mut columns, "b_mix", "dimensionless", b_mix);
        }

        "eos.rachford_rice_binary" => {
            let (z1, k1, k2) = (
                take(&inputs, "z1")?,
                take(&inputs, "K1")?,
                take(&inputs, "K2")?,
            );
            let mut beta = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::rachford_rice_binary(z1[i], k1[i], k2[i]),
                    &mut warnings,
                )?;
                beta.push(r.beta);
            }
            push_values(&mut columns, "beta", "dimensionless", beta);
        }

        "eos.pr_molar_volume" => {
            let (z, t, p) = (
                take(&inputs, "z")?,
                take(&inputs, "T")?,
                take(&inputs, "P")?,
            );
            let mut v = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::pr_molar_volume(z[i], kelvins(t[i]), pascals(p[i])),
                    &mut warnings,
                )?;
                v.push(r.v.value);
            }
            push_values(&mut columns, "v", "m**3/mol", v);
        }

        "eos.ideal_gas_cp" => {
            let (a, b, c, d, e, t) = (
                take(&inputs, "cp_a")?,
                take(&inputs, "cp_b")?,
                take(&inputs, "cp_c")?,
                take(&inputs, "cp_d")?,
                take(&inputs, "cp_e")?,
                take(&inputs, "T")?,
            );
            let mut cp = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::ideal_gas_cp(a[i], b[i], c[i], d[i], e[i], kelvins(t[i])),
                    &mut warnings,
                )?;
                cp.push(r.cp.value);
            }
            push_values(&mut columns, "cp", "J/(mol*K)", cp);
        }

        "eos.pr_mass_density" => {
            let (m, v) = (take(&inputs, "M")?, take(&inputs, "v")?);
            let mut rho = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::pr_mass_density(kilograms_per_mole(m[i]), cubic_meters_per_mole(v[i])),
                    &mut warnings,
                )?;
                rho.push(r.rho.value);
            }
            push_values(&mut columns, "rho", "kg/m**3", rho);
        }

        "eos.pr_peneloux_shift" => {
            let (omega, tc, pc) = (
                take(&inputs, "omega")?,
                take(&inputs, "Tc")?,
                take(&inputs, "Pc")?,
            );
            let mut c = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::pr_peneloux_shift(omega[i], kelvins(tc[i]), pascals(pc[i])),
                    &mut warnings,
                )?;
                c.push(r.c.value);
            }
            push_values(&mut columns, "c", "m**3/mol", c);
        }

        "eos.srk_peneloux_shift" => {
            let (omega, tc, pc) = (
                take(&inputs, "omega")?,
                take(&inputs, "Tc")?,
                take(&inputs, "Pc")?,
            );
            let mut c = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::srk_peneloux_shift(omega[i], kelvins(tc[i]), pascals(pc[i])),
                    &mut warnings,
                )?;
                c.push(r.c.value);
            }
            push_values(&mut columns, "c", "m**3/mol", c);
        }

        "eos.heat_of_vaporization" => {
            let (c0, c1, c2, c3, tc, t) = (
                take(&inputs, "c0")?,
                take(&inputs, "c1")?,
                take(&inputs, "c2")?,
                take(&inputs, "c3")?,
                take(&inputs, "Tc")?,
                take(&inputs, "T")?,
            );
            let mut hov = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::heat_of_vaporization(
                        c0[i],
                        c1[i],
                        c2[i],
                        c3[i],
                        kelvins(tc[i]),
                        kelvins(t[i]),
                    ),
                    &mut warnings,
                )?;
                hov.push(r.hov.value);
            }
            push_values(&mut columns, "hov", "J/mol", hov);
        }

        "eos.liquid_heat_capacity" => {
            let (c0, c1, c2, c3, c4, t) = (
                take(&inputs, "c0")?,
                take(&inputs, "c1")?,
                take(&inputs, "c2")?,
                take(&inputs, "c3")?,
                take(&inputs, "c4")?,
                take(&inputs, "T")?,
            );
            let mut cp = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::liquid_heat_capacity(c0[i], c1[i], c2[i], c3[i], c4[i], kelvins(t[i])),
                    &mut warnings,
                )?;
                cp.push(r.cp.value);
            }
            push_values(&mut columns, "cp", "J/(mol*K)", cp);
        }

        other => {
            return Err(pyo3::exceptions::PyNotImplementedError::new_err(format!(
                "no batch arm for `{other}`"
            )));
        }
    }

    Ok(PyBatchResult { columns, warnings })
}
