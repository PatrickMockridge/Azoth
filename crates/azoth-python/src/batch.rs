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
    cubic_meters_per_mole, cubic_meters_per_second, joules_per_mole_kelvin, kelvin_intervals,
    kelvins, kilograms_per_cubic_meter, kilograms_per_mole, kilograms_per_second,
    kilograms_per_square_meter_second, meters, meters_per_second, moles_per_cubic_meter,
    pascal_seconds, pascals, seconds, square_meters, watts_per_meter_kelvin,
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

        "eos.beggs_brill_segment" => {
            let (gas_flow, liquid_flow) =
                (take(&inputs, "gas_flow")?, take(&inputs, "liquid_flow")?);
            let (gas_density, liquid_density) = (
                take(&inputs, "gas_density")?,
                take(&inputs, "liquid_density")?,
            );
            let (gas_viscosity, liquid_viscosity) = (
                take(&inputs, "gas_viscosity")?,
                take(&inputs, "liquid_viscosity")?,
            );
            let (temperature, pressure) =
                (take(&inputs, "temperature")?, take(&inputs, "pressure")?);
            let (diameter, length) = (take(&inputs, "diameter")?, take(&inputs, "length")?);
            let (elevation, roughness) = (take(&inputs, "elevation")?, take(&inputs, "roughness")?);
            let mut regime = Vec::with_capacity(n);
            let mut holdup = Vec::with_capacity(n);
            let mut v_sg = Vec::with_capacity(n);
            let mut v_sl = Vec::with_capacity(n);
            let mut v_mix = Vec::with_capacity(n);
            let mut rho_ns = Vec::with_capacity(n);
            let mut mu_ns = Vec::with_capacity(n);
            let mut reynolds = Vec::with_capacity(n);
            let mut friction = Vec::with_capacity(n);
            let mut friction_two_phase = Vec::with_capacity(n);
            let mut hydrostatic = Vec::with_capacity(n);
            let mut friction_drop = Vec::with_capacity(n);
            let mut total = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::beggs_brill_segment::beggs_brill_segment(
                        cubic_meters_per_second(gas_flow[i]),
                        cubic_meters_per_second(liquid_flow[i]),
                        kilograms_per_cubic_meter(gas_density[i]),
                        kilograms_per_cubic_meter(liquid_density[i]),
                        pascal_seconds(gas_viscosity[i]),
                        pascal_seconds(liquid_viscosity[i]),
                        kelvins(temperature[i]),
                        pascals(pressure[i]),
                        meters(diameter[i]),
                        meters(length[i]),
                        meters(elevation[i]),
                        meters(roughness[i]),
                    ),
                    &mut warnings,
                )?;
                regime.push(Some(r.regime.as_str().to_string()));
                holdup.push(r.liquid_holdup);
                v_sg.push(r.superficial_gas_velocity.value);
                v_sl.push(r.superficial_liquid_velocity.value);
                v_mix.push(r.superficial_mixture_velocity.value);
                rho_ns.push(r.no_slip_density.value);
                mu_ns.push(r.no_slip_viscosity.value);
                reynolds.push(r.reynolds);
                friction.push(r.friction_factor);
                friction_two_phase.push(r.two_phase_friction_factor);
                hydrostatic.push(r.hydrostatic_drop.value);
                friction_drop.push(r.friction_drop.value);
                total.push(r.total_drop.value);
            }
            push_labels(&mut columns, "regime", regime);
            push_values(&mut columns, "liquid_holdup", "dimensionless", holdup);
            push_values(&mut columns, "superficial_gas_velocity", "m/s", v_sg);
            push_values(&mut columns, "superficial_liquid_velocity", "m/s", v_sl);
            push_values(&mut columns, "superficial_mixture_velocity", "m/s", v_mix);
            push_values(&mut columns, "no_slip_density", "kg/m**3", rho_ns);
            push_values(&mut columns, "no_slip_viscosity", "Pa*s", mu_ns);
            push_values(&mut columns, "reynolds", "dimensionless", reynolds);
            push_values(&mut columns, "friction_factor", "dimensionless", friction);
            push_values(
                &mut columns,
                "two_phase_friction_factor",
                "dimensionless",
                friction_two_phase,
            );
            push_values(&mut columns, "hydrostatic_drop", "Pa", hydrostatic);
            push_values(&mut columns, "friction_drop", "Pa", friction_drop);
            push_values(&mut columns, "total_drop", "Pa", total);
        }

        "eos.de_waard_milliams_corrosion" => {
            let (temperature, co2, ph, h2s, velocity) = (
                take(&inputs, "temperature")?,
                take(&inputs, "co2_partial_pressure")?,
                take(&inputs, "ph")?,
                take(&inputs, "h2s_partial_pressure")?,
                take(&inputs, "flow_velocity")?,
            );
            let (inhibitor, glycol) = (
                take(&inputs, "inhibitor_efficiency")?,
                take(&inputs, "glycol_fraction")?,
            );
            let (sulfur, life) = (
                take(&inputs, "sulfur_deposition_rate")?,
                take(&inputs, "design_life")?,
            );
            let mut baseline = Vec::with_capacity(n);
            let mut ph_factor = Vec::with_capacity(n);
            let mut scale_factor = Vec::with_capacity(n);
            let mut glycol_factor = Vec::with_capacity(n);
            let mut flow_factor = Vec::with_capacity(n);
            let mut inhibitor_factor = Vec::with_capacity(n);
            let mut corrected = Vec::with_capacity(n);
            let mut sulfur_rate = Vec::with_capacity(n);
            let mut fes_rate = Vec::with_capacity(n);
            let mut co2_corrected = Vec::with_capacity(n);
            let mut total = Vec::with_capacity(n);
            let mut feco3_si = Vec::with_capacity(n);
            let mut feco3_forming = Vec::with_capacity(n);
            let mut sour = Vec::with_capacity(n);
            let mut severity = Vec::with_capacity(n);
            let mut morphology = Vec::with_capacity(n);
            let mut allowance = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::de_waard_milliams_corrosion::de_waard_milliams_corrosion(
                        kelvins(temperature[i]),
                        pascals(co2[i]),
                        ph[i],
                        pascals(h2s[i]),
                        meters_per_second(velocity[i]),
                        inhibitor[i],
                        glycol[i],
                        kilograms_per_square_meter_second(sulfur[i]),
                        seconds(life[i]),
                    ),
                    &mut warnings,
                )?;
                baseline.push(r.baseline_rate.value);
                ph_factor.push(r.ph_correction_factor);
                scale_factor.push(r.scale_correction_factor);
                glycol_factor.push(r.glycol_correction_factor);
                flow_factor.push(r.flow_correction_factor);
                inhibitor_factor.push(r.inhibitor_factor);
                corrected.push(r.corrected_rate.value);
                sulfur_rate.push(r.sulfur_rate.value);
                fes_rate.push(r.fes_rate.value);
                co2_corrected.push(r.co2_corrected_rate.value);
                total.push(r.total_corrosion_rate.value);
                feco3_si.push(r.feco3_saturation_index);
                feco3_forming.push(flag(r.feco3_scale_forming));
                sour.push(flag(r.is_sour_service));
                severity.push(Some(r.corrosion_severity.as_str().to_string()));
                morphology.push(Some(r.fes_morphology.as_str().to_string()));
                allowance.push(r.corrosion_allowance.value);
            }
            push_values(&mut columns, "baseline_rate", "m/s", baseline);
            push_values(
                &mut columns,
                "ph_correction_factor",
                "dimensionless",
                ph_factor,
            );
            push_values(
                &mut columns,
                "scale_correction_factor",
                "dimensionless",
                scale_factor,
            );
            push_values(
                &mut columns,
                "glycol_correction_factor",
                "dimensionless",
                glycol_factor,
            );
            push_values(
                &mut columns,
                "flow_correction_factor",
                "dimensionless",
                flow_factor,
            );
            push_values(
                &mut columns,
                "inhibitor_factor",
                "dimensionless",
                inhibitor_factor,
            );
            push_values(&mut columns, "corrected_rate", "m/s", corrected);
            push_values(&mut columns, "sulfur_rate", "m/s", sulfur_rate);
            push_values(&mut columns, "fes_rate", "m/s", fes_rate);
            push_values(&mut columns, "co2_corrected_rate", "m/s", co2_corrected);
            push_values(&mut columns, "total_corrosion_rate", "m/s", total);
            push_values(
                &mut columns,
                "feco3_saturation_index",
                "dimensionless",
                feco3_si,
            );
            push_values(
                &mut columns,
                "feco3_scale_forming",
                "dimensionless",
                feco3_forming,
            );
            push_values(&mut columns, "is_sour_service", "dimensionless", sour);
            push_labels(&mut columns, "corrosion_severity", severity);
            push_labels(&mut columns, "fes_morphology", morphology);
            push_values(&mut columns, "corrosion_allowance", "m", allowance);
        }

        "eos.electrolyte_scale_calculator" => {
            let (calcium, barium, strontium) = (
                take(&inputs, "calcium")?,
                take(&inputs, "barium")?,
                take(&inputs, "strontium")?,
            );
            let (magnesium, sodium, potassium) = (
                take(&inputs, "magnesium")?,
                take(&inputs, "sodium")?,
                take(&inputs, "potassium")?,
            );
            let (iron, chloride, sulphate) = (
                take(&inputs, "iron")?,
                take(&inputs, "chloride")?,
                take(&inputs, "sulphate")?,
            );
            let (bicarbonate, carbonate) =
                (take(&inputs, "bicarbonate")?, take(&inputs, "carbonate")?);
            let (temperature, pressure, ph) = (
                take(&inputs, "temperature")?,
                take(&inputs, "pressure")?,
                take(&inputs, "ph")?,
            );
            let mut ionic_strength = Vec::with_capacity(n);
            let mut debye_hueckel_a = Vec::with_capacity(n);
            let mut divalent = Vec::with_capacity(n);
            let mut monovalent = Vec::with_capacity(n);
            let mut si_caco3 = Vec::with_capacity(n);
            let mut si_baso4 = Vec::with_capacity(n);
            let mut si_caso4 = Vec::with_capacity(n);
            let mut si_srso4 = Vec::with_capacity(n);
            let mut si_brucite = Vec::with_capacity(n);
            let mut hydroxide = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::electrolyte_scale_calculator::electrolyte_scale_calculator(
                        kilograms_per_cubic_meter(calcium[i]),
                        kilograms_per_cubic_meter(barium[i]),
                        kilograms_per_cubic_meter(strontium[i]),
                        kilograms_per_cubic_meter(magnesium[i]),
                        kilograms_per_cubic_meter(sodium[i]),
                        kilograms_per_cubic_meter(potassium[i]),
                        kilograms_per_cubic_meter(iron[i]),
                        kilograms_per_cubic_meter(chloride[i]),
                        kilograms_per_cubic_meter(sulphate[i]),
                        kilograms_per_cubic_meter(bicarbonate[i]),
                        kilograms_per_cubic_meter(carbonate[i]),
                        kelvins(temperature[i]),
                        pascals(pressure[i]),
                        ph[i],
                    ),
                    &mut warnings,
                )?;
                ionic_strength.push(r.ionic_strength.value);
                debye_hueckel_a.push(r.debye_hueckel_a);
                divalent.push(r.divalent_activity_coefficient);
                monovalent.push(r.monovalent_activity_coefficient);
                si_caco3.push(r.si_caco3);
                si_baso4.push(r.si_baso4);
                si_caso4.push(r.si_caso4);
                si_srso4.push(r.si_srso4);
                si_brucite.push(r.si_brucite);
                hydroxide.push(r.hydroxide_activity);
            }
            push_values(&mut columns, "ionic_strength", "mol/kg", ionic_strength);
            push_values(
                &mut columns,
                "debye_hueckel_a",
                "dimensionless",
                debye_hueckel_a,
            );
            push_values(
                &mut columns,
                "divalent_activity_coefficient",
                "dimensionless",
                divalent,
            );
            push_values(
                &mut columns,
                "monovalent_activity_coefficient",
                "dimensionless",
                monovalent,
            );
            push_values(&mut columns, "si_caco3", "dimensionless", si_caco3);
            push_values(&mut columns, "si_baso4", "dimensionless", si_baso4);
            push_values(&mut columns, "si_caso4", "dimensionless", si_caso4);
            push_values(&mut columns, "si_srso4", "dimensionless", si_srso4);
            push_values(&mut columns, "si_brucite", "dimensionless", si_brucite);
            push_values(
                &mut columns,
                "hydroxide_activity",
                "dimensionless",
                hydroxide,
            );
        }

        "eos.scale_kinetics" => {
            let (saturation_index, nucleation_constant, nucleation_shape_factor) = (
                take(&inputs, "saturation_index")?,
                take(&inputs, "nucleation_constant")?,
                take(&inputs, "nucleation_shape_factor")?,
            );
            let (surface_rate_constant, surface_reaction_order) = (
                take(&inputs, "surface_rate_constant")?,
                take(&inputs, "surface_reaction_order")?,
            );
            let (mass_transfer_coefficient, bulk, equilibrium) = (
                take(&inputs, "mass_transfer_coefficient")?,
                take(&inputs, "bulk_concentration")?,
                take(&inputs, "equilibrium_concentration")?,
            );
            let mut supersaturation = Vec::with_capacity(n);
            let mut induction = Vec::with_capacity(n);
            let mut surface = Vec::with_capacity(n);
            let mut transport_limited = Vec::with_capacity(n);
            let mut effective = Vec::with_capacity(n);
            let mut regime = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::scale_kinetics::scale_kinetics(
                        saturation_index[i],
                        nucleation_constant[i],
                        nucleation_shape_factor[i],
                        surface_rate_constant[i],
                        surface_reaction_order[i],
                        mass_transfer_coefficient[i],
                        moles_per_cubic_meter(bulk[i]),
                        moles_per_cubic_meter(equilibrium[i]),
                    ),
                    &mut warnings,
                )?;
                supersaturation.push(r.supersaturation_ratio);
                induction.push(r.induction_time.value);
                surface.push(r.surface_limited_growth.value);
                transport_limited.push(r.transport_limited_growth.value);
                effective.push(r.effective_growth_rate.value);
                regime.push(Some(r.limiting_regime));
            }
            push_values(
                &mut columns,
                "supersaturation_ratio",
                "dimensionless",
                supersaturation,
            );
            push_values(&mut columns, "induction_time", "s", induction);
            push_values(&mut columns, "surface_limited_growth", "m/s", surface);
            push_values(
                &mut columns,
                "transport_limited_growth",
                "m/s",
                transport_limited,
            );
            push_values(&mut columns, "effective_growth_rate", "m/s", effective);
            push_labels(&mut columns, "limiting_regime", regime);
        }

        "eos.iapws_saturation_pressure" => {
            let t_values = take(&inputs, "T")?;
            let mut p_sat = Vec::with_capacity(n);
            for value in &t_values {
                // The batch carries SI magnitudes, so the kelvin conversion the scalar
                // boundary does is not repeated here - the column already is one. An element
                // outside the Region 4 curve raises rather than being dropped: the spec
                // declares both ends `error`, because there the question is a different one.
                let r = element(
                    py,
                    eos::iapws_saturation_pressure(kelvins(*value)),
                    &mut warnings,
                )?;
                p_sat.push(r.p_sat.value);
            }
            push_values(&mut columns, "p_sat", "Pa", p_sat);
        }

        "eos.nitric_sulfuric_acid_vapor_pressure" => {
            let t_values = take(&inputs, "T")?;
            let mut p_water = Vec::with_capacity(n);
            let mut p_nitric_acid = Vec::with_capacity(n);
            let mut p_sulfuric_acid = Vec::with_capacity(n);
            for value in &t_values {
                // The batch carries SI magnitudes, so the kelvin conversion the scalar
                // boundary does is not repeated here - the column already is one.
                let r = element(
                    py,
                    eos::nitric_sulfuric_acid_vapor_pressure(kelvins(*value)),
                    &mut warnings,
                )?;
                p_water.push(r.p_water.value);
                p_nitric_acid.push(r.p_nitric_acid.value);
                p_sulfuric_acid.push(r.p_sulfuric_acid.value);
            }
            push_values(&mut columns, "p_water", "Pa", p_water);
            push_values(&mut columns, "p_nitric_acid", "Pa", p_nitric_acid);
            push_values(&mut columns, "p_sulfuric_acid", "Pa", p_sulfuric_acid);
        }

        "eos.pr_lee_kesler_alpha" => {
            let (omega, tr) = (take(&inputs, "omega")?, take(&inputs, "Tr")?);
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::pr_lee_kesler_alpha(omega[i], tr[i]), &mut warnings)?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.matcop_alpha" => {
            let (mc1, mc2, mc3, tr) = (
                take(&inputs, "mc1")?,
                take(&inputs, "mc2")?,
                take(&inputs, "mc3")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::matcop_alpha(mc1[i], mc2[i], mc3[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.matcop_pr_alpha" => {
            let (omega, mc1, mc2, mc3, tr) = (
                take(&inputs, "omega")?,
                take(&inputs, "mc1")?,
                take(&inputs, "mc2")?,
                take(&inputs, "mc3")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::matcop_pr_alpha(omega[i], mc1[i], mc2[i], mc3[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.matcop_prumr_alpha" => {
            let (omega, mc1, mc2, mc3, tr) = (
                take(&inputs, "omega")?,
                take(&inputs, "mc1")?,
                take(&inputs, "mc2")?,
                take(&inputs, "mc3")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::matcop_prumr_alpha(omega[i], mc1[i], mc2[i], mc3[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.matcop_prumr_new_alpha" => {
            let (omega, mc1, mc2, mc3, mc4, mc5, tr) = (
                take(&inputs, "omega")?,
                take(&inputs, "mc1")?,
                take(&inputs, "mc2")?,
                take(&inputs, "mc3")?,
                take(&inputs, "mc4")?,
                take(&inputs, "mc5")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::matcop_prumr_new_alpha(
                        omega[i], mc1[i], mc2[i], mc3[i], mc4[i], mc5[i], tr[i],
                    ),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.matcop5_prumr_alpha" => {
            let (omega, mc1, mc2, mc3, mc4, mc5, tr) = (
                take(&inputs, "omega")?,
                take(&inputs, "mc1")?,
                take(&inputs, "mc2")?,
                take(&inputs, "mc3")?,
                take(&inputs, "mc4")?,
                take(&inputs, "mc5")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::matcop5_prumr_alpha(
                        omega[i], mc1[i], mc2[i], mc3[i], mc4[i], mc5[i], tr[i],
                    ),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.mollerup_alpha" => {
            let (p1, p2, p3, tr) = (
                take(&inputs, "p1")?,
                take(&inputs, "p2")?,
                take(&inputs, "p3")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::mollerup_alpha(p1[i], p2[i], p3[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.pr_danesh_alpha" => {
            let (omega, tr) = (take(&inputs, "omega")?, take(&inputs, "Tr")?);
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::pr_danesh_alpha(omega[i], tr[i]), &mut warnings)?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.pr_delft1998_alpha" => {
            let (omega, tr) = (take(&inputs, "omega")?, take(&inputs, "Tr")?);
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::pr_delft1998_alpha(omega[i], tr[i]), &mut warnings)?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.pr_gassem2001_alpha" => {
            let (omega, tr) = (take(&inputs, "omega")?, take(&inputs, "Tr")?);
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::pr_gassem2001_alpha(omega[i], tr[i]), &mut warnings)?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
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

        "eos.tbp_fraction_properties" => {
            let (molar_mass, density) = (take(&inputs, "molar_mass")?, take(&inputs, "density")?);
            let mut tc = Vec::with_capacity(n);
            let mut pc = Vec::with_capacity(n);
            let mut boiling_temperature = Vec::with_capacity(n);
            let mut acentric_factor = Vec::with_capacity(n);
            let mut attraction_exponent = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::tbp_fraction_properties(molar_mass[i], density[i]),
                    &mut warnings,
                )?;
                tc.push(r.tc.value);
                pc.push(r.pc.value);
                boiling_temperature.push(r.boiling_temperature.value);
                acentric_factor.push(r.acentric_factor);
                attraction_exponent.push(r.attraction_exponent);
            }
            push_values(&mut columns, "tc", "K", tc);
            push_values(&mut columns, "pc", "Pa", pc);
            push_values(
                &mut columns,
                "boiling_temperature",
                "K",
                boiling_temperature,
            );
            push_values(
                &mut columns,
                "acentric_factor",
                "dimensionless",
                acentric_factor,
            );
            push_values(
                &mut columns,
                "attraction_exponent",
                "dimensionless",
                attraction_exponent,
            );
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

        "eos.twucoon_alpha" => {
            let (omega, tr) = (take(&inputs, "omega")?, take(&inputs, "Tr")?);
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::twucoon_alpha(omega[i], tr[i]), &mut warnings)?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.twucoon_param_alpha" => {
            let (a, b, c, tr) = (
                take(&inputs, "a")?,
                take(&inputs, "b")?,
                take(&inputs, "c")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::twucoon_param_alpha(a[i], b[i], c[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.twucoon_statoil_alpha" => {
            let (a, b, c, tr) = (
                take(&inputs, "a")?,
                take(&inputs, "b")?,
                take(&inputs, "c")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::twucoon_statoil_alpha(a[i], b[i], c[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.schwartzentruber_alpha" => {
            let (omega, p1, p2, p3, tr) = (
                take(&inputs, "omega")?,
                take(&inputs, "p1")?,
                take(&inputs, "p2")?,
                take(&inputs, "p3")?,
                take(&inputs, "Tr")?,
            );
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::schwartzentruber_alpha(omega[i], p1[i], p2[i], p3[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.soreide_whitson_alpha" => {
            let (salinity, tr) = (take(&inputs, "salinity")?, take(&inputs, "Tr")?);
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::soreide_whitson_alpha(salinity[i], tr[i]),
                    &mut warnings,
                )?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
        }

        "eos.umrpr_alpha" => {
            let (omega, tr) = (take(&inputs, "omega")?, take(&inputs, "Tr")?);
            let mut alpha = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(py, eos::umrpr_alpha(omega[i], tr[i]), &mut warnings)?;
                alpha.push(r.alpha);
            }
            push_values(&mut columns, "alpha", "dimensionless", alpha);
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
                beta.push(r.vapour_fraction);
            }
            push_values(&mut columns, "vapour_fraction", "dimensionless", beta);
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
            // The Rackett compressibility is an **optional** input, so a sweep that does
            // not carry the column is a sweep over the fallback correlation - which is
            // NeqSim's own reading of a zero column, and not a missing value.
            let z_ra = take_optional(&inputs, "z_ra");
            let mut c = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::pr_peneloux_shift(
                        omega[i],
                        kelvins(tc[i]),
                        pascals(pc[i]),
                        z_ra.as_ref().map(|column| column[i]),
                    ),
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
            let z_ra = take_optional(&inputs, "z_ra");
            let mut c = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::srk_peneloux_shift(
                        omega[i],
                        kelvins(tc[i]),
                        pascals(pc[i]),
                        z_ra.as_ref().map(|values| values[i]),
                    ),
                    &mut warnings,
                )?;
                c.push(r.c.value);
            }
            push_values(&mut columns, "c", "m**3/mol", c);
        }

        "eos.costald_molar_volume" => {
            let (omega, tc, vc, m, rho, t) = (
                take(&inputs, "omega")?,
                take(&inputs, "Tc")?,
                take(&inputs, "Vc")?,
                take(&inputs, "M")?,
                take(&inputs, "rho_normal")?,
                take(&inputs, "T")?,
            );
            let mut v = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::costald_molar_volume(
                        omega[i],
                        kelvins(tc[i]),
                        cubic_meters_per_mole(vc[i]),
                        kilograms_per_mole(m[i]),
                        kilograms_per_cubic_meter(rho[i]),
                        kelvins(t[i]),
                    ),
                    &mut warnings,
                )?;
                v.push(r.v.value);
            }
            push_values(&mut columns, "v", "m**3/mol", v);
        }

        "eos.chung_viscosity" => {
            let (omega, tc, vc, m, dipole, kappa, t, v) = (
                take(&inputs, "omega")?,
                take(&inputs, "Tc")?,
                take(&inputs, "Vc")?,
                take(&inputs, "M")?,
                take(&inputs, "dipole")?,
                take(&inputs, "kappa")?,
                take(&inputs, "T")?,
                take(&inputs, "V")?,
            );
            let mut mu = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::chung_viscosity(
                        omega[i],
                        kelvins(tc[i]),
                        cubic_meters_per_mole(vc[i]),
                        kilograms_per_mole(m[i]),
                        dipole[i],
                        kappa[i],
                        kelvins(t[i]),
                        cubic_meters_per_mole(v[i]),
                    ),
                    &mut warnings,
                )?;
                mu.push(r.mu.value);
            }
            push_values(&mut columns, "mu", "Pa*s", mu);
        }

        "eos.chung_conductivity" => {
            let (cv0, m, omega, tc, vc, dipole, kappa, t) = (
                take(&inputs, "Cv0")?,
                take(&inputs, "M")?,
                take(&inputs, "omega")?,
                take(&inputs, "Tc")?,
                take(&inputs, "Vc")?,
                take(&inputs, "dipole")?,
                take(&inputs, "kappa")?,
                take(&inputs, "T")?,
            );
            let mut k = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::chung_conductivity(
                        joules_per_mole_kelvin(cv0[i]),
                        kilograms_per_mole(m[i]),
                        omega[i],
                        kelvins(tc[i]),
                        cubic_meters_per_mole(vc[i]),
                        dipole[i],
                        kappa[i],
                        kelvins(t[i]),
                    ),
                    &mut warnings,
                )?;
                k.push(r.k.value);
            }
            push_values(&mut columns, "k", "W/(m*K)", k);
        }

        "eos.tyn_calus_diffusivity" => {
            let (va, vb, t, eta) = (
                take(&inputs, "VA")?,
                take(&inputs, "VB")?,
                take(&inputs, "T")?,
                take(&inputs, "eta")?,
            );
            let mut d = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::tyn_calus_diffusivity(
                        cubic_meters_per_mole(va[i]),
                        cubic_meters_per_mole(vb[i]),
                        kelvins(t[i]),
                        pascal_seconds(eta[i]),
                    ),
                    &mut warnings,
                )?;
                d.push(r.d.value);
            }
            push_values(&mut columns, "d", "m**2/s", d);
        }

        "eos.chapman_enskog_diffusivity" => {
            let (ma, mb, sigma, eps, temperature, p) = (
                take(&inputs, "MA")?,
                take(&inputs, "MB")?,
                take(&inputs, "sigma")?,
                take(&inputs, "eps")?,
                take(&inputs, "T")?,
                take(&inputs, "P")?,
            );
            let mut d = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::chapman_enskog_diffusivity(
                        kilograms_per_mole(ma[i]),
                        kilograms_per_mole(mb[i]),
                        sigma[i],
                        kelvins(eps[i]),
                        kelvins(temperature[i]),
                        pascals(p[i]),
                    ),
                    &mut warnings,
                )?;
                d.push(r.d.value);
            }
            push_values(&mut columns, "d", "m**2/s", d);
        }

        "eos.fuller_schettler_giddings_diffusivity" => {
            let (ma, mb, va, vb, t, p) = (
                take(&inputs, "MA")?,
                take(&inputs, "MB")?,
                take(&inputs, "VA")?,
                take(&inputs, "VB")?,
                take(&inputs, "T")?,
                take(&inputs, "P")?,
            );
            let mut d = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::fuller_schettler_giddings_diffusivity(
                        kilograms_per_mole(ma[i]),
                        kilograms_per_mole(mb[i]),
                        cubic_meters_per_mole(va[i]),
                        cubic_meters_per_mole(vb[i]),
                        kelvins(t[i]),
                        pascals(p[i]),
                    ),
                    &mut warnings,
                )?;
                d.push(r.d.value);
            }
            push_values(&mut columns, "d", "m**2/s", d);
        }

        "eos.wilke_chang_diffusivity" => {
            let (phi, m, t, eta, va) = (
                take(&inputs, "phi")?,
                take(&inputs, "M")?,
                take(&inputs, "T")?,
                take(&inputs, "eta")?,
                take(&inputs, "VA")?,
            );
            let mut d = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::wilke_chang_diffusivity(
                        phi[i],
                        kilograms_per_mole(m[i]),
                        kelvins(t[i]),
                        pascal_seconds(eta[i]),
                        cubic_meters_per_mole(va[i]),
                    ),
                    &mut warnings,
                )?;
                d.push(r.d.value);
            }
            push_values(&mut columns, "d", "m**2/s", d);
        }

        "eos.co2_water_diffusivity" => {
            let t = take(&inputs, "T")?;
            let mut d = Vec::with_capacity(n);
            for &ti in &t {
                let r = element(py, eos::co2_water_diffusivity(kelvins(ti)), &mut warnings)?;
                d.push(r.d.value);
            }
            push_values(&mut columns, "d", "m**2/s", d);
        }

        "eos.parachor_surface_tension" => {
            let (parachor, rho_l, rho_v, m) = (
                take(&inputs, "parachor")?,
                take(&inputs, "rho_l")?,
                take(&inputs, "rho_v")?,
                take(&inputs, "M")?,
            );
            let mut sigma = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::parachor_surface_tension(
                        parachor[i],
                        kilograms_per_cubic_meter(rho_l[i]),
                        kilograms_per_cubic_meter(rho_v[i]),
                        kilograms_per_mole(m[i]),
                    ),
                    &mut warnings,
                )?;
                sigma.push(r.sigma.value);
            }
            push_values(&mut columns, "sigma", "N/m", sigma);
        }

        "eos.firozabadi_ramley_interface_tension" => {
            let (t, pseudo, rho_1, rho_2) = (
                take(&inputs, "T")?,
                take(&inputs, "pseudo_critical_temperature")?,
                take(&inputs, "rho_1")?,
                take(&inputs, "rho_2")?,
            );
            let mut sigma = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::firozabadi_ramley_interface_tension(
                        kelvins(t[i]),
                        kelvins(pseudo[i]),
                        kilograms_per_cubic_meter(rho_1[i]),
                        kilograms_per_cubic_meter(rho_2[i]),
                    ),
                    &mut warnings,
                )?;
                sigma.push(r.sigma.value);
            }
            push_values(&mut columns, "sigma", "N/m", sigma);
        }

        "eos.rackett_molar_volume" => {
            let (omega, tc, pc, t) = (
                take(&inputs, "omega")?,
                take(&inputs, "Tc")?,
                take(&inputs, "Pc")?,
                take(&inputs, "T")?,
            );
            let mut v = Vec::with_capacity(n);
            for i in 0..n {
                let r = element(
                    py,
                    eos::rackett_molar_volume(
                        omega[i],
                        kelvins(tc[i]),
                        pascals(pc[i]),
                        kelvins(t[i]),
                    ),
                    &mut warnings,
                )?;
                v.push(r.v.value);
            }
            push_values(&mut columns, "v", "m**3/mol", v);
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
