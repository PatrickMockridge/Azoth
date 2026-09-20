//! `eos.freezing_point` - the temperature at which a pure substance's solid appears.
//!
//! ```text
//! solve  (g_fluid(T, P) - g_solid(T, P)) / (R T) = 0   for T
//! ```
//!
//! NeqSim's `FreezingPointTemperatureFlash`, on the route that uses a solid reference
//! equation. The other route - a tabulated solid built from COMP.csv's melting point, heat
//! of fusion and densities - is carried rather than ported: it is a second solid model, and
//! the substance this one exists for is para-hydrogen, whose solid `thermo/util/solid/` gives
//! an equation.
//!
//! # The calibration, which is not in the solid
//!
//! NeqSim's `ParaHydrogenSolidHelmholtzEquation` is the **raw** thesis equation, and the
//! reference it is reported on is fixed by its *system*:
//! `SystemLeachmanEos.createCalibratedParaHydrogenSolidEquation` shifts it to the
//! para-Leachman liquid at the triple point,
//!
//! ```text
//! gibbsShift   = g_liquid(T_tp, P_tp) - g_raw_solid(T_tp, P_tp)
//! entropyShift = s_liquid(T_tp, P_tp) - s_raw_solid(T_tp, P_tp) - H_fus / T_tp
//! ```
//!
//! so that solid and liquid meet there - which is what a triple point *is*, and what the raw
//! equation alone does not do. This module applies it, because it is the thing that compares
//! the two phases; `eos.parahydrogen_solid_phase` stays raw and says so.
//!
//! # Which root of the fluid
//!
//! The flash re-initialises its fluid phase as a gas below the triple-point pressure and a
//! liquid above, and takes the root that follows. That is what [`FluidRoot`] reproduces: a
//! pure substance's isotherm crosses a pressure three times, and the two roots are different
//! states at the same temperature and pressure.

use azoth_core::units::{Pressure, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::leachman::{self, HydrogenType};
use crate::model_gen;
use crate::parahydrogen_solid;
use crate::results::FreezingPointResult;

/// NeqSim's `ThermodynamicConstantsInterface.R`, which its residual divides by. **Not** the
/// equation's own `R`: the two reference equations here carry their own rounded constants, and
/// this is the third.
pub const R: f64 = 8.314_462_1;

/// The tolerance on the dimensionless residual, and on the bracket.
const RESIDUAL_TOLERANCE: f64 = 1.0e-10;
const BRACKET_TOLERANCE: f64 = 1.0e-12;

/// How far the search may walk, in kelvin, and how many steps it may take.
const MINIMUM_TEMPERATURE: f64 = 0.5;
const MAXIMUM_TEMPERATURE: f64 = 5000.0;
const SPAN_GROWTH: f64 = 1.8;
const MAXIMUM_BRACKET_STEPS: usize = 50;
const MAXIMUM_SOLVER_STEPS: usize = 100;

/// The solid's reference shift at the triple point, in J/mol and J/(mol.K).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Calibration {
    /// Added to the raw solid's Gibbs energy, and carried with `T - T_tp` by the entropy
    /// shift.
    pub gibbs_shift: f64,
    /// What the raw solid's entropy is short of the liquid's at the triple point.
    pub entropy_shift: f64,
}

/// Which root of the fluid's isotherm the equilibrium is with.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum FluidRoot {
    /// Below the triple-point pressure.
    Gas,
    /// At and above it.
    Liquid,
}

/// The calibration NeqSim's `SystemLeachmanEos` applies, computed the way it computes it.
///
/// The liquid is `para`-hydrogen's dense root at the triple point - the fluid the shift is
/// measured against - and the solid is this crate's raw equation at the same state. The two
/// constants are measurements of those two models rather than tunables: at the triple point
/// they come out `-3.122683` J/mol and a positive entropy shift whose size is fixed by the
/// fusion enthalpy.
#[must_use]
pub fn calibration() -> Calibration {
    let t_tp = parahydrogen_solid::TRIPLE_POINT_TEMPERATURE;
    let p_tp = parahydrogen_solid::TRIPLE_POINT_PRESSURE;

    let rho = leachman::solve_density_dense(t_tp, p_tp, HydrogenType::Para)
        .expect("the para-hydrogen liquid root exists at the triple point");
    let liquid = leachman::properties(t_tp, rho, HydrogenType::Para);
    let raw_solid = parahydrogen_solid::properties(t_tp, p_tp);

    Calibration {
        gibbs_shift: liquid.g - raw_solid.g,
        entropy_shift: liquid.s
            - raw_solid.s
            - parahydrogen_solid::TRIPLE_POINT_ENTHALPY_OF_FUSION / t_tp,
    }
}

/// The solid's molar Gibbs energy at a state, on the reference the liquid is on.
#[must_use]
pub fn calibrated_solid_gibbs(t: f64, p: f64, calibration: &Calibration) -> f64 {
    let raw = parahydrogen_solid::properties(t, p);
    // `g = a + Pv` gains the shift directly, and the entropy shift's own temperature
    // dependence comes with it: a shift of `-S (T - T_tp) + G` in the Helmholtz energy is
    // `-S (T - T_tp) + G` in the Gibbs energy too, since the shift carries no volume.
    raw.g - calibration.entropy_shift * (t - parahydrogen_solid::TRIPLE_POINT_TEMPERATURE)
        + calibration.gibbs_shift
}

/// The dimensionless residual the solve drives to zero, and the one NeqSim's
/// `calculateLogEquilibriumResidual` returns on this route.
///
/// # Errors
/// * [`AzothError::OutOfRange`] if the fluid's dense root does not exist at `t`.
pub fn residual(t: f64, p: f64, root: FluidRoot, calibration: &Calibration) -> Result<f64> {
    let rho = match root {
        FluidRoot::Gas => leachman::solve_density(t, p, HydrogenType::Para),
        FluidRoot::Liquid => {
            leachman::solve_density_dense(t, p, HydrogenType::Para).ok_or_else(|| {
                AzothError::out_of_range(
                    "P",
                    p,
                    "the freezing search asks for the fluid's dense root, and this pressure \
                     has none at this temperature within the range the equation is fitted to",
                )
            })?
        }
    };
    let fluid = leachman::properties(t, rho, HydrogenType::Para);
    let solid = calibrated_solid_gibbs(t, p, calibration);
    Ok((fluid.g - solid) / (R * t))
}

/// The freezing-point temperature of para-hydrogen at a pressure.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if no component is named, or one is named that has no solid
///   equation here.
/// * [`AzothError::OutOfRange`] if `P` is not positive.
/// * [`AzothError::SolverNotConverged`] if the residual will not bracket or the bracket
///   collapses without reaching the tolerance.
pub fn freezing_point(components: &[String], p: Pressure) -> Result<FreezingPointResult> {
    let spec = &model_gen::FREEZING_POINT_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    if components.len() != 1 {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "{} entries; a freezing point is a pure substance's, and this crate has a \
                 solid equation for para-hydrogen alone",
                components.len()
            ),
        ));
    }
    if !matches!(
        components[0].trim().to_lowercase().as_str(),
        "para-hydrogen"
    ) {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "`{}`; `thermo/util/solid/` gives an equation for para-hydrogen (and argon, \
                 which is not a hydrogen phase), so that is the name this takes",
                components[0]
            ),
        ));
    }

    let algorithm = crate::algorithm_of(spec)?;
    let p_pa = p.value;
    let root = if p_pa < parahydrogen_solid::TRIPLE_POINT_PRESSURE {
        FluidRoot::Gas
    } else {
        FluidRoot::Liquid
    };
    let calibration = calibration();

    // The search starts where the algorithm says and expands a span either side of it, which
    // is NeqSim's own walk: the answer is a property of the substance and the pressure, so
    // the start is a path and not a state.
    let start = algorithm.initial_temperature.unwrap_or(14.0);
    let start_residual = residual(start, p_pa, root, &calibration)?;
    if start_residual.abs() < RESIDUAL_TOLERANCE {
        return Ok(FreezingPointResult {
            temperature: kelvins(start),
            iterations: 1,
            residual: start_residual,
            warnings,
        });
    }

    let mut span = (start * 0.01).max(0.1);
    let (mut lo, mut lo_residual) = (start, start_residual);
    // The bracket's upper side is only read through its residual's sign, which the walk
    // refreshes below; it starts at the same point as the lower one because a bracket needs
    // two sides and either can move first.
    let (mut hi, mut hi_residual) = (start, start_residual);
    let mut iterations = 1;
    let mut bracketed = false;
    for _ in 0..MAXIMUM_BRACKET_STEPS {
        let candidate_lo = (start - span).max(MINIMUM_TEMPERATURE);
        let candidate_hi = (start + span).min(MAXIMUM_TEMPERATURE);
        // A trial the equation cannot evaluate is skipped rather than fatal, which is
        // NeqSim's `evaluateResidualIfValid`: near the triple point one side of a span can
        // ask for a dense root that does not exist.
        if let Ok(value) = residual(candidate_lo, p_pa, root, &calibration) {
            lo = candidate_lo;
            lo_residual = value;
            iterations += 1;
        }
        if let Ok(value) = residual(candidate_hi, p_pa, root, &calibration) {
            hi = candidate_hi;
            hi_residual = value;
            iterations += 1;
        }
        if lo_residual * hi_residual <= 0.0 {
            bracketed = true;
            break;
        }
        span *= SPAN_GROWTH;
    }
    if !bracketed {
        return Err(AzothError::SolverNotConverged {
            iterations,
            residual: start_residual,
            tolerance: algorithm.tolerance,
        });
    }

    for _ in 0..MAXIMUM_SOLVER_STEPS {
        let mid = 0.5 * (lo + hi);
        let value = residual(mid, p_pa, root, &calibration)?;
        iterations += 1;
        if value.abs() < RESIDUAL_TOLERANCE || hi - lo < BRACKET_TOLERANCE {
            return Ok(FreezingPointResult {
                temperature: kelvins(mid),
                iterations,
                residual: value,
                warnings,
            });
        }
        if lo_residual * value <= 0.0 {
            // The upper side moves; its residual is not stored because the next step
            // recomputes one at the midpoint and compares against `lo_residual` alone.
            hi = mid;
        } else {
            lo = mid;
            lo_residual = value;
        }
    }

    Err(AzothError::SolverNotConverged {
        iterations,
        residual: 0.5 * (lo + hi),
        tolerance: algorithm.tolerance,
    })
}
