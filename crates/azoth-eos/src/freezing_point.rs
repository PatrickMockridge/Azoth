//! `eos.freezing_point` - the temperature at which a fluid's solid-forming substance appears.
//!
//! NeqSim's `FreezingPointTemperatureFlash`, both of the routes it has, chosen by an
//! `instanceof` upstream and by the candidate's name here:
//!
//! - **[`SolidRoute::Helmholtz`]** - the substance has a solid reference equation
//!   (`thermo/util/solid/`), and the residual is the molar Gibbs difference:
//!
//!   ```text
//!   solve  (g_fluid(T, P) - g_solid(T, P)) / (R T) = 0   for T
//!   ```
//!
//! - **[`SolidRoute::Tabulated`]** - every other substance, against the tabulated solid built
//!   from COMP.csv's melting point, heat of fusion, heat capacities and densities. Its residual
//!   is [`tabulated_residual`]'s multiphase appearance condition, which is a log-sum-exp over
//!   the fluid's own phases rather than a difference of two Gibbs energies.
//!
//! Both are wrapped by the same [`solve`], which is NeqSim's own walk: a span either side of a
//! starting temperature, expanded until it brackets a sign change, then bisected.
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

use azoth_core::units::{Pressure, kelvins, pascals};
use azoth_core::{AzothError, Result, apply_checks};

use crate::leachman::{self, HydrogenType};
use crate::mixture::Mixture;
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

/// Which of NeqSim's two solid routes a substance is solved on.
///
/// NeqSim dispatches on the *system*: `FreezingPointTemperatureFlash` asks whether its solid
/// phase is a `PhaseSolidHelmholtzEos`, which `SystemLeachmanEos` and `SystemSolidHelmholtzEos`
/// build and every other system does not. The set of substances with a Helmholtz solid here is
/// the one `thermo/util/solid/` gives an equation for, so this is derived from the name rather
/// than declared - and it agrees with NeqSim's because NeqSim's own systems are what carry those
/// equations.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SolidRoute {
    /// A solid reference equation: the calibrated para-hydrogen one, whose residual is the
    /// molar Gibbs difference.
    Helmholtz,
    /// The tabulated solid: `ComponentSolid.fugcoef2` over the databank's melting point, heat
    /// of fusion, heat capacities and density correlations, and the residual is the multiphase
    /// appearance condition over the fluid's own phases.
    Tabulated,
}

/// The route a named substance is solved on.
#[must_use]
pub fn route_for(solid: &str) -> SolidRoute {
    if matches!(normalised(solid).as_str(), "para-hydrogen") {
        SolidRoute::Helmholtz
    } else {
        SolidRoute::Tabulated
    }
}

/// A name as the databank spells it: trimmed and lower-cased.
fn normalised(name: &str) -> String {
    name.trim().to_lowercase()
}

/// The dimensionless residual the tabulated route drives to zero, for one candidate.
///
/// NeqSim's `calculateLogEquilibriumResidual`, outside its Helmholtz branch:
///
/// ```text
/// residual = ln z_k - ln( sum_p beta_p phi_solid_k / phi_kp )
/// ```
///
/// **It is the multiphase appearance condition, and the log-sum-exp is why it is written that
/// way.** A fluid can be several phases at once and the solid meets all of them, so the sum is
/// over the phases the flash found, weighted by their amounts. NeqSim takes the maximum
/// logarithm out before exponentiating, which is not a rearrangement for its own sake: the
/// terms are ratios of fugacity coefficients across twenty orders of magnitude, and a direct
/// sum underflows to zero and makes the logarithm infinite.
///
/// # Errors
/// * [`AzothError`] from the flash, from the solid coefficient, or from a candidate the fluid
///   does not contain.
fn tabulated_residual(
    mixture: &Mixture,
    z: &[f64],
    index: usize,
    solid: &str,
    eos: &str,
    t: f64,
    p: f64,
) -> Result<f64> {
    // A `NaN` fraction is as much a refusal as a zero one, and `<=` alone would let it past.
    if z[index] <= 0.0 || !z[index].is_finite() {
        return Err(AzothError::invalid_input(
            "solid",
            format!(
                "`{solid}` has no overall mole fraction, and a solid cannot appear from a \
                 substance the feed does not have"
            ),
        ));
    }
    let flash = crate::pt_flash::pt_flash(mixture, kelvins(t), pascals(p), z)?;
    // **Methane never freezes in NeqSim.** `ComponentSolid.fugcoef` opens with
    // `if (componentName.equals("methane")) return 1e30;`, and that is the entry point the
    // solid phase calls - so a methane candidate's solid is infinitely volatile, its residual
    // has no sign change, and the search refuses rather than reporting a freezing point. This
    // is that guard.
    let phi_solid = if normalised(solid) == "methane" {
        crate::tp_solid_flash::METHANE_NEVER_FREEZES
    } else {
        crate::tp_solid_flash::tabulated_solid_fugacity(
            mixture,
            index,
            solid,
            kelvins(t),
            pascals(p),
            eos,
        )?
    };
    let ln_phi_solid = phi_solid.ln();

    // The fluid's phases, as `(amount, ln phi of this component)`. A two-phase flash has two;
    // a single-phase one has one, carrying the whole feed.
    let mut phases: Vec<(f64, f64)> = Vec::new();
    match flash.phase {
        crate::results::Phase::TwoPhase => {
            let beta = flash.beta.unwrap_or(0.0);
            phases.push((1.0 - beta, flash.ln_phi_liquid[index]));
            phases.push((beta, flash.ln_phi_vapour[index]));
        }
        crate::results::Phase::AllLiquid => phases.push((1.0, flash.ln_phi_liquid[index])),
        _ => phases.push((1.0, flash.ln_phi_vapour[index])),
    }

    let terms: Vec<f64> = phases
        .iter()
        .filter(|(amount, _)| *amount > 0.0)
        .map(|(amount, ln_phi)| amount.ln() + ln_phi_solid - ln_phi)
        .collect();
    let maximum = terms.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !maximum.is_finite() {
        return Err(AzothError::out_of_range(
            "P",
            p,
            "no phase of the flashed fluid has a finite fugacity contribution at this state",
        ));
    }
    let scaled: f64 = terms.iter().map(|term| (term - maximum).exp()).sum();
    Ok(z[index].ln() - maximum - scaled.ln())
}

/// Where a solve starts, and how far it may walk.
#[derive(Debug, Clone, Copy)]
struct Search {
    start: f64,
}

/// The temperature a residual crosses zero at, by NeqSim's own walk.
///
/// The search starts where the algorithm says and expands a span either side of it: the answer
/// is a property of the substance and the pressure, so the start is a path and not a state. A
/// trial the equation cannot evaluate is skipped rather than fatal, which is NeqSim's
/// `evaluateResidualIfValid` - near a triple point one side of a span can ask for a root that
/// does not exist.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if no span brackets a sign change, or the bracket
///   collapses without reaching the tolerance.
fn solve<F>(residual: F, search: Search, tolerance: f64) -> Result<(f64, f64, u32)>
where
    F: Fn(f64) -> Result<f64>,
{
    let start = search.start;
    let start_residual = residual(start)?;
    if start_residual.abs() < RESIDUAL_TOLERANCE {
        return Ok((start, start_residual, 1));
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
        if let Ok(value) = residual(candidate_lo) {
            lo = candidate_lo;
            lo_residual = value;
            iterations += 1;
        }
        if let Ok(value) = residual(candidate_hi) {
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
            tolerance,
        });
    }

    for _ in 0..MAXIMUM_SOLVER_STEPS {
        let mid = 0.5 * (lo + hi);
        let value = residual(mid)?;
        iterations += 1;
        if value.abs() < RESIDUAL_TOLERANCE || hi - lo < BRACKET_TOLERANCE {
            return Ok((mid, value, iterations));
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
        tolerance,
    })
}

/// The freezing-point temperature of a fluid's solid-forming substance.
///
/// Two routes, as NeqSim has two: para-hydrogen is solved against its calibrated solid
/// reference equation, and every other substance against the tabulated solid built from its
/// melting point, heat of fusion and densities.
///
/// **A fluid may have more than one candidate and the highest freezing point is the answer.**
/// NeqSim loops the components the caller enabled a solid check for and keeps the maximum,
/// because a fluid at a temperature where any one of its substances freezes has frozen; the
/// result names the one that set it. `solid` is that set, and naming one substance is the
/// single-candidate case.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` does not match `components`, if `solid` is not one of
///   them, if `P` is not positive, or if the tabulated route's candidate carries no melt data.
/// * [`AzothError::OutOfRange`] from a trial's state.
/// * [`AzothError::SolverNotConverged`] if no candidate's residual brackets a sign change.
pub fn freezing_point(
    components: &[String],
    z: &[f64],
    solid: &str,
    p: Pressure,
) -> Result<FreezingPointResult> {
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

    if components.len() != z.len() {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "{} components and {} mole fractions; a candidate's own fraction is what \
                 decides whether it can appear",
                components.len(),
                z.len()
            ),
        ));
    }
    let names: Vec<String> = components.iter().map(|name| normalised(name)).collect();
    let index = names
        .iter()
        .position(|name| *name == normalised(solid))
        .ok_or_else(|| {
            AzothError::invalid_input(
                "solid",
                format!(
                    "`{solid}` is not one of the fluid's components, so it is not a substance \
                     whose freezing point this can solve for"
                ),
            )
        })?;

    let algorithm = crate::algorithm_of(spec)?;
    let search = Search {
        start: algorithm.initial_temperature.unwrap_or(14.0),
    };
    let p_pa = p.value;
    let eos = "srk";

    let (temperature, residual, iterations) = match route_for(solid) {
        SolidRoute::Helmholtz => {
            let root = if p_pa < parahydrogen_solid::TRIPLE_POINT_PRESSURE {
                FluidRoot::Gas
            } else {
                FluidRoot::Liquid
            };
            let calibration = calibration();
            solve(
                |t| residual(t, p_pa, root, &calibration),
                search,
                algorithm.tolerance,
            )?
        }
        SolidRoute::Tabulated => {
            let borrowed: Vec<&str> = components.iter().map(String::as_str).collect();
            let (mixture, _) = crate::databank::mixture_of(&borrowed, crate::Cubic::Srk, None)?;
            solve(
                |t| tabulated_residual(&mixture, z, index, solid, eos, t, p_pa),
                search,
                algorithm.tolerance,
            )?
        }
    };

    Ok(FreezingPointResult {
        temperature: kelvins(temperature),
        component: components[index].clone(),
        iterations,
        residual,
        warnings,
    })
}
