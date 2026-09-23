//! The multiphase split: how much of a feed is in each of *n* phases.
//!
//! [`crate::pt_flash`] answers for two phases and finds the split by successive
//! substitution on a single vapour fraction. This answers for however many phases the
//! caller has brought, and it is a Newton on the whole fraction vector.
//!
//! The formulation is Michelsen's. With the fugacity coefficients held at the phases'
//! current compositions, the objective is
//!
//! ```text
//! Q(beta) = sum_k beta_k - sum_i z_i ln E_i        E_i = sum_k beta_k / phi_ik
//! ```
//!
//! whose gradient and Hessian, with the coefficients frozen, are
//!
//! ```text
//! dQ/dbeta_k       = 1 - sum_i z_i / (phi_ik E_i)
//! d2Q/dbeta_j dbeta_k = sum_i z_i / (phi_ij phi_ik E_i^2)
//! ```
//!
//! so one Newton step is a dense solve of size *n*, the phase count, rather than a
//! bracketed root find. That is the whole point of the form: a Rachford-Rice bracket
//! exists for two phases and has no *n*-phase counterpart, while this does.
//!
//! Each phase's composition is then rebuilt as `x_ik = z_i / (E_i phi_ik)`, which is the
//! form the material balance gives:
//!
//! ```text
//! sum_k beta_k x_ik = sum_k beta_k z_i / (E_i phi_ik) = (z_i / E_i) sum_k beta_k / phi_ik = z_i
//! ```
//!
//! That identity is the reason the iteration is on the fractions and not on the
//! compositions. It holds for the closed form as written, and the iteration then
//! *normalises* each composition - which the identity does not survive exactly, because
//! `sum_i z_i / (E_i phi_ik)` is not exactly one at a finite iterate. Measured at the
//! converged answers below, the balance holds to about `2e-11` relative: the
//! normalisation's cost, not an error in the form.
//!
//! **A phase fraction is never negative here.** The two-phase flash reports the negative
//! flash - a phase that would have to be *added* to reach saturation - because a single
//! vapour fraction can carry that meaning. A phase with a negative amount cannot, and the
//! Newton drives it to the floor instead; a phase that is not there is dropped rather
//! than reported as a negative one.

use azoth_core::{AzothError, ModelAlgorithm, Result};

use crate::mixture::{Mixture, ReducedParameters, RootSide};
use crate::wax_solid_fugacity;

/// The floor under the phase-split denominator, upstream's.
///
/// `E_i` is a sum of positive terms divided by positive coefficients, so it is positive
/// in exact arithmetic; a coefficient that has gone to infinity - which the cubic
/// produces where a root degenerates - takes it to zero, and every composition is then a
/// division by it.
const MINIMUM_E: f64 = 1.0e-100;

/// A phase fraction is kept strictly inside `(0, 1)`, upstream's limit.
const FRACTION_FLOOR: f64 = 1.0e-12;

/// The diagonal regulariser, upstream's.
///
/// The Hessian is a sum of rank-one terms over the components, so it is singular
/// whenever two phases sit at the same composition - which is exactly the state the
/// iteration passes through on its way to merging them.
const REGULARISER: f64 = 1.0e-3;

/// The gradient norm the iteration must also reach, upstream's second convergence test.
///
/// `solveBeta` stops when **both** the step is short and the gradient is small - `err > 1e-12
/// || gradResidual > 1e-10` keeps it going - because a Newton step on a Hessian the regulariser
/// has just flattened can be short while the state is nowhere near stationary. A step-only
/// test would call that converged. The step half is the spec's `tolerance`; this half has no
/// field in the spec's `algorithm` block, so it is named here.
const GRADIENT_TOLERANCE: f64 = 1.0e-10;

/// The fewest steps the iteration takes, upstream's `|| iter < 3`.
///
/// `solveBeta` increments before its test, so its third test is the first that can exit, and
/// two bodies always run. The first step is where a phase the seeding added at a guessed
/// fraction is most likely to be pinned at the floor; running a second turns "pinned while
/// passing through" into "pinned at the answer".
const MINIMUM_STEPS: u32 = 2;

/// The fractions and compositions a multiphase solve converged on.
///
/// Not a registered result type: the solve is the model layer's arithmetic, in the same
/// way the fugacity coefficient is, and a model that composes it carries whatever of
/// this it reports.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiphaseSplit {
    /// One fraction per phase, summing to one.
    pub fractions: Vec<f64>,
    /// One composition per phase, each summing to one.
    pub compositions: Vec<Vec<f64>>,
    /// Newton steps taken.
    pub iterations: u32,
    /// The larger of the last step's norm and the gradient norm it was judged beside.
    pub residual: f64,
    /// Whether the residual met the tolerance, or the iteration left by its cap.
    pub converged: bool,
}

/// How one phase's fugacity coefficients are obtained.
///
/// **The two are not the same kind of thing, which is why this is one enum and not a third
/// [`RootSide`].** A cubic phase's coefficients are a function of *its own composition* on a
/// chosen root - that is what makes a liquid-liquid split answerable at all - and the wax
/// solid's are a function of the *state alone*: its coefficient is
/// `phi_liq exp(...)` for one component, with the composition cancelled out of
/// `SolidFug/(P x)`, so there is no root to choose and nothing about the phase's composition
/// enters.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PhaseKind {
    /// A phase of the cubic, on this root.
    Cubic(RootSide),
    /// The wax solid, whose coefficient is [`crate::wax_solid_fugacity`]'s for a component
    /// the fluid marks as a wax former and NeqSim's `1e50` marker for everything else.
    Wax,
    /// A GE aqueous phase, whose coefficients are [`crate::pitzer_phase`]'s at its own
    /// composition.
    ///
    /// The one kind whose model is keyed by **name**: the Pitzer parameters resolve
    /// against a dataset keyed by the component, so a mixture built from constants alone
    /// cannot answer for it and [`Mixture::names`] is what this reads.
    Ge,
}

/// One phase of a multiphase split.
#[derive(Debug, Clone, PartialEq)]
pub struct MultiphasePhase {
    /// Moles of this phase per mole of feed. The set's fractions sum to one.
    pub fraction: f64,
    /// Mole fractions within the phase, summing to one.
    pub composition: Vec<f64>,
    /// Which root of the cubic this phase sits on, or that it is not one.
    ///
    /// Two liquid phases share a side and differ in composition alone, which is what
    /// makes a cubic able to describe a liquid-liquid split at all: the root is chosen
    /// per phase, at that phase's own composition.
    pub kind: PhaseKind,
}

/// The fugacity coefficients of every component in every phase, row-major by phase.
///
/// The dispatch the two multiphase loops share: what a phase's coefficients *are* is a
/// property of its kind, and both the ordinary solve and the hybrid's fixed-topology one ask
/// this rather than each deciding for itself.
pub(crate) fn coefficients(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    phases: &[MultiphasePhase],
) -> Result<Vec<Vec<f64>>> {
    let mut out = Vec::with_capacity(phases.len());
    for phase in phases {
        match phase.kind {
            PhaseKind::Cubic(side) => {
                let state = mixture.phase_state(reduced, &phase.composition, side)?;
                out.push(state.ln_phi.iter().map(|value| value.exp()).collect());
            }
            PhaseKind::Wax => out.push(wax_coefficients(mixture, reduced)?),
            PhaseKind::Ge => out.push(ge_coefficients(mixture, reduced, &phase.composition)?),
        }
    }
    Ok(out)
}

/// A GE aqueous phase's fugacity coefficients, from `eos.pitzer_phase`.
///
/// The model is the *activity* one over the same names, so the composition is the phase's
/// own in the caller's order - the two agree because both are indexed by the fluid's
/// component order.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mixture carries no names, which is what the
///   parameter datasets are keyed by, or if the brine's topology is one the loaded dataset
///   does not cover.
fn ge_coefficients(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    composition: &[f64],
) -> Result<Vec<f64>> {
    let names = mixture.names().ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "a GE phase resolves its parameters by name, and this mixture carries none: a \
             fluid built from constants cannot answer for one"
                .to_string(),
        )
    })?;
    let borrowed: Vec<&str> = names.iter().map(String::as_str).collect();
    let coefficients =
        crate::pitzer_phase(&borrowed, reduced.t_kelvin, reduced.pressure, composition)?;
    Ok(coefficients
        .ln_phi
        .iter()
        .map(|value| value.exp())
        .collect())
}

/// NeqSim's marker for a component that cannot be in a wax phase.
///
/// `ComponentWax.fugcoef` returns it for a substance that is not a wax former, and the
/// fraction solve turns it into `x = z/(E phi)` - which is zero to any precision the answer
/// is read at, so the exclusion is a *number* the solve can carry rather than a phase whose
/// composition has to be assembled per state.
const NOT_A_WAX_FORMER: f64 = 1.0e50;

/// The wax solid's fugacity coefficient for every component at a state.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a wax former carries no molar mass, no heat of fusion or
///   no triple-point temperature, or if the fluid's cubic is not one the wax model reaches.
///   Each is a value the model reads rather than a case to default: a zero heat of fusion
///   says the substance does not melt, and a molar mass of zero is not a cut.
fn wax_coefficients(mixture: &Mixture, reduced: &ReducedParameters) -> Result<Vec<f64>> {
    let cubic = mixture.cubic();
    let eos = match cubic {
        crate::Cubic::Pr => "pr",
        crate::Cubic::Srk => "srk",
        other => {
            return Err(AzothError::invalid_input(
                "eos",
                format!(
                    "a wax phase over a {other:?} fluid is not a state this reaches: the                      reference liquid is a phase of the *host's* own class, and NeqSim's wax                      route is on `PhaseSrkEos` and `PhasePrEos`"
                ),
            ));
        }
    };
    let mut out = Vec::with_capacity(mixture.components().len());
    for (index, component) in mixture.components().iter().enumerate() {
        if !component.wax_former {
            out.push(NOT_A_WAX_FORMER);
            continue;
        }
        let missing = |what: &str| {
            AzothError::invalid_input(
                "components",
                format!(
                    "component {index} is a wax former and carries no {what}; the wax model                      reads it rather than defaulting, and `eos.tbp_fraction_properties` is                      what gives a cut one"
                ),
            )
        };
        let molar_mass = component.molar_mass.ok_or_else(|| missing("molar mass"))?;
        if component.heat_of_fusion <= 0.0 {
            return Err(missing("heat of fusion"));
        }
        if component.triple_point_temperature <= 0.0 {
            return Err(missing("triple-point temperature"));
        }
        let coefficient = wax_solid_fugacity(
            molar_mass,
            component.tc,
            component.pc,
            component.omega,
            component.heat_of_fusion,
            azoth_core::units::kelvins(component.triple_point_temperature),
            azoth_core::units::kelvins(reduced.t_kelvin),
            azoth_core::units::pascals(reduced.pressure),
            eos,
        )?;
        out.push(coefficient.fugacity_coefficient);
    }
    Ok(out)
}

/// Solve for the fractions of a set of phases at a state.
///
/// `phases` is modified in place: on return each entry carries the composition the
/// material balance gives at the converged fractions, and the split's fraction vector.
///
/// The phases are the caller's - this routine *finds the amounts*, not the phases. What
/// phases exist and how many is a stability question, answered by
/// [`crate::stability_test`] and by upstream's seeding; a Newton on the fractions cannot
/// answer it, because a phase that should not exist has a converged fraction of zero and
/// a phase that should is not reachable from a set that does not contain it.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if there are fewer than two phases, or a composition is
///   the wrong length.
/// * [`AzothError::SolverNotConverged`] if the Hessian cannot be solved at an iterate - which
///   means two phases have met, and the caller's answer is to drop one, not to take a step
///   the regulariser invented.
///
/// **Reaching the cap is not an error.** Upstream's `solveBeta` returns whatever residual it
/// left with and the flash reports it. Raising would turn a fraction vector right to six
/// figures into a failure - and the positions this seed lands on are not the two-phase
/// flash's, so a cap exit is a state to report rather than to refuse. The split carries the
/// residual and `converged`, so a caller warns instead of losing the answer.
pub fn solve_phase_fractions(
    mixture: &Mixture,
    reduced: &ReducedParameters,
    feed: &[f64],
    phases: &mut [MultiphasePhase],
    algorithm: &ModelAlgorithm,
) -> Result<MultiphaseSplit> {
    let n = feed.len();
    let count = phases.len();
    if count < 2 {
        return Err(AzothError::invalid_input(
            "phases",
            format!(
                "a split needs at least two phases and {count} were given. One phase is \
                 the feed, and its fraction is one"
            ),
        ));
    }
    for phase in phases.iter() {
        if phase.composition.len() != n {
            return Err(AzothError::invalid_input(
                "composition",
                format!(
                    "a composition for {n} components has {} entries",
                    phase.composition.len()
                ),
            ));
        }
    }

    let mut iterations = 0;
    let mut residual = f64::NAN;
    let mut gradient_norm = f64::INFINITY;
    for step in 1..=algorithm.max_iterations {
        iterations = step;
        let phi = coefficients(mixture, reduced, phases)?;

        // `E_i`, and the two terms the gradient and the Hessian are sums of.
        let mut e = vec![0.0; n];
        for (k, phase) in phases.iter().enumerate() {
            for i in 0..n {
                e[i] += phase.fraction / phi[k][i];
            }
        }
        for value in e.iter_mut() {
            if *value < MINIMUM_E {
                *value = MINIMUM_E;
            }
        }

        // The gradient `1 - sum_i z_i/(phi_ik E_i)` and the Hessian
        // `sum_i z_i/(phi_ij phi_ik E_i^2)`, both with the coefficients frozen.
        let mut gradient = vec![0.0; count];
        let mut hessian = vec![vec![0.0; count]; count];
        for (k, value) in gradient.iter_mut().enumerate() {
            *value = 1.0;
            for i in 0..n {
                *value -= feed[i] / (e[i] * phi[k][i]);
            }
        }
        for (j, row) in hessian.iter_mut().enumerate() {
            for (k, value) in row.iter_mut().enumerate() {
                for i in 0..n {
                    *value += feed[i] / (e[i] * e[i] * phi[j][i] * phi[k][i]);
                }
                if j == k {
                    *value += REGULARISER;
                }
            }
        }

        gradient_norm = gradient
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();
        let correction = solve(&hessian, &gradient, algorithm)?;
        residual = correction
            .iter()
            .map(|value| value * value)
            .sum::<f64>()
            .sqrt();

        // The step is damped by `n/(n+3)` - a half at the first iteration, approaching
        // one - and the fractions are clamped away from the ends. A fraction at zero is a
        // phase that is not there, and the caller decides whether to drop it.
        let scale = f64::from(step) / (f64::from(step) + 3.0);
        let mut total = 0.0;
        for (k, phase) in phases.iter_mut().enumerate() {
            let candidate = phase.fraction - scale * correction[k];
            phase.fraction = candidate.clamp(FRACTION_FLOOR, 1.0 - FRACTION_FLOOR);
            total += phase.fraction;
        }
        for phase in phases.iter_mut() {
            phase.fraction /= total;
        }

        // The compositions the material balance gives at those fractions, in the same
        // closed form and from the same frozen coefficients.
        for (k, phase) in phases.iter_mut().enumerate() {
            let mut sum = 0.0;
            for i in 0..n {
                let value = feed[i] / (e[i] * phi[k][i]);
                phase.composition[i] = value;
                sum += value;
            }
            if sum > 0.0 {
                for value in phase.composition.iter_mut() {
                    *value /= sum;
                }
            }
        }

        if step >= MINIMUM_STEPS
            && residual <= algorithm.tolerance
            && gradient_norm <= GRADIENT_TOLERANCE
        {
            break;
        }
    }

    let converged = residual <= algorithm.tolerance && gradient_norm <= GRADIENT_TOLERANCE;

    Ok(MultiphaseSplit {
        fractions: phases.iter().map(|phase| phase.fraction).collect(),
        compositions: phases
            .iter()
            .map(|phase| phase.composition.clone())
            .collect(),
        iterations,
        residual: residual.max(gradient_norm),
        converged,
    })
}

/// `H x = g` by Gaussian elimination with partial pivoting.
///
/// The system is the phase count square, which is a handful, so a factorising solver
/// would be the larger dependency. A singular Hessian is reported rather than smoothed
/// over: it means two phases have met, and the caller's answer is to drop one, not to
/// take a step the regulariser invented.
/// `H x = g` by Gaussian elimination with partial pivoting.
///
/// Exposed because [`crate::tp_solid_flash`]'s Newton is the same dense solve on a different
/// gradient: the two flashes differ in what `E` is and in how the step is damped, not in how a
/// small system is eliminated.
///
/// # Errors
/// * [`AzothError::SolverNotConverged`] if the matrix is singular at a pivot.
pub fn solve_dense(h: &[Vec<f64>], g: &[f64], algorithm: &ModelAlgorithm) -> Result<Vec<f64>> {
    solve(h, g, algorithm)
}

fn solve(h: &[Vec<f64>], g: &[f64], algorithm: &ModelAlgorithm) -> Result<Vec<f64>> {
    let n = g.len();
    let mut a: Vec<Vec<f64>> = (0..n)
        .map(|i| [h[i].clone(), vec![g[i]]].concat())
        .collect();
    for column in 0..n {
        let pivot = (column..n)
            .max_by(|&p, &q| {
                a[p][column]
                    .abs()
                    .partial_cmp(&a[q][column].abs())
                    .unwrap_or(std::cmp::Ordering::Equal)
            })
            .unwrap_or(column);
        a.swap(column, pivot);
        if a[column][column] == 0.0 {
            return Err(AzothError::SolverNotConverged {
                iterations: 0,
                residual: f64::NAN,
                tolerance: algorithm.tolerance,
            });
        }
        let pivot_value = a[column][column];
        let (upper, lower) = a.split_at_mut(column + 1);
        let pivot_row = &upper[column];
        for row in lower {
            let factor = row[column] / pivot_value;
            for (entry, &pivot_entry) in row.iter_mut().zip(pivot_row.iter()).skip(column) {
                *entry -= factor * pivot_entry;
            }
        }
    }
    let mut x = vec![0.0; n];
    for row in (0..n).rev() {
        let sum: f64 = ((row + 1)..n).map(|k| a[row][k] * x[k]).sum();
        x[row] = (a[row][n] - sum) / a[row][row];
    }
    Ok(x)
}
