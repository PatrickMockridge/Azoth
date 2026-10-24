//! `GibbsReactor`'s solve: the Lagrange-multiplier Newton iteration over the element balances.
//!
//! **This is the class's solver and not a better one.** It is a constrained minimisation of
//! the Gibbs energy written as a root find on
//!
//! ```text
//! F_i = Gf0_i(T) + RT ln phi_i + RT ln y_i + RT ln(P_bara) - sum_j lambda_j a_ij
//! ```
//!
//! for the species rows, with one element-balance row per *active* element, and the whole
//! thing driven by a damped Newton step whose composition update is `alpha` and whose
//! multipliers take the full step.
//!
//! **The measurements the port is held to are in
//! `validation/neqsim/captures/process_gibbs_reactor.tsv`**, which prints the path and not
//! only the destination: the iteration count, the final error, the element balances on the
//! class's own seven element names, the seven multipliers, the objective vector, the
//! fugacity coefficients, and the Gibbs energy history per iteration.
//!
//! Four behaviours are the class's and are reproduced rather than tidied:
//!
//! * **The convergence test is on the undamped step.** `deltaX` is `-J^-1 F`, and
//!   `deltaXNorm` is its norm *before* `alpha` is applied, while the update takes
//!   `alpha = dampingComposition` (0.05 by default) of it. So the test and the step disagree
//!   about scale, and the tolerance of `1e-3` is loose against the step actually taken.
//! * **`minIterations = 100` is a floor and it bites.** Three of the capture's six rows stop
//!   at exactly 100 with errors of `8.1e-7`, `5.8e-15` and `5.7e-5` - the first iteration at
//!   which `iteration >= minIterations` holds. A port that stopped on the tolerance alone
//!   would stop earlier.
//! * **Reaching `maxIterations` returns success.** The class sets `converged` from the
//!   tolerance, so a run that hits the cap reports `hasConverged() == false` - but the
//!   method returns `true` and the caller reads what is there. Both are carried.
//! * **A component the database does not carry is still a variable.** It enters the Jacobian
//!   with a composition row and no element entry, its objective value is absent (and reads as
//!   zero), and its update is skipped - which is how `testComponentNotInDatabaseMolesUnchanged`
//!   gets its name.
//!
//! **The step controls and the regularisation are ported, and the capture measures what each
//! does.** All three flags are `false` by default and none is reached by a row of
//! `process_gibbs_reactor.tsv`; `validation/neqsim/captures/process_gibbs_algorithms.tsv` drives
//! them from the class's own algorithm test. The Armijo arm accepts its maximum on that test's
//! state and only the undamped row backtracks; the adaptive limit answers `1.0, 1.0, 0.0193,
//! ...`; and the augmentation fires on the regularized methane row, whose first condition number
//! is `2.9e12` against the class's `1e10` threshold.
//!
//! **The condition number is recorded on every pass after the first.** The class's disjunction
//! is `useRegularization || !conditionNumberHistory.isEmpty() || iteration <= 5`, so the history
//! is the first five passes and then the rest of the run - `getConditionNumberHistory()` is a
//! diagnostic of the whole trajectory, not a sample of it.
//!
//! **This iteration is ill-conditioned, and that bounds how closely a port can match it.** Two
//! things compound, both the class's own: the `MIN_MOLES = 1e-6` floor puts an `RT/n` entry of
//! about `1e6` in the Jacobian against entries near `1`, so its 2-norm condition number is
//! `6.1e6` on the ammonia case; and the residual is a cancelling sum, where hydrogen's `F` of
//! `6.5e-5` is built from terms summing to `-0.345`, so a relative disagreement of `e` between
//! two implementations in the terms is `5300e` in `F`. One iteration from the same state agrees
//! to `1e-13` on the major species and `1.6e-10` on the Jacobian; over the hundred damped steps
//! the class takes, that is enough to separate two runs. The tests say which rows close and
//! which do not, and why.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};
use azoth_eos::mixture::{Mixture, RootSide};
use azoth_eos::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use azoth_reactions::linalg::{condition_number_2, solve_least_squares, solve_lu};

use crate::reactor::gibbs_database::{GibbsDatabase, GibbsSpecies, R_KJ};

/// `MIN_MOLES`, the class's floor for a mole number in a logarithm.
const MIN_MOLES: f64 = 1e-6;

/// `MIN_JACOBIAN_MOLES`, the same floor where the Jacobian divides by it.
const MIN_JACOBIAN_MOLES: f64 = 1e-6;

/// `ELEMENT_ZERO_THRESHOLD`, the class's test for "this element is not present".
const ELEMENT_ZERO_THRESHOLD: f64 = 1e-6;

/// The floor `performIterationUpdate` puts under a composition after the damped step.
const MIN_COMPOSITION: f64 = 1e-15;

/// `MAX_STEP_LIMIT`, `Math.log(5.0)`: the largest fractional mole change one adaptive step
/// allows a significant component.
const MAX_STEP_LIMIT: f64 = 1.609_437_912_434_100_4;
/// `calculateAdaptiveAlpha`'s floor under the step it returns.
const MIN_ADAPTIVE_ALPHA: f64 = 1e-4;
/// `armijoLineSearch`'s floor when every backtrack fails.
const MIN_ARMIJO_ALPHA: f64 = 1e-6;
/// `evaluateTotalGibbsEnergy`'s `1e-30` floor under a mole number.
const GIBBS_ZERO_MOLES: f64 = 1e-30;

/// The class's controls, as `GibbsReactor`'s own fields default them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct GibbsSettings {
    /// Whether the energy balance moves the temperature inside the loop.
    pub adiabatic: bool,
    /// `dampingComposition`, the fraction of the Newton step the compositions take.
    pub damping: f64,
    /// `maxIterations`.
    pub max_iterations: u32,
    /// `convergenceTolerance`, applied to the undamped step's norm.
    pub tolerance: f64,
    /// `minIterations`, the floor below which the tolerance is not consulted.
    pub min_iterations: u32,
    /// `useArmijoLineSearch`, the backtracking line search on the total Gibbs energy.
    pub armijo_line_search: bool,
    /// `armijoC1`, the sufficient-decrease constant.
    pub armijo_c1: f64,
    /// `armijoRho`, the backtracking contraction.
    pub armijo_rho: f64,
    /// `armijoMaxBacktracks`.
    pub armijo_max_backtracks: u32,
    /// `useAdaptiveStepSize`, the NASA CEA-style step limiter.
    pub adaptive_step_size: bool,
    /// `useRegularization`, the Tikhonov diagonal augmentation.
    pub regularization: bool,
    /// `regularizationThreshold`, the condition number above which it fires.
    pub regularization_threshold: f64,
    /// `regularizationTau`, the diagonal addition before scaling.
    pub regularization_tau: f64,
}

impl Default for GibbsSettings {
    fn default() -> Self {
        Self {
            adiabatic: false,
            damping: 0.05,
            max_iterations: 5000,
            tolerance: 1e-3,
            min_iterations: 100,
            armijo_line_search: false,
            armijo_c1: 1e-4,
            armijo_rho: 0.5,
            armijo_max_backtracks: 20,
            adaptive_step_size: false,
            regularization: false,
            regularization_threshold: 1e10,
            regularization_tau: 1e-6,
        }
    }
}

/// What the solve reached, with the trace a fixed point needs to be checkable.
#[derive(Debug, Clone, PartialEq)]
pub struct GibbsState {
    /// The outlet moles per component, in the fluid's order.
    pub moles: Vec<f64>,
    /// The outlet temperature, K: the inlet's unless the energy mode is adiabatic.
    pub temperature: f64,
    /// The seven Lagrange multipliers, on the class's element names.
    pub lambdas: [f64; 7],
    /// `hasConverged()`, which is the *tolerance* test and not the return value.
    pub converged: bool,
    /// `getActualIterations()`.
    pub iterations: u32,
    /// `getFinalConvergenceError()`, the last undamped step norm.
    pub final_error: f64,
    /// The element balance on the inlet, mol/s.
    pub element_in: [f64; 7],
    /// The element balance on the outlet.
    pub element_out: [f64; 7],
    /// `element_out - element_in`.
    pub element_diff: [f64; 7],
    /// The total Gibbs energy at the top of each iteration, kJ/mol.
    pub gibbs_history: Vec<f64>,
    /// The objective value per component, `None` where the database carries no species.
    pub objective: Vec<Option<f64>>,
    /// Which components are variables - neither inert nor excluded by the feed.
    pub variables: Vec<usize>,
    /// Which element indices are active, in the class's own order.
    pub active_elements: Vec<usize>,
    /// The Jacobian as the last pass computed it: species rows and columns first, then the
    /// active element balances. Carried because it is the only way to compare one iteration of
    /// this port against one of the class's - `getJacobianMatrix()` is public there, and a
    /// trajectory that separates says nothing about which of the two maps moved.
    pub jacobian: Vec<Vec<f64>>,
    /// `getConditionNumberHistory()`, one 2-norm condition number per pass it was recorded on:
    /// the first five always, and every pass when regularization is on.
    pub condition_number_history: Vec<f64>,
}

/// Minimise the Gibbs energy of a fluid at a temperature and pressure.
///
/// `inlet_moles` is the feed's moles per component and `mixture` must be the fluid those
/// components name, with `ideal_gas` its heat-capacity polynomials - which the class reads
/// per component through `getCpA()`..`getCpD()` and which
/// [`azoth_eos::databank::mixture_of`] returns beside the mixture.
///
/// **The pressure is taken in the unit the class's objective writes it in.** The residual
/// carries `RT ln(P/1 bara)`, so `pressure_bara` is the numeric bar value and not a pascal
/// quantity - a caller that passed pascals would shift every species row equally and converge
/// to the same composition with a different `F`, which is the failure the parameter's name is
/// there to prevent.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a component count does not match the mixture, if the
///   total moles are not positive, if a mole number or fugacity coefficient is not positive -
///   the class would take the logarithm of a non-positive number and carry a NaN through the
///   whole solve - or if the Newton system is singular, which is the one case the class
///   answers through a pseudo-inverse this port does not carry.
pub fn solve(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    database: &GibbsDatabase,
    inlet_moles: &[f64],
    temperature: ThermodynamicTemperature,
    pressure_bara: f64,
    settings: &GibbsSettings,
) -> Result<GibbsState> {
    let n = mixture.len();
    let names = mixture.names().ok_or_else(|| {
        AzothError::invalid_input(
            "mixture",
            "the solve looks its species up by name, and this mixture carries none",
        )
    })?;
    if inlet_moles.len() != n {
        return Err(AzothError::invalid_input(
            "inlet_moles",
            format!(
                "a feed for {n} components has {} entries",
                inlet_moles.len()
            ),
        ));
    }
    if !pressure_bara.is_finite() || pressure_bara <= 0.0 {
        return Err(AzothError::invalid_input(
            "pressure",
            format!(
                "the objective takes ln(P), so P must be positive, and it is {pressure_bara} bara"
            ),
        ));
    }

    // The species each component resolves to, or `None` for a component outside the database.
    let species: Vec<Option<&GibbsSpecies>> = names.iter().map(|name| database.get(name)).collect();
    // The component's own `A`..`D`, which `calculateCorrectedHeatCapacityCoeffs` subtracts the
    // elements from. Only the first four of the five are read, as the class reads them.
    let cp: Vec<[f64; 4]> = (0..n)
        .map(|i| {
            [
                ideal_gas.cp_a[i],
                ideal_gas.cp_b[i],
                ideal_gas.cp_c[i],
                ideal_gas.cp_d[i],
            ]
        })
        .collect();

    // **The element balance on the *feed*, before any floor is applied.** The class reads
    // `inlet_mole` for this and `outlet_mole` for the other side, and both exclusions below
    // are decided from this one.
    let element_in = element_balance(&species, inlet_moles);

    // `determineFeedExcludedComponents`: a species that needs an element the feed has none of
    // cannot form, so it is frozen at its feed amount and dropped from the matrix.
    let excluded: Vec<bool> = species
        .iter()
        .map(|species| match species {
            None => false,
            Some(species) => species
                .balance_elements()
                .iter()
                .enumerate()
                .any(|(j, &count)| {
                    count.abs() > ELEMENT_ZERO_THRESHOLD
                        && element_in[j].abs() <= ELEMENT_ZERO_THRESHOLD
                }),
        })
        .collect();

    // `performGibbsMinimization` and `enforceMinimumConcentrations` both floor a database
    // component at `MIN_MOLES`; the two are the same operation applied twice, and a
    // feed-excluded species keeps its feed amount exactly.
    let mut moles: Vec<f64> = (0..n)
        .map(|i| match species[i] {
            Some(_) if !excluded[i] => inlet_moles[i].max(MIN_MOLES),
            _ => inlet_moles[i],
        })
        .collect();

    // `variableComponents`: everything the database carries and the feed does not exclude -
    // **including a component the database does not carry**, which is the class's own
    // handling of a substance it has no species for.
    let variables: Vec<usize> = (0..n).filter(|&i| !excluded[i]).collect();

    let active_elements = active_elements(&species, &element_in);
    let mut lambdas = [0.0f64; 7];
    let mut temperature_value = temperature.value;
    let mut gibbs_history = Vec::new();
    let mut condition_history: Vec<f64> = Vec::new();
    let mut objective: Vec<Option<f64>> = vec![None; n];
    let mut converged = false;
    let mut iterations = 0u32;
    let mut final_error = f64::MAX;
    let mut enthalpy_old = 0.0;
    let mut last_jacobian: Vec<Vec<f64>> = Vec::new();

    for iteration in 1..=settings.max_iterations {
        iterations = iteration;

        // `calculateObjectiveFunctionValues`, then the fugacity coefficients it needs.
        objective = objective_values(
            mixture,
            database,
            names,
            &species,
            &cp,
            &moles,
            temperature_value,
            pressure_bara,
            &lambdas,
        )?;

        // `performNewtonRaphsonIteration`: the Jacobian, the objective vector, and the solve.
        let total: f64 = moles.iter().sum();
        let (ln_phi, d_ln_phi_dn) =
            phase_derivatives(mixture, &moles, temperature_value, pressure_bara)?;
        let mut jacobian = jacobian(
            mixture,
            database,
            &species,
            &moles,
            &variables,
            &active_elements,
            total,
            temperature_value,
            &ln_phi,
            &d_ln_phi_dn,
        );
        let vector = objective_vector(
            &objective,
            &variables,
            &active_elements,
            &element_in,
            &moles,
            &species,
        );
        let negated: Vec<f64> = vector.iter().map(|v| -v).collect();
        // **The class's `solveNewtonSystem` catches the LU failure and retries through
        // `SimpleMatrix.pseudoInverse`**, and `solve_least_squares` is that inverse's own
        // construction - its doc says so, and it is the same one-sided Jacobi decomposition.
        // The reachable state is a feed with fewer species than active elements:
        // `steam_methane_hot` carries methane and water alone, so two species meet three
        // balance rows, and the capture's row converges there.
        let mut delta = match solve_lu(&jacobian, &negated) {
            Ok(delta) if delta.iter().all(|value| value.is_finite()) => delta,
            _ => solve_least_squares(&jacobian, &negated)?,
        };

        // **`applyRegularization`, called every pass while `iteration <= 5`** even with the
        // flag off, because the class records the condition number there; and every pass when
        // it is on. When the augmentation fires the condition drops, and the class then
        // re-solves the Newton system against the mutated Jacobian - its own condition at
        // `:3016-3025`, reproduced as written.
        // **`|| !condition_history.is_empty()` is the class's own disjunction**: once the first
        // number is recorded, every later pass records one too, so the history is the first
        // five passes followed by the rest of the run.
        if settings.regularization || !condition_history.is_empty() || iteration <= 5 {
            let number = apply_regularization(&settings, &mut jacobian, variables.len())?;
            condition_history.push(number);
            if settings.regularization && number < settings.regularization_threshold {
                delta = match solve_lu(&jacobian, &negated) {
                    Ok(delta) if delta.iter().all(|value| value.is_finite()) => delta,
                    _ => solve_least_squares(&jacobian, &negated)?,
                };
            }
        }

        // The state carries the matrix `getJacobianMatrix()` would answer with, which is the
        // *regularized* one when the augmentation fired - the class mutates its field in place.
        last_jacobian.clone_from(&jacobian);
        // The class's `columnScaleFactors` are the identity, so `unscaleNewtonStep` is one.
        let delta_norm = delta.iter().map(|v| v * v).sum::<f64>().sqrt();

        // The Gibbs energy at the top of the iteration, before any update.
        let gibbs = mixture_gibbs_energy(&species, &moles, temperature_value, &cp);
        gibbs_history.push(gibbs);

        if settings.adiabatic {
            let enthalpy = mixture_enthalpy(&species, &moles, temperature_value, &cp);
            if iteration == 1 {
                enthalpy_old = enthalpy;
            } else {
                let dh = enthalpy - enthalpy_old;
                enthalpy_old = enthalpy;
                let cp_total = total
                    * total_heat_capacity(
                        mixture,
                        ideal_gas,
                        &moles,
                        temperature_value,
                        pressure_bara,
                    )?;
                let t_out = temperature_value - dh * 1000.0 / cp_total;
                if (t_out - temperature_value).abs() > 1000.0 {
                    return Err(AzothError::invalid_input(
                        "temperature",
                        format!(
                            "the adiabatic step moves the temperature by {} K in one iteration, \
                             and the class refuses beyond 1000 - reduce the damping",
                            (t_out - temperature_value).abs()
                        ),
                    ));
                }
                temperature_value = t_out;
            }
        }

        // **The convergence test, on the undamped step, and the class's floor under it.**
        if (delta_norm < settings.tolerance && iteration >= settings.min_iterations)
            || iteration == settings.max_iterations
        {
            converged = delta_norm < settings.tolerance;
            final_error = delta_norm;
            break;
        }

        // **The step size: adaptive first, then Armijo, then the fixed damping.** The class's
        // own order, with the adaptive branch replacing the damping rather than bounding it.
        let mut alpha = settings.damping;
        if settings.adaptive_step_size {
            alpha = adaptive_alpha(&moles, &variables, &delta);
        }
        if settings.armijo_line_search {
            let current_g = total_gibbs_energy(
                mixture,
                &species,
                &moles,
                temperature_value,
                pressure_bara,
                &cp,
            )?;
            alpha = armijo_alpha(
                mixture,
                &species,
                &cp,
                &moles,
                temperature_value,
                pressure_bara,
                &objective,
                &variables,
                &delta,
                alpha,
                current_g,
                &settings,
            )?;
        }

        // `performIterationUpdate`: the compositions take `alpha`, the multipliers take all.
        for (slot, &i) in variables.iter().enumerate() {
            moles[i] = (moles[i] + delta[slot] * alpha).max(MIN_COMPOSITION);
        }
        for (slot, &element) in active_elements.iter().enumerate() {
            lambdas[element] += delta[variables.len() + slot];
        }
        delta.clear();
    }

    let element_out = element_balance(&species, &moles);
    let mut element_diff = [0.0; 7];
    for j in 0..7 {
        element_diff[j] = element_out[j] - element_in[j];
    }

    Ok(GibbsState {
        moles,
        temperature: temperature_value,
        lambdas,
        converged,
        iterations,
        final_error,
        element_in,
        element_out,
        element_diff,
        gibbs_history,
        objective,
        variables,
        active_elements,
        jacobian: last_jacobian,
        condition_number_history: condition_history,
    })
}

/// `applyRegularization`: the 2-norm condition number, and the diagonal augmentation.
///
/// **The augmentation fires only above the threshold**, and `tau` is raised to
/// `maxDiag * 1e-8` when that is larger than the configured value. The first
/// `variables` diagonal entries are the composition block; the element rows are untouched.
fn apply_regularization(
    settings: &GibbsSettings,
    jacobian: &mut [Vec<f64>],
    variables: usize,
) -> Result<f64> {
    let mut number = condition_number_2(jacobian)?;
    if settings.regularization && number > settings.regularization_threshold {
        let max_diag = (0..variables)
            .map(|i| jacobian[i][i].abs())
            .fold(0.0_f64, f64::max);
        let mut tau = settings.regularization_tau;
        if max_diag > 0.0 {
            tau = tau.max(max_diag * 1e-8);
        }
        for i in 0..variables {
            jacobian[i][i] += tau;
        }
        number = condition_number_2(jacobian)?;
    }
    Ok(number)
}

/// `evaluateTotalGibbsEnergy`: the total Gibbs energy **with the fugacity coefficients**, kJ/mol.
///
/// **A second Gibbs energy in this file, and not the history's.** The history pushes
/// `calculateMixtureGibbsEnergy` - the formation sum alone - while Armijo compares this one,
/// which adds `RT ln phi`, `RT ln y` and `RT ln P`. The class keeps both and they are not
/// interchangeable.
fn total_gibbs_energy(
    mixture: &Mixture,
    species: &[Option<&GibbsSpecies>],
    moles: &[f64],
    temperature: f64,
    pressure_bara: f64,
    cp: &[[f64; 4]],
) -> Result<f64> {
    let rt = R_KJ * temperature;
    let total: f64 = moles.iter().sum();
    let (ln_phi, _) = phase_derivatives(mixture, moles, temperature, pressure_bara)?;
    let ln_p = pressure_bara.ln();
    let mut total_g = 0.0;
    for (i, entry) in species.iter().enumerate() {
        let Some(species) = entry else {
            continue;
        };
        let ni = moles[i];
        if ni < GIBBS_ZERO_MOLES {
            continue;
        }
        let y = ni / total;
        total_g += ni
            * (species.gibbs_energy(temperature, cp[i]) + rt * ln_phi[i] + rt * y.ln() + rt * ln_p);
    }
    Ok(total_g)
}

/// The Armijo slope: `sum_i F_i * delta_i` over the composition variables.
///
/// **The class sums only the composition entries**, not the multiplier rows, and reads `F` from
/// the objective map at the current point.
fn directional_derivative(objective: &[Option<f64>], variables: &[usize], delta: &[f64]) -> f64 {
    variables
        .iter()
        .enumerate()
        .map(|(slot, &i)| objective[i].unwrap_or(0.0) * delta[slot])
        .sum()
}

/// `calculateAdaptiveAlpha`: the NASA CEA-style step limit, from a requested `1.0`.
///
/// **The requested step is `1.0` and not `dampingComposition`**, which is what the class passes
/// at its call site; the fixed damping is what this replaces.
fn adaptive_alpha(moles: &[f64], variables: &[usize], delta: &[f64]) -> f64 {
    let mut alpha = 1.0_f64;
    let total: f64 = moles.iter().sum();
    let significance = total * 1e-3;
    for (slot, &i) in variables.iter().enumerate() {
        let current = moles[i];
        if current < significance {
            continue;
        }
        let abs_delta = delta[slot].abs();
        if abs_delta < 1e-30 {
            continue;
        }
        alpha = alpha.min(MAX_STEP_LIMIT * current / abs_delta);
    }
    for (slot, &i) in variables.iter().enumerate() {
        let current = moles[i];
        let step = delta[slot];
        if step < 0.0 && current > MIN_COMPOSITION {
            alpha = alpha.min(0.9 * current / (-step));
        }
    }
    alpha.max(MIN_ADAPTIVE_ALPHA)
}

/// `armijoLineSearch`: the backtracking step that satisfies `G(trial) <= G + c1*alpha*slope`.
///
/// **The trial needs the fugacity coefficients at the trial composition** - the class calls
/// `updateSystemWithNewCompositions(); system.init(3)` before `evaluateTotalGibbsEnergy` - so
/// each backtrack is one `phase_derivatives` call. The class's save/restore is a copy here: the
/// trial is its own vector and `moles` is untouched.
#[allow(clippy::too_many_arguments)] // the class's own argument list, less its state
fn armijo_alpha(
    mixture: &Mixture,
    species: &[Option<&GibbsSpecies>],
    cp: &[[f64; 4]],
    moles: &[f64],
    temperature: f64,
    pressure_bara: f64,
    objective: &[Option<f64>],
    variables: &[usize],
    delta: &[f64],
    alpha_max: f64,
    current_g: f64,
    settings: &GibbsSettings,
) -> Result<f64> {
    let slope = directional_derivative(objective, variables, delta);
    if slope >= 0.0 {
        return Ok(alpha_max);
    }
    let mut alpha = alpha_max;
    for _ in 0..settings.armijo_max_backtracks {
        let mut trial = moles.to_vec();
        for (slot, &i) in variables.iter().enumerate() {
            trial[i] = (moles[i] + alpha * delta[slot]).max(MIN_COMPOSITION);
        }
        let trial_g = total_gibbs_energy(mixture, species, &trial, temperature, pressure_bara, cp)?;
        if trial_g <= current_g + settings.armijo_c1 * alpha * slope {
            return Ok(alpha);
        }
        alpha *= settings.armijo_rho;
    }
    Ok(alpha.max(MIN_ARMIJO_ALPHA))
}

/// `calculateObjectiveFunctionValues`' own loop, skipping a component with no species.
#[allow(clippy::too_many_arguments)] // The class's own argument list, less its state.
fn objective_values(
    mixture: &Mixture,
    database: &GibbsDatabase,
    names: &[String],
    species: &[Option<&GibbsSpecies>],
    cp: &[[f64; 4]],
    moles: &[f64],
    temperature: f64,
    pressure_bara: f64,
    lambdas: &[f64; 7],
) -> Result<Vec<Option<f64>>> {
    let _ = database;
    let rt = R_KJ * temperature;
    let total: f64 = moles.iter().sum();
    let (ln_phi, _) = phase_derivatives(mixture, moles, temperature, pressure_bara)?;

    let mut values = vec![None; moles.len()];
    for i in 0..moles.len() {
        let Some(species) = species[i] else {
            continue;
        };
        if total <= 0.0 {
            return Err(AzothError::invalid_input(
                "moles",
                "the objective takes ln(y), and the fluid carries no moles",
            ));
        }
        let species = &species;
        let y = moles[i] / total;
        if y <= 0.0 {
            return Err(AzothError::invalid_input(
                "moles",
                format!(
                    "{:?} has {y} as its mole fraction, and the objective takes its logarithm",
                    names[i]
                ),
            ));
        }
        let phi = ln_phi[i].exp();
        if !phi.is_finite() || phi <= 0.0 {
            return Err(AzothError::invalid_input(
                "fugacity coefficient",
                format!("{:?} has a fugacity coefficient of {phi}", names[i]),
            ));
        }
        let mut lagrange = 0.0;
        let elements = species.balance_elements();
        for (j, lambda) in lambdas.iter().enumerate() {
            lagrange += lambda * elements[j];
        }
        values[i] = Some(
            species.gibbs_energy(temperature, cp[i])
                + rt * ln_phi[i]
                + rt * y.ln()
                + rt * pressure_bara.ln()
                - lagrange,
        );
    }
    Ok(values)
}

/// `getObjectiveVectorForVariables`: the species rows, then the active element balances.
fn objective_vector(
    objective: &[Option<f64>],
    variables: &[usize],
    active_elements: &[usize],
    element_in: &[f64; 7],
    moles: &[f64],
    species: &[Option<&GibbsSpecies>],
) -> Vec<f64> {
    // The element rows read `elementMoleBalanceDiff`, which the class keeps from the last
    // `updateSystemWithNewCompositions` - so it is the balance of the *current* moles.
    let element_out = element_balance(species, moles);
    let mut vector = Vec::with_capacity(variables.len() + active_elements.len());
    for &i in variables {
        vector.push(objective[i].unwrap_or(0.0));
    }
    for &element in active_elements {
        vector.push(element_out[element] - element_in[element]);
    }
    vector
}

/// `calculateJacobian`, in the class's own block order.
#[allow(clippy::too_many_arguments)] // The class's own argument list, less its state.
fn jacobian(
    mixture: &Mixture,
    database: &GibbsDatabase,
    species: &[Option<&GibbsSpecies>],
    moles: &[f64],
    variables: &[usize],
    active_elements: &[usize],
    total: f64,
    temperature: f64,
    ln_phi: &[f64],
    d_ln_phi_dn: &[Vec<f64>],
) -> Vec<Vec<f64>> {
    let _ = (mixture, database, temperature, ln_phi);
    let nv = variables.len();
    let na = active_elements.len();
    let size = nv + na;
    let mut matrix = vec![vec![0.0; size]; size];
    let rt = R_KJ * temperature;

    for (i, &component_i) in variables.iter().enumerate() {
        let ni = moles[component_i].max(MIN_JACOBIAN_MOLES);
        for (j, &component_j) in variables.iter().enumerate() {
            let dfugdn = d_ln_phi_dn
                .get(component_i)
                .and_then(|row| row.get(component_j))
                .copied()
                .unwrap_or(0.0);
            let composition = if i == j { 1.0 / ni } else { 0.0 };
            matrix[i][j] = rt * (composition - 1.0 / total + dfugdn);
        }
        if let Some(species) = species[component_i] {
            let elements = species.balance_elements();
            for (k, &element) in active_elements.iter().enumerate() {
                matrix[i][nv + k] = -elements[element];
            }
        }
    }

    for (i, &element) in active_elements.iter().enumerate() {
        for (j, &component_j) in variables.iter().enumerate() {
            if let Some(species) = species[component_j] {
                matrix[nv + i][j] = species.balance_elements()[element];
            }
        }
        // The multiplier columns of a balance row are zero, left as they were allocated.
    }
    matrix
}

/// `findActiveElements`: an element whose feed balance is non-zero *and* which some species
/// carries.
///
/// **The threshold is `ELEMENT_ZERO_THRESHOLD` here and `1e-10` in the class's other copy.**
/// `getActiveElementIndices` runs the same test with a literal `1e-10` for the coefficient,
/// and `performIterationUpdate` uses *that* one to size its half of the update vector. The
/// two can disagree only for a coefficient strictly between `1e-10` and `1e-6`, and every
/// element count this database carries is a whole number - so on this table they cannot, and
/// the port carries one function. Recorded rather than unified away: another table would make
/// the two disagree, and the class would then refuse its own delta vector as the wrong length.
fn active_elements(species: &[Option<&GibbsSpecies>], element_in: &[f64; 7]) -> Vec<usize> {
    (0..7)
        .filter(|&j| element_in[j].abs() > ELEMENT_ZERO_THRESHOLD)
        .filter(|&j| {
            species
                .iter()
                .flatten()
                .any(|species| species.balance_elements()[j].abs() > ELEMENT_ZERO_THRESHOLD)
        })
        .collect()
}

/// `calculateElementMoleBalance`, over the class's seven names.
fn element_balance(species: &[Option<&GibbsSpecies>], moles: &[f64]) -> [f64; 7] {
    let mut balance = [0.0; 7];
    for (i, entry) in species.iter().enumerate() {
        if let Some(species) = entry {
            let elements = species.balance_elements();
            for (j, value) in balance.iter_mut().enumerate() {
                *value += elements[j] * moles[i];
            }
        }
    }
    balance
}

/// `calculateMixtureGibbsEnergy`, kJ/mol: the total, not the molar one.
fn mixture_gibbs_energy(
    species: &[Option<&GibbsSpecies>],
    moles: &[f64],
    temperature: f64,
    cp: &[[f64; 4]],
) -> f64 {
    let mut total = 0.0;
    for (i, entry) in species.iter().enumerate() {
        if let Some(species) = entry {
            total += moles[i] * species.gibbs_energy(temperature, cp[i]);
        }
    }
    total
}

/// `calculateMixtureEnthalpy`, kJ/mol.
fn mixture_enthalpy(
    species: &[Option<&GibbsSpecies>],
    moles: &[f64],
    temperature: f64,
    cp: &[[f64; 4]],
) -> f64 {
    let mut total = 0.0;
    for (i, entry) in species.iter().enumerate() {
        if let Some(species) = entry {
            total += moles[i] * species.enthalpy(temperature, cp[i]);
        }
    }
    total
}

/// The fluid's total heat capacity in J/K, which is what `getCp("J/K")` returns.
///
/// The class's is extensive - it is the whole system's, not one mole's - and it is the *real*
/// capacity, carrying the departure, because it comes from the initialised fluid rather than
/// from the ideal-gas polynomial.
fn total_heat_capacity(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    moles: &[f64],
    temperature: f64,
    pressure_bara: f64,
) -> Result<f64> {
    let total: f64 = moles.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "moles",
            "a heat capacity needs moles to be per",
        ));
    }
    let x: Vec<f64> = moles.iter().map(|m| m / total).collect();
    let reduced = mixture.reduced_parameters(
        azoth_core::units::kelvins(temperature),
        pascals_of(pressure_bara),
    )?;
    let state = mixture.phase_state(&reduced, &x, RootSide::Vapour)?;
    let molar = molar_enthalpy_entropy(
        mixture,
        ideal_gas,
        azoth_core::units::kelvins(temperature),
        pascals_of(pressure_bara),
        &x,
        state.z,
    )?;
    Ok(molar.cp.value)
}

/// The fugacity coefficients and the composition derivative at a mole vector.
///
/// **The class reads the *gas* phase's coefficients** - `getPhase(0)` - and the capture records
/// the root: `phase0_type=gas` on all six rows, and `feed_phases=1` on all six, so the fluid is
/// one phase and the root is the largest admissible one.
fn phase_derivatives(
    mixture: &Mixture,
    moles: &[f64],
    temperature: f64,
    pressure_bara: f64,
) -> Result<(Vec<f64>, Vec<Vec<f64>>)> {
    let total: f64 = moles.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "moles",
            "a composition needs a positive total",
        ));
    }
    let x: Vec<f64> = moles.iter().map(|m| m / total).collect();
    let reduced = mixture.reduced_parameters(
        azoth_core::units::kelvins(temperature),
        pascals_of(pressure_bara),
    )?;
    let state = mixture.phase_state(&reduced, &x, RootSide::Vapour)?;
    let derivatives = mixture.phase_derivatives(&reduced, &x, state.z)?;
    Ok((derivatives.ln_phi, derivatives.d_ln_phi_dn))
}

/// The bar value the class's objective carries, as the pascal quantity azoth's surfaces take.
fn pascals_of(pressure_bara: f64) -> Pressure {
    azoth_core::units::pascals(pressure_bara * 1.0e5)
}
