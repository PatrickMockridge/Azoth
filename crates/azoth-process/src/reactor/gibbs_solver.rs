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
//! **Two branches are refused rather than ported**, because the default path never takes them
//! and a branch nothing reaches is a claim nobody can check: the Armijo backtracking line
//! search and the Tikhonov regularisation, both `false` by default. `applyRegularization` is
//! *called* every iteration while `iteration <= 5` even with the flag off - it computes a
//! 2-norm condition number and stores it in a diagnostic history - and that is refused too,
//! since it consumes no state and modifies the Jacobian only when the flag is set.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result};
use azoth_eos::mixture::{Mixture, RootSide};
use azoth_eos::molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
use azoth_reactions::linalg::solve_lu;

use crate::reactor::gibbs_database::{GibbsDatabase, GibbsSpecies, R_KJ};

/// `MIN_MOLES`, the class's floor for a mole number in a logarithm.
const MIN_MOLES: f64 = 1e-6;

/// `MIN_JACOBIAN_MOLES`, the same floor where the Jacobian divides by it.
const MIN_JACOBIAN_MOLES: f64 = 1e-6;

/// `ELEMENT_ZERO_THRESHOLD`, the class's test for "this element is not present".
const ELEMENT_ZERO_THRESHOLD: f64 = 1e-6;

/// The floor `performIterationUpdate` puts under a composition after the damped step.
const MIN_COMPOSITION: f64 = 1e-15;

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
}

impl Default for GibbsSettings {
    fn default() -> Self {
        Self {
            adiabatic: false,
            damping: 0.05,
            max_iterations: 5000,
            tolerance: 1e-3,
            min_iterations: 100,
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
    let mut objective: Vec<Option<f64>> = vec![None; n];
    let mut converged = false;
    let mut iterations = 0u32;
    let mut final_error = f64::MAX;
    let mut enthalpy_old = 0.0;

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
        let jacobian = jacobian(
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
        let mut delta = match solve_lu(&jacobian, &negated) {
            Ok(delta) => delta,
            // **The class's fallback here is a pseudo-inverse and this port does not carry it.**
            // `solveNewtonSystem` catches the LU failure and retries through
            // `SimpleMatrix.pseudoInverse`, and the reachable case is a feed with fewer species
            // than elements - `steam_methane_hot` carries methane and water alone, so two
            // species meet three active balance rows and the system is over-determined. That
            // row is in the capture and converges there, so the fallback is load-bearing rather
            // than defensive; a port without it refuses the row rather than answering it
            // differently. `azoth-reactions`' `linalg` has `row_space_basis` and
            // `project_onto_null_space`, which is where a pseudo-inverse would be built.
            Err(error) => {
                return Err(AzothError::invalid_input(
                    "jacobian",
                    format!(
                        "the Newton system is singular at iteration {iteration} ({error}), and \
                         the class's fallback for that is `SimpleMatrix.pseudoInverse`, which \
                         this port does not carry - a feed with fewer species than active \
                         elements reaches it, and `steam_methane_hot` in the capture is one"
                    ),
                ));
            }
        };

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

        // `performIterationUpdate`: the compositions take `alpha`, the multipliers take all.
        let alpha = settings.damping;
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
    })
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
