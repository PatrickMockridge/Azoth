//! `unit_ops.gibbs_reactor` - the equilibrium composition of a fluid at its own state.
//!
//! ```text
//! minimise   G(n) = sum_i n_i [ Gf0_i(T) + RT ln phi_i + RT ln y_i + RT ln(P/1 bara) ]
//! subject to sum_i a_ij n_i = b_j     one row per active element
//! ```
//!
//! Spec: `specs/models/process/gibbs_reactor.toml`.
//!
//! **The class carries its own solver and its own data, and neither is P10's.** `GibbsReactor`
//! never calls `ChemicalEquilibrium`: it is a Lagrange-multiplier Newton iteration over the
//! element balances, against a species database whose formation properties are not the
//! databank's. A port composed from `reactions.chemical_equilibrium` would answer with different
//! numbers than the machine it ports, which is why the tier is a database, a solver and this.
//!
//! **The equilibrium temperature is the feed's, and there is no setter.** The class reads
//! `system.getTemperature()` from the fluid it is handed; the palette entry used to declare a
//! `temperature` parameter, and it does not have one.
//!
//! **Two of the class's behaviours are the reasons this tier exists and are in the solver**:
//! `ensureProductComponentsExist` is *not* one of them, because `useAllDatabaseSpecies` is
//! `false` by default and the variable set is the feed's own species; and a component the
//! database does not carry is still a variable, which is what makes the class's own
//! `testComponentNotInDatabaseMolesUnchanged` pass - see [`crate::reactor::gibbs_solver`].

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, Result};

use crate::reactor::gibbs_database::GibbsDatabase;
use crate::reactor::gibbs_solver::{GibbsSettings, GibbsState, solve};
use crate::stream::Stream;

/// How the reactor's temperature moves while the composition is solved.
///
/// The class's `EnergyMode`. **`Isothermal` is its default**, unlike the plug-flow reactor's -
/// so a Gibbs reactor with no stated mode holds its feed's temperature.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EnergyMode {
    /// The temperature is the feed's and does not move.
    Isothermal,
    /// The temperature moves with the reaction enthalpy inside the loop.
    Adiabatic,
}

impl EnergyMode {
    /// The mode a declaration's name selects, where any name but `adiabatic` is isothermal -
    /// the class's `setEnergyMode(String)` refuses an unknown name, and a declaration cannot.
    #[must_use]
    pub fn named(name: &str) -> Self {
        if name.eq_ignore_ascii_case("adiabatic") {
            Self::Adiabatic
        } else {
            Self::Isothermal
        }
    }
}

/// The class's controls, as its own fields default them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReactorSetup {
    /// Which energy mode the loop carries.
    pub energy_mode: EnergyMode,
    /// `dampingComposition`, the fraction of the Newton step the compositions take.
    pub damping_composition: f64,
    /// `maxIterations`.
    pub max_iterations: u32,
    /// `convergenceTolerance`, applied to the undamped step's norm.
    pub convergence_tolerance: f64,
    /// `minIterations`, the floor below which the tolerance is not consulted.
    pub min_iterations: u32,
}

impl Default for ReactorSetup {
    fn default() -> Self {
        Self {
            energy_mode: EnergyMode::Isothermal,
            damping_composition: 0.05,
            max_iterations: 5000,
            convergence_tolerance: 1e-3,
            min_iterations: 100,
        }
    }
}

/// What the solve reached, beyond the outlet stream: the trace and the balances.
///
/// **Carried because a Gibbs solve has a fixed point instead of a formula.** The outlet
/// composition alone cannot say whether the port reproduced the class - the capture's six rows
/// show the iteration count, the final error, the multipliers and the element balances, and
/// those are what a divergence shows up in first.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactorNumbers {
    /// `hasConverged()`, which is the tolerance test and not the solver's return value.
    pub converged: bool,
    /// `getActualIterations()`.
    pub iterations: u32,
    /// `getFinalConvergenceError()`, the last undamped step norm.
    pub final_error: f64,
    /// The seven Lagrange multipliers, on the class's own element names.
    pub lagrange_multipliers: [f64; 7],
    /// The element balance's inlet side, mol/s, on the class's seven element names.
    pub element_balance_in: [f64; 7],
    /// The outlet side.
    pub element_balance_out: [f64; 7],
    /// `outlet - inlet`, which the class refuses to call converged beyond a tolerance of its own.
    pub element_balance_difference: [f64; 7],
    /// The total Gibbs energy at the top of each iteration, kJ/mol.
    pub gibbs_energy_history: Vec<f64>,
    /// The outlet temperature, K.
    pub temperature: f64,
}

/// Bring a feed to its Gibbs equilibrium at its own temperature and pressure.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the mixture cannot be resolved, if the feed carries no
///   moles, or if the solve meets a state it cannot take a logarithm of.
/// * Whatever [`solve`] returns, including its refusal of a singular Newton system - the one
///   case the class answers through a pseudo-inverse this port does not carry.
pub fn gibbs_reactor(feed: &Stream, setup: &ReactorSetup) -> Result<(Stream, ReactorNumbers)> {
    if feed.n <= 0.0 {
        return Err(AzothError::invalid_input(
            "feed_n",
            "a reactor needs a positive molar flow",
        ));
    }
    let (mixture, ideal_gas) = feed.mixture()?;
    let database = GibbsDatabase::shipped()?;
    let inlet_moles: Vec<f64> = feed.z.iter().map(|z| z * feed.n).collect();

    let settings = GibbsSettings {
        adiabatic: setup.energy_mode == EnergyMode::Adiabatic,
        damping: setup.damping_composition,
        max_iterations: setup.max_iterations,
        tolerance: setup.convergence_tolerance,
        min_iterations: setup.min_iterations,
    };
    let state = solve(
        &mixture,
        &ideal_gas,
        &database,
        &inlet_moles,
        feed.t,
        feed.p.value / 1.0e5,
        &settings,
    )?;

    let outlet = outlet_stream(feed, &state)?;
    let numbers = ReactorNumbers {
        converged: state.converged,
        iterations: state.iterations,
        final_error: state.final_error,
        lagrange_multipliers: state.lambdas,
        element_balance_in: state.element_in,
        element_balance_out: state.element_out,
        element_balance_difference: state.element_diff,
        gibbs_energy_history: state.gibbs_history,
        temperature: state.temperature,
    };
    Ok((outlet, numbers))
}

/// The outlet stream: the same substances, the solved amounts, at the solved temperature.
///
/// **The enthalpy is the outlet state's, computed rather than carried.** The class's outlet
/// stream is run after the solve, so its `h` is a flash at the outlet state - which is what
/// [`Stream::from_pt`] evaluates.
fn outlet_stream(feed: &Stream, state: &GibbsState) -> Result<Stream> {
    let total: f64 = state.moles.iter().sum();
    if total <= 0.0 {
        return Err(AzothError::invalid_input(
            "moles",
            "the solve left no moles, so there is no outlet state",
        ));
    }
    let z: Vec<f64> = state.moles.iter().map(|moles| moles / total).collect();
    Stream::from_pt(
        feed.components.clone(),
        z,
        total,
        pascals(feed.p.value),
        kelvins(state.temperature),
    )
}
