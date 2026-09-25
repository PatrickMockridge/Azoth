//! The reactions namespace's extension surface.
//!
//! One function per registered `reactions.*` calc, thin: extract the SI magnitude,
//! call the kernel, wrap the result. The spec-driven behaviour lives in
//! `azoth-reactions`; what is here is the boundary.
//!
//! **The reaction's name crosses unresolved**, as a string. The Rust side reads the
//! table itself - which source, which row, which coefficients - so the Python caller
//! never handles a coefficient and the two languages cannot disagree about which row
//! answered.

use pyo3::prelude::*;

use azoth_reactions::chemical_equilibrium::ConcentrationBasis;
use azoth_reactions::databank::ReactionDataSource;
use azoth_reactions::reactive_phase_equilibrium::ReactionSeed;

use crate::errors::to_pyerr;
use crate::results::{
    PyEquilibriumConstantResult, PyKineticRateLawResult, PyKineticsResult,
    PyReactiveHybridEosGeFlashResult, PyReactivePhFlashResult, PyReactivePhaseEquilibriumResult,
    PyReactiveTpFlashResult,
};

/// One reaction's equilibrium constant, its derivative and its heat of reaction.
///
/// `source` is parsed rather than matched here, so a source this library does not carry
/// is refused by the same code path the Rust tests exercise.
#[pyfunction]
#[pyo3(signature = (source, reaction, T))]
#[pyo3(text_signature = "(source, reaction, T)")]
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn equilibrium_constant(
    py: Python<'_>,
    source: &str,
    reaction: &str,
    T: f64,
) -> PyResult<PyEquilibriumConstantResult> {
    let parsed: ReactionDataSource = source.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::equilibrium_constant::equilibrium_constant(
        parsed,
        reaction,
        azoth_core::units::kelvins(T),
    )
    .map(|r| PyEquilibriumConstantResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The standard-state reference potentials of a fluid's reactive components.
///
/// `components` crosses as names and unresolved, for the same reason the reaction name
/// does: the Rust side reads the tables, chooses the basis and propagates, so no
/// stoichiometric coefficient and no rank decision reaches Python.
#[pyfunction]
#[pyo3(signature = (components, source, T))]
#[pyo3(text_signature = "(components, source, T)")]
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn reference_potentials(
    py: Python<'_>,
    components: Vec<String>,
    source: &str,
    T: f64,
) -> PyResult<crate::results::PyReferencePotentialsResult> {
    let parsed: ReactionDataSource = source.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::reference_potentials::reference_potentials(
        &components,
        parsed,
        azoth_core::units::kelvins(T),
    )
    .map(|r| crate::results::PyReferencePotentialsResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The reactive equilibrium composition of a phase.
///
/// **The matrix crosses nested**, so the parameter list is exactly the spec's declared
/// inputs - a flat vector plus a row count would be two parameters where the spec has one,
/// and the generated stub would disagree with this signature.
///
/// `whole_system` is NeqSim's `getNumberOfPhases() == 1`, which decides whether the solve
/// corrects its conservation coupling to the phase's own element amounts.
#[pyfunction]
#[pyo3(signature = (a_matrix, b, whole_system, moles, chem_ref, log_activity, T, max_iterations, tolerance, concentration_basis, solvent_weight, solvent_mask, phase_moles))]
#[pyo3(
    text_signature = "(a_matrix, b, whole_system, moles, chem_ref, log_activity, T, max_iterations, tolerance, concentration_basis, solvent_weight, solvent_mask, phase_moles)"
)]
#[allow(non_snake_case)] // `T` is the symbol in the chemistry
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are thirteen
pub fn chemical_equilibrium(
    py: Python<'_>,
    a_matrix: Vec<Vec<f64>>,
    b: Vec<f64>,
    whole_system: bool,
    moles: Vec<f64>,
    chem_ref: Vec<f64>,
    log_activity: Vec<f64>,
    T: f64,
    max_iterations: u32,
    tolerance: f64,
    concentration_basis: &str,
    solvent_weight: f64,
    solvent_mask: Vec<f64>,
    phase_moles: f64,
) -> PyResult<crate::results::PyChemicalEquilibriumResult> {
    let parsed_basis: ConcentrationBasis =
        concentration_basis.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::chemical_equilibrium::chemical_equilibrium(
        &a_matrix,
        &b,
        whole_system,
        &moles,
        &chem_ref,
        &log_activity,
        T,
        max_iterations,
        tolerance,
        parsed_basis,
        solvent_weight,
        &solvent_mask,
        phase_moles,
    )
    .map(|r| crate::results::PyChemicalEquilibriumResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The phase's reactive equilibrium, as the operation the facade is.
///
/// **`skipped` is part of the result and not an error.** A phase that is neither aqueous
/// nor liquid nor oil has no phase for a water-based equilibrium to be solved in, and the
/// composition comes back as it went in - so the flag is what separates that from a solve
/// that ran and did not converge.
#[pyfunction]
#[pyo3(
    signature = (components, source, phase, moles, phase_charge, phase_moles, whole_system, log_activity, T, max_iterations, tolerance, concentration_basis, seed)
)]
#[pyo3(
    text_signature = "(components, source, phase, moles, phase_charge, phase_moles, whole_system, log_activity, T, max_iterations, tolerance, concentration_basis, seed)"
)]
#[allow(non_snake_case)] // `T` is the symbol in the chemistry
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are thirteen
pub fn reactive_phase_equilibrium(
    py: Python<'_>,
    components: Vec<String>,
    source: &str,
    phase: &str,
    moles: Vec<f64>,
    phase_charge: f64,
    phase_moles: f64,
    whole_system: bool,
    log_activity: Vec<f64>,
    T: f64,
    max_iterations: u32,
    tolerance: f64,
    concentration_basis: &str,
    seed: &str,
) -> PyResult<PyReactivePhaseEquilibriumResult> {
    let parsed: ReactionDataSource = source.parse().map_err(|e| to_pyerr(py, e))?;
    let parsed_basis: ConcentrationBasis =
        concentration_basis.parse().map_err(|e| to_pyerr(py, e))?;
    let parsed_seed: ReactionSeed = seed.parse().map_err(|e| to_pyerr(py, e))?;
    azoth_reactions::reactive_phase_equilibrium::reactive_phase_equilibrium(
        &components,
        parsed,
        phase,
        &moles,
        phase_charge,
        phase_moles,
        whole_system,
        &log_activity,
        T,
        max_iterations,
        tolerance,
        parsed_seed,
        parsed_basis,
    )
    .map(|r| PyReactivePhaseEquilibriumResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The reactive flash: simultaneous chemical and phase equilibrium at fixed `T` and `P`.
///
/// `moles` is the **overall component amounts and not a composition**, which is what the
/// driver's `getOverallMoles` reads and what its frozen element inventory is built from.
/// A charged component is refused: the RAND solve's ionic branch is not ported.
#[pyfunction]
#[pyo3(signature = (components, T, P, moles, max_phases, cubic = "srk"))]
#[pyo3(text_signature = "(components, T, P, moles, max_phases, cubic = \"srk\")")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the flash's own name
pub fn reactive_tp_flash(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    moles: Vec<f64>,
    max_phases: f64,
    cubic: &str,
) -> PyResult<PyReactiveTpFlashResult> {
    azoth_reactions::reactive_tp_flash::reactive_tp_flash(
        &components,
        cubic.parse().map_err(pyo3::exceptions::PyValueError::new_err)?,
        T,
        P,
        &moles,
        max_phases as usize,
    )
    .map(|r| PyReactiveTpFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The coupled reactive hybrid flash: chemistry and a fixed gas-oil-brine topology at once.
///
/// `components` crosses **unresolved**: the two EoS roles' constants, the seeding's classes and
/// the brine's ion mask all resolve on the Rust side from the same databank, so the two
/// languages cannot disagree about which substance is which. `moles` is the feed's own mole
/// numbers, and the inventory the loop's later passes are solved at is not it - a reaction
/// changes the number of species moles, and `coupled_moles` is what came back.
#[pyfunction]
#[pyo3(signature = (components, cubic, T, P, moles))]
#[pyo3(text_signature = "(components, cubic, T, P, moles)")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the flash's own name
pub fn reactive_hybrid_eos_ge_flash(
    py: Python<'_>,
    components: Vec<String>,
    cubic: String,
    T: f64,
    P: f64,
    moles: Vec<f64>,
) -> PyResult<PyReactiveHybridEosGeFlashResult> {
    azoth_reactions::reactive_hybrid_eos_ge_flash::reactive_hybrid_eos_ge_flash(
        &components,
        cubic
            .parse()
            .map_err(pyo3::exceptions::PyValueError::new_err)?,
        T,
        P,
        &moles,
    )
    .map(|r| PyReactiveHybridEosGeFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The reactive PH flash: the temperature a specified enthalpy asks for.
///
/// `enthalpy` is the **thermochemical** specification - the fluid's sensible enthalpy plus the
/// formation inventory - and `T` is where the secant search starts.
#[pyfunction]
#[pyo3(signature = (components, T, P, moles, enthalpy, max_phases, cubic = "srk"))]
#[pyo3(text_signature = "(components, T, P, moles, enthalpy, max_phases, cubic = \"srk\")")]
#[allow(non_snake_case)] // `T` and `P` are the symbols in the flash's own name
pub fn reactive_ph_flash(
    py: Python<'_>,
    components: Vec<String>,
    T: f64,
    P: f64,
    moles: Vec<f64>,
    enthalpy: f64,
    max_phases: f64,
    cubic: &str,
) -> PyResult<PyReactivePhFlashResult> {
    azoth_reactions::reactive_ph_flash::reactive_ph_flash(
        &components,
        cubic.parse().map_err(pyo3::exceptions::PyValueError::new_err)?,
        T,
        P,
        &moles,
        enthalpy,
        max_phases as usize,
    )
    .map(|r| PyReactivePhFlashResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// A reaction's kinetic rate factor, by the law its selector names.
///
/// `law` is parsed by the kernel rather than matched here, so a branch this library does
/// not carry is refused by the same code path the Rust tests exercise.
#[pyfunction]
#[pyo3(signature = (law, T, reference_rate, activation_energy, reference_temperature))]
#[pyo3(text_signature = "(law, T, reference_rate, activation_energy, reference_temperature)")]
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn kinetic_rate_law(
    py: Python<'_>,
    law: &str,
    T: f64,
    reference_rate: f64,
    activation_energy: f64,
    reference_temperature: f64,
) -> PyResult<PyKineticRateLawResult> {
    azoth_reactions::kinetic_rate_law::kinetic_rate_law(
        law,
        T,
        reference_rate,
        activation_energy,
        reference_temperature,
    )
    .map(|r| PyKineticRateLawResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}

/// The Krishna-Standart mass-transfer rate matrix of one phase's reactions.
///
/// The reactions cross as **concatenated name lists** rather than a matrix over the
/// phase's species, because a reaction's own name order decides which sibling's
/// `phiInfinite` it keeps - and the coefficient a matrix would encode is not enough to
/// recover that order.
#[pyfunction]
#[pyo3(signature = (
    components,
    reaction_components,
    reaction_lengths,
    reaction_coefficients,
    rate_factors,
    equilibrium_constants,
    fractions,
    molar_masses,
    density,
    inter_fractions,
    inter_density,
    diffusion
))]
#[pyo3(
    text_signature = "(components, reaction_components, reaction_lengths, \
                        reaction_coefficients, rate_factors, equilibrium_constants, fractions, \
                        molar_masses, density, inter_fractions, inter_density, diffusion)"
)]
#[allow(clippy::too_many_arguments)] // the class's own parameter list, one input per fact
pub fn kinetics(
    py: Python<'_>,
    components: Vec<String>,
    reaction_components: Vec<String>,
    reaction_lengths: Vec<f64>,
    reaction_coefficients: Vec<f64>,
    rate_factors: Vec<f64>,
    equilibrium_constants: Vec<f64>,
    fractions: Vec<f64>,
    molar_masses: Vec<f64>,
    density: f64,
    inter_fractions: Vec<f64>,
    inter_density: f64,
    diffusion: Vec<f64>,
) -> PyResult<PyKineticsResult> {
    let names: Vec<&str> = components.iter().map(String::as_str).collect();
    let reaction_names: Vec<&str> = reaction_components.iter().map(String::as_str).collect();
    azoth_reactions::kinetics::kinetics(
        &names,
        &reaction_names,
        &reaction_lengths,
        &reaction_coefficients,
        &rate_factors,
        &equilibrium_constants,
        &fractions,
        &molar_masses,
        density,
        &inter_fractions,
        inter_density,
        &diffusion,
    )
    .map(|r| PyKineticsResult::from(&r))
    .map_err(|e| to_pyerr(py, e))
}
