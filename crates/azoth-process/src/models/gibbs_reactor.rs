//! `process.gibbs_reactor` - the registry surface over the kernel.
//!
//! The kernel is [`crate::kernels::gibbs_reactor`], which takes a [`Stream`]; this is the same
//! arithmetic as an id the registry addresses, with the declared input list and the checks the
//! spec carries. **Two surfaces, one arithmetic**, for the reason `models`' own module doc gives:
//! the executor calls kernels and a case calls models, and neither wraps the other's types.

use azoth_core::units::{MolarEnergy, Pressure, ThermodynamicTemperature, kelvins, pascals};
use azoth_core::{CalcResult, Result, Warning, apply_checks};
use serde::Serialize;

use crate::executor::json::{scalar, warnings as wire_warnings};
use crate::kernels::gibbs_reactor::{
    EnergyMode, ReactorNumbers, ReactorSetup, gibbs_reactor as solve_equilibrium,
};
use crate::model_gen;
use crate::stream::Stream;

/// Result of `process.gibbs_reactor`.
#[derive(Debug, Clone, PartialEq, Serialize)]
pub struct GibbsReactorResult {
    /// Product molar flow, mol/s.
    pub product_n: f64,
    /// Product composition, over the feed's own species order.
    pub product_z: Vec<f64>,
    /// Product pressure.
    #[serde(serialize_with = "scalar")]
    pub product_p: Pressure,
    /// Product temperature.
    #[serde(serialize_with = "scalar")]
    pub product_t: ThermodynamicTemperature,
    /// Product molar enthalpy.
    #[serde(serialize_with = "scalar")]
    pub product_h: MolarEnergy,
    /// `hasConverged()`, which a run that reached the cap reports as false.
    pub converged: bool,
    /// The pass the loop stopped on.
    pub iterations: f64,
    /// The last undamped step norm.
    pub final_error: f64,
    /// The element Lagrange multipliers, J/mol, on the class's own seven element names.
    pub lagrange_multipliers: Vec<f64>,
    /// The outlet element balance less the inlet's, mol/s.
    pub element_balance_difference: Vec<f64>,
    /// The total Gibbs energy at the top of every iteration, W.
    pub gibbs_energy_history: Vec<f64>,
    /// Caveats.
    #[serde(serialize_with = "wire_warnings")]
    pub warnings: Vec<Warning>,
}

impl GibbsReactorResult {
    /// The result of one kernel call, from the solve's own two answers.
    ///
    /// **The warnings are the caller's, because the two callers have different ones.** A case
    /// runs the spec's checks through [`apply_checks`]; a flowsheet's equivalent is the
    /// *checker's*, which reports through the envelope's diagnostics rather than through a result.
    #[must_use]
    pub fn of(outlet: &Stream, numbers: &ReactorNumbers, warnings: Vec<Warning>) -> Self {
        Self {
            product_n: outlet.n,
            product_z: outlet.z.clone(),
            product_p: pascals(outlet.p.value),
            product_t: kelvins(outlet.t.value),
            product_h: outlet.h,
            converged: numbers.converged,
            iterations: f64::from(numbers.iterations),
            final_error: numbers.final_error,
            lagrange_multipliers: numbers
                .lagrange_multipliers
                .iter()
                .map(|v| v * 1000.0)
                .collect(),
            element_balance_difference: numbers.element_balance_difference.to_vec(),
            gibbs_energy_history: numbers
                .gibbs_energy_history
                .iter()
                .map(|v| v * 1000.0)
                .collect(),
            warnings,
        }
    }
}

impl CalcResult for GibbsReactorResult {
    const CALC_ID: &'static str = "process.gibbs_reactor";
    const FIELDS: &'static [&'static str] = &[
        "product_n",
        "product_z",
        "product_p",
        "product_t",
        "product_h",
        "converged",
        "iterations",
        "final_error",
        "lagrange_multipliers",
        "element_balance_difference",
        "gibbs_energy_history",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// Bring a feed to its Gibbs equilibrium at its own temperature and pressure.
///
/// **The class's own arithmetic is in kJ/mol and two of its outputs cross in other units.** The
/// multipliers are kJ/mol there and J/mol here, and the Gibbs energy history is a sum of `mol/s`
/// against `kJ/mol` - a rate of energy, reported in watts. The spec records both conversions.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] on a feed the spec's ranges refuse, on a fluid that
///   cannot be resolved, or on the singular Newton system the class answers with a pseudo-inverse
///   this port does not carry.
#[allow(clippy::too_many_arguments)] // one parameter per declared input, and there are ten
pub fn gibbs_reactor(
    components: &[String],
    feed_n: f64,
    feed_z: &[f64],
    feed_p: Pressure,
    feed_t: ThermodynamicTemperature,
    energy_mode: &str,
    damping_composition: f64,
    max_iterations: f64,
    convergence_tolerance: f64,
    min_iterations: f64,
) -> Result<GibbsReactorResult> {
    let spec = &model_gen::GIBBS_REACTOR_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "feed_n" => Some(feed_n),
            "damping_composition" => Some(damping_composition),
            "max_iterations" => Some(max_iterations),
            "convergence_tolerance" => Some(convergence_tolerance),
            _ => None,
        },
        &mut warnings,
    )?;

    let feed = Stream::from_pt(components.to_vec(), feed_z.to_vec(), feed_n, feed_p, feed_t)?;
    let setup = ReactorSetup {
        energy_mode: EnergyMode::named(energy_mode),
        damping_composition,
        max_iterations: max_iterations as u32,
        convergence_tolerance,
        min_iterations: min_iterations as u32,
    };
    let (outlet, numbers) = solve_equilibrium(&feed, &setup)?;

    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "iterations" => Some(f64::from(numbers.iterations)),
            "final_error" => Some(numbers.final_error),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(GibbsReactorResult::of(&outlet, &numbers, warnings))
}
