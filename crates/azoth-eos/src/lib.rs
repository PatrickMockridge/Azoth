//! Equations of state.
//!
//! The cubic equation of state, decomposed into constitutive coefficients rather than
//! into one calculation per equation: each coefficient has its own spec and worked
//! example under `specs/calcs/eos/`, and the models that compose them - procedures
//! rather than equations, with a scheme, a tolerance and an iteration cap - are under
//! `specs/models/eos/`. `tools/gen_models.py` generates those into [`model_gen`], the
//! calculations into `spec_gen`.
//!
//! Everything here works in reduced variables, which is why no `uom` type appears in
//! this crate: the dimensional conversion belongs at the boundary, in one place.

use azoth_core::{AzothError, ModelAlgorithm, ModelSpec, Result};

/// The algorithm a model's spec fixes, which its `kind` says it has.
///
/// The schema requires an `algorithm` block for a `procedure` and forbids one for a
/// `direct` model, so a missing one where a procedure is expected means the generated
/// table and the schema disagree - a generator defect rather than a caller error.
/// Returned rather than panicked, per this crate's no-panic rule.
///
/// Here rather than in each of the models that need it, and `pub` because it is part
/// of the model layer's contract rather than of any one model.
pub fn algorithm_of(spec: &ModelSpec) -> Result<&'static ModelAlgorithm> {
    spec.algorithm.ok_or_else(|| AzothError::InvalidInput {
        field: "algorithm".to_string(),
        reason: format!(
            "model `{}` is declared `{}` but its spec carries no algorithm, so there is \
             no procedure to run",
            spec.id, spec.kind
        ),
    })
}

pub mod alpha_term;
pub mod antoine_vapor_pressure;
pub mod bubble_pressure;
pub mod card;
pub mod chung_viscosity;
pub mod costald_molar_volume;
pub mod critical_point;
pub mod cubic;
pub mod databank;
pub mod dew_pressure;
pub mod flash_property;
pub mod heat_of_vaporization;
pub mod ideal_gas_cp;
pub mod liquid_heat_capacity;
pub mod mixture;
pub mod model_gen;
pub mod molar_enthalpy_entropy;
pub mod ph_flash;
pub mod phase_boundary;
pub mod pr78_kappa;
pub mod pr_alpha_ab;
pub mod pr_departure;
pub mod pr_kappa;
pub mod pr_mass_density;
pub mod pr_molar_volume;
pub mod pr_peneloux_shift;
pub mod pr_z_factor;
pub mod prsv_kappa;
pub mod ps_flash;
pub mod pt_flash;
pub mod pure_saturation;
pub mod rachford_rice_binary;
pub mod rackett_molar_volume;
pub mod results;
pub mod rk_alpha_ab;
pub mod rk_departure;
pub mod spec_gen;
pub mod srk_alpha_ab;
pub mod srk_departure;
pub mod srk_kappa;
pub mod srk_peneloux_shift;
pub mod srk_z_factor;
pub mod stability_test;
pub mod twu_kappa;
pub mod vdw1f_mix_binary;
pub mod wilke_viscosity;

pub use alpha_term::Alpha;
pub use antoine_vapor_pressure::{AntoineForm, antoine_vapor_pressure, form_from_type};
pub use bubble_pressure::bubble_pressure;
pub use chung_viscosity::chung_viscosity;
pub use costald_molar_volume::costald_molar_volume;
pub use critical_point::critical_point;
pub use cubic::Cubic;
pub use dew_pressure::dew_pressure;
pub use heat_of_vaporization::heat_of_vaporization;
pub use ideal_gas_cp::ideal_gas_cp;
pub use liquid_heat_capacity::liquid_heat_capacity;
pub use mixture::{Component, Mixture, PhaseState, ReducedParameters, RootSide};
pub use molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
pub use ph_flash::ph_flash;
pub use phase_boundary::{Incipient, PhaseBoundary, phase_boundary_pressure};
pub use pr_alpha_ab::{OMEGA_A, OMEGA_B, pr_alpha_ab};
pub use pr_departure::pr_departure;
pub use pr_kappa::pr_kappa;
pub use pr_mass_density::pr_mass_density;
pub use pr_molar_volume::{MOLAR_GAS_CONSTANT, pr_molar_volume};
pub use pr_peneloux_shift::pr_peneloux_shift;
pub use pr_z_factor::pr_z_factor;
pub use pr78_kappa::pr78_kappa;
pub use prsv_kappa::prsv_kappa;
pub use ps_flash::ps_flash;
pub use pt_flash::pt_flash;
pub use pure_saturation::pure_saturation;
pub use rachford_rice_binary::rachford_rice_binary;
pub use rackett_molar_volume::rackett_molar_volume;
pub use results::{
    AntoineVaporPressureResult, BubblePressureResult, ChungViscosityResult,
    CostaldMolarVolumeResult, CriticalPointResult, DewPressureResult, HeatOfVaporizationResult,
    IdealGasCpResult, LiquidHeatCapacityResult, MolarEnthalpyEntropyResult, Phase, Pr78KappaResult,
    PrAlphaAbResult, PrDepartureResult, PrKappaResult, PrMassDensityResult, PrMolarVolumeResult,
    PrPenelouxShiftResult, PrZFactorResult, PrsvKappaResult, PtFlashResult, PureSaturationResult,
    RachfordRiceBinaryResult, RackettMolarVolumeResult, RkAlphaAbResult, RkDepartureResult,
    RootStructure, SrkAlphaAbResult, SrkDepartureResult, SrkKappaResult, SrkPenelouxShiftResult,
    SrkZFactorResult, StabilityTestResult, StabilityVerdict, TwuKappaResult, Vdw1fMixBinaryResult,
    WilkeViscosityResult,
};
pub use rk_alpha_ab::rk_alpha_ab;
pub use rk_departure::rk_departure;
pub use srk_alpha_ab::srk_alpha_ab;
pub use srk_departure::srk_departure;
pub use srk_kappa::srk_kappa;
pub use srk_peneloux_shift::srk_peneloux_shift;
pub use srk_z_factor::srk_z_factor;
pub use stability_test::stability_test;
pub use twu_kappa::twu_kappa;
pub use vdw1f_mix_binary::vdw1f_mix_binary;
pub use wilke_viscosity::wilke_viscosity;
