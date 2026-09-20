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
pub mod ammonia;
pub mod ammonia_phase;
pub mod antoine_vapor_pressure;
pub mod argon_solid;
pub mod argon_solid_phase;
pub mod association;
pub mod bubble_pressure;
pub mod bubble_temperature;
pub mod bwrs;
pub mod bwrs_mixture;
pub mod bwrs_phase;
pub mod capillary_dew_point;
pub mod card;
pub mod chung_conductivity;
pub mod chung_viscosity;
pub mod co2_phase;
pub mod co2_water_diffusivity;
pub mod costald_molar_volume;
pub mod critical_point;
pub mod cubic;
pub mod databank;
pub mod dew_pressure;
pub mod dew_temperature;
pub mod dual;
pub mod electrolyte;
pub mod eos_cg;
mod eos_cg_data;
pub mod eos_cg_phase;
pub mod flash_iteration;
pub mod flash_newton;
pub mod flash_property;
pub mod ge_nrtl_flash;
pub mod ge_nrtl_phase;
pub mod ge_phase;
pub mod ge_unifac_phase;
pub mod ge_uniquac_phase;
pub mod ge_van_laar_acid_phase;
pub mod ge_wilson_phase;
pub mod gerg2008;
mod gerg2008_data;
pub mod gerg2008_phase;
pub mod hayduk_minhas_diffusivity;
pub mod heat_of_vaporization;
pub mod helium_phase;
pub mod henry;
pub mod hv_ge;
pub mod hydrogen_phase;
pub mod iapws_if97;
pub mod ideal_gas_cp;
pub mod leachman;
pub mod liquid_heat_capacity;
pub mod mason_saxena_conductivity;
pub mod matcop5_prumr_alpha;
pub mod matcop_alpha;
pub mod matcop_pr_alpha;
pub mod matcop_prumr_alpha;
pub mod matcop_prumr_new_alpha;
pub mod mixing_rule;
pub mod mixture;
pub mod model_gen;
pub mod molar_enthalpy_entropy;
pub mod mollerup_alpha;
pub mod multiphase;
pub mod nitric_sulfuric_acid_vapor_pressure;
pub mod nrtl_activity_coefficients;
pub mod parachor_surface_tension;
pub mod parahydrogen_solid;
pub mod parahydrogen_solid_phase;
pub mod pcsaft;
pub mod pcsaft_phase;
pub mod pcsaft_rahmat_phase;
pub mod ph_flash;
pub mod phase_boundary;
pub mod pitzer_catalog;
pub mod pitzer_phase;
pub mod pr78_kappa;
pub mod pr_alpha_ab;
pub mod pr_danesh_alpha;
pub mod pr_delft1998_alpha;
pub mod pr_departure;
pub mod pr_gassem2001_alpha;
pub mod pr_kappa;
pub mod pr_lee_kesler_alpha;
pub mod pr_mass_density;
pub mod pr_molar_volume;
pub mod pr_peneloux_shift;
pub mod pr_z_factor;
pub mod prsv_kappa;
pub mod ps_flash;
pub mod pt_flash;
pub mod pt_phase_envelope;
pub mod pu_flash;
pub mod pure_saturation;
pub mod pv_flash;
pub mod pv_reflux_flash;
pub mod pvf_flash;
pub mod rachford_rice;
pub mod rachford_rice_binary;
pub mod rackett_molar_volume;
pub mod results;
pub mod rk_alpha_ab;
pub mod rk_departure;
pub mod saft_vr_mie;
pub mod saft_vr_mie_phase;
pub mod saturation_temperature;
pub mod schwartzentruber_alpha;
pub mod siddiqi_lucas_diffusivity;
pub mod soreide_whitson_alpha;
pub mod spanwagner;
pub mod spec_gen;
pub mod srk_alpha_ab;
pub mod srk_cpa_phase;
pub mod srk_departure;
// The CPA phase state's shared arithmetic, private: the two models that use it are
// `srk_cpa_phase` and `pr_cpa_phase`, which name the family and the result type.
mod cpa_phase;
pub mod pr_cpa_phase;
pub mod srk_kappa;
pub mod srk_peneloux_shift;
pub mod srk_z_factor;
pub mod stability_test;
pub mod th_flash;
pub mod thermal_conductivity;
pub mod tp_flash_saft;
pub mod tp_multiflash;
pub mod ts_flash;
pub mod tu_flash;
pub mod tv_flash;
pub mod tv_fraction_flash;
pub mod twu_kappa;
pub mod twucoon_alpha;
pub mod twucoon_param_alpha;
pub mod twucoon_statoil_alpha;
pub mod tyn_calus_diffusivity;
pub mod umr_cpa_phase;
pub mod umrpr_alpha;
pub mod unifac_activity_coefficients;
pub mod unifac_psrk_activity_coefficients;
pub mod unifac_umrpru_activity_coefficients;
pub mod uniquac_activity_coefficients;
pub mod van_laar_acid_activity_coefficients;
pub mod vdw1f_mix_binary;
pub mod vega;
pub mod vh_flash;
pub mod viscosity;
pub mod volume_solve;
pub mod vs_flash;
pub mod vu_flash;
pub mod vu_flash_single_comp;
pub mod water_phase;
pub mod wilke_chang_diffusivity;
pub mod wilke_viscosity;
pub mod wilson_activity_coefficients;

pub use alpha_term::Alpha;
pub use ammonia_phase::ammonia_phase;
pub use antoine_vapor_pressure::{AntoineForm, antoine_vapor_pressure, form_from_type};
pub use argon_solid_phase::argon_solid_phase;
pub use bubble_pressure::bubble_pressure;
pub use bubble_temperature::bubble_temperature;
pub use bwrs_phase::bwrs_phase;
pub use capillary_dew_point::capillary_dew_point;
pub use chung_conductivity::chung_conductivity;
pub use chung_viscosity::chung_viscosity;
pub use co2_phase::co2_phase;
pub use co2_water_diffusivity::co2_water_diffusivity;
pub use costald_molar_volume::costald_molar_volume;
pub use critical_point::critical_point;
pub use cubic::Cubic;
pub use dew_pressure::dew_pressure;
pub use dew_temperature::dew_temperature;
pub use eos_cg_phase::eos_cg_phase;
pub use ge_nrtl_phase::ge_nrtl_phase;
pub use ge_unifac_phase::ge_unifac_phase;
pub use ge_uniquac_phase::ge_uniquac_phase;
pub use ge_van_laar_acid_phase::ge_van_laar_acid_phase;
pub use ge_wilson_phase::ge_wilson_phase;
pub use gerg2008_phase::gerg2008_phase;
pub use hayduk_minhas_diffusivity::{HaydukMinhasForm, hayduk_minhas_diffusivity};
pub use heat_of_vaporization::heat_of_vaporization;
pub use helium_phase::helium_phase;
pub use hydrogen_phase::hydrogen_phase;
pub use ideal_gas_cp::ideal_gas_cp;
pub use liquid_heat_capacity::liquid_heat_capacity;
pub use mason_saxena_conductivity::mason_saxena_conductivity;
pub use matcop_alpha::matcop_alpha;
pub use matcop_pr_alpha::matcop_pr_alpha;
pub use matcop_prumr_alpha::matcop_prumr_alpha;
pub use matcop_prumr_new_alpha::matcop_prumr_new_alpha;
pub use matcop5_prumr_alpha::matcop5_prumr_alpha;
pub use mixing_rule::{MixingRule, SoreideWhitsonRole};
pub use mixture::{Component, Mixture, PhaseState, ReducedParameters, RootSide};
pub use molar_enthalpy_entropy::{IdealGasModel, molar_enthalpy_entropy};
pub use mollerup_alpha::mollerup_alpha;
pub use nitric_sulfuric_acid_vapor_pressure::nitric_sulfuric_acid_vapor_pressure;
pub use nrtl_activity_coefficients::nrtl_activity_coefficients;
pub use parachor_surface_tension::parachor_surface_tension;
pub use parahydrogen_solid_phase::parahydrogen_solid_phase;
pub use pcsaft_rahmat_phase::pcsaft_rahmat_phase;
pub use ph_flash::ph_flash;
pub use phase_boundary::{Incipient, PhaseBoundary, phase_boundary_pressure};
pub use pr_alpha_ab::{OMEGA_A, OMEGA_B, pr_alpha_ab};
pub use pr_cpa_phase::{pr_cpa_phase, pr_cpa_phase_of};
pub use pr_danesh_alpha::pr_danesh_alpha;
pub use pr_delft1998_alpha::pr_delft1998_alpha;
pub use pr_departure::pr_departure;
pub use pr_gassem2001_alpha::pr_gassem2001_alpha;
pub use pr_kappa::pr_kappa;
pub use pr_lee_kesler_alpha::pr_lee_kesler_alpha;
pub use pr_mass_density::pr_mass_density;
pub use pr_molar_volume::{MOLAR_GAS_CONSTANT, pr_molar_volume};
pub use pr_peneloux_shift::pr_peneloux_shift;
pub use pr_z_factor::pr_z_factor;
pub use pr78_kappa::pr78_kappa;
pub use prsv_kappa::prsv_kappa;
pub use ps_flash::ps_flash;
pub use pt_flash::pt_flash;
pub use pt_phase_envelope::pt_phase_envelope;
pub use pu_flash::pu_flash;
pub use pure_saturation::pure_saturation;
pub use pv_flash::pv_flash;
pub use rachford_rice_binary::rachford_rice_binary;
pub use rackett_molar_volume::rackett_molar_volume;
pub use results::{
    AmmoniaPhaseResult, AntoineVaporPressureResult, ArgonSolidPhaseResult, BubblePressureResult,
    BubbleTemperatureResult, BwrsPhaseResult, ChungConductivityResult, ChungViscosityResult,
    Co2PhaseResult, Co2WaterDiffusivityResult, CostaldMolarVolumeResult, CriticalPointResult,
    DewPressureResult, DewTemperatureResult, EosCgPhaseResult, GeNrtlFlashResult,
    GeNrtlPhaseResult, GeUnifacPhaseResult, GeUniquacPhaseResult, GeVanLaarAcidPhaseResult,
    GeWilsonPhaseResult, Gerg2008PhaseResult, HaydukMinhasDiffusivityResult,
    HeatOfVaporizationResult, HeliumPhaseResult, HydrogenPhaseResult, IdealGasCpResult,
    LiquidHeatCapacityResult, MasonSaxenaConductivityResult, Matcop5PrumrAlphaResult,
    MatcopAlphaResult, MatcopPrAlphaResult, MatcopPrumrAlphaResult, MatcopPrumrNewAlphaResult,
    MolarEnthalpyEntropyResult, MollerupAlphaResult, NitricSulfuricAcidVaporPressureResult,
    NrtlActivityCoefficientsResult, ParachorSurfaceTensionResult, ParahydrogenSolidPhaseResult,
    Phase, Pr78KappaResult, PrAlphaAbResult, PrDaneshAlphaResult, PrDelft1998AlphaResult,
    PrDepartureResult, PrGassem2001AlphaResult, PrKappaResult, PrLeeKeslerAlphaResult,
    PrMassDensityResult, PrMolarVolumeResult, PrPenelouxShiftResult, PrZFactorResult,
    PrsvKappaResult, PsFlashResult, PtFlashResult, PtPhaseEnvelopeResult, PuFlashResult,
    PureSaturationResult, PvRefluxFlashResult, PvfFlashResult, RachfordRiceBinaryResult,
    RachfordRiceResult, RackettMolarVolumeResult, RkAlphaAbResult, RkDepartureResult,
    RootStructure, SchwartzentruberAlphaResult, SiddiqiLucasDiffusivityResult,
    SoreideWhitsonAlphaResult, SrkAlphaAbResult, SrkDepartureResult, SrkKappaResult,
    SrkPenelouxShiftResult, SrkZFactorResult, StabilityTestResult, StabilityVerdict, ThFlashResult,
    ThermalConductivityResult, TsFlashResult, TuFlashResult, TwuKappaResult, TwucoonAlphaResult,
    TwucoonParamAlphaResult, TwucoonStatoilAlphaResult, TynCalusDiffusivityResult,
    UmrCpaPhaseResult, UmrprAlphaResult, UnifacActivityCoefficientsResult,
    UnifacPsrkActivityCoefficientsResult, UnifacUmrpruActivityCoefficientsResult,
    UniquacActivityCoefficientsResult, VanLaarAcidActivityCoefficientsResult, Vdw1fMixBinaryResult,
    VhFlashResult, ViscosityResult, VuFlashResult, VuFlashSingleCompResult, WaterPhaseResult,
    WilkeChangDiffusivityResult, WilkeViscosityResult, WilsonActivityCoefficientsResult,
};
pub use rk_alpha_ab::rk_alpha_ab;
pub use rk_departure::rk_departure;
pub use schwartzentruber_alpha::schwartzentruber_alpha;
pub use siddiqi_lucas_diffusivity::{SiddiqiLucasForm, siddiqi_lucas_diffusivity};
pub use soreide_whitson_alpha::soreide_whitson_alpha;
pub use srk_alpha_ab::srk_alpha_ab;
pub use srk_cpa_phase::{srk_cpa_phase, srk_cpa_phase_of};
pub use srk_departure::srk_departure;
pub use srk_kappa::srk_kappa;
pub use srk_peneloux_shift::srk_peneloux_shift;
pub use srk_z_factor::srk_z_factor;
pub use stability_test::stability_test;
pub use th_flash::th_flash;
pub use thermal_conductivity::thermal_conductivity;
pub use tp_multiflash::tp_multiflash;
pub use ts_flash::ts_flash;
pub use tu_flash::tu_flash;
pub use tv_flash::tv_flash;
pub use twu_kappa::twu_kappa;
pub use twucoon_alpha::twucoon_alpha;
pub use twucoon_param_alpha::twucoon_param_alpha;
pub use twucoon_statoil_alpha::twucoon_statoil_alpha;
pub use tyn_calus_diffusivity::tyn_calus_diffusivity;
pub use umr_cpa_phase::umr_cpa_phase;
pub use umrpr_alpha::umrpr_alpha;
pub use unifac_activity_coefficients::unifac_activity_coefficients;
pub use unifac_psrk_activity_coefficients::unifac_psrk_activity_coefficients;
pub use unifac_umrpru_activity_coefficients::unifac_umrpru_activity_coefficients;
pub use uniquac_activity_coefficients::uniquac_activity_coefficients;
pub use van_laar_acid_activity_coefficients::van_laar_acid_activity_coefficients;
pub use vdw1f_mix_binary::vdw1f_mix_binary;
pub use viscosity::viscosity;
pub use vu_flash::vu_flash;
pub use water_phase::water_phase;
pub use wilke_chang_diffusivity::wilke_chang_diffusivity;
pub use wilke_viscosity::wilke_viscosity;
pub use wilson_activity_coefficients::wilson_activity_coefficients;
