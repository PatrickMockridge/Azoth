//! `azoth._core` - the compiled core, exposed to Python.
//!
//! # What this module is and is not
//!
//! It is a numerical core with a private name. The public Python API is
//! `azoth.hydraulics`, which decides whether to call this module or the pure
//! Python reference implementation and presents identical results either way.
//!
//! Two consequences worth stating, because they look like omissions:
//!
//! * **Arguments are SI magnitudes, not quantities.** Unit handling happens once,
//!   in Python, so the two backends cannot disagree about what a number is in.
//!   See `hydraulics.rs` and `thermal.rs`.
//! * **Exception classes are imported from `azoth.core.errors`, not defined
//!   here.** Both backends then raise the *same* class objects, so
//!   ``except OutOfRangeError`` works regardless of which one answered. See
//!   `errors.rs`.
//!
//! # Introspection
//!
//! `warning_codes`, `unit_names`, `solver_kinds`, `result_fields`, `calc_ids` and
//! `version` exist so the test suite can assert cross-language agreement without
//! parsing Rust source. They are the mechanism behind the claims that the two
//! implementations share a warning vocabulary, a unit vocabulary, a solver
//! vocabulary and a result shape.
//!
//! `data_files`, `fittings_rows` and `fluid_rows` do the same job for the data the
//! calcs are built from: the claim that both languages read the same tables is
//! asserted by comparing parsed rows and embedded bytes. See `data.rs`.

use pyo3::prelude::*;
use pyo3::types::PyModule;

mod batch;
mod data;
mod eos;
mod errors;
mod hydraulics;
mod overlay;
mod process;
mod reactions;
mod results;
mod standards;
mod thermal;

use results::{
    PyChokedFlowAreaResult, PyColebrookResult, PyConductionPlaneWallResult, PyControlValveCvResult,
    PyCriticalPointResult, PyDarcyWeisbachResult, PyFreezingPointResult, PyHaalandResult,
    PyHydrateFormationPressureResult, PyHydrateFormationTemperatureResult, PyHydrateFractionResult,
    PyKComponent, PyKFactorsResult, PyOrificeFlowResult, PyPhFlashResult, PyPr78KappaResult,
    PyPrAlphaAbResult, PyPrDepartureResult, PyPrKappaResult, PyPrMassDensityResult,
    PyPrMolarVolumeResult, PyPrZFactorResult, PyPrsvKappaResult, PyPsFlashResult,
    PyPumpPowerResult, PyPureSaturationResult, PyQty, PyRachfordRiceBinaryResult,
    PyReynoldsNumberResult, PyRkAlphaAbResult, PyRkDepartureResult, PySaltPrecipitationResult,
    PyScaleSaturationRatioResult, PySolidFugacityResult, PySrkAlphaAbResult, PySrkDepartureResult,
    PySrkKappaResult, PySrkZFactorResult, PySwameeJainResult, PyTbpFractionPropertiesResult,
    PyTpMultiflashWaxResult, PyTwuKappaResult, PyVdw1fMixBinaryResult, PyWarning,
    PyWaxSolidFugacityResult,
};

#[pymodule]
fn _core(py: Python<'_>, m: &Bound<'_, PyModule>) -> PyResult<()> {
    m.add(
        "__doc__",
        "Compiled core for azoth. Use azoth.hydraulics instead.",
    )?;

    // Transport and result types.
    m.add_class::<PyQty>()?;
    m.add_class::<PyWarning>()?;
    m.add_class::<PyKComponent>()?;
    m.add_class::<PyReynoldsNumberResult>()?;
    m.add_class::<PyColebrookResult>()?;
    m.add_class::<PySwameeJainResult>()?;
    m.add_class::<PyHaalandResult>()?;
    m.add_class::<PyConductionPlaneWallResult>()?;
    m.add_class::<PyPrKappaResult>()?;
    m.add_class::<PyPrAlphaAbResult>()?;
    m.add_class::<PyPrZFactorResult>()?;
    m.add_class::<PyPr78KappaResult>()?;
    m.add_class::<PyPrsvKappaResult>()?;
    m.add_class::<PyPrDepartureResult>()?;
    m.add_class::<PySrkKappaResult>()?;
    m.add_class::<PySrkAlphaAbResult>()?;
    m.add_class::<PySrkZFactorResult>()?;
    m.add_class::<PySrkDepartureResult>()?;
    m.add_class::<PyTwuKappaResult>()?;
    m.add_class::<PyRkAlphaAbResult>()?;
    m.add_class::<PyRkDepartureResult>()?;
    m.add_class::<PyVdw1fMixBinaryResult>()?;
    m.add_class::<PyRachfordRiceBinaryResult>()?;
    m.add_class::<PyPrMolarVolumeResult>()?;
    m.add_class::<PyPrMassDensityResult>()?;
    m.add_class::<PyPureSaturationResult>()?;
    m.add_class::<PyFreezingPointResult>()?;
    m.add_class::<PyHydrateFormationTemperatureResult>()?;
    m.add_class::<PyHydrateFractionResult>()?;
    m.add_class::<PyHydrateFormationPressureResult>()?;
    m.add_class::<PyTbpFractionPropertiesResult>()?;
    m.add_class::<PyWaxSolidFugacityResult>()?;
    m.add_class::<PyTpMultiflashWaxResult>()?;
    m.add_class::<PyScaleSaturationRatioResult>()?;
    m.add_class::<PySolidFugacityResult>()?;
    m.add_class::<PySaltPrecipitationResult>()?;
    m.add_class::<PyPhFlashResult>()?;
    m.add_class::<PyPsFlashResult>()?;
    m.add_class::<PyCriticalPointResult>()?;
    // An associating mixture's parameters, which every model whose Python side takes a
    // `Mixture` requires. A record rather than a defaulted argument: a caller that forgets
    // must get a `TypeError`, not a plausible answer for a different fluid.
    m.add_class::<eos::PyAssociationSpec>()?;
    m.add_class::<PyPumpPowerResult>()?;
    m.add_class::<PyOrificeFlowResult>()?;
    m.add_class::<PyControlValveCvResult>()?;
    m.add_class::<PyChokedFlowAreaResult>()?;

    // Data transport, for the cross-language data comparison.
    m.add_class::<batch::PyBatchColumn>()?;
    m.add_class::<batch::PyBatchResult>()?;
    m.add_class::<data::PyDataFile>()?;
    m.add_class::<data::PyFittingRow>()?;
    m.add_class::<data::PyFluidRow>()?;
    m.add_class::<data::PyComponentRow>()?;
    m.add_class::<data::PyKijRow>()?;
    m.add_class::<PyKFactorsResult>()?;
    m.add_class::<PyDarcyWeisbachResult>()?;

    // Calculations.
    m.add_function(wrap_pyfunction!(hydraulics::reynolds_number, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::friction_factor_colebrook, m)?)?;
    m.add_function(wrap_pyfunction!(
        hydraulics::friction_factor_swamee_jain,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(hydraulics::friction_factor_haaland, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::crane_k_factors, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::darcy_weisbach, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::packing_hydraulics, m)?)?;

    m.add_function(wrap_pyfunction!(hydraulics::pump_power, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::orifice_flow, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::control_valve_cv, m)?)?;
    m.add_function(wrap_pyfunction!(hydraulics::choked_flow_area, m)?)?;

    // Thermal calculations.
    m.add_function(wrap_pyfunction!(thermal::conduction_plane_wall, m)?)?;

    // Equations of state.
    m.add_function(wrap_pyfunction!(eos::pr_kappa, m)?)?;
    m.add_function(wrap_pyfunction!(eos::effective_diffusion, m)?)?;
    m.add_function(wrap_pyfunction!(reactions::equilibrium_constant, m)?)?;
    m.add_function(wrap_pyfunction!(reactions::reference_potentials, m)?)?;
    m.add_function(wrap_pyfunction!(reactions::chemical_equilibrium, m)?)?;
    m.add_function(wrap_pyfunction!(reactions::reactive_phase_equilibrium, m)?)?;
    m.add_function(wrap_pyfunction!(
        reactions::reactive_hybrid_eos_ge_flash,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(reactions::reactive_tp_flash, m)?)?;
    m.add_function(wrap_pyfunction!(reactions::reactive_ph_flash, m)?)?;
    m.add_function(wrap_pyfunction!(reactions::kinetic_rate_law, m)?)?;
    m.add_function(wrap_pyfunction!(reactions::kinetics, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_lee_kesler_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::matcop5_prumr_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::matcop_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::matcop_pr_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::matcop_prumr_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::matcop_prumr_new_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::mollerup_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_danesh_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_delft1998_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_gassem2001_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_alpha_ab, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_z_factor, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr78_kappa, m)?)?;
    m.add_function(wrap_pyfunction!(eos::prsv_kappa, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_departure, m)?)?;
    m.add_function(wrap_pyfunction!(eos::srk_kappa, m)?)?;
    m.add_function(wrap_pyfunction!(eos::srk_alpha_ab, m)?)?;
    m.add_function(wrap_pyfunction!(eos::srk_z_factor, m)?)?;
    m.add_function(wrap_pyfunction!(eos::srk_departure, m)?)?;
    m.add_function(wrap_pyfunction!(eos::rk_alpha_ab, m)?)?;
    m.add_function(wrap_pyfunction!(eos::rk_departure, m)?)?;
    m.add_function(wrap_pyfunction!(eos::twu_kappa, m)?)?;
    m.add_function(wrap_pyfunction!(eos::twucoon_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::twucoon_param_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::twucoon_statoil_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::vdw1f_mix_binary, m)?)?;
    m.add_function(wrap_pyfunction!(eos::rachford_rice, m)?)?;
    m.add_function(wrap_pyfunction!(eos::rachford_rice_binary, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_molar_volume, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_mass_density, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_peneloux_shift, m)?)?;
    m.add_function(wrap_pyfunction!(eos::srk_peneloux_shift, m)?)?;
    m.add_function(wrap_pyfunction!(eos::heat_of_vaporization, m)?)?;
    m.add_function(wrap_pyfunction!(eos::liquid_heat_capacity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::antoine_vapor_pressure, m)?)?;
    m.add_function(wrap_pyfunction!(
        eos::nitric_sulfuric_acid_vapor_pressure,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(eos::rackett_molar_volume, m)?)?;
    m.add_function(wrap_pyfunction!(eos::costald_molar_volume, m)?)?;
    m.add_function(wrap_pyfunction!(eos::chung_viscosity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::chapman_enskog_diffusivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::liquid_viscosity_pure, m)?)?;
    m.add_function(wrap_pyfunction!(eos::liquid_conductivity_polynom, m)?)?;
    m.add_function(wrap_pyfunction!(eos::phase_transport, m)?)?;
    m.add_function(wrap_pyfunction!(eos::chung_conductivity, m)?)?;

    // Models: the same shape, a different spec tree and generator.
    m.add_function(wrap_pyfunction!(eos::pure_saturation, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pt_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pt_phase_envelope, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ph_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ps_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tv_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tv_fraction_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pv_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::th_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ts_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tu_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pu_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pv_reflux_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pvf_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::vh_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::vs_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::vu_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::vu_flash_single_comp, m)?)?;
    m.add_function(wrap_pyfunction!(eos::stability_test, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tp_multiflash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::bubble_pressure, m)?)?;
    m.add_function(wrap_pyfunction!(eos::bubble_temperature, m)?)?;
    m.add_function(wrap_pyfunction!(eos::critical_point, m)?)?;
    m.add_function(wrap_pyfunction!(eos::dew_pressure, m)?)?;
    m.add_function(wrap_pyfunction!(eos::capillary_dew_point, m)?)?;
    m.add_function(wrap_pyfunction!(eos::dew_temperature, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ideal_gas_cp, m)?)?;
    m.add_function(wrap_pyfunction!(eos::molar_enthalpy_entropy, m)?)?;
    m.add_function(wrap_pyfunction!(eos::bwrs_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::srk_cpa_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pcsaft_rahmat_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::saft_vr_mie_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tp_flash_saft, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pr_cpa_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::umr_cpa_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::soreide_whitson_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::furst_electrolyte_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::furst_electrolyte_mod2004_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ammonia_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::co2_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::helium_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::freezing_point, m)?)?;
    m.add_function(wrap_pyfunction!(
        eos::fuller_schettler_giddings_diffusivity,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(eos::hydrate_formation_temperature, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hydrate_formation_pressure, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tbp_fraction_properties, m)?)?;
    m.add_function(wrap_pyfunction!(eos::wax_solid_fugacity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tp_multiflash_wax, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tp_solid_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::scale_saturation_ratio, m)?)?;
    m.add_function(wrap_pyfunction!(eos::solid_fugacity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::salt_precipitation, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hydrate_equilibrium_line, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hydrate_inhibitor_concentration, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hydrate_inhibitor_wt, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hydrate_fraction, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hydrogen_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hybrid_eos_ge_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::iapws_henry_law, m)?)?;
    m.add_function(wrap_pyfunction!(eos::water_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::argon_solid_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::parahydrogen_solid_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::eos_cg_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::gerg2008_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::viscosity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::aqueous_viscosity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::thermal_conductivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::wilke_viscosity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::mason_saxena_conductivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::nrtl_activity_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ge_nrtl_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ge_nrtl_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ge_flash, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ge_unifac_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ge_uniquac_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ge_van_laar_acid_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::ge_wilson_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::pitzer_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::kent_eisenberg_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::desmukh_mather_phase, m)?)?;
    m.add_function(wrap_pyfunction!(eos::umrpr_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::unifac_activity_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(eos::unifac_psrk_activity_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(
        eos::unifac_umrpru_activity_coefficients,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(
        eos::van_laar_acid_activity_coefficients,
        m
    )?)?;
    m.add_function(wrap_pyfunction!(eos::uniquac_activity_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(eos::wilson_activity_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(eos::tyn_calus_diffusivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::wilke_chang_diffusivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::hayduk_minhas_diffusivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::schwartzentruber_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::soreide_whitson_alpha, m)?)?;
    m.add_function(wrap_pyfunction!(eos::siddiqi_lucas_diffusivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::co2_water_diffusivity, m)?)?;
    m.add_function(wrap_pyfunction!(eos::parachor_mixture_surface_tension, m)?)?;
    m.add_function(wrap_pyfunction!(eos::parachor_surface_tension, m)?)?;

    // Introspection.
    m.add_function(wrap_pyfunction!(batch::batch_run, m)?)?;
    m.add_function(wrap_pyfunction!(data::data_files, m)?)?;
    m.add_function(wrap_pyfunction!(data::fittings_rows, m)?)?;
    m.add_function(wrap_pyfunction!(data::fluid_rows, m)?)?;
    m.add_function(wrap_pyfunction!(data::component_rows, m)?)?;
    m.add_function(wrap_pyfunction!(data::kij_rows, m)?)?;
    m.add_function(wrap_pyfunction!(results::warning_codes, m)?)?;
    m.add_class::<overlay::PyOverlay>()?;
    m.add_function(wrap_pyfunction!(overlay::overlay, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::card_overlay, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::card_coefficients, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_entry_row, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_component_rows, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_kij_rows, m)?)?;
    m.add_function(wrap_pyfunction!(overlay::overlay_cpa_kij_rows, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_names, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_dimensions, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_si_factor, m)?)?;
    m.add_function(wrap_pyfunction!(results::unit_slots, m)?)?;
    m.add_function(wrap_pyfunction!(results::solver_kinds, m)?)?;
    m.add_function(wrap_pyfunction!(eos::model_ids, m)?)?;
    m.add_function(wrap_pyfunction!(eos::model_schemes, m)?)?;
    m.add_function(wrap_pyfunction!(eos::model_kind, m)?)?;
    m.add_function(wrap_pyfunction!(results::result_fields, m)?)?;
    m.add_function(wrap_pyfunction!(results::calc_ids, m)?)?;
    m.add_function(wrap_pyfunction!(results::version, m)?)?;

    // The process layer: a stream value and the unit-operation kernels, plus the
    // flowsheet checker. A binding, not a second implementation.
    m.add_class::<process::PyStream>()?;
    m.add_function(wrap_pyfunction!(process::splitter_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::mixer_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::separator_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::throttling_valve_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::heat_exchanger_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::pump, m)?)?;
    m.add_function(wrap_pyfunction!(process::pump_stream, m)?)?;
    m.add_function(wrap_pyfunction!(process::heater, m)?)?;
    m.add_function(wrap_pyfunction!(process::cooler, m)?)?;
    m.add_function(wrap_pyfunction!(process::filter, m)?)?;
    m.add_function(wrap_pyfunction!(process::compressor, m)?)?;
    m.add_function(wrap_pyfunction!(process::expander, m)?)?;
    m.add_function(wrap_pyfunction!(process::pipe, m)?)?;
    m.add_function(wrap_pyfunction!(process::mixer, m)?)?;
    m.add_function(wrap_pyfunction!(process::separator, m)?)?;
    m.add_function(wrap_pyfunction!(process::gas_scrubber, m)?)?;
    m.add_function(wrap_pyfunction!(process::component_splitter, m)?)?;
    m.add_function(wrap_pyfunction!(process::shortcut_distillation_column, m)?)?;
    m.add_function(wrap_pyfunction!(process::distillation_column, m)?)?;
    m.add_function(wrap_pyfunction!(process::absorption_column, m)?)?;
    m.add_function(wrap_pyfunction!(process::stripping_column, m)?)?;
    m.add_function(wrap_pyfunction!(process::packed_column, m)?)?;
    m.add_function(wrap_pyfunction!(process::rate_based_packed_column, m)?)?;
    m.add_function(wrap_pyfunction!(process::heat_exchanger, m)?)?;
    m.add_function(wrap_pyfunction!(process::throttling_valve, m)?)?;
    m.add_function(wrap_pyfunction!(process::splitter, m)?)?;
    m.add_function(wrap_pyfunction!(process::manifold, m)?)?;
    m.add_function(wrap_pyfunction!(process::tank, m)?)?;
    m.add_function(wrap_pyfunction!(process::three_phase_separator, m)?)?;
    m.add_function(wrap_pyfunction!(process::ejector, m)?)?;
    m.add_function(wrap_pyfunction!(process::flare, m)?)?;
    m.add_function(wrap_pyfunction!(process::stirred_tank_reactor, m)?)?;
    m.add_function(wrap_pyfunction!(process::plug_flow_reactor, m)?)?;
    m.add_function(wrap_pyfunction!(process::gibbs_reactor, m)?)?;
    m.add_function(wrap_pyfunction!(process::validate_flowsheet, m)?)?;
    m.add_function(wrap_pyfunction!(process::run_flowsheet, m)?)?;
    m.add_function(wrap_pyfunction!(process::catalogue, m)?)?;
    m.add_class::<process::PySession>()?;

    // The standards namespace.
    m.add_function(wrap_pyfunction!(standards::iso6976, m)?)?;

    // Re-export the Python exception classes so both backends raise the same
    // objects rather than two lookalike hierarchies.
    errors::register(py, m)?;

    Ok(())
}
