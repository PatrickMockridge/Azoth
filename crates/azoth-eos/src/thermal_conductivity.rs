//! `eos.thermal_conductivity` - the liquid thermal conductivity from the Pedersen (PFCT)
//! corresponding-states correlation.
//!
//! Spec: `specs/models/eos/thermal_conductivity.toml`. A faithful port of NeqSim's
//! `PFCTConductivityMethodMod86.calcConductivity`: a methane SRK reference flash, the
//! Tc/Pc/molar-mass mixing rule, the dilute-gas methane viscosity, the ideal-gas heat
//! capacity, the reference-conductivity correlation and the corresponding-states scaling.

use azoth_core::units::{Pressure, ThermodynamicTemperature, watts_per_meter_kelvin};
use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::Mixture;
use crate::model_gen;
use crate::molar_enthalpy_entropy::IdealGasModel;
use crate::results::ThermalConductivityResult;
use crate::viscosity::{
    METHANE_M, METHANE_PC, METHANE_TC, NEQSIM_GAS_CONSTANT, methane_srk_density,
};

/// One atmosphere, Pa - the ideal-gas density's reference pressure.
const ATM: f64 = 101_325.0;

/// The reference critical molar density of `calcConductivity`, mol/dm³.
const CRIT_MOL_DENS: f64 = 10.152_119_7;

/// Methane's ideal-gas heat-capacity coefficients, J/(mol·K**n) - the reference `Cp0`.
const METHANE_CP_A: f64 = 37.978_352;
const METHANE_CP_B: f64 = -0.074_618_15;
const METHANE_CP_C: f64 = 0.000_301_881;
const METHANE_CP_D: f64 = -2.83e-7;
const METHANE_CP_E: f64 = 9.070_574e-11;

/// The conductivity `GVcoef` polynomial, from Friend, Ely and Ingham (1989).
const COND_GV: [f64; 9] = [
    -2.147_621e5,
    2.190_461e5,
    -8.618_097e4,
    1.496_099e4,
    -4.730_660e2,
    -2.331_178e2,
    3.778_439e1,
    -2.320_481,
    5.311_764e-2,
];
const COND_REF_A: f64 = -0.252_762_92;
const COND_REF_B: f64 = 0.334_328_59;
const COND_REF_C: f64 = 1.12;
const COND_REF_F: f64 = 168.0;
const COND_REF_J: [f64; 7] = [
    -7.040_363_399_07,
    12.319_512_908,
    -8.852_597_993_3e2,
    72.835_897_919,
    0.744_214_629_02,
    -2.970_691_454_0,
    2.220_975_850_1e3,
];
const COND_REF_K: [f64; 7] = [
    -8.551_09, 12.553_9, -1_020.85, 238.394, 1.315_63, -72.575_9, 1_411.6,
];

/// The dilute-gas viscosity `GVcoef` and `J` coefficients, reused by the conductivity for
/// the reference viscosity.
const VISC_GV: [f64; 9] = [
    -2.090_975e5,
    2.647_269e5,
    -1.472_818e5,
    4.716_740e4,
    -9.491_872e3,
    1.219_979e3,
    -9.627_993e1,
    4.274_152,
    -8.141_531e-2,
];
const VISC_REF_J: [f64; 7] = [
    -1.035_060_586e1,
    1.757_159_967_1e1,
    -3.019_391_865_6e3,
    1.887_301_159_4e2,
    4.290_360_948_8e-2,
    1.452_902_344_4e2,
    6.127_681_870_6e3,
];

/// Methane's ideal-gas heat capacity, J/(mol·K).
fn methane_cp0(t: f64) -> f64 {
    METHANE_CP_A
        + METHANE_CP_B * t
        + METHANE_CP_C * t.powi(2)
        + METHANE_CP_D * t.powi(3)
        + METHANE_CP_E * t.powi(4)
}

/// The dilute-gas methane viscosity at `t`, `getRefComponentViscosity`. Ideal-gas density
/// at one atmosphere, so the pressure argument NeqSim passes is unused.
fn dilute_gas_viscosity(t: f64) -> f64 {
    let mol_dens = ATM / NEQSIM_GAS_CONSTANT / t / 1e3;
    let red = (mol_dens - 10.15) / 10.15;
    let visc_ref_o = VISC_GV[0] * t.powi(-1)
        + VISC_GV[1] * t.powf(-2.0 / 3.0)
        + VISC_GV[2] * t.powf(-1.0 / 3.0)
        + VISC_GV[3]
        + VISC_GV[4] * t.powf(1.0 / 3.0)
        + VISC_GV[5] * t.powf(2.0 / 3.0)
        + VISC_GV[6] * t
        + VISC_GV[7] * t.powf(4.0 / 3.0)
        + VISC_GV[8] * t.powf(5.0 / 3.0);
    let temp_1 = mol_dens.powf(0.1) * (VISC_REF_J[1] + VISC_REF_J[2] / t.powf(1.5));
    let temp_2 =
        red * mol_dens.powf(0.5) * (VISC_REF_J[4] + VISC_REF_J[5] / t + VISC_REF_J[6] / t.powi(2));
    let temp_3 = (temp_1 + temp_2).exp();
    let htan = (t - 90.69).tanh();
    let vis_ref_e = (htan + 1.0) / 2.0;
    let visc_ref_2 = vis_ref_e * (VISC_REF_J[0] + VISC_REF_J[3] / t).exp() * (temp_3 - 1.0);
    let visc_ref_2 = if visc_ref_2.is_nan() { 0.0 } else { visc_ref_2 };
    (visc_ref_o + visc_ref_2) / 1.0e7
}

/// The methane reference conductivity, `getRefComponentConductivity`.
fn reference_conductivity(t: f64, p_pa: f64) -> Result<f64> {
    let rho = methane_srk_density(t, p_pa)?;
    let mol_dens_molar = rho / METHANE_M * 1e-3;
    let red = (mol_dens_molar - 10.15) / 10.15;
    let mol_dens = rho * 1e-3;

    let visc_ref_o = COND_GV[0] * t.powi(-1)
        + COND_GV[1] * t.powf(-2.0 / 3.0)
        + COND_GV[2] * t.powf(-1.0 / 3.0)
        + COND_GV[3]
        + COND_GV[4] * t.powf(1.0 / 3.0)
        + COND_GV[5] * t.powf(2.0 / 3.0)
        + COND_GV[6] * t
        + COND_GV[7] * t.powf(4.0 / 3.0)
        + COND_GV[8] * t.powf(5.0 / 3.0);
    let visc_ref_1 =
        (COND_REF_A + COND_REF_B * (COND_REF_C - (t / COND_REF_F).ln()).powi(2)) * mol_dens;

    let temp_1 = mol_dens.powf(0.1) * (COND_REF_J[1] + COND_REF_J[2] / t.powf(1.5));
    let temp_2 =
        red * mol_dens.powf(0.5) * (COND_REF_J[4] + COND_REF_J[5] / t + COND_REF_J[6] / t.powi(2));
    let temp_3 = (temp_1 + temp_2).exp();
    let htan = (t - 90.69).tanh();
    let cond_e = (htan + 1.0) / 2.0;
    let visc_ref_2 = cond_e * (COND_REF_J[0] + COND_REF_J[3] / t).exp() * (temp_3 - 1.0);
    let visc_ref_2 = if visc_ref_2.is_nan() { 0.0 } else { visc_ref_2 };

    let temp_4 = mol_dens.powf(0.1) * (COND_REF_K[1] + COND_REF_K[2] / t.powf(1.5));
    let temp_5 =
        red * mol_dens.powf(0.5) * (COND_REF_K[4] + COND_REF_K[5] / t + COND_REF_K[6] / t.powi(2));
    let temp_6 = (temp_4 + temp_5).exp();
    let cond_g = (1.0 - htan) / 2.0;
    let visc_ref_3 = cond_g * (COND_REF_K[0] + COND_REF_K[3] / t).exp() * (temp_6 - 1.0);
    let visc_ref_3 = if visc_ref_3.is_nan() { 0.0 } else { visc_ref_3 };

    let mut ref_cond = (visc_ref_o + visc_ref_1 + visc_ref_2 + visc_ref_3) * 1e-3;
    if t > 400.0 {
        let correction = (1.0 + (-9.0e-4) * (t - 400.0)).max(0.70);
        ref_cond *= correction;
    }
    Ok(ref_cond)
}

/// The liquid thermal conductivity of a mixture, from the Pedersen (PFCT) correlation.
///
/// `mixture` carries `Tc`, `Pc` and molar mass; `ideal_gas` the per-component heat-capacity
/// coefficients; `z` the mole fractions. The reference is pure methane, so the mixture's
/// cubic and alpha are not consulted.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is not a composition of the mixture's length.
/// * [`AzothError::PropertyUnavailable`] if a component carries no molar mass.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
pub fn thermal_conductivity(
    mixture: &Mixture,
    ideal_gas: &IdealGasModel,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<ThermalConductivityResult> {
    let spec = &model_gen::THERMAL_CONDUCTIVITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let components = mixture.components();
    let n = components.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a mixture of {n} components needs {n} mole fractions, got {}",
                z.len()
            ),
        ));
    }
    let molar_mass: Vec<f64> = components
        .iter()
        .map(|c| {
            c.molar_mass.ok_or_else(|| {
                AzothError::property_unavailable(
                    "component",
                    "molar mass",
                    "a card-added component needs its own molar mass",
                )
            })
        })
        .collect::<Result<_>>()?;

    // The Tc/Pc/molar-mass mixing rule, identical to the viscosity's.
    let mut temp_tc1 = 0.0;
    let mut temp_tc2 = 0.0;
    let mut temp_pc1 = 0.0;
    let mut temp_pc2 = 0.0;
    let mut mw_temp = 0.0;
    let mut mm_temp = 0.0;
    for i in 0..n {
        let tc_i = components[i].tc.value;
        let pc_i = components[i].pc.value;
        for j in 0..n {
            let tc_j = components[j].tc.value;
            let pc_j = components[j].pc.value;
            let temp_var = z[i]
                * z[j]
                * ((tc_i / pc_i).powf(1.0 / 3.0) + (tc_j / pc_j).powf(1.0 / 3.0)).powi(3);
            temp_tc1 += temp_var * (tc_i * tc_j).sqrt();
            temp_tc2 += temp_var;
            temp_pc1 += temp_var * (tc_i * tc_j).sqrt();
            temp_pc2 += temp_var;
        }
        mw_temp += z[i] * molar_mass[i].powi(2);
        mm_temp += z[i] * molar_mass[i];
    }
    let pc_mix = 8.0 * temp_pc1 / (temp_pc2 * temp_pc2);
    let tc_mix = temp_tc1 / temp_tc2;
    let m_mix =
        (mm_temp + 1.304e-4 * ((mw_temp / mm_temp).powf(2.303) - mm_temp.powf(2.303))) * 1e3;

    // Reference flash at the scaled state, then the per-component alpha mixing rule.
    let t_o_ref = t.value * METHANE_TC / tc_mix;
    let p_o_ref = p.value * METHANE_PC / pc_mix;
    let rho = methane_srk_density(t_o_ref, p_o_ref)?;
    let red_dens = (rho / METHANE_M * 1e-3) / CRIT_MOL_DENS;
    let alpha = |m: f64| 1.0 + 6.004e-4 * red_dens.powf(2.043) * (m * 1e3).powf(1.086);
    let mut alfa_mix = 0.0;
    for i in 0..n {
        for j in 0..n {
            alfa_mix += z[i] * z[j] * (alpha(molar_mass[i]) * alpha(molar_mass[j])).sqrt();
        }
    }
    let alfa0 = 1.0 + 6.004e-4 * red_dens.powf(2.043) * (METHANE_M * 1e3).powf(1.086);
    if alfa_mix < 1e-10 {
        return Ok(ThermalConductivityResult {
            k: watts_per_meter_kelvin(0.0),
            warnings,
        });
    }
    let t0 = t_o_ref * alfa0 / alfa_mix;
    let p0 = p_o_ref * alfa0 / alfa_mix;

    // The dilute-gas reference conductivity.
    let nstar_ref = dilute_gas_viscosity(t0);
    let f_func =
        1.0 + 0.053_432 * red_dens - 0.030_182 * red_dens.powi(2) - 0.029_725 * red_dens.powi(3);
    let cond_int_ref =
        1.186_53 * nstar_ref * (methane_cp0(t0) - 2.5 * NEQSIM_GAS_CONSTANT) * f_func / METHANE_M;

    // The mixture dilute-gas viscosity (`calcMixLPViscosity`), with its own ideal-gas
    // density and the alpha ratio inverted.
    let red_dens_2 = ATM / NEQSIM_GAS_CONSTANT / t.value / 1e3 / 10.15;
    let alpha_2 = |m: f64| 1.0 + 6.004e-4 * red_dens_2.powf(2.043) * (m * 1e3).powf(1.086);
    let mut alfa_mix_2 = 0.0;
    for i in 0..n {
        for j in 0..n {
            alfa_mix_2 += z[i] * z[j] * (alpha_2(molar_mass[i]) * alpha_2(molar_mass[j])).sqrt();
        }
    }
    let alfa0_2 = 1.0 + 6.004e-4 * red_dens_2.powf(2.043) * (METHANE_M * 1e3).powf(1.086);
    let t0_b = t.value * METHANE_TC / tc_mix * alfa_mix_2 / alfa0_2;
    let p0_b = ATM * METHANE_PC / pc_mix * alfa_mix_2 / alfa0_2;
    let _ = p0_b; // the dilute-gas viscosity does not read pressure
    let ref_visc_2 = dilute_gas_viscosity(t0_b);
    let nstar_mix = ref_visc_2
        * (tc_mix / METHANE_TC).powf(-1.0 / 6.0)
        * (pc_mix / METHANE_PC).powf(2.0 / 3.0)
        * (m_mix / (METHANE_M * 1e3)).powf(0.5)
        * alfa_mix_2
        / alfa0_2;

    let cp_id_mix: f64 = (0..n)
        .map(|i| {
            z[i] * (ideal_gas.cp_a[i]
                + ideal_gas.cp_b[i] * t.value
                + ideal_gas.cp_c[i] * t.value.powi(2)
                + ideal_gas.cp_d[i] * t.value.powi(3)
                + ideal_gas.cp_e[i] * t.value.powi(4))
        })
        .sum();
    let cond_int_mix =
        1.186_53 * nstar_mix * (cp_id_mix - 2.5 * NEQSIM_GAS_CONSTANT) * f_func / (m_mix / 1e3);

    let ref_conductivity = reference_conductivity(t0, p0)?;
    let conductivity = (tc_mix / METHANE_TC).powf(-1.0 / 6.0)
        * (pc_mix / METHANE_PC).powf(2.0 / 3.0)
        * (m_mix / (METHANE_M * 1e3)).powf(-0.5)
        * alfa_mix
        / alfa0
        * (ref_conductivity - cond_int_ref)
        + cond_int_mix;

    apply_checks(
        spec.derived_checks(),
        |name| (name == "k").then_some(conductivity),
        &mut warnings,
    )?;

    Ok(ThermalConductivityResult {
        k: watts_per_meter_kelvin(conductivity),
        warnings,
    })
}
