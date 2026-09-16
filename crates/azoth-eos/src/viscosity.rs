//! `eos.viscosity` - the liquid viscosity from the Pedersen (PFCT) heavy-oil
//! corresponding-states correlation.
//!
//! Spec: `specs/models/eos/viscosity.toml`. This is a faithful port of NeqSim's
//! `PFCTViscosityMethodHeavyOil.calcViscosity`: a methane SRK reference flash, a
//! Tc/Pc/molar-mass mixing rule, the `getRefComponentViscosity` correlation and the
//! corresponding-states scaling.
//!
//! The reference is pure methane, solved with the SRK pieces this crate already
//! registers ([`crate::srk_kappa`], [`crate::srk_alpha_ab`], [`crate::srk_z_factor`],
//! [`crate::srk_departure`]). The correlation is exponentially sensitive to the
//! reference density, so the reference flash is composed here rather than approximated.

use azoth_core::units::{Pressure, ThermodynamicTemperature, pascal_seconds};
use azoth_core::{AzothError, Result, apply_checks};

use crate::mixture::Mixture;
use crate::model_gen;
use crate::results::ViscosityResult;
use crate::{srk_alpha_ab, srk_departure, srk_kappa, srk_z_factor};

/// NeqSim's gas constant, J/(mol·K). The reference density `rho = M*P/(Z*R*T)` is the
/// one place `R` does not cancel out of the cubic, and NeqSim's `getDensity()` uses its
/// own (truncated) value, so matching the oracle means using it rather than
/// [`crate::pr_molar_volume::MOLAR_GAS_CONSTANT`].
pub(crate) const NEQSIM_GAS_CONSTANT: f64 = 8.314_462_1;

/// Methane's critical temperature, K - the reference component's.
pub(crate) const METHANE_TC: f64 = 190.56;
/// Methane's critical pressure, Pa.
pub(crate) const METHANE_PC: f64 = 4.599e6;
/// Methane's molar mass, kg/mol.
pub(crate) const METHANE_M: f64 = 0.016043;
/// Methane's acentric factor.
pub(crate) const METHANE_OMEGA: f64 = 0.0115;

/// The `GVcoef` polynomial of `PFCTViscosityMethodHeavyOil`, in the order it is read.
const GVCOEF: [f64; 9] = [
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
const VIS_REF_A: f64 = 1.696_985_927;
const VIS_REF_B: f64 = -0.133_372_346;
const VIS_REF_C: f64 = 1.4;
const VIS_REF_F: f64 = 168.0;
const VISC_REF_J: [f64; 7] = [
    -1.035_060_586e1,
    1.757_159_967_1e1,
    -3.019_391_865_6e3,
    1.887_301_159_4e2,
    4.290_360_948_8e-2,
    1.452_902_344_4e2,
    6.127_681_870_6e3,
];
const VISC_REF_K: [f64; 7] = [
    -9.746_02, 18.083_4, -4_126.66, 44.605_5, 0.976_544, 81.813_4, 15_649.9,
];
/// The reference critical molar density, mol/dm³.
const CRIT_MOL_DENS: f64 = 10.15;

/// Methane's SRK density in kg/m³, the stable phase's, at `(t, p)`.
///
/// The stable root is the one with the lower `ln_phi` - lower fugacity is lower Gibbs,
/// which is what NeqSim's `getLowestGibbsEnergyPhase()` selects. `srk_z_factor`
/// already discards the middle (unstable) root, so only `z_min` and `z_max` remain.
pub(crate) fn methane_srk_density(t: f64, p_pa: f64) -> Result<f64> {
    let kappa = srk_kappa(METHANE_OMEGA)?.kappa;
    let tr = t / METHANE_TC;
    let pr = p_pa / METHANE_PC;
    let ab = srk_alpha_ab(kappa, tr, pr)?;
    let z = srk_z_factor(ab.a_reduced, ab.b_reduced)?;
    let phi_min = srk_departure(ab.a_reduced, ab.b_reduced, z.z_min, kappa, tr)?.ln_phi;
    let phi_max = srk_departure(ab.a_reduced, ab.b_reduced, z.z_max, kappa, tr)?.ln_phi;
    let z_stable = if phi_min <= phi_max { z.z_min } else { z.z_max };
    Ok(METHANE_M * p_pa / (z_stable * NEQSIM_GAS_CONSTANT * t))
}

/// The methane reference viscosity, `PFCTViscosityMethodHeavyOil.getRefComponentViscosity`.
///
/// `t0` and `p0` are the alpha-corrected reference temperature (K) and pressure (Pa).
fn get_ref_viscosity(t0: f64, p0: f64) -> Result<f64> {
    let rho = methane_srk_density(t0, p0)?;
    let mol_dens_molar = rho / METHANE_M * 1e-3; // mol/dm^3
    let red_mol_dens = (mol_dens_molar - CRIT_MOL_DENS) / CRIT_MOL_DENS; // subtraction here
    let mol_dens = rho * 1e-3; // g/cm^3

    let visc_ref_o = GVCOEF[0] * t0.powf(-1.0)
        + GVCOEF[1] * t0.powf(-2.0 / 3.0)
        + GVCOEF[2] * t0.powf(-1.0 / 3.0)
        + GVCOEF[3]
        + GVCOEF[4] * t0.powf(1.0 / 3.0)
        + GVCOEF[5] * t0.powf(2.0 / 3.0)
        + GVCOEF[6] * t0
        + GVCOEF[7] * t0.powf(4.0 / 3.0)
        + GVCOEF[8] * t0.powf(5.0 / 3.0);
    let visc_ref_1 =
        (VIS_REF_A + VIS_REF_B * (VIS_REF_C - (t0 / VIS_REF_F).ln()).powi(2)) * mol_dens;

    let temp_1 = mol_dens.powf(0.1) * (VISC_REF_J[1] + VISC_REF_J[2] / t0.powf(1.5));
    let temp_2 = red_mol_dens
        * mol_dens.powf(0.5)
        * (VISC_REF_J[4] + VISC_REF_J[5] / t0 + VISC_REF_J[6] / t0.powi(2));
    let temp_3 = (temp_1 + temp_2).exp();

    let htan = (t0 - 90.69).tanh();
    let vis_ref_e = (htan + 1.0) / 2.0;
    let visc_ref_2 = vis_ref_e * (VISC_REF_J[0] + VISC_REF_J[3] / t0).exp() * (temp_3 - 1.0);
    let visc_ref_2 = if visc_ref_2.is_nan() { 0.0 } else { visc_ref_2 };

    let temp_4 = mol_dens.powf(0.1) * (VISC_REF_K[1] + VISC_REF_K[2] / t0.powf(1.5));
    let temp_5 = red_mol_dens
        * mol_dens.powf(0.5)
        * (VISC_REF_K[4] + VISC_REF_K[5] / t0 + VISC_REF_K[6] / t0.powi(2));
    let temp_6 = (temp_4 + temp_5).exp();
    let vis_ref_g = (1.0 - htan) / 2.0;
    let visc_ref_3 = vis_ref_g * (VISC_REF_K[0] + VISC_REF_K[3] / t0).exp() * (temp_6 - 1.0);
    let visc_ref_3 = if visc_ref_3.is_nan() { 0.0 } else { visc_ref_3 };

    Ok((visc_ref_o + visc_ref_1 + visc_ref_2 + visc_ref_3) / 1.0e7)
}

/// The corresponding-states scaling, `calculateCorrespondingStatesViscosity`, with the
/// four correction factors at their default `1.0`.
#[allow(clippy::too_many_arguments)]
fn csp(
    reference_viscosity: f64,
    tc_mix: f64,
    pc_mix: f64,
    m_mix: f64,
    alfa_mix: f64,
    alfa0: f64,
) -> f64 {
    reference_viscosity
        * (tc_mix / METHANE_TC).powf(-1.0 / 6.0)
        * (pc_mix / METHANE_PC).powf(2.0 / 3.0)
        * (m_mix / (METHANE_M * 1e3)).powf(0.5)
        * (alfa_mix / alfa0)
}

/// The liquid viscosity of a mixture, from the Pedersen (PFCT) heavy-oil correlation.
///
/// `mixture` carries the components' `Tc`, `Pc` and molar mass; `z` is the mole
/// fraction, checked rather than renormalised. The reference flash is pure methane,
/// so the mixture's cubic and alpha are not consulted.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `z` is not a composition of the mixture's length.
/// * [`AzothError::PropertyUnavailable`] if a component carries no molar mass - a
///   card-added substance the databank could not complete.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
pub fn viscosity(
    mixture: &Mixture,
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<ViscosityResult> {
    let spec = &model_gen::VISCOSITY_SPEC;
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

    // The Tc/Pc/molar-mass mixing rule, in g/mol for the mass (NeqSim's getMolarMass
    // is kg/mol, so the sums are in kg/mol and the `* 1e3` recovers g/mol).
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
    if temp_tc2 < 1e-10 {
        return Ok(ViscosityResult {
            mu: pascal_seconds(0.0),
            warnings,
        });
    }
    let pc_mix = 8.0 * temp_pc1 / (temp_pc2 * temp_pc2);
    let tc_mix = temp_tc1 / temp_tc2;
    let m_mix =
        (mm_temp + 1.304e-4 * ((mw_temp / mm_temp).powf(2.303) - mm_temp.powf(2.303))) * 1e3;

    // Reference flash at the scaled state, then the alpha correction.
    let t_scaled = t.value * METHANE_TC / tc_mix;
    let p_scaled = p.value * METHANE_PC / pc_mix;
    let rho = methane_srk_density(t_scaled, p_scaled)?;
    let mol_dens_molar = rho / METHANE_M * 1e-3;
    let red_dens = mol_dens_molar / CRIT_MOL_DENS; // division, not the subtraction below
    let alfa_mix = 1.0 + 7.378e-3 * red_dens.powf(1.847) * m_mix.powf(0.5173);
    let alfa0 = 1.0 + 7.378e-3 * red_dens.powf(1.847) * (METHANE_M * 1e3).powf(0.5173);
    let t0 = t_scaled * alfa0 / alfa_mix;
    let p0 = p_scaled * alfa0 / alfa_mix;

    let mu = if t0 < 75.0 {
        // The heavy-oil branch, above the freeze switch of the reference correlation.
        let mol_m = if mw_temp / mm_temp / mm_temp <= 1.5 {
            mm_temp * 1e3
        } else {
            mm_temp * (mw_temp / mm_temp / (1.5 * mm_temp)).powf(0.5) * 1e3
        };
        let sign = if t.value > 564.49 { -1.0 } else { 1.0 };
        let termm = -0.079_55 - sign * 0.011_01 * mol_m - 371.8 / t.value + 6.215 * mol_m / t.value;
        let ho_viscosity = 10.0f64.powf(termm) * 1e-3;
        let ho_viscosity = ho_viscosity + ho_viscosity * 0.008 * (p.value / 1e5 - 1.0);
        if t0 < 65.0 {
            ho_viscosity
        } else {
            let ref_viscosity = get_ref_viscosity(t0, p0)?;
            let lo_viscosity = csp(ref_viscosity, tc_mix, pc_mix, m_mix, alfa_mix, alfa0);
            lo_viscosity * (1.0 - (75.0 - t0) / 10.0) + ho_viscosity * (75.0 - t0) / 10.0
        }
    } else {
        let ref_viscosity = get_ref_viscosity(t0, p0)?;
        csp(ref_viscosity, tc_mix, pc_mix, m_mix, alfa_mix, alfa0)
    };

    apply_checks(
        spec.derived_checks(),
        |name| (name == "mu").then_some(mu),
        &mut warnings,
    )?;

    Ok(ViscosityResult {
        mu: pascal_seconds(mu),
        warnings,
    })
}
