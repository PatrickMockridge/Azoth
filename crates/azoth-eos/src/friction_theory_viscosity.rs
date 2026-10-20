//! `eos.friction_theory_viscosity` - the friction-theory viscosity of Quinones-Cisneros
//! and Firoozabadi, over a Chung baseline and Wilke mixing.
//!
//! Spec: `specs/models/eos/friction_theory_viscosity.toml`, which records
//! NeqSim's `FrictionTheoryViscosityMethod.calcViscosity` and its helpers.

use azoth_core::units::{
    DynamicViscosity, MolarMass, MolarVolume, Pressure, ThermodynamicTemperature, pascal_seconds,
};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::FrictionTheoryViscosityResult;

/// NeqSim's gas constant, `ThermodynamicConstantsInterface.R`.
const R: f64 = 8.3144621;
/// The floor `MIN_VISCOSITY_PA_S` puts under the answer.
const MIN_VISCOSITY_PA_S: f64 = 1.0e-6;
/// Micropoise to Pa*s.
const MICRO_POISE_TO_PA_S: f64 = 1.0e-7;

/// The reduced friction constants, one set per cubic family.
struct FrictionConstants {
    kapac: f64,
    kaprc: f64,
    kaprrc: f64,
    kapa: [[f64; 3]; 3],
    kapr: [[f64; 3]; 3],
    kaprr: f64,
}

/// `setEosDependentConstants`: PR for `PhasePrEos`, the SRK set for every other phase.
fn constants(peng_robinson: bool) -> FrictionConstants {
    if peng_robinson {
        FrictionConstants {
            kapac: -0.140464,
            kaprc: 1.19902e-2,
            kaprrc: 8.55115e-4,
            kapa: [
                [-4.89197e-2, 0.270572, -4.48111e-2],
                [0.270572, -1.10473e-4, 4.08972e-5],
                [-4.48111e-2, 4.08972e-5, -5.79765e-9],
            ],
            kapr: [
                [-0.357875, 0.637572, -7.9024e-2],
                [0.637572, -6.02128e-5, 3.72408e-5],
                [-7.9024e-2, 3.72408e-5, -5.65610e-9],
            ],
            kaprr: 1.37290e-8,
        }
    } else {
        FrictionConstants {
            kapac: -0.165302,
            kaprc: 6.99574e-3,
            kaprrc: 1.26358e-3,
            kapa: [
                [-0.114804, 0.246622, -3.94638e-2],
                [0.246622, -1.15648e-4, 4.18863e-5],
                [-3.94638e-2, 4.18863e-5, -5.91999e-9],
            ],
            kapr: [
                [-0.315903, 0.566713, -7.29995e-2],
                [0.566713, -1.0086e-4, 5.17459e-5],
                [-7.29995e-2, 5.17459e-5, -5.68708e-9],
            ],
            kaprr: 1.35994e-8,
        }
    }
}

fn reduced_kapr(phi: f64, gamma: f64, c: &FrictionConstants) -> f64 {
    c.kaprc
        + c.kapr[0][0] * (gamma - 1.0)
        + (c.kapr[1][0] + c.kapr[1][1] * phi) * (gamma - 1.0).exp_m1()
        + (c.kapr[2][0] + c.kapr[2][1] * phi + c.kapr[2][2] * phi * phi)
            * (2.0 * gamma - 2.0).exp_m1()
}

fn reduced_kaprr(phi: f64, gamma: f64, c: &FrictionConstants) -> f64 {
    c.kaprrc + c.kaprr * phi * (2.0 * gamma).exp_m1() * (gamma - 1.0).powi(2)
}

fn reduced_kapa(phi: f64, gamma: f64, c: &FrictionConstants) -> f64 {
    c.kapac
        + c.kapa[0][0] * (gamma - 1.0)
        + (c.kapa[1][0] + c.kapa[1][1] * phi) * (gamma - 1.0).exp_m1()
        + (c.kapa[2][0] + c.kapa[2][1] * phi + c.kapa[2][2] * phi * phi)
            * (2.0 * gamma - 2.0).exp_m1()
}

/// The friction-theory viscosity of a mixture, from its components' constants and the
/// EOS's two pressure contributions.
///
/// The pure-component baseline is the class's own Chung form, the dilute mixture is
/// Wilke's rule, and the excess friction is
/// `kappa_r P_rep + kappa_a P_att + kappa_rr P_rep^2`, all in micropoise before the final
/// conversion to Pa*s. `p_repulsive` and `p_attractive` are the caller's EOS values; this
/// id does not build one.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] if the vectors differ in length.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     cubic_meters_per_mole, kilograms_per_mole, kelvins, pascal_seconds, pascals,
/// };
/// use azoth_eos::friction_theory_viscosity::friction_theory_viscosity;
///
/// // The SRK oil phase of methane/n-butane at 300 K and 20 bar.
/// let r = friction_theory_viscosity(
///     &[kelvins(190.56), kelvins(425.12)],
///     &[pascals(4.599e6), pascals(3.796e6)],
///     &[cubic_meters_per_mole(9.9e-5), cubic_meters_per_mole(2.55e-4)],
///     &[0.0115, 0.2002],
///     &[kilograms_per_mole(0.016043), kilograms_per_mole(0.058123)],
///     &[pascal_seconds(1.5293e-5), pascal_seconds(2.5768e-5)],
///     &[0.0, 0.0],
///     &[0.0, 0.0],
///     &[0.09542082642782022, 0.9045791735721797],
///     kelvins(300.0),
///     pascals(829.1753495002592),
///     pascals(-809.1753495002592),
///     false,
/// )?;
/// assert!((r.mu.value - 1.383523806147524e-4).abs() < 1e-16);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` is the symbol in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn friction_theory_viscosity(
    tc: &[ThermodynamicTemperature],
    pc: &[Pressure],
    critical_volume: &[MolarVolume],
    omega: &[f64],
    molar_mass: &[MolarMass],
    critical_viscosity: &[DynamicViscosity],
    viscosity_correction_factor: &[f64],
    dipole: &[f64],
    x: &[f64],
    T: ThermodynamicTemperature,
    p_repulsive: Pressure,
    p_attractive: Pressure,
    peng_robinson: bool,
) -> Result<FrictionTheoryViscosityResult> {
    let spec = &model_gen::FRICTION_THEORY_VISCOSITY_SPEC;
    let mut warnings = Vec::new();

    let count = x.len();
    if tc.len() != count
        || pc.len() != count
        || critical_volume.len() != count
        || omega.len() != count
        || molar_mass.len() != count
        || critical_viscosity.len() != count
        || viscosity_correction_factor.len() != count
        || dipole.len() != count
    {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "{} mole fraction(s) against {} critical temperature(s), {} critical \
                 pressure(s), {} critical volume(s), {} acentric factor(s), {} molar mass(es), \
                 {} critical viscosity(ies), {} correction factor(s) and {} dipole moment(s): \
                 every one is one entry per component",
                count,
                tc.len(),
                pc.len(),
                critical_volume.len(),
                omega.len(),
                molar_mass.len(),
                critical_viscosity.len(),
                viscosity_correction_factor.len(),
                dipole.len()
            ),
        ));
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "tc" => tc.iter().map(|v| v.value).reduce(f64::min),
            "pc" => pc.iter().map(|v| v.value).reduce(f64::min),
            "critical_volume" => critical_volume.iter().map(|v| v.value).reduce(f64::min),
            "molar_mass" => molar_mass.iter().map(|v| v.value).reduce(f64::min),
            "critical_viscosity" => critical_viscosity.iter().map(|v| v.value).reduce(f64::min),
            "x" => x.iter().copied().reduce(f64::min),
            "p_repulsive" => Some(p_repulsive.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The Chung pure-component baseline, in micropoise.
    let mut pure = Vec::with_capacity(count);
    for i in 0..count {
        let vc_cm3 = critical_volume[i].value * 1.0e6;
        let relative_dipole = 131.3 * dipole[i] / (vc_cm3 * tc[i].value).sqrt();
        let fc = 1.0 - 0.2756 * omega[i]
            + 0.059035 * relative_dipole.powi(4)
            + viscosity_correction_factor[i];
        let t_star = 1.2593 * T.value / tc[i].value;
        let var_last =
            -6.435e-4 * t_star.powf(0.14874) * (18.0323 * t_star.powf(-0.76830) - 7.27371).sin();
        let omega_visc = 1.16145 / t_star.powf(0.14874)
            + 0.52487 / (0.77320 * t_star).exp()
            + 2.16178 / (2.43787 * t_star).exp()
            + var_last;
        pure.push(
            40.785 * (molar_mass[i].value * 1000.0 * T.value).sqrt()
                / (vc_cm3.powf(2.0 / 3.0) * omega_visc)
                * fc,
        );
    }

    // Wilke's mixing rule, in micropoise.
    let molar_mass_pow: Vec<f64> = molar_mass.iter().map(|mass| mass.value.powf(0.3)).collect();
    let mass_sum: f64 = x
        .iter()
        .zip(&molar_mass_pow)
        .map(|(fraction, pow)| fraction / pow)
        .sum();
    let mut visk0 = 0.0;
    for i in 0..count {
        let mut denominator = 0.0;
        for j in 0..count {
            let phi = (1.0
                + (pure[i] / pure[j]).sqrt()
                    * (molar_mass[j].value / molar_mass[i].value).powf(0.25))
            .powi(2)
                / (8.0 * (1.0 + molar_mass[i].value / molar_mass[j].value)).sqrt();
            denominator += x[j] * phi;
        }
        visk0 += x[i] * pure[i] / denominator;
    }

    // The friction term, in micropoise.
    let c = constants(peng_robinson);
    let mut kaprmx = 0.0;
    let mut kapamx = 0.0;
    let mut kaprrmx = 0.0;
    for i in 0..count {
        let nci = critical_viscosity[i].value * 1.0e7;
        let pc_bar = pc[i].value / 1.0e5;
        let phi = 1.0e1 * R * tc[i].value / pc_bar;
        let gamma = tc[i].value / T.value;
        let zi = x[i] / (molar_mass_pow[i] * mass_sum);
        kaprmx += zi * nci * reduced_kapr(phi, gamma, &c) / pc_bar;
        kapamx += zi * nci * reduced_kapa(phi, gamma, &c) / pc_bar;
        kaprrmx += zi * nci * reduced_kaprr(phi, gamma, &c) / (pc_bar * pc_bar);
    }
    let visk1 = kaprmx * p_repulsive.value
        + kapamx * p_attractive.value
        + kaprrmx * p_repulsive.value.powi(2);

    let mu = ((visk0 + visk1) * MICRO_POISE_TO_PA_S).max(MIN_VISCOSITY_PA_S);

    apply_checks(
        spec.derived_checks(),
        |name| (name == "mu").then_some(mu),
        &mut warnings,
    )?;

    Ok(FrictionTheoryViscosityResult {
        mu: pascal_seconds(mu),
        pure_component_viscosity: pure
            .iter()
            .map(|value| value * MICRO_POISE_TO_PA_S)
            .collect(),
        dilute_viscosity: pascal_seconds(visk0 * MICRO_POISE_TO_PA_S),
        friction_viscosity: pascal_seconds(visk1 * MICRO_POISE_TO_PA_S),
        warnings,
    })
}
