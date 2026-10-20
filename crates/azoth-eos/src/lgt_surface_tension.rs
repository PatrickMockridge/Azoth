//! `eos.lgt_surface_tension` - the linear-gradient-theory surface tension.
//!
//! Spec: `specs/models/eos/lgt_surface_tension.toml`, which records NeqSim's
//! `LGTSurfaceTension.calcSurfaceTension`.

use azoth_core::units::{Pressure, ThermodynamicTemperature, newtons_per_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::Mixture;
use crate::model_gen;
use crate::pt_flash::pt_flash;
use crate::results::{LgtSurfaceTensionResult, Phase};

/// NeqSim's gas constant, `ThermodynamicConstantsInterface.R`.
const R: f64 = 8.314_462_1;
/// The class's fixed number of steps along the reference component's density.
const STEPS: usize = 200;

/// The linear-gradient-theory interface tension of a flashed fluid.
///
/// The density path between the two phases is a straight line in molar density over
/// [`STEPS`] steps; at each step the T,V chemical potential
/// `mu_i/(RT) = ln(rho_i R T) + dA_R/dn_i` is evaluated, and the tension is the sum of
/// `sqrt(2 kappa mu_times_rho)` over the path. **The path crosses the cubic's own loop**,
/// where the pressure is negative and a pressure-based potential does not exist; the T,V
/// form does.
///
/// `influence_parameters` and `pair_interactions` are the class's own per-component data,
/// which it derives from `getSurfaceTenisionInfluenceParameter` and the component names;
/// they are inputs here so the model does not need the names.
///
/// # Errors
/// * [`AzothError::InvalidInput`] for a single-phase feed, a zero phase density, an
///   absent component, or a step whose integrand is not positive.
/// * [`AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::Cubic;
/// use azoth_eos::databank;
/// use azoth_eos::lgt_surface_tension::lgt_surface_tension;
///
/// let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], Cubic::Srk, None)?;
/// let r = lgt_surface_tension(
///     &mixture,
///     &[9.634812734820117e-21, 2.7571945977326497e-19],
///     1.0,
///     &[vec![0.0, 0.0], vec![0.0, 0.0]],
///     kelvins(300.0),
///     pascals(2.0e6),
///     &[0.5, 0.5],
/// )?;
/// assert!((r.sigma.value - 0.010011889708138906).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn lgt_surface_tension(
    mixture: &Mixture,
    influence_parameters: &[f64],
    reference_component: f64,
    pair_interactions: &[Vec<f64>],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<LgtSurfaceTensionResult> {
    let spec = &model_gen::LGT_SURFACE_TENSION_SPEC;
    let mut warnings = Vec::new();

    let count = mixture.len();
    if influence_parameters.len() != count
        || pair_interactions.len() != count
        || pair_interactions.iter().any(|row| row.len() != count)
        || z.len() != count
    {
        return Err(AzothError::invalid_input(
            "influence_parameters",
            format!(
                "{count} component(s) against {} influence parameter(s), {} interaction \
                 row(s) and {} mole fraction(s): every one is one entry per component",
                influence_parameters.len(),
                pair_interactions.len(),
                z.len()
            ),
        ));
    }
    if !reference_component.is_finite()
        || reference_component < 0.0
        || reference_component.fract() != 0.0
        || reference_component as usize >= count
    {
        return Err(AzothError::invalid_input(
            "reference_component",
            format!(
                "`{reference_component}` is not a whole index in 0..{count}: it is the \
                 component the profile integrates along"
            ),
        ));
    }
    let reference = reference_component as usize;

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t.value),
            "P" => Some(p.value),
            "influence_parameters" => influence_parameters.iter().copied().reduce(f64::min),
            "z" => z.iter().copied().reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    let flash = pt_flash(mixture, t, p, z)?;
    if flash.phase != Phase::TwoPhase {
        return Err(AzothError::invalid_input(
            "z",
            "the feed is single phase at this state, so there is no interface to have a \
             tension; a gradient-theory profile runs between two phases",
        ));
    }

    let rt = R * t.value;
    let vapour = &flash.y;
    let liquid = &flash.x;
    if flash.z_vapour <= 0.0 || flash.z_liquid <= 0.0 {
        return Err(AzothError::invalid_input(
            "z",
            "a phase's compressibility root is not positive, so it has no density",
        ));
    }
    // The class's stored densities are molar densities, mol/m³.
    let rho_vapour: Vec<f64> = vapour
        .iter()
        .map(|fraction| fraction * p.value / (flash.z_vapour * rt))
        .collect();
    let rho_liquid: Vec<f64> = liquid
        .iter()
        .map(|fraction| fraction * p.value / (flash.z_liquid * rt))
        .collect();
    for fraction in vapour.iter().chain(liquid) {
        if *fraction <= 0.0 {
            return Err(AzothError::invalid_input(
                "z",
                "a component is absent from one of the phases, so its chemical potential \
                 would be `ln(0)`; the class stops at the resulting NaN and this port \
                 refuses the partial integral",
            ));
        }
    }

    let reduced = mixture.reduced_parameters(t, p)?;
    let gradient_vapour = mixture.helmholtz_gradient(&reduced, vapour, flash.z_vapour)?;
    let mu_equi: Vec<f64> = (0..count)
        .map(|i| rt * ((rho_vapour[i] * rt).ln() + gradient_vapour[i]))
        .collect();

    // The absolute attraction and covolume at the flashed state, from the reduced pair.
    let a_abs: Vec<f64> = reduced
        .a
        .iter()
        .map(|value| value * rt * rt / p.value)
        .collect();
    let b_abs: Vec<f64> = reduced.b.iter().map(|value| value * rt / p.value).collect();
    let cubic = mixture.cubic();
    let (delta1, delta2) = (cubic.delta1(), cubic.delta2());

    let delta_rho: Vec<f64> = (0..count)
        .map(|i| (rho_liquid[i] - rho_vapour[i]) / STEPS as f64)
        .collect();
    if delta_rho[reference] == 0.0 {
        return Err(AzothError::invalid_input(
            "reference_component",
            "the reference component's density is the same in both phases, so the profile \
             has no integration variable",
        ));
    }

    let mut sigma = 0.0;
    for step in 1..STEPS {
        let fraction = step as f64 / STEPS as f64;
        let rho: Vec<f64> = (0..count)
            .map(|i| rho_vapour[i] + (rho_liquid[i] - rho_vapour[i]) * fraction)
            .collect();
        let total: f64 = rho.iter().sum();
        let x: Vec<f64> = rho.iter().map(|value| value / total).collect();
        let volume = 1.0 / total;

        let mut b_mix = 0.0;
        let mut a_mix = 0.0;
        for i in 0..count {
            b_mix += x[i] * b_abs[i];
            for j in 0..count {
                let kij = reduced.kij[i * count + j];
                a_mix += x[i] * x[j] * (1.0 - kij) * (a_abs[i] * a_abs[j]).sqrt();
            }
        }
        let p_step =
            rt / (volume - b_mix) - a_mix / ((volume + delta1 * b_mix) * (volume + delta2 * b_mix));
        let compressibility = p_step * volume / rt;
        // **The state pressure can be negative, and `reduced_parameters` refuses that.**
        // The gradient needs only the reduced `A_i` and `B_i` at *this* state, and those
        // are the flashed state's scaled by `P_step/P`: `A_i = a_i P/(RT)^2` is linear in
        // the pressure. Scaling the pair leaves `kij` and everything the gradient does not
        // read untouched, and keeps the primitive available exactly where the class's own
        // T,V potential is.
        let ratio = p_step / p.value;
        let mut reduced_step = reduced.clone();
        for value in &mut reduced_step.a {
            *value *= ratio;
        }
        for value in &mut reduced_step.b {
            *value *= ratio;
        }
        let gradient = mixture.helmholtz_gradient(&reduced_step, &x, compressibility)?;
        let mu: Vec<f64> = (0..count)
            .map(|i| rt * ((rho[i] * rt).ln() + gradient[i]))
            .collect();

        let mut kappa = 0.0;
        for i in 0..count {
            let ki = if i == reference {
                1.0
            } else {
                delta_rho[i] / delta_rho[reference]
            };
            for k in 0..count {
                let kk = if k == reference {
                    1.0
                } else {
                    delta_rho[k] / delta_rho[reference]
                };
                kappa += (influence_parameters[i] * influence_parameters[k]).sqrt()
                    * ki
                    * kk
                    * (1.0 - pair_interactions[i][k]);
            }
        }
        let mu_times_rho: f64 = (0..count)
            .map(|i| rho[i] * (mu[i] - mu_equi[i]))
            .sum::<f64>()
            - (p_step - p.value);
        if kappa <= 0.0 || mu_times_rho <= 0.0 {
            return Err(AzothError::invalid_input(
                "z",
                format!(
                    "the gradient-theory integrand is not positive at step {step} \
                     (kappa = {kappa}, mu_times_rho = {mu_times_rho}); the class takes its \
                     square root and stops at the NaN, and this port refuses a partial \
                     integral"
                ),
            ));
        }
        sigma += (2.0 * kappa * mu_times_rho).sqrt() * delta_rho[reference];
    }

    Ok(LgtSurfaceTensionResult {
        sigma: newtons_per_meter(sigma.abs()),
        warnings,
    })
}
