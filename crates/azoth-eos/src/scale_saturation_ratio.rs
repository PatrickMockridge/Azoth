//! `eos.scale_saturation_ratio` - one salt's saturation ratio in a brine.
//!
//! ```text
//! m_i = x_i / (x_water M_water)                    [mol/kg water]
//! Ksp = exp(A/T + B + C ln T + D T + E/T**2) exp(-dV (P - P0)/(R T))
//! SR  = (gamma1 m1)**stoc1 (gamma2 m2)**stoc2 a_water**waterstoc / Ksp
//! ```
//!
//! Spec: `specs/calcs/eos/scale_saturation_ratio.toml`, which carries the three name
//! overrides, the clamps and the measurement that pins both worked cases.
//!
//! NeqSim's `CheckScalePotential`, one row of its `compsalt` walk. The operation reports one
//! ratio per salt whose two ions the phase has; this is that row, so porting the operation is
//! looping over this rather than restating it.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank;
use crate::results::ScaleSaturationRatioResult;
use crate::spec_gen;

/// The gas constant in the units `Vdelta` and the pressure term are stated in: cm³ bar per
/// (mol K).
///
/// **Not `8.3144621`.** `Vdelta` is cm³/mol and the pressure in NeqSim's own code is bar, so
/// the term is dimensionless as written; the SI constant against a cm³ volume is `1e4` out.
const R_CM3_BAR: f64 = 83.1446;

/// The pressure the volume term is referred to, in bar - one atmosphere.
const REFERENCE_PRESSURE_BAR: f64 = 1.01325;

/// Below this molality the ion activity product is taken as zero.
const MIN_MOLALITY: f64 = 1.0e-30;

/// The clamp on `ln SR`, NeqSim's own: `exp(69)` is about `1e30`.
const LN_SR_CLAMP: f64 = 69.0;

/// One salt's saturation ratio in an aqueous phase.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the salt is not in the table, if `FeS` is asked for
///   without a hydrogen-ion molality, or if water carries no molar mass.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive, or if the salt's solubility
///   product has no value.
#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn scale_saturation_ratio(
    salt: &str,
    x1: f64,
    x2: f64,
    x_water: f64,
    gamma1: f64,
    gamma2: f64,
    water_activity: f64,
    h3o_molality: Option<f64>,
    T: ThermodynamicTemperature,
    P: Pressure,
) -> Result<ScaleSaturationRatioResult> {
    let spec = &spec_gen::SCALE_SATURATION_RATIO_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let record = databank::salt(salt).ok_or_else(|| {
        AzothError::invalid_input(
            "salt",
            format!("`{salt}` is not a row of `compsalt`, so it has no solubility product"),
        )
    })?;
    let t = T.value;
    let p_bar = P.value / 1.0e5;

    // The row's own correlation, which three names then override.
    let mut ksp = (record.ksp[0] / t
        + record.ksp[1]
        + record.ksp[2] * t.ln()
        + record.ksp[3] * t
        + record.ksp[4] / (t * t))
        .exp();
    match salt {
        "NaCl" => ksp = 92.78 - 0.407 * t + 0.000747 * t * t,
        "CaCO3" => {
            // Plummer & Busenberg (1982), for calcite.
            let log10_ksp = -171.9065 - 0.077993 * t + 2839.319 / t + 71.595 * t.log10();
            ksp = 10.0_f64.powf(log10_ksp);
        }
        "FeCO3" => {
            // Greenberg & Tomson (1992).
            let log10_ksp = -59.3498 - 0.041377 * t - 2.1963 / t + 24.5724 * t.log10();
            ksp = 10.0_f64.powf(log10_ksp);
        }
        "FeS" => {
            // **The one salt whose product carries an ion molality**, and NeqSim skips the
            // salt outright where the phase has no `H3O+` - which a scalar calc cannot do, so
            // the caller says which it is by stating the molality or by not calling.
            let hydrogen = h3o_molality.ok_or_else(|| {
                AzothError::invalid_input(
                    "h3o_molality",
                    "`FeS`'s solubility product is multiplied by the hydrogen-ion molality, \
                     and NeqSim skips the salt altogether where the phase has no `H3O+`; this \
                     calc has no phase to look in, so the caller states it"
                        .to_string(),
                )
            })?;
            ksp *= hydrogen;
        }
        _ => {}
    }

    // The pressure term, guarded exactly as NeqSim guards it.
    if record.volume_delta.abs() > 1.0e-10 && p_bar > 1.013 {
        let correction = -record.volume_delta * (p_bar - REFERENCE_PRESSURE_BAR) / (R_CM3_BAR * t);
        ksp *= correction.clamp(-LN_SR_CLAMP, LN_SR_CLAMP).exp();
    }
    if !ksp.is_finite() || ksp.partial_cmp(&0.0) != Some(std::cmp::Ordering::Greater) {
        return Err(AzothError::out_of_range(
            "Ksp",
            ksp,
            format!("the solubility product of {salt} has no value at {t} K and {p_bar} bar"),
        ));
    }

    // The molalities, per kilogram of water: NeqSim reads the phase's own water component's
    // molar mass, which is the table's.
    let molar_mass_water = databank::entry("water", None)?.molar_mass.ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "water carries no molar mass, and a molality is per kilogram of it".to_string(),
        )
    })?;
    let m1 = x1 / (x_water * molar_mass_water);
    let m2 = x2 / (x_water * molar_mass_water);

    // A positive finite test written out rather than as `!(x > 0.0)`, which clippy reads as
    // the double negative it is: a `NaN` is neither, and it is not a coefficient.
    let positive = |value: f64| value.is_finite() && value > 0.0;
    // **A ratio too far below saturation to matter is a zero**, and the guard is before the
    // logarithm rather than after it.
    let (saturation_ratio, ion_activity_product) =
        if m1 < MIN_MOLALITY || m2 < MIN_MOLALITY || !positive(gamma1) || !positive(gamma2) {
            (0.0, 0.0)
        } else {
            let mut ln_iap = record.cation_stoichiometry * (gamma1 * m1).ln()
                + record.anion_stoichiometry * (gamma2 * m2).ln();
            let mut iap = (gamma1 * m1).powf(record.cation_stoichiometry)
                * (gamma2 * m2).powf(record.anion_stoichiometry);
            if record.water_stoichiometry > 0.0 {
                if !positive(water_activity) {
                    return Ok(ScaleSaturationRatioResult {
                        saturation_ratio: 0.0,
                        ion_activity_product: 0.0,
                        solubility_product: ksp,
                        warnings,
                    });
                }
                ln_iap += record.water_stoichiometry * water_activity.ln();
                iap *= water_activity.powf(record.water_stoichiometry);
            }
            let ln_sr = ln_iap - ksp.ln();
            let ratio = if ln_sr < -LN_SR_CLAMP {
                0.0
            } else if ln_sr > LN_SR_CLAMP {
                LN_SR_CLAMP.exp()
            } else {
                ln_sr.exp()
            };
            (ratio, iap)
        };

    Ok(ScaleSaturationRatioResult {
        saturation_ratio,
        ion_activity_product,
        solubility_product: ksp,
        warnings,
    })
}
