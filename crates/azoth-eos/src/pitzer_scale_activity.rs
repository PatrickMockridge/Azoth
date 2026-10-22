//! `eos.pitzer_scale_activity` - NeqSim's `PitzerScaleActivityModel`.
//!
//! A binary Pitzer model over the PHREEQC `pitzer.dat` parameter sets: mean activity
//! coefficients from the six-term temperature expression, and trace ions mapped at an
//! equivalent ionic strength by the MacInnes convention.

use azoth_core::units::ThermodynamicTemperature;
use azoth_core::{AzothError, Result, apply_checks};

use crate::pitzer_catalog::{Family, find};
use crate::results::PitzerScaleActivityResult;
use crate::spec_gen;

/// Pitzer's universal `b`, the class's `PITZER_B`.
const PITZER_B: f64 = 1.2;
/// The floor on an ionic strength under a square root.
const MIN_IONIC_STRENGTH: f64 = 1.0e-12;
/// The temperature the six-term expressions are written about.
const REFERENCE_TEMPERATURE_K: f64 = 298.15;

/// The nine ions the model carries, in the class's order.
const IONS: [&str; 9] = [
    "sodium",
    "chloride",
    "calcium",
    "barium",
    "strontium",
    "ferrous",
    "sulphate",
    "carbonate",
    "bicarbonate",
];

/// One binary salt: its two ions and Pitzer's stoichiometric constants.
struct SaltSpec {
    cation: &'static str,
    anion: &'static str,
    cation_stoichiometry: f64,
    anion_stoichiometry: f64,
    cation_charge: f64,
    anion_charge: f64,
    alpha1: f64,
    alpha2: f64,
}

fn salt_spec(name: &str) -> Option<SaltSpec> {
    let (cation, anion, ns, na, zs, za, alpha1, alpha2) = match name {
        "nacl" => ("Na+", "Cl-", 1.0, 1.0, 1.0, 1.0, 2.0, 0.0),
        "cacl2" => ("Ca++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 12.0),
        "bacl2" => ("Ba++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 0.0),
        "srcl2" => ("Sr++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 0.0),
        "fecl2" => ("Fe++", "Cl-", 1.0, 2.0, 2.0, 1.0, 2.0, 0.0),
        "na2so4" => ("Na+", "SO4--", 2.0, 1.0, 1.0, 2.0, 2.0, 0.0),
        "na2co3" => ("Na+", "CO3--", 2.0, 1.0, 1.0, 2.0, 2.0, 0.0),
        "nahco3" => ("Na+", "HCO3-", 1.0, 1.0, 1.0, 1.0, 2.0, 12.0),
        _ => return None,
    };
    Some(SaltSpec {
        cation,
        anion,
        cation_stoichiometry: ns,
        anion_stoichiometry: na,
        cation_charge: zs,
        anion_charge: za,
        alpha1,
        alpha2,
    })
}

/// The binary Pitzer activity coefficients of a scale-forming brine.
///
/// `ionic_strength_molal` configures the trace-ion mapping; `salt_molality` is the mean
/// coefficient's own molality, from which the class computes its own ionic strength.
///
/// # Errors
/// * [`AzothError::InvalidInput`] for a salt or ion name outside the model's set.
/// * [`AzothError::OutOfRange`] if a temperature or strength is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::kelvins;
/// use azoth_eos::pitzer_scale_activity::pitzer_scale_activity;
///
/// let r = pitzer_scale_activity(
///     kelvins(298.15),
///     1.0,
///     0.5107140231833494,
///     "nacl",
///     1.0,
///     "sulphate",
///     "barium",
///     "sulphate",
/// )?;
/// assert!((r.mean_activity_coefficient - 0.6566052005348894).abs() < 1e-12);
/// assert!((r.activity_coefficient_product - 0.01298921509374087).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn pitzer_scale_activity(
    temperature: ThermodynamicTemperature,
    ionic_strength_molal: f64,
    debye_huckel_a: f64,
    salt: &str,
    salt_molality: f64,
    ion: &str,
    cation: &str,
    anion: &str,
) -> Result<PitzerScaleActivityResult> {
    let spec = &spec_gen::PITZER_SCALE_ACTIVITY_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "ionic_strength_molal" => Some(ionic_strength_molal),
            "debye_huckel_a" => Some(debye_huckel_a),
            "salt_molality" => Some(salt_molality),
            _ => None,
        },
        &mut warnings,
    )?;

    let selected = salt_spec(salt).ok_or_else(|| {
        AzothError::invalid_input(
            "salt",
            format!("`{salt}` is not one of the model's eight binary salts"),
        )
    })?;
    for (name, role) in [(ion, "ion"), (cation, "cation"), (anion, "anion")] {
        if !IONS.contains(&name) {
            return Err(AzothError::invalid_input(
                role,
                format!("`{name}` is not one of the model's nine ions: {IONS:?}"),
            ));
        }
    }

    let t = temperature.value;
    let mean = mean_activity_coefficient(&selected, salt_molality, t, debye_huckel_a);
    let gammas = ion_coefficients(t, debye_huckel_a, ionic_strength_molal);
    let gamma_of = |name: &str| gammas[IONS.iter().position(|value| *value == name).unwrap()];
    let ion_value = gamma_of(ion);
    let product = gamma_of(cation) * gamma_of(anion);

    Ok(PitzerScaleActivityResult {
        mean_activity_coefficient: mean,
        ion_activity_coefficient: ion_value,
        activity_coefficient_product: product,
        gamma_sodium: gammas[0],
        gamma_chloride: gammas[1],
        gamma_calcium: gammas[2],
        gamma_barium: gammas[3],
        gamma_strontium: gammas[4],
        gamma_ferrous: gammas[5],
        gamma_sulphate: gammas[6],
        gamma_carbonate: gammas[7],
        gamma_bicarbonate: gammas[8],
        warnings,
    })
}

/// The binary mean molal activity coefficient at a salt molality.
fn mean_activity_coefficient(
    salt: &SaltSpec,
    molality: f64,
    temperature_k: f64,
    debye_huckel_a: f64,
) -> f64 {
    if molality == 0.0 {
        return 1.0;
    }
    let ionic_strength = 0.5
        * molality
        * (salt.cation_stoichiometry * salt.cation_charge * salt.cation_charge
            + salt.anion_stoichiometry * salt.anion_charge * salt.anion_charge);
    let sqrt_i = ionic_strength.max(MIN_IONIC_STRENGTH).sqrt();
    let a_phi = debye_huckel_a * std::f64::consts::LN_10 / 3.0;
    let f_gamma = -a_phi
        * (sqrt_i / (1.0 + PITZER_B * sqrt_i) + 2.0 / PITZER_B * (1.0 + PITZER_B * sqrt_i).ln());

    let beta0 = temperature_parameter(family(salt, Family::B0), temperature_k);
    let beta1 = temperature_parameter(family(salt, Family::B1), temperature_k);
    let beta2 = temperature_parameter(family(salt, Family::B2), temperature_k);
    let mut b_gamma = 2.0 * beta0 + beta_gamma_term(beta1, salt.alpha1, ionic_strength);
    if salt.alpha2 > 0.0 && beta2.abs() > 0.0 {
        b_gamma += beta_gamma_term(beta2, salt.alpha2, ionic_strength);
    }

    let nu = salt.cation_stoichiometry + salt.anion_stoichiometry;
    let b_factor = 2.0 * salt.cation_stoichiometry * salt.anion_stoichiometry / nu;
    let c_factor = 3.0 * (salt.cation_stoichiometry * salt.anion_stoichiometry).powf(1.5) / nu;
    let c_phi = temperature_parameter(family(salt, Family::C0), temperature_k);
    let ln_gamma = (salt.cation_charge * salt.anion_charge).abs() * f_gamma
        + molality * b_factor * b_gamma
        + molality * molality * c_factor * c_phi;
    ln_gamma.exp()
}

/// The salt's six-term set for one family, or zeros where the catalogue carries none.
fn family(salt: &SaltSpec, family: Family) -> [f64; 6] {
    find(family, &[salt.cation, salt.anion]).unwrap_or([0.0; 6])
}

fn beta_gamma_term(beta: f64, alpha: f64, ionic_strength: f64) -> f64 {
    if alpha <= 0.0 || beta == 0.0 {
        return 0.0;
    }
    let strength = ionic_strength.max(MIN_IONIC_STRENGTH);
    let x = alpha * strength.sqrt();
    let numerator = 1.0 - (1.0 + x - 0.5 * x * x) * (-x).exp();
    2.0 * beta * numerator / (alpha * alpha * strength)
}

/// The PHREEQC six-term temperature expression about 298.15 K.
fn temperature_parameter(coefficients: [f64; 6], temperature_k: f64) -> f64 {
    let tr = REFERENCE_TEMPERATURE_K;
    coefficients[0]
        + coefficients[1] * (1.0 / temperature_k - 1.0 / tr)
        + coefficients[2] * (temperature_k / tr).ln()
        + coefficients[3] * (temperature_k - tr)
        + coefficients[4] * (temperature_k * temperature_k - tr * tr)
        + coefficients[5] * (1.0 / (temperature_k * temperature_k) - 1.0 / (tr * tr))
}

/// The mean coefficient as if the configured strength were a salt's own molality.
fn mean_at_equivalent_strength(name: &str, ionic_strength: f64, t: f64, a: f64) -> f64 {
    let salt = salt_spec(name).expect("the model's own salts");
    let per_molality = 0.5
        * (salt.cation_stoichiometry * salt.cation_charge * salt.cation_charge
            + salt.anion_stoichiometry * salt.anion_charge * salt.anion_charge);
    mean_activity_coefficient(&salt, ionic_strength / per_molality, t, a)
}

/// The nine ions' coefficients at the configured strength, by the MacInnes mapping.
fn ion_coefficients(t: f64, a: f64, ionic_strength: f64) -> [f64; 9] {
    let strength = ionic_strength.max(MIN_IONIC_STRENGTH);
    let gamma_nacl = mean_at_equivalent_strength("nacl", strength, t, a);
    let divalent = |salt: &str| {
        let mean = mean_at_equivalent_strength(salt, strength, t, a);
        mean.powi(3) / (gamma_nacl * gamma_nacl)
    };
    let gamma_nahco3 = mean_at_equivalent_strength("nahco3", strength, t, a);
    [
        gamma_nacl,
        gamma_nacl,
        divalent("cacl2"),
        divalent("bacl2"),
        divalent("srcl2"),
        divalent("fecl2"),
        divalent("na2so4"),
        divalent("na2co3"),
        gamma_nahco3 * gamma_nahco3 / gamma_nacl,
    ]
}
