//! `eos.electrolyte_scale_calculator` - NeqSim's `ElectrolyteScaleCalculator`.
//!
//! An ionic strength from the water's ions, Davies activity coefficients over a
//! temperature-adjusted Debye-Huckel A, and five saturation indices - calcite, barite,
//! gypsum, celestite and brucite - over the class's own temperature correlations.

use azoth_core::units::{MassDensity, Pressure, ThermodynamicTemperature, moles_per_kilogram};
use azoth_core::{Result, apply_checks};

use crate::results::ElectrolyteScaleCalculatorResult;
use crate::spec_gen;

/// The class's molar masses, in g/mol, as written upstream.
pub(crate) const MM_CALCIUM: f64 = 40.08;
pub(crate) const MM_BARIUM: f64 = 137.33;
pub(crate) const MM_STRONTIUM: f64 = 87.62;
pub(crate) const MM_MAGNESIUM: f64 = 24.31;
pub(crate) const MM_SODIUM: f64 = 22.99;
pub(crate) const MM_POTASSIUM: f64 = 39.10;
pub(crate) const MM_IRON: f64 = 55.85;
pub(crate) const MM_CHLORIDE: f64 = 35.45;
pub(crate) const MM_SULPHATE: f64 = 96.06;
pub(crate) const MM_BICARBONATE: f64 = 61.02;
pub(crate) const MM_CARBONATE: f64 = 60.01;

/// The class's own `2.302585093`, which is `ln 10` rounded to nine decimals. Kept as the
/// literal upstream writes, because the port reproduces its arithmetic.
#[allow(clippy::approx_constant)]
const LN_10_UPSTREAM: f64 = 2.302585093;

/// The class's inputs in its own units: mg/L, degrees Celsius and bara.
pub(crate) struct ElectrolyteInputs {
    pub calcium_mg_l: f64,
    pub barium_mg_l: f64,
    pub strontium_mg_l: f64,
    pub magnesium_mg_l: f64,
    pub sodium_mg_l: f64,
    pub potassium_mg_l: f64,
    pub iron_mg_l: f64,
    pub chloride_mg_l: f64,
    pub sulphate_mg_l: f64,
    pub bicarbonate_mg_l: f64,
    pub carbonate_mg_l: f64,
    pub temperature_c: f64,
    pub pressure_bara: f64,
    pub ph: f64,
}

/// `calculate()`'s state in the class's own units: mol/kg and bara.
pub(crate) struct ElectrolyteState {
    pub ionic_strength: f64,
    pub debye_hueckel_a: f64,
    pub gamma_divalent: f64,
    pub gamma_monovalent: f64,
    pub si_caco3: f64,
    pub si_baso4: f64,
    pub si_caso4: f64,
    pub si_srso4: f64,
    pub si_brucite: f64,
    pub hydroxide_activity: f64,
}

/// The activity-corrected saturation indices of a brine.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
/// use azoth_eos::electrolyte_scale_calculator::electrolyte_scale_calculator;
///
/// let r = electrolyte_scale_calculator(
///     kilograms_per_cubic_meter(0.42),   // calcium 420 mg/L
///     kilograms_per_cubic_meter(0.00005),
///     kilograms_per_cubic_meter(0.008),
///     kilograms_per_cubic_meter(1.28),
///     kilograms_per_cubic_meter(10.8),
///     kilograms_per_cubic_meter(0.4),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(19.4),
///     kilograms_per_cubic_meter(2.7),
///     kilograms_per_cubic_meter(0.142),
///     kilograms_per_cubic_meter(0.0),
///     kelvins(333.15),
///     pascals(5.0e6),
///     6.5,
/// )?;
/// assert!((r.si_caco3 - -0.8060791957491737).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn electrolyte_scale_calculator(
    calcium: MassDensity,
    barium: MassDensity,
    strontium: MassDensity,
    magnesium: MassDensity,
    sodium: MassDensity,
    potassium: MassDensity,
    iron: MassDensity,
    chloride: MassDensity,
    sulphate: MassDensity,
    bicarbonate: MassDensity,
    carbonate: MassDensity,
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
    ph: f64,
) -> Result<ElectrolyteScaleCalculatorResult> {
    let spec = &spec_gen::ELECTROLYTE_SCALE_CALCULATOR_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "calcium" => Some(calcium.value),
            "barium" => Some(barium.value),
            "strontium" => Some(strontium.value),
            "magnesium" => Some(magnesium.value),
            "sodium" => Some(sodium.value),
            "potassium" => Some(potassium.value),
            "iron" => Some(iron.value),
            "chloride" => Some(chloride.value),
            "sulphate" => Some(sulphate.value),
            "bicarbonate" => Some(bicarbonate.value),
            "carbonate" => Some(carbonate.value),
            "ph" => Some(ph),
            _ => None,
        },
        &mut warnings,
    )?;

    let inputs = ElectrolyteInputs {
        calcium_mg_l: calcium.value * 1000.0,
        barium_mg_l: barium.value * 1000.0,
        strontium_mg_l: strontium.value * 1000.0,
        magnesium_mg_l: magnesium.value * 1000.0,
        sodium_mg_l: sodium.value * 1000.0,
        potassium_mg_l: potassium.value * 1000.0,
        iron_mg_l: iron.value * 1000.0,
        chloride_mg_l: chloride.value * 1000.0,
        sulphate_mg_l: sulphate.value * 1000.0,
        bicarbonate_mg_l: bicarbonate.value * 1000.0,
        carbonate_mg_l: carbonate.value * 1000.0,
        temperature_c: temperature.value - 273.15,
        pressure_bara: pressure.value / 1.0e5,
        ph,
    };
    let state = state(&inputs);

    Ok(ElectrolyteScaleCalculatorResult {
        ionic_strength: moles_per_kilogram(state.ionic_strength),
        debye_hueckel_a: state.debye_hueckel_a,
        divalent_activity_coefficient: state.gamma_divalent,
        monovalent_activity_coefficient: state.gamma_monovalent,
        si_caco3: state.si_caco3,
        si_baso4: state.si_baso4,
        si_caso4: state.si_caso4,
        si_srso4: state.si_srso4,
        si_brucite: state.si_brucite,
        hydroxide_activity: state.hydroxide_activity,
        warnings,
    })
}

/// `calculate()` in the class's own units, which the brine-mixing model calls.
pub(crate) fn state(inputs: &ElectrolyteInputs) -> ElectrolyteState {
    let m_ca = inputs.calcium_mg_l / MM_CALCIUM / 1000.0;
    let m_ba = inputs.barium_mg_l / MM_BARIUM / 1000.0;
    let m_sr = inputs.strontium_mg_l / MM_STRONTIUM / 1000.0;
    let m_mg = inputs.magnesium_mg_l / MM_MAGNESIUM / 1000.0;
    let m_na = inputs.sodium_mg_l / MM_SODIUM / 1000.0;
    let m_k = inputs.potassium_mg_l / MM_POTASSIUM / 1000.0;
    let m_fe = inputs.iron_mg_l / MM_IRON / 1000.0;
    let m_cl = inputs.chloride_mg_l / MM_CHLORIDE / 1000.0;
    let m_so4 = inputs.sulphate_mg_l / MM_SULPHATE / 1000.0;
    let m_hco3 = inputs.bicarbonate_mg_l / MM_BICARBONATE / 1000.0;
    let m_co3 = inputs.carbonate_mg_l / MM_CARBONATE / 1000.0;

    let ionic_strength = 0.5
        * (m_ca * 4.0
            + m_ba * 4.0
            + m_sr * 4.0
            + m_mg * 4.0
            + m_na
            + m_k
            + m_fe * 4.0
            + m_cl
            + m_so4 * 4.0
            + m_hco3
            + m_co3 * 4.0);

    let mut debye_hueckel_a = 0.5092 + 8.5e-4 * (inputs.temperature_c - 25.0);
    if debye_hueckel_a < 0.5 {
        debye_hueckel_a = 0.5;
    }

    let gamma_divalent = 10.0f64.powf(davies_log_gamma(2, ionic_strength, debye_hueckel_a));
    let gamma_monovalent = 10.0f64.powf(davies_log_gamma(1, ionic_strength, debye_hueckel_a));

    let temperature_k = inputs.temperature_c + 273.15;
    let pk2 = 10.329 + 0.013 * (temperature_k - 298.15) / 25.0;
    let m_co3_eff = m_hco3 * 10.0f64.powf(inputs.ph - pk2);

    let p_ksp_caco3 = p_ksp(inputs, temperature_k, caco3_ksp(temperature_k), -58.4);
    let p_ksp_baso4 = p_ksp(inputs, temperature_k, baso4_ksp(temperature_k), -46.4);
    let p_ksp_caso4 = p_ksp(inputs, temperature_k, caso4_ksp(temperature_k), -33.0);
    let p_ksp_srso4 = p_ksp(inputs, temperature_k, srso4_ksp(temperature_k), -47.0);

    let si_caco3 = saturation_index(
        gamma_divalent * m_ca.max(1e-15),
        gamma_divalent * (m_co3_eff + m_co3).max(1e-15),
        p_ksp_caco3,
    );
    let si_baso4 = saturation_index(
        gamma_divalent * m_ba.max(1e-15),
        gamma_divalent * m_so4.max(1e-15),
        p_ksp_baso4,
    );
    let si_caso4 = saturation_index(
        gamma_divalent * m_ca.max(1e-15),
        gamma_divalent * m_so4.max(1e-15),
        p_ksp_caso4,
    );
    let si_srso4 = saturation_index(
        gamma_divalent * m_sr.max(1e-15),
        gamma_divalent * m_so4.max(1e-15),
        p_ksp_srso4,
    );

    let hydroxide_activity = 10.0f64.powf(-p_kw(temperature_k) + inputs.ph);
    let p_ksp_brucite = p_ksp(inputs, temperature_k, brucite_ksp(temperature_k), -20.0);
    let si_brucite = saturation_index(
        gamma_divalent * m_mg.max(1e-15),
        hydroxide_activity * hydroxide_activity,
        p_ksp_brucite,
    );

    ElectrolyteState {
        ionic_strength,
        debye_hueckel_a,
        gamma_divalent,
        gamma_monovalent,
        si_caco3,
        si_baso4,
        si_caso4,
        si_srso4,
        si_brucite,
        hydroxide_activity,
    }
}

fn davies_log_gamma(charge: i32, ionic_strength: f64, debye_hueckel_a: f64) -> f64 {
    if ionic_strength <= 0.0 {
        return 0.0;
    }
    let sqrt_i = ionic_strength.sqrt();
    let z = f64::from(charge);
    -debye_hueckel_a * z * z * (sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_strength)
}

fn saturation_index(activity_m: f64, activity_x: f64, p_ksp: f64) -> f64 {
    if activity_m <= 0.0 || activity_x <= 0.0 {
        return -10.0;
    }
    (activity_m * activity_x).log10() + p_ksp
}

/// `-log10(Ksp)` after the class's pressure correction, which its `pKsp` helper applies.
fn p_ksp(inputs: &ElectrolyteInputs, temperature_k: f64, ksp0: f64, volume_delta: f64) -> f64 {
    let corrected = if inputs.pressure_bara <= 1.013 || volume_delta.abs() < 1.0e-10 {
        ksp0
    } else {
        let r_cm3_bar = 83.1446;
        let delta_p = inputs.pressure_bara - 1.01325;
        ksp0 * (-volume_delta * delta_p / (r_cm3_bar * temperature_k)).exp()
    };
    -corrected.log10()
}

fn caco3_ksp(temperature_k: f64) -> f64 {
    let log10_ksp = -171.9065 - 0.077993 * temperature_k
        + 2839.319 / temperature_k
        + 71.595 * temperature_k.log10();
    10.0f64.powf(log10_ksp)
}

fn baso4_ksp(temperature_k: f64) -> f64 {
    let ln_ksp = 37588.0 / temperature_k - 747.610 + temperature_k.ln() * 119.280
        - temperature_k * 0.162830
        - 2880000.0 / (temperature_k * temperature_k);
    ln_ksp.exp()
}

fn caso4_ksp(temperature_k: f64) -> f64 {
    let ln_ksp = -26309.9 / temperature_k
        + 815.978
        + temperature_k.ln() * -138.361
        + temperature_k * 0.167863
        + 18.6143 / (temperature_k * temperature_k);
    ln_ksp.exp()
}

fn srso4_ksp(temperature_k: f64) -> f64 {
    let ln_ksp = -10452.9 / temperature_k + 19.790;
    ln_ksp.exp()
}

fn brucite_ksp(temperature_k: f64) -> f64 {
    let log10_ksp25 = -11.16;
    let delta_h_kj = -113.4;
    let gas_constant_kj = 0.0083145;
    let log10_ksp = log10_ksp25
        - (delta_h_kj / (LN_10_UPSTREAM * gas_constant_kj)) * (1.0 / temperature_k - 1.0 / 298.15);
    10.0f64.powf(log10_ksp)
}

fn p_kw(temperature_k: f64) -> f64 {
    4470.99 / temperature_k - 6.0875 + 0.01706 * temperature_k
}
