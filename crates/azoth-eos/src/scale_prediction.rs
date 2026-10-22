//! `eos.scale_prediction` - NeqSim's `ScalePredictionCalculator`.
//!
//! The water-analysis route to mineral-scale saturation: an ionic strength from the TDS, a
//! Davies activity coefficient from a Debye-Huckel A built on water's density and dielectric,
//! six ion-pair corrections, and five saturation indices over the class's own Ksp set.

use azoth_core::units::{MassDensity, Pressure, ThermodynamicTemperature, moles_per_cubic_meter};
use azoth_core::{Result, apply_checks};

use crate::results::ScalePredictionResult;
use crate::spec_gen;

/// The class's molar masses, in mg/mol, and its two sodium constants.
pub(crate) const MW_CALCIUM: f64 = 40078.0;
pub(crate) const MW_BARIUM: f64 = 137327.0;
pub(crate) const MW_STRONTIUM: f64 = 87620.0;
pub(crate) const MW_IRON: f64 = 55845.0;
pub(crate) const MW_MAGNESIUM: f64 = 24305.0;
pub(crate) const MW_SODIUM: f64 = 22990.0;
/// The getter's sodium molar mass, which differs from [`MW_SODIUM`] in the class.
pub(crate) const MW_SODIUM_GETTER: f64 = 22989.77;
pub(crate) const MW_BICARBONATE: f64 = 61017.0;
pub(crate) const MW_SULPHATE: f64 = 96060.0;

/// The class's inputs in its own units: mg/L, degrees Celsius, bara and bar.
pub(crate) struct ScaleInputs {
    pub calcium_mg_l: f64,
    pub barium_mg_l: f64,
    pub strontium_mg_l: f64,
    pub iron_mg_l: f64,
    pub magnesium_mg_l: f64,
    pub sodium_mg_l: f64,
    pub bicarbonate_mg_l: f64,
    pub sulphate_mg_l: f64,
    pub tds_mg_l: f64,
    pub temperature_c: f64,
    pub pressure_bara: f64,
    pub co2_partial_pressure_bar: f64,
    pub ph: f64,
    pub auto_ph: bool,
    pub second_order_pressure_correction: bool,
}

/// The class's state, in its own units: mol/L, bar and (mol/L)^2.
pub(crate) struct ScaleState {
    pub ph_used: f64,
    pub ionic_strength_mol_l: f64,
    pub debye_huckel_a: f64,
    pub gamma_divalent: f64,
    pub gamma_monovalent: f64,
    pub free_calcium: f64,
    pub free_sulphate: f64,
    pub free_carbonate: f64,
    pub free_bicarbonate: f64,
    pub total_barium: f64,
    pub total_strontium: f64,
    pub total_iron: f64,
    pub total_sodium: f64,
    pub tds_mg_l: f64,
    pub ksp: [f64; 5],
    pub si: [f64; 5],
}

/// The saturation indices of a water analysis.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     kelvins, kilograms_per_cubic_meter, pascals,
/// };
/// use azoth_eos::scale_prediction::scale_prediction;
///
/// let r = scale_prediction(
///     kilograms_per_cubic_meter(0.4),   // calcium 400 mg/L
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.15),  // bicarbonate 150 mg/L
///     kilograms_per_cubic_meter(0.01),  // sulphate 10 mg/L
///     kilograms_per_cubic_meter(35.0),  // TDS 35000 mg/L
///     kelvins(353.15),
///     pascals(1.0e7),
///     pascals(2.0e5),
///     6.5,
///     false,
///     false,
/// )?;
/// assert!((r.si_caco3.expect("calcite has both ions") - -0.48139366009468343).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn scale_prediction(
    calcium: MassDensity,
    barium: MassDensity,
    strontium: MassDensity,
    iron: MassDensity,
    magnesium: MassDensity,
    sodium: MassDensity,
    bicarbonate: MassDensity,
    sulphate: MassDensity,
    total_dissolved_solids: MassDensity,
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
    co2_partial_pressure: Pressure,
    ph: f64,
    auto_ph: bool,
    second_order_pressure_correction: bool,
) -> Result<ScalePredictionResult> {
    let spec = &spec_gen::SCALE_PREDICTION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "co2_partial_pressure" => Some(co2_partial_pressure.value),
            "calcium" => Some(calcium.value),
            "barium" => Some(barium.value),
            "strontium" => Some(strontium.value),
            "iron" => Some(iron.value),
            "magnesium" => Some(magnesium.value),
            "sodium" => Some(sodium.value),
            "bicarbonate" => Some(bicarbonate.value),
            "sulphate" => Some(sulphate.value),
            "total_dissolved_solids" => Some(total_dissolved_solids.value),
            "ph" => Some(ph),
            _ => None,
        },
        &mut warnings,
    )?;

    let inputs = ScaleInputs {
        calcium_mg_l: calcium.value * 1000.0,
        barium_mg_l: barium.value * 1000.0,
        strontium_mg_l: strontium.value * 1000.0,
        iron_mg_l: iron.value * 1000.0,
        magnesium_mg_l: magnesium.value * 1000.0,
        sodium_mg_l: sodium.value * 1000.0,
        bicarbonate_mg_l: bicarbonate.value * 1000.0,
        sulphate_mg_l: sulphate.value * 1000.0,
        tds_mg_l: total_dissolved_solids.value * 1000.0,
        temperature_c: temperature.value - 273.15,
        pressure_bara: pressure.value / 1.0e5,
        co2_partial_pressure_bar: co2_partial_pressure.value / 1.0e5,
        ph,
        auto_ph,
        second_order_pressure_correction,
    };
    let state = state(&inputs);

    let si = state.si;
    let options = |index: usize| (!si[index].is_nan()).then_some(si[index]);
    let risk_count = si
        .iter()
        .filter(|value| value.is_finite() && **value > 0.0)
        .count() as u32;

    Ok(ScalePredictionResult {
        ionic_strength: moles_per_cubic_meter(state.ionic_strength_mol_l * 1000.0),
        debye_huckel_a: state.debye_huckel_a,
        divalent_activity_coefficient: state.gamma_divalent,
        monovalent_activity_coefficient: state.gamma_monovalent,
        ph_used: state.ph_used,
        free_calcium: moles_per_cubic_meter(state.free_calcium * 1000.0),
        free_sulphate: moles_per_cubic_meter(state.free_sulphate * 1000.0),
        free_carbonate: moles_per_cubic_meter(state.free_carbonate * 1000.0),
        free_bicarbonate: moles_per_cubic_meter(state.free_bicarbonate * 1000.0),
        total_barium: moles_per_cubic_meter(state.total_barium * 1000.0),
        total_strontium: moles_per_cubic_meter(state.total_strontium * 1000.0),
        total_iron: moles_per_cubic_meter(state.total_iron * 1000.0),
        total_sodium: moles_per_cubic_meter(state.total_sodium * 1000.0),
        ksp_calcite: state.ksp[0],
        ksp_barite: state.ksp[1],
        ksp_celestite: state.ksp[2],
        ksp_anhydrite: state.ksp[3],
        ksp_siderite: state.ksp[4],
        log10_ksp_calcite: state.ksp[0].log10(),
        log10_ksp_barite: state.ksp[1].log10(),
        log10_ksp_celestite: state.ksp[2].log10(),
        log10_ksp_anhydrite: state.ksp[3].log10(),
        log10_ksp_siderite: state.ksp[4].log10(),
        si_caco3: options(0),
        si_baso4: options(1),
        si_srso4: options(2),
        si_caso4: options(3),
        si_feco3: options(4),
        has_scaling_risk: risk_count > 0,
        scale_risk_count: risk_count,
        warnings,
    })
}

/// `calculate()` in the class's own units, which the later scale models call.
pub(crate) fn state(inputs: &ScaleInputs) -> ScaleState {
    let ionic_strength = inputs.tds_mg_l / 40000.0;
    let temperature_k = inputs.temperature_c + 273.15;
    let debye_huckel_a = debye_huckel_a(temperature_k);
    let gamma_divalent = activity_coefficient(2, ionic_strength, debye_huckel_a);
    let gamma_monovalent = activity_coefficient(1, ionic_strength, debye_huckel_a);

    let mut ph = inputs.ph;
    if inputs.auto_ph && inputs.co2_partial_pressure_bar > 0.0 && inputs.bicarbonate_mg_l > 0.0 {
        ph = estimate_ph(inputs, temperature_k);
    }

    let c_ca = inputs.calcium_mg_l / MW_CALCIUM;
    let c_ba = inputs.barium_mg_l / MW_BARIUM;
    let c_sr = inputs.strontium_mg_l / MW_STRONTIUM;
    let c_fe = inputs.iron_mg_l / MW_IRON;
    let c_mg = inputs.magnesium_mg_l / MW_MAGNESIUM;
    let c_na = inputs.sodium_mg_l / MW_SODIUM;
    let c_hco3 = inputs.bicarbonate_mg_l / MW_BICARBONATE;
    let c_so4 = inputs.sulphate_mg_l / MW_SULPHATE;
    let c_co3 = carbonate_from_bicarbonate(c_hco3, temperature_k, ph);

    let free = ion_pairing(
        c_ca,
        c_mg,
        c_na,
        c_so4,
        c_co3,
        c_hco3,
        temperature_k,
        gamma_divalent,
        gamma_monovalent,
    );
    let (free_ca, free_so4, free_co3, free_hco3) = free;

    let ksp = [
        calcite_ksp(temperature_k, inputs),
        barite_ksp(temperature_k, inputs),
        celestite_ksp(temperature_k, inputs),
        anhydrite_ksp(temperature_k, inputs),
        siderite_ksp(temperature_k, inputs),
    ];
    let iap_caco3 = (free_ca * gamma_divalent) * (free_co3 * gamma_divalent);
    let mut si = [f64::NAN; 5];
    si[0] = saturation_index(iap_caco3, ksp[0]);
    si[1] = if inputs.barium_mg_l > 0.0 && inputs.sulphate_mg_l > 0.0 {
        saturation_index(
            (c_ba * gamma_divalent) * (free_so4 * gamma_divalent),
            ksp[1],
        )
    } else {
        f64::NAN
    };
    si[2] = if inputs.strontium_mg_l > 0.0 && inputs.sulphate_mg_l > 0.0 {
        saturation_index(
            (c_sr * gamma_divalent) * (free_so4 * gamma_divalent),
            ksp[2],
        )
    } else {
        f64::NAN
    };
    si[3] = if inputs.sulphate_mg_l > 0.0 {
        saturation_index(
            (free_ca * gamma_divalent) * (free_so4 * gamma_divalent),
            ksp[3],
        )
    } else {
        f64::NAN
    };
    si[4] = if inputs.iron_mg_l > 0.0 {
        saturation_index(
            (c_fe * gamma_divalent) * (free_co3 * gamma_divalent),
            ksp[4],
        )
    } else {
        f64::NAN
    };

    ScaleState {
        ph_used: ph,
        ionic_strength_mol_l: ionic_strength,
        debye_huckel_a,
        gamma_divalent,
        gamma_monovalent,
        free_calcium: free_ca,
        free_sulphate: free_so4,
        free_carbonate: free_co3,
        free_bicarbonate: free_hco3,
        total_barium: c_ba,
        total_strontium: c_sr,
        total_iron: c_fe,
        total_sodium: inputs.sodium_mg_l / MW_SODIUM_GETTER,
        tds_mg_l: inputs.tds_mg_l,
        ksp,
        si,
    }
}

fn saturation_index(iap: f64, ksp: f64) -> f64 {
    if iap > 0.0 && ksp > 0.0 {
        (iap / ksp).log10()
    } else {
        f64::NAN
    }
}

/// The class's Debye-Huckel A: water density (Kell 1975) and dielectric (Archer-Wang 1990),
/// with the class's silent floors at 0.85 g/cm3 and 20.
fn debye_huckel_a(temperature_k: f64) -> f64 {
    let tc = temperature_k - 273.15;
    let mut rho = 0.99983 + 5.0948e-5 * tc - 7.5722e-6 * tc * tc + 3.8907e-8 * tc * tc * tc
        - 1.2e-10 * tc * tc * tc * tc;
    if rho < 0.85 {
        rho = 0.85;
    }
    let mut eps = 87.740 - 0.40008 * tc + 9.398e-4 * tc * tc - 1.410e-6 * tc * tc * tc;
    if eps < 20.0 {
        eps = 20.0;
    }
    1.8246e6 * rho.sqrt() / (eps * temperature_k).powf(1.5)
}

/// The Davies equation with Debye-Huckel A in log10 units.
fn activity_coefficient(charge: i32, ionic_strength: f64, debye_huckel_a: f64) -> f64 {
    let sqrt_i = ionic_strength.sqrt();
    let log_gamma = -debye_huckel_a
        * f64::from(charge * charge)
        * (sqrt_i / (1.0 + sqrt_i) - 0.3 * ionic_strength);
    10.0f64.powf(log_gamma)
}

/// `[CO3--]` from bicarbonate and pH: `K2 [HCO3-] / [H+]`.
fn carbonate_from_bicarbonate(c_hco3: f64, temperature_k: f64, ph: f64) -> f64 {
    let k2 = k2_carbonate(temperature_k);
    let h = 10.0f64.powf(-ph);
    if h > 0.0 { k2 * c_hco3 / h } else { 0.0 }
}

/// The class's pH estimate. **Its Henry constant is the class's as shipped** - about 4.9e48 at
/// 25 C - so auto pH lands near -45; see the spec.
fn estimate_ph(inputs: &ScaleInputs, temperature_k: f64) -> f64 {
    let kh = henry_co2(temperature_k);
    let k1 = k1_carbonate(temperature_k);
    let co2_aq = kh * inputs.co2_partial_pressure_bar;
    let c_hco3 = inputs.bicarbonate_mg_l / MW_BICARBONATE;
    if c_hco3 > 0.0 && co2_aq > 0.0 {
        -((k1 * co2_aq / c_hco3).log10())
    } else {
        7.0
    }
}

/// Five successive-substitution sweeps over six complexes, as the class runs them.
#[allow(clippy::too_many_arguments)] // The class's own list of ions.
fn ion_pairing(
    c_ca: f64,
    c_mg: f64,
    c_na: f64,
    c_so4: f64,
    c_co3: f64,
    c_hco3: f64,
    temperature_k: f64,
    gamma_divalent: f64,
    gamma_monovalent: f64,
) -> (f64, f64, f64, f64) {
    let correction = |log_k_25: f64, enthalpy: f64| {
        10.0f64.powf(log_k_25 + enthalpy / (2.303 * 8.314) * (1.0 / 298.15 - 1.0 / temperature_k))
    };
    let k_caso4 = correction(2.31, 6100.0);
    let k_mgso4 = correction(2.37, 5600.0);
    let k_naso4 = correction(0.70, 3000.0);
    let k_cahco3 = correction(1.11, 3600.0);
    let k_mghco3 = correction(1.16, 2900.0);
    let k_caco3 = correction(3.22, 14000.0);

    let mut free_ca = c_ca;
    let mut free_so4 = c_so4;
    let mut free_co3 = c_co3;
    let mut free_hco3 = c_hco3;
    let g2 = gamma_divalent * gamma_divalent;
    let gm2 = gamma_monovalent * gamma_monovalent;

    // The class computes six complex amounts here and never reads them; the sweeps below are
    // the whole of what it uses.
    for _ in 0..5 {
        free_ca = c_ca
            / (1.0
                + k_caso4 * g2 * free_so4 * g2
                + k_cahco3 * g2 * free_hco3 * gamma_monovalent / gamma_monovalent
                + k_caco3 * g2 * free_co3 * g2);
        if free_ca < 0.0 {
            free_ca = c_ca * 0.01;
        }
        let free_mg = {
            let value =
                c_mg / (1.0 + k_mgso4 * g2 * free_so4 * g2 + k_mghco3 * g2 * free_hco3 * gm2 / gm2);
            if value < 0.0 { c_mg * 0.01 } else { value }
        };
        let free_na = {
            let value =
                c_na / (1.0 + k_naso4 * gamma_monovalent * free_so4 * g2 / gamma_monovalent);
            if value < 0.0 { c_na * 0.01 } else { value }
        };
        free_so4 = c_so4
            / (1.0
                + k_caso4 * free_ca * g2 * g2
                + k_mgso4 * free_mg * g2 * g2
                + k_naso4 * free_na * gamma_monovalent * g2 / gamma_monovalent);
        if free_so4 < 0.0 {
            free_so4 = c_so4 * 0.01;
        }
        free_co3 = c_co3 / (1.0 + k_caco3 * free_ca * g2 * g2);
        if free_co3 < 0.0 {
            free_co3 = c_co3 * 0.01;
        }
        free_hco3 = c_hco3
            / (1.0
                + k_cahco3 * free_ca * g2 * gamma_monovalent / gamma_monovalent
                + k_mghco3 * free_mg * g2);
        if free_hco3 < 0.0 {
            free_hco3 = c_hco3 * 0.01;
        }
    }
    (free_ca, free_so4, free_co3, free_hco3)
}

/// The class's pressure correction, including the optional compressibility term.
fn pressure_correction(
    inputs: &ScaleInputs,
    ksp0: f64,
    temperature_k: f64,
    volume_delta: f64,
    kappa_delta: f64,
) -> f64 {
    if inputs.pressure_bara <= 1.013 || volume_delta.abs() < 1.0e-10 {
        return ksp0;
    }
    let r_cm3_bar = 83.1446;
    let delta_p = inputs.pressure_bara - 1.01325;
    let mut ln_correction = -volume_delta * delta_p / (r_cm3_bar * temperature_k);
    if inputs.second_order_pressure_correction && kappa_delta.abs() > 1.0e-12 {
        ln_correction += 0.5 * kappa_delta * delta_p * delta_p / (r_cm3_bar * temperature_k);
    }
    ksp0 * ln_correction.exp()
}

fn calcite_ksp(temperature_k: f64, inputs: &ScaleInputs) -> f64 {
    let log_ksp = -171.9065 - 0.077993 * temperature_k
        + 2839.319 / temperature_k
        + 71.595 * temperature_k.log10();
    pressure_correction(inputs, 10.0f64.powf(log_ksp), temperature_k, -58.4, -8.0e-3)
}

fn barite_ksp(temperature_k: f64, inputs: &ScaleInputs) -> f64 {
    let log_ksp = 136.035 - 7680.41 / temperature_k - 48.595 * temperature_k.log10();
    pressure_correction(inputs, 10.0f64.powf(log_ksp), temperature_k, -46.4, -6.0e-3)
}

fn celestite_ksp(temperature_k: f64, inputs: &ScaleInputs) -> f64 {
    let ln_ksp = -10452.9 / temperature_k + 19.790;
    pressure_correction(inputs, ln_ksp.exp(), temperature_k, -47.0, -6.5e-3)
}

fn anhydrite_ksp(temperature_k: f64, inputs: &ScaleInputs) -> f64 {
    let ln_ksp = -19966.8 / temperature_k + 454.860 - 69.840 * temperature_k.ln();
    pressure_correction(inputs, ln_ksp.exp(), temperature_k, -52.4, -7.0e-3)
}

fn siderite_ksp(temperature_k: f64, inputs: &ScaleInputs) -> f64 {
    let log_ksp = -59.3498 - 0.041377 * temperature_k - 2.1963 / temperature_k
        + 24.5724 * temperature_k.log10();
    pressure_correction(inputs, 10.0f64.powf(log_ksp), temperature_k, -52.9, -7.5e-3)
}

fn k1_carbonate(temperature_k: f64) -> f64 {
    let log_k1 = -356.3094 - 0.06091964 * temperature_k
        + 21834.37 / temperature_k
        + 126.8339 * temperature_k.log10()
        - 1684915.0 / (temperature_k * temperature_k);
    10.0f64.powf(log_k1)
}

fn k2_carbonate(temperature_k: f64) -> f64 {
    let log_k2 = -107.8871 - 0.03252849 * temperature_k
        + 5151.79 / temperature_k
        + 38.92561 * temperature_k.log10()
        - 563713.9 / (temperature_k * temperature_k);
    10.0f64.powf(log_k2)
}

/// The class's CO2 Henry constant, **as shipped**: its coefficients put `log10(KH)` at about
/// 48.7 at 25 C, so `[CO2(aq)]` is enormous and auto pH is about -45.
fn henry_co2(temperature_k: f64) -> f64 {
    let log_kh = -6.8346 + 1684.88 / temperature_k + 21.6215 * temperature_k.log10()
        - 0.012174 * temperature_k;
    10.0f64.powf(log_kh)
}
