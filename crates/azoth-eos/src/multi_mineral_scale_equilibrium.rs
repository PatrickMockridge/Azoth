//! `eos.multi_mineral_scale_equilibrium` - NeqSim's `MultiMineralScaleEquilibrium`.
//!
//! The coupled solid/liquid equilibrium of five minerals over one brine, by the class's greedy
//! coordinate descent: the single most-violating mineral is advanced by a relaxed step until no
//! extent exceeds the threshold, with under-saturated solids allowed to redissolve.

use azoth_core::units::{
    MassDensity, MolarConcentration, Pressure, ThermodynamicTemperature, Volume, kilograms,
    kilograms_per_cubic_meter, moles_per_cubic_meter,
};
use azoth_core::{Result, Warning, WarningCode, apply_checks};

use crate::model_gen;
use crate::pitzer_scale_activity::ion_coefficients;
use crate::results::{
    ActivityCoefficientRow, MineralEquilibriumRow, MultiMineralScaleEquilibriumResult,
};
use crate::scale_prediction::{ScaleInputs, state};

/// The class's free-ion floor, `EPS`.
const EPS: f64 = 1.0e-18;
/// The class's Debye-Huckel B, `DH_B`, at about 25 C.
const DH_B: f64 = 0.3283;
/// The class's molar masses, in g/mol.
const MW_BASO4: f64 = 233.39;
const MW_SRSO4: f64 = 183.68;
const MW_CASO4: f64 = 136.14;
const MW_CACO3: f64 = 100.09;
const MW_FECO3: f64 = 115.85;

/// The class's free-concentration vector layout.
const CA: usize = 0;
const BA: usize = 1;
const SR: usize = 2;
const FE: usize = 3;
const SO4: usize = 4;
const CO3: usize = 5;

/// The index of each ion in `pitzer_scale_activity`'s nine, for the Pitzer activity product.
const PITZER_CA: usize = 2;
const PITZER_BA: usize = 3;
const PITZER_SR: usize = 4;
const PITZER_FE: usize = 5;
const PITZER_SO4: usize = 6;
const PITZER_CO3: usize = 7;

/// The coupled multi-mineral solid/liquid equilibrium of a water analysis.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] for an activity-model name outside the three.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{
///     cubic_meters, kilograms_per_cubic_meter, kelvins, pascals,
/// };
/// use azoth_eos::multi_mineral_scale_equilibrium::multi_mineral_scale_equilibrium;
///
/// let r = multi_mineral_scale_equilibrium(
///     kilograms_per_cubic_meter(0.0), kilograms_per_cubic_meter(1.0),
///     kilograms_per_cubic_meter(0.0), kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0), kilograms_per_cubic_meter(0.0),
///     kilograms_per_cubic_meter(0.0), kilograms_per_cubic_meter(0.7),
///     kilograms_per_cubic_meter(0.02), kelvins(298.15), pascals(1.013e5),
///     pascals(0.0), 7.0, false, false, "davies", None, None, None,
///     Some(cubic_meters(0.001)), None, None, None,
/// )?;
/// assert!((r.total_mass.value - 1.6973864119588732).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn multi_mineral_scale_equilibrium(
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
    activity_model: &str,
    bdot_ion_size: Option<f64>,
    bdot_b: Option<f64>,
    pitzer_ionic_strength_molal: Option<f64>,
    water_volume: Option<Volume>,
    max_iterations: Option<usize>,
    relaxation: Option<f64>,
    convergence_threshold: Option<MolarConcentration>,
) -> Result<MultiMineralScaleEquilibriumResult> {
    let spec = &model_gen::MULTI_MINERAL_SCALE_EQUILIBRIUM_SPEC;
    let mut warnings = Vec::new();
    let threshold_mol_l = convergence_threshold.map_or(1.0e-15, |value| value.value * 1000.0);
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "calcium" => Some(calcium.value),
            "barium" => Some(barium.value),
            "strontium" => Some(strontium.value),
            "iron" => Some(iron.value),
            "magnesium" => Some(magnesium.value),
            "sodium" => Some(sodium.value),
            "bicarbonate" => Some(bicarbonate.value),
            "sulphate" => Some(sulphate.value),
            "total_dissolved_solids" => Some(total_dissolved_solids.value),
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "co2_partial_pressure" => Some(co2_partial_pressure.value),
            "ph" => Some(ph),
            "bdot_ion_size" => bdot_ion_size,
            "bdot_b" => bdot_b,
            "pitzer_ionic_strength_molal" => pitzer_ionic_strength_molal,
            "water_volume" => water_volume.map(|value| value.value),
            "relaxation" => relaxation,
            "convergence_threshold" => convergence_threshold.map(|value| value.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let model = match activity_model {
        "davies" | "bdot" | "pitzer_binary" => activity_model,
        other => {
            return Err(azoth_core::AzothError::invalid_input(
                "activity_model",
                format!("`{other}` is not one of davies, bdot, pitzer_binary"),
            ));
        }
    };
    let ion_size = bdot_ion_size.unwrap_or(5.0);
    let b_dot = bdot_b.unwrap_or(0.041);
    let max_iterations = max_iterations.unwrap_or(200_000);
    let relaxation = relaxation.unwrap_or(0.5);
    let litres = water_volume.map_or(1.0, |volume| volume.value * 1000.0);

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
    let temperature_k = temperature.value;

    let (gamma2, activity_coefficients, ion_products) = match model {
        "pitzer_binary" => {
            let sodium = state.total_sodium;
            let strength = pitzer_ionic_strength_molal.unwrap_or_else(|| {
                if sodium > 0.0 {
                    sodium
                } else {
                    (state.tds_mg_l / 58440.0).max(0.0)
                }
            });
            let gammas = ion_coefficients(temperature_k, state.debye_huckel_a, strength);
            let rows = [
                ("Na+", 0),
                ("Cl-", 1),
                ("Ca++", PITZER_CA),
                ("Ba++", PITZER_BA),
                ("Sr++", PITZER_SR),
                ("Fe++", PITZER_FE),
                ("SO4--", PITZER_SO4),
                ("CO3--", PITZER_CO3),
                ("HCO3-", 8),
            ]
            .iter()
            .map(|(ion, index)| ActivityCoefficientRow {
                ion: (*ion).to_string(),
                activity_coefficient: gammas[*index],
            })
            .collect();
            (gammas[PITZER_CA], rows, Some(gammas))
        }
        "bdot" => {
            let strength = state.ionic_strength_mol_l.max(0.0);
            let sqrt_i = strength.sqrt();
            let log_gamma = -state.debye_huckel_a * 4.0 * sqrt_i / (1.0 + ion_size * DH_B * sqrt_i)
                + b_dot * strength;
            let gamma = 10.0f64.powf(log_gamma);
            (
                gamma,
                vec![ActivityCoefficientRow {
                    ion: "divalent".to_string(),
                    activity_coefficient: gamma,
                }],
                None,
            )
        }
        _ => (
            state.gamma_divalent,
            vec![ActivityCoefficientRow {
                ion: "divalent".to_string(),
                activity_coefficient: state.gamma_divalent,
            }],
            None,
        ),
    };

    let product = |cation: usize, anion: usize| match &ion_products {
        Some(gammas) => {
            let of = |index: usize| match index {
                CA => gammas[PITZER_CA],
                BA => gammas[PITZER_BA],
                SR => gammas[PITZER_SR],
                FE => gammas[PITZER_FE],
                SO4 => gammas[PITZER_SO4],
                _ => gammas[PITZER_CO3],
            };
            of(cation) * of(anion)
        }
        None => gamma2 * gamma2,
    };

    let mut free = [0.0; 6];
    free[CA] = state.free_calcium.max(0.0);
    free[BA] = state.total_barium.max(0.0);
    free[SR] = state.total_strontium.max(0.0);
    free[FE] = state.total_iron.max(0.0);
    free[SO4] = state.free_sulphate.max(0.0);
    free[CO3] = state.free_carbonate.max(0.0);
    let initial_free = free;

    let mut minerals: Vec<Mineral> = Vec::new();
    let register = |free: &[f64; 6],
                    minerals: &mut Vec<Mineral>,
                    cation: usize,
                    anion: usize,
                    ksp: f64,
                    molar_mass: f64,
                    name: &str,
                    common: &str| {
        // A NaN Ksp is refused with the non-positive ones, as the class's `!(ksp > 0.0)` does.
        if free[cation] <= EPS || free[anion] <= EPS || ksp.is_nan() || ksp <= 0.0 {
            return;
        }
        minerals.push(Mineral {
            cation,
            anion,
            ksp,
            molar_mass,
            name: name.to_string(),
            common_name: common.to_string(),
            activity_coefficient_product: product(cation, anion),
            initial_saturation_index: 0.0,
            precipitated_mol_l: 0.0,
        });
    };
    register(
        &free,
        &mut minerals,
        BA,
        SO4,
        state.ksp[1],
        MW_BASO4,
        "BaSO4",
        "barite",
    );
    register(
        &free,
        &mut minerals,
        SR,
        SO4,
        state.ksp[2],
        MW_SRSO4,
        "SrSO4",
        "celestite",
    );
    register(
        &free,
        &mut minerals,
        CA,
        SO4,
        state.ksp[3],
        MW_CASO4,
        "CaSO4",
        "anhydrite",
    );
    register(
        &free,
        &mut minerals,
        CA,
        CO3,
        state.ksp[0],
        MW_CACO3,
        "CaCO3",
        "calcite",
    );
    register(
        &free,
        &mut minerals,
        FE,
        CO3,
        state.ksp[4],
        MW_FECO3,
        "FeCO3",
        "siderite",
    );

    for mineral in &mut minerals {
        mineral.initial_saturation_index = saturation_index(
            free[mineral.cation],
            free[mineral.anion],
            mineral.activity_coefficient_product,
            mineral.ksp,
        );
    }

    let mut iterations = 0usize;
    let mut stopped_before_limit = false;
    for _ in 0..max_iterations {
        let mut best: Option<usize> = None;
        let mut best_abs = 0.0;
        let mut best_dxi = 0.0;
        for (index, mineral) in minerals.iter().enumerate() {
            let cation = free[mineral.cation];
            let anion = free[mineral.anion];
            let k_prime = mineral.ksp / mineral.activity_coefficient_product;
            let b = cation + anion;
            let c = cation * anion - k_prime;
            let disc = (b * b - 4.0 * c).max(0.0);
            let xi_eq = 0.5 * (b - disc.sqrt());
            let dxi = if xi_eq > 0.0 {
                let cap = cation.min(anion) - EPS;
                xi_eq.min(cap).max(0.0)
            } else if xi_eq < 0.0 && mineral.precipitated_mol_l > EPS {
                xi_eq.max(-mineral.precipitated_mol_l)
            } else {
                0.0
            };
            if dxi.abs() > best_abs {
                best_abs = dxi.abs();
                best_dxi = dxi;
                best = Some(index);
            }
        }
        let Some(best) = best else {
            stopped_before_limit = true;
            break;
        };
        if best_abs < threshold_mol_l {
            stopped_before_limit = true;
            break;
        }
        let step = relaxation * best_dxi;
        let mineral = &mut minerals[best];
        free[mineral.cation] -= step;
        free[mineral.anion] -= step;
        mineral.precipitated_mol_l += step;
        if mineral.precipitated_mol_l < 0.0 {
            mineral.precipitated_mol_l = 0.0;
        }
        if free[mineral.cation] < 0.0 {
            free[mineral.cation] = 0.0;
        }
        if free[mineral.anion] < 0.0 {
            free[mineral.anion] = 0.0;
        }
        iterations += 1;
    }
    let iteration_limit_reached = !stopped_before_limit && iterations >= max_iterations;
    if iteration_limit_reached {
        warnings.push(Warning::new(
            WarningCode::SolverNotConverged,
            "the coordinate descent reached its iteration cap; the answer is the last iterate",
        ));
    }

    let mut rows = Vec::with_capacity(minerals.len());
    let mut total_mass_mg_l = 0.0;
    let mut violation = 0.0f64;
    for mineral in &minerals {
        let mass_mg_l = mineral.precipitated_mol_l * mineral.molar_mass * 1000.0;
        let final_si = saturation_index(
            free[mineral.cation],
            free[mineral.anion],
            mineral.activity_coefficient_product,
            mineral.ksp,
        );
        total_mass_mg_l += mass_mg_l;
        let limiting = (free[mineral.cation] <= 10.0 * EPS || free[mineral.anion] <= 10.0 * EPS)
            && mineral.precipitated_mol_l > EPS;
        rows.push(MineralEquilibriumRow {
            name: mineral.name.clone(),
            common_name: mineral.common_name.clone(),
            initial_saturation_index: mineral.initial_saturation_index,
            final_saturation_index: final_si,
            activity_coefficient_product: mineral.activity_coefficient_product,
            precipitated_amount: moles_per_cubic_meter(mineral.precipitated_mol_l * 1000.0),
            precipitated_mass: kilograms_per_cubic_meter(mass_mg_l * 1.0e-3),
            // mg/L is g/m3, so the per-cubic-metre mass in kg is `mass * 1e-3`, and the
            // volume it is over is `litres * 1e-3` m3.
            total_precipitated_mass: kilograms(mass_mg_l * 1.0e-6 * litres),
            limiting_ion_exhausted: limiting,
        });
        let residual = if mineral.precipitated_mol_l > EPS {
            final_si.abs()
        } else {
            final_si.max(0.0)
        };
        let residual = if residual.is_finite() {
            residual
        } else {
            f64::INFINITY
        };
        violation = violation.max(residual);
    }

    let mut precipitated_by_ion = [0.0; 6];
    for mineral in &minerals {
        precipitated_by_ion[mineral.cation] += mineral.precipitated_mol_l;
        precipitated_by_ion[mineral.anion] += mineral.precipitated_mol_l;
    }
    let ion_balance = (0..6)
        .map(|ion| (initial_free[ion] - free[ion] - precipitated_by_ion[ion]).abs())
        .fold(0.0, f64::max);

    Ok(MultiMineralScaleEquilibriumResult {
        minerals: rows,
        total_mass: kilograms_per_cubic_meter(total_mass_mg_l * 1.0e-3),
        activity_coefficients,
        divalent_activity_coefficient: gamma2,
        iterations: u32::try_from(iterations).unwrap_or(u32::MAX),
        iteration_limit_reached,
        maximum_complementarity_violation: violation,
        maximum_ion_balance_residual: moles_per_cubic_meter(ion_balance * 1000.0),
        residual_calcium: moles_per_cubic_meter(free[CA] * 1000.0),
        residual_barium: moles_per_cubic_meter(free[BA] * 1000.0),
        residual_strontium: moles_per_cubic_meter(free[SR] * 1000.0),
        residual_ferrous: moles_per_cubic_meter(free[FE] * 1000.0),
        residual_sulphate: moles_per_cubic_meter(free[SO4] * 1000.0),
        residual_carbonate: moles_per_cubic_meter(free[CO3] * 1000.0),
        warnings,
    })
}

/// One registered mineral and its ledger.
struct Mineral {
    cation: usize,
    anion: usize,
    ksp: f64,
    molar_mass: f64,
    name: String,
    common_name: String,
    activity_coefficient_product: f64,
    initial_saturation_index: f64,
    precipitated_mol_l: f64,
}

/// `SI = log10(product max(c_cat, EPS) max(c_an, EPS) / Ksp)`, NaN where either side is not positive.
fn saturation_index(c_cation: f64, c_anion: f64, product: f64, ksp: f64) -> f64 {
    let iap = product * c_cation.max(EPS) * c_anion.max(EPS);
    if iap <= 0.0 || ksp <= 0.0 {
        f64::NAN
    } else {
        (iap / ksp).log10()
    }
}
