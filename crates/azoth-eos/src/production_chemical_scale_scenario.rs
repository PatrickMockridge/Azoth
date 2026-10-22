//! `eos.production_chemical_scale_scenario` - NeqSim's `ProductionChemicalScaleScenario`.
//!
//! A produced water before and after chemical dosing: pH adjusters as acid/base equivalents
//! on a closed carbonate balance, H2S scavengers by capacity, scale inhibitors as kinetic
//! controls, and the coupled mineral equilibrium at both states.

use azoth_core::units::{
    MassDensity, Pressure, ThermodynamicTemperature, kilograms_per_cubic_meter,
    moles_per_cubic_meter,
};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::multi_mineral_scale_equilibrium::multi_mineral_scale_equilibrium;
use crate::results::{
    ChemicalWarningRow, MultiMineralScaleEquilibriumResult, ProductionChemicalScaleScenarioResult,
    ScenarioChangeRow, ScenarioIndexRow,
};

/// The class's constants: bicarbonate and H2S molar masses, and the 25 C water ion product.
const MW_HCO3: f64 = 61.016;
const MW_H2S: f64 = 34.08;
const WATER_KW_25C: f64 = 1.0e-14;

/// One dosed chemical, as the class's `ProductionChemical` record carries it.
#[derive(Debug, Clone, PartialEq)]
pub struct ProductionChemicalDose {
    /// The functional type: `acid`, `ph_adjuster`, `h2s_scavenger`, `scale_inhibitor` or `other`.
    pub chemical_type: String,
    /// The chemical's name, which its warning code carries.
    pub name: String,
    /// The active ingredient family.
    pub active_ingredient: String,
    /// The active ingredient's weight percent in the neat product.
    pub active_wt_pct: f64,
    /// The neat-product dose.
    pub dosage: MassDensity,
    /// The signed acid/base capacity per kg active.
    pub alkalinity_capacity: Option<f64>,
    /// Kilograms of H2S removed per kilogram active.
    pub h2s_capacity: Option<f64>,
}

fn warn(
    warnings: &mut Vec<ChemicalWarningRow>,
    code: impl Into<String>,
    message: impl Into<String>,
) {
    warnings.push(ChemicalWarningRow {
        code: code.into(),
        message: message.into(),
    });
}

fn key(active: &str) -> String {
    active.to_lowercase().replace(['-', ' '], "")
}

fn base_equivalents(active: &str, active_mg_l: f64) -> Option<f64> {
    let key = key(active);
    let grams_per_litre = active_mg_l / 1000.0;
    if key.contains("naoh") || key.contains("sodiumhydroxide") || key.contains("caustic") {
        return Some(grams_per_litre / 39.997);
    }
    if key.contains("na2co3") || key.contains("sodiumcarbonate") || key.contains("sodaash") {
        return Some(2.0 * grams_per_litre / 105.9888);
    }
    if key.contains("mdea") || key.contains("methyldiethanolamine") {
        return Some(grams_per_litre / 119.163);
    }
    if key == "mea" || key.contains("monoethanolamine") {
        return Some(grams_per_litre / 61.084);
    }
    if key.contains("hcl") || key.contains("hydrochloricacid") {
        return Some(-grams_per_litre / 36.461);
    }
    None
}

fn scavenger_capacity_mg_l(active: &str, active_mg_l: f64) -> Option<f64> {
    let key = key(active);
    let (molar_mass, mol_h2s_per_mol) = if key.contains("meatriazine") {
        (189.26, 1.0)
    } else if key.contains("mmatriazine") {
        (129.21, 1.0)
    } else if key.contains("iron") || key == "fe" {
        (55.845, 1.0)
    } else if key.contains("glyoxal") {
        (58.04, 1.0 / 1.5)
    } else {
        return None;
    };
    Some(active_mg_l / molar_mass * mol_h2s_per_mol * MW_H2S)
}

fn carbonate_k1(temperature_k: f64) -> f64 {
    let log_k1 = -356.3094 - 0.06091964 * temperature_k
        + 21834.37 / temperature_k
        + 126.8339 * temperature_k.log10()
        - 1684915.0 / (temperature_k * temperature_k);
    10.0f64.powf(log_k1)
}

fn carbonate_k2(temperature_k: f64) -> f64 {
    let log_k2 = -107.8871 - 0.03252849 * temperature_k
        + 5151.79 / temperature_k
        + 38.92561 * temperature_k.log10()
        - 563713.9 / (temperature_k * temperature_k);
    10.0f64.powf(log_k2)
}

fn carbonate_fractions(ph: f64, k1: f64, k2: f64) -> [f64; 3] {
    let h = 10.0f64.powf(-ph);
    let denominator = h * h + k1 * h + k1 * k2;
    [
        h * h / denominator,
        k1 * h / denominator,
        k1 * k2 / denominator,
    ]
}

fn alkalinity(total_carbon: f64, ph: f64, k1: f64, k2: f64) -> f64 {
    let alpha = carbonate_fractions(ph, k1, k2);
    let h = 10.0f64.powf(-ph);
    let oh = WATER_KW_25C / h;
    total_carbon * (alpha[1] + 2.0 * alpha[2]) + oh - h
}

/// The treatment-aware scale scenario of a produced water.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_cubic_meter, pascals};
/// use azoth_eos::production_chemical_scale_scenario::{
///     ProductionChemicalDose, production_chemical_scale_scenario,
/// };
///
/// let caustic = ProductionChemicalDose {
///     chemical_type: "ph_adjuster".to_owned(),
///     name: "caustic".to_owned(),
///     active_ingredient: "NaOH".to_owned(),
///     active_wt_pct: 50.0,
///     dosage: kilograms_per_cubic_meter(0.1),
///     alkalinity_capacity: None,
///     h2s_capacity: None,
/// };
/// let r = production_chemical_scale_scenario(
///     kilograms_per_cubic_meter(2.0),
///     kilograms_per_cubic_meter(0.05),
///     kilograms_per_cubic_meter(0.1),
///     kilograms_per_cubic_meter(0.005),
///     kilograms_per_cubic_meter(0.5),
///     kilograms_per_cubic_meter(20.0),
///     kilograms_per_cubic_meter(0.5),
///     kilograms_per_cubic_meter(0.8),
///     kilograms_per_cubic_meter(80.0),
///     6.5,
///     kilograms_per_cubic_meter(0.02),
///     kelvins(333.15),
///     pascals(1.0e7),
///     "davies",
///     None,
///     &[caustic],
/// )?;
/// assert!((r.treated_ph - 6.612286718254484).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn production_chemical_scale_scenario(
    calcium: MassDensity,
    barium: MassDensity,
    strontium: MassDensity,
    iron: MassDensity,
    magnesium: MassDensity,
    sodium: MassDensity,
    sulphate: MassDensity,
    bicarbonate: MassDensity,
    total_dissolved_solids: MassDensity,
    ph: f64,
    dissolved_h2s: MassDensity,
    temperature: ThermodynamicTemperature,
    pressure: Pressure,
    activity_model: &str,
    chemicals: &[ProductionChemicalDose],
    pitzer_ionic_strength_molal: Option<f64>,
) -> Result<ProductionChemicalScaleScenarioResult> {
    let spec = &model_gen::PRODUCTION_CHEMICAL_SCALE_SCENARIO_SPEC;
    let mut caveats = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "temperature" => Some(temperature.value),
            "pressure" => Some(pressure.value),
            "ph" => Some(ph),
            _ => None,
        },
        &mut caveats,
    )?;

    let calcium_mg_l = calcium.value * 1000.0;
    let bicarbonate_mg_l = bicarbonate.value * 1000.0;
    let dissolved_h2s_mg_l = dissolved_h2s.value * 1000.0;

    let mut equivalents_mol_l = 0.0;
    let mut capacity_mg_l = 0.0;
    let mut scale_inhibitor = false;
    let mut model_warnings = Vec::new();

    for chemical in chemicals {
        let active_mg_l =
            chemical.dosage.value.max(0.0) * 1000.0 * chemical.active_wt_pct.max(0.0) / 100.0;
        match chemical.chemical_type.as_str() {
            "acid" | "ph_adjuster" => {
                let equivalents = chemical.alkalinity_capacity.map_or_else(
                    || base_equivalents(&chemical.active_ingredient, active_mg_l),
                    |capacity| Some(capacity * active_mg_l * 1.0e-6),
                );
                match equivalents {
                    Some(value) => equivalents_mol_l += value,
                    None => warn(
                        &mut model_warnings,
                        format!("unsupported_ph_adjuster_{}", chemical.name),
                        format!(
                            "No acid/base-equivalent model for active ingredient {}",
                            chemical.active_ingredient
                        ),
                    ),
                }
            }
            "h2s_scavenger" => {
                let capacity = chemical.h2s_capacity.map_or_else(
                    || scavenger_capacity_mg_l(&chemical.active_ingredient, active_mg_l),
                    |value| Some(value * active_mg_l),
                );
                match capacity {
                    Some(value) => {
                        capacity_mg_l += value;
                        if key(&chemical.active_ingredient).contains("triazine") {
                            warn(
                                &mut model_warnings,
                                "triazine_spent_product",
                                "Check dithiazine/amorphous spent-product deposition; scavenger removal is not mineral equilibrium",
                            );
                        }
                    }
                    None => warn(
                        &mut model_warnings,
                        format!("unsupported_h2s_scavenger_{}", chemical.name),
                        format!(
                            "No H2S stoichiometric model for active ingredient {}",
                            chemical.active_ingredient
                        ),
                    ),
                }
            }
            "scale_inhibitor" => scale_inhibitor = true,
            _ => warn(
                &mut model_warnings,
                format!("no_equilibrium_effect_{}", chemical.name),
                "Chemical retained in audit trail; no supported stoichiometric effect on this aqueous scale scenario",
            ),
        }
    }

    let (treated_ph, treated_bicarbonate_mol_l) = apply_base_dose(
        ph,
        bicarbonate_mg_l,
        equivalents_mol_l,
        temperature.value,
        &mut model_warnings,
    );
    let treated_bicarbonate_mg_l = treated_bicarbonate_mol_l * MW_HCO3 * 1000.0;
    let residual_h2s_mg_l = (dissolved_h2s_mg_l - capacity_mg_l).max(0.0);

    let baseline = equilibrium(
        calcium_mg_l,
        barium.value * 1000.0,
        strontium.value * 1000.0,
        iron.value * 1000.0,
        magnesium.value * 1000.0,
        sodium.value * 1000.0,
        sulphate.value * 1000.0,
        bicarbonate_mg_l,
        total_dissolved_solids.value * 1000.0,
        ph,
        temperature.value,
        pressure.value,
        activity_model,
        pitzer_ionic_strength_molal,
    )?;
    let treated = equilibrium(
        calcium_mg_l,
        barium.value * 1000.0,
        strontium.value * 1000.0,
        iron.value * 1000.0,
        magnesium.value * 1000.0,
        sodium.value * 1000.0,
        sulphate.value * 1000.0,
        treated_bicarbonate_mg_l,
        total_dissolved_solids.value * 1000.0,
        treated_ph,
        temperature.value,
        pressure.value,
        activity_model,
        pitzer_ionic_strength_molal,
    )?;

    let baseline_si: Vec<ScenarioIndexRow> = baseline
        .minerals
        .iter()
        .map(|mineral| ScenarioIndexRow {
            mineral: mineral.name.clone(),
            saturation_index: Some(mineral.initial_saturation_index),
        })
        .collect();
    let treated_si: Vec<ScenarioIndexRow> = treated
        .minerals
        .iter()
        .map(|mineral| ScenarioIndexRow {
            mineral: mineral.name.clone(),
            saturation_index: Some(mineral.initial_saturation_index),
        })
        .collect();
    let si_change: Vec<ScenarioChangeRow> = baseline_si
        .iter()
        .map(|row| {
            let after = treated_si
                .iter()
                .find(|other| other.mineral == row.mineral)
                .and_then(|other| other.saturation_index);
            ScenarioChangeRow {
                mineral: row.mineral.clone(),
                change: match (row.saturation_index, after) {
                    (Some(before), Some(after)) => Some(after - before),
                    _ => None,
                },
            }
        })
        .collect();

    let delta_calcite = si_change
        .iter()
        .find(|row| row.mineral == "CaCO3")
        .and_then(|row| row.change);
    if scale_inhibitor {
        warn(
            &mut model_warnings,
            "inhibitor_is_kinetic",
            "Scale inhibitor does not change thermodynamic SI; evaluate MIC and growth control separately with ScaleControlAssessor",
        );
    }
    if let Some(delta) = delta_calcite {
        if delta > 0.25 {
            warn(
                &mut model_warnings,
                "chemical_induced_calcite_risk",
                format!("Base-equivalent dosing increased calcite SI by {delta:.2}"),
            );
        }
    }
    if capacity_mg_l > 0.0 && residual_h2s_mg_l > 0.0 {
        warn(
            &mut model_warnings,
            "h2s_scavenger_under_capacity",
            "Calculated dissolved-H2S load exceeds the practical active scavenger capacity",
        );
    }

    Ok(ProductionChemicalScaleScenarioResult {
        treated_ph,
        treated_bicarbonate: kilograms_per_cubic_meter(treated_bicarbonate_mg_l * 1e-3),
        base_dose_equivalents: moles_per_cubic_meter(equivalents_mol_l * 1000.0),
        scavenger_capacity: kilograms_per_cubic_meter(capacity_mg_l * 1e-3),
        residual_dissolved_h2s: kilograms_per_cubic_meter(residual_h2s_mg_l * 1e-3),
        baseline_si,
        treated_si,
        si_change,
        model_warnings,
        warnings: caveats,
    })
}

fn apply_base_dose(
    ph: f64,
    bicarbonate_mg_l: f64,
    equivalents_mol_l: f64,
    temperature_k: f64,
    warnings: &mut Vec<ChemicalWarningRow>,
) -> (f64, f64) {
    if bicarbonate_mg_l <= 0.0 || equivalents_mol_l == 0.0 {
        return (ph, bicarbonate_mg_l.max(0.0) / (MW_HCO3 * 1000.0));
    }
    let k1 = carbonate_k1(temperature_k);
    let k2 = carbonate_k2(temperature_k);
    let initial_hco3 = bicarbonate_mg_l / (MW_HCO3 * 1000.0);
    let initial_alpha1 = carbonate_fractions(ph, k1, k2)[1];
    if initial_alpha1 < 1.0e-12 {
        warn(
            warnings,
            "carbonate_balance_ill_conditioned",
            "Initial pH/HCO3 combination cannot define total inorganic carbon",
        );
        return (ph, initial_hco3);
    }
    let total_carbon = initial_hco3 / initial_alpha1;
    let target = alkalinity(total_carbon, ph, k1, k2) + equivalents_mol_l;

    let (mut low, mut high) = (2.0, 14.0);
    for _ in 0..100 {
        let mid = 0.5 * (low + high);
        if alkalinity(total_carbon, mid, k1, k2) < target {
            low = mid;
        } else {
            high = mid;
        }
    }
    let new_ph = 0.5 * (low + high);
    if new_ph > 13.99 {
        warn(
            warnings,
            "ph_solver_upper_bound",
            "Dose exceeds the closed-carbonate model pH 14 search bound",
        );
    }
    let new_hco3 = total_carbon * carbonate_fractions(new_ph, k1, k2)[1];
    (new_ph, new_hco3)
}

#[allow(clippy::too_many_arguments)]
fn equilibrium(
    calcium_mg_l: f64,
    barium_mg_l: f64,
    strontium_mg_l: f64,
    iron_mg_l: f64,
    magnesium_mg_l: f64,
    sodium_mg_l: f64,
    sulphate_mg_l: f64,
    bicarbonate_mg_l: f64,
    tds_mg_l: f64,
    ph: f64,
    temperature_k: f64,
    pressure_pa: f64,
    activity_model: &str,
    pitzer_ionic_strength_molal: Option<f64>,
) -> Result<MultiMineralScaleEquilibriumResult> {
    multi_mineral_scale_equilibrium(
        kilograms_per_cubic_meter(calcium_mg_l * 1e-3),
        kilograms_per_cubic_meter(barium_mg_l * 1e-3),
        kilograms_per_cubic_meter(strontium_mg_l * 1e-3),
        kilograms_per_cubic_meter(iron_mg_l * 1e-3),
        kilograms_per_cubic_meter(magnesium_mg_l * 1e-3),
        kilograms_per_cubic_meter(sodium_mg_l * 1e-3),
        kilograms_per_cubic_meter(bicarbonate_mg_l * 1e-3),
        kilograms_per_cubic_meter(sulphate_mg_l * 1e-3),
        kilograms_per_cubic_meter(tds_mg_l * 1e-3),
        azoth_core::units::kelvins(temperature_k),
        azoth_core::units::pascals(pressure_pa),
        azoth_core::units::pascals(0.0),
        ph,
        false,
        false,
        activity_model,
        None,
        None,
        pitzer_ionic_strength_molal,
        None,
        None,
        None,
        None,
    )
}
