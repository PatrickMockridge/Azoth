//! `eos.interface_properties` - NeqSim's phase-pair dispatch for interfacial tension.
//!
//! Spec: `specs/models/eos/interface_properties.toml`, which records
//! `InterfaceProperties.setInterfacialTensionModel`'s six arms, the string selector and
//! `getSurfaceTension`'s phase-type dispatch.

use azoth_core::units::{
    MassDensity, MolarMass, Pressure, ThermodynamicTemperature, kelvins, kilograms_per_cubic_meter,
    kilograms_per_mole, newtons_per_meter,
};
use azoth_core::{AzothError, Result, Warning, apply_checks};

use crate::cdft_surface_tension::cdft_surface_tension;
use crate::firozabadi_ramley_interface_tension::firozabadi_ramley_interface_tension;
use crate::gt_surface_tension::gt_surface_tension;
use crate::gt_surface_tension_simple::gt_surface_tension_simple;
use crate::hydrate_inhibitor_wt::{PhaseLabel, label};
use crate::lgt_surface_tension::lgt_surface_tension;
use crate::model_gen;
use crate::parachor_mixture_surface_tension::parachor_mixture_surface_tension;
use crate::results::InterfacePropertiesResult;
use crate::tp_multiflash::tp_multiflash;
use crate::{Cubic, Mixture};

/// NeqSim's gas constant.
const R: f64 = 8.314_462_1;
/// The class the gas-liquid and gas-aqueous slots install in model 0.
const PARACHOR: &str = "ParachorSurfaceTension";
/// The class the liquid-liquid slot installs in model 0, and every slot's Firozabadi arm.
const FIROZABADI: &str = "FirozabadiRamleyInterfaceTension";

/// One phase as the dispatch's calculators read it.
///
/// **Two volumes, because the class reads two.** `getMolarVolume()` is the untranslated
/// cubic root, which the profile calculators and CDFT take; `getPhysicalProperties()
/// .getDensity()` is the Peneloux-translated one, which the Parachor and Firozabadi
/// calculators divide their molar mass by.
struct PhaseState {
    x: Vec<f64>,
    z_factor: f64,
    molar_volume: f64,
    translated_density: f64,
    molar_mass: f64,
}

impl PhaseState {
    fn mass_density(&self) -> MassDensity {
        kilograms_per_cubic_meter(self.translated_density * self.molar_mass)
    }

    fn molar_volume(&self) -> f64 {
        self.molar_volume
    }

    fn molar_mass(&self) -> MolarMass {
        kilograms_per_mole(self.molar_mass)
    }
}

/// The interface tension of a flashed phase pair by NeqSim's dispatch.
///
/// **The pair decides the slot**: GAS/OIL takes the gas-liquid slot, GAS/AQUEOUS the
/// gas-aqueous slot, and every other pair - including a reversed oil/gas one - the
/// liquid-liquid slot. Models 3 and 5 are aliases of 0; the string selector's three
/// calculators can be asked for directly, and models 1, 2 and 4 install their classes.
///
/// The inputs are the union of the arms': the feed, the fluid, the influence parameters,
/// reference component and interaction matrix the profile calculators take, and the
/// selection. An arm reads only what its calculator uses.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] for a feed that is not two phases, an
///   unknown selection, a component with no parachor or molar mass, or vectors that are
///   not one entry per component.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, pascals};
/// use azoth_eos::Cubic;
/// use azoth_eos::databank;
/// use azoth_eos::interface_properties::interface_properties;
///
/// // The capture's gas/oil pair: methane/n-butane at 300 K and 20 bar, model 0.
/// let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], Cubic::Srk, None)?;
/// let r = interface_properties(
///     &mixture,
///     "model_0",
///     &[9.634812734820117e-21, 2.7571945977326497e-19],
///     1.0,
///     &[vec![0.0, 0.0], vec![0.0, 0.0]],
///     kelvins(300.0),
///     pascals(2.0e6),
///     &[0.5, 0.5],
/// )?;
/// assert!((r.sigma.value - 0.009424889197268286).abs() < 1e-6);
/// assert_eq!(r.selected_model, "ParachorSurfaceTension");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn interface_properties(
    mixture: &Mixture,
    selection: &str,
    influence_parameters: &[f64],
    reference_component: f64,
    pair_interactions: &[Vec<f64>],
    t: ThermodynamicTemperature,
    p: Pressure,
    z: &[f64],
) -> Result<InterfacePropertiesResult> {
    let spec = &model_gen::INTERFACE_PROPERTIES_SPEC;
    let mut warnings: Vec<Warning> = Vec::new();

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
            format!("`{reference_component}` is not a whole index in 0..{count}"),
        ));
    }
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

    let flash = tp_multiflash(mixture, t, p, z)?;
    if flash.phase_count != 2 {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "the feed splits into {} phase(s) at this state, so there is no pair to \
                 dispatch; the class's `getSurfaceTension` is called with two phases",
                flash.phase_count
            ),
        ));
    }
    let rt = R * t.value;
    // NeqSim's two cubic families carry different Peneloux constants and the same
    // `RACKETZ`; the RK form has no translation at all.
    let shift: Vec<f64> = match mixture.cubic() {
        Cubic::Srk => mixture
            .components()
            .iter()
            .map(|component| component.volume_shift_srk)
            .collect(),
        Cubic::Pr => mixture
            .components()
            .iter()
            .map(|component| component.volume_shift)
            .collect(),
        // RK, and every other form, carries no translation.
        _ => vec![0.0; count],
    };
    let mut phases: Vec<PhaseState> = (0..2)
        .map(|index| {
            let x = flash.x[index].to_vec();
            let mut molar_mass = 0.0;
            for (i, fraction) in x.iter().enumerate() {
                let mass = mixture.components()[i].molar_mass.ok_or_else(|| {
                    AzothError::invalid_input(
                        "components",
                        format!(
                            "component {i} carries no molar mass, and the mass densities the \
                             parachor and Firozabadi slots take need one"
                        ),
                    )
                })?;
                molar_mass += fraction * mass;
            }
            let molar_volume = flash.z_factor[index] * rt / p.value;
            let translated: f64 = x.iter().zip(&shift).map(|(f, c)| f * c).sum();
            Ok(PhaseState {
                x,
                z_factor: flash.z_factor[index],
                molar_volume,
                translated_density: 1.0 / (molar_volume - translated),
                molar_mass,
            })
        })
        .collect::<Result<Vec<_>>>()?;
    // NeqSim's system orders its phases lightest first, which is the order
    // `getSurfaceTension` hands the slot; the flash's own order is denser first.
    phases.sort_by(|a, b| b.molar_volume.total_cmp(&a.molar_volume));

    let reduced = mixture.reduced_parameters(t, p)?;
    let labels = [
        label(mixture, &reduced, &phases[0].x, phases[0].z_factor)?,
        label(mixture, &reduced, &phases[1].x, phases[1].z_factor)?,
    ];
    let slot = match (labels[0], labels[1]) {
        (PhaseLabel::Gas, PhaseLabel::Oil) => "gas_oil",
        (PhaseLabel::Gas, PhaseLabel::Aqueous) => "gas_aqueous",
        _ => "liquid_liquid",
    };

    let (sigma, class, arm_warnings) = match selection {
        "model_0" | "model_3" | "model_5" => {
            if slot == "liquid_liquid" {
                firozabadi_arm(mixture, &phases, t)?
            } else {
                parachor_arm(mixture, &phases)?
            }
        }
        "parachor" => parachor_arm(mixture, &phases)?,
        "firozabadi_ramley" => firozabadi_arm(mixture, &phases, t)?,
        "model_1" => {
            if slot == "gas_oil" {
                let result = gt_surface_tension(mixture, influence_parameters, t, p, z)?;
                (result.sigma.value, "GTSurfaceTension", result.warnings)
            } else {
                let result = gt_surface_tension_simple(
                    mixture,
                    influence_parameters,
                    reference_component,
                    pair_interactions,
                    t,
                    p,
                    z,
                )?;
                (
                    result.sigma.value,
                    "GTSurfaceTensionSimple",
                    result.warnings,
                )
            }
        }
        "model_2" => {
            let result = lgt_surface_tension(
                mixture,
                influence_parameters,
                reference_component,
                pair_interactions,
                t,
                p,
                z,
            )?;
            (result.sigma.value, "LGTSurfaceTension", result.warnings)
        }
        "model_4" => match slot {
            "gas_oil" => {
                let result = gt_surface_tension_simple(
                    mixture,
                    influence_parameters,
                    reference_component,
                    pair_interactions,
                    t,
                    p,
                    z,
                )?;
                (
                    result.sigma.value,
                    "GTSurfaceTensionSimple",
                    result.warnings,
                )
            }
            "gas_aqueous" => {
                let result = lgt_surface_tension(
                    mixture,
                    influence_parameters,
                    reference_component,
                    pair_interactions,
                    t,
                    p,
                    z,
                )?;
                (result.sigma.value, "LGTSurfaceTension", result.warnings)
            }
            _ => parachor_arm(mixture, &phases)?,
        },
        "cdft" => {
            let result = cdft_surface_tension(
                mixture,
                t,
                &phases[0].x,
                &phases[1].x,
                phases[0].molar_volume(),
                phases[1].molar_volume(),
            )?;
            (result.sigma.value, "CDFTSurfaceTension", result.warnings)
        }
        other => {
            return Err(AzothError::invalid_input(
                "selection",
                format!(
                    "`{other}` is not one of `model_0`, `model_3`, `model_5`, `parachor`, \
                     `firozabadi_ramley`, `model_1`, `model_2`, `model_4` or `cdft`"
                ),
            ));
        }
    };
    warnings.extend(arm_warnings);

    Ok(InterfacePropertiesResult {
        sigma: newtons_per_meter(sigma),
        selected_model: class.to_string(),
        warnings,
    })
}

/// The Parachor slot: the databank's component parachors over both phases' mass densities.
fn parachor_arm(
    mixture: &Mixture,
    phases: &[PhaseState],
) -> Result<(f64, &'static str, Vec<Warning>)> {
    let mut parachors = Vec::with_capacity(mixture.len());
    for (index, component) in mixture.components().iter().enumerate() {
        if component.parachor <= 0.0 {
            return Err(AzothError::invalid_input(
                "components",
                format!(
                    "component {index} carries no parachor, so the slot's sum would count it \
                     as zero; the class reads `getParachor()` and would answer a number built \
                     on a missing coefficient"
                ),
            ));
        }
        parachors.push(component.parachor);
    }
    let result = parachor_mixture_surface_tension(
        &parachors,
        phases[0].mass_density(),
        phases[0].molar_mass(),
        &phases[0].x,
        phases[1].mass_density(),
        phases[1].molar_mass(),
        &phases[1].x,
    )?;
    Ok((result.sigma.value, PARACHOR, result.warnings))
}

/// The Firozabadi-Ramley slot: the pseudo-critical temperature is the first phase's.
fn firozabadi_arm(
    mixture: &Mixture,
    phases: &[PhaseState],
    t: ThermodynamicTemperature,
) -> Result<(f64, &'static str, Vec<Warning>)> {
    let mut pseudo_critical = 0.0;
    for (index, fraction) in phases[0].x.iter().enumerate() {
        pseudo_critical += fraction * mixture.components()[index].tc.value;
    }
    let result = firozabadi_ramley_interface_tension(
        t,
        kelvins(pseudo_critical),
        phases[0].mass_density(),
        phases[1].mass_density(),
    )?;
    Ok((result.sigma.value, FIROZABADI, result.warnings))
}
