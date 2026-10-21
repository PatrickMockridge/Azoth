//! `eos.interface_properties` - NeqSim's phase-pair dispatch for interfacial tension.
//!
//! Spec: `specs/models/eos/interface_properties.toml`, which records
//! `InterfaceProperties.setInterfacialTensionModel`'s model-0 arm and
//! `getSurfaceTension`'s phase-type dispatch.

use azoth_core::units::{MassDensity, MolarMass, ThermodynamicTemperature, kelvins};
use azoth_core::{AzothError, Result, apply_checks};

use crate::firozabadi_ramley_interface_tension::firozabadi_ramley_interface_tension;
use crate::model_gen;
use crate::parachor_mixture_surface_tension::parachor_mixture_surface_tension;
use crate::results::InterfacePropertiesResult;

/// The class the gas-liquid and gas-aqueous slots install in model 0.
const PARACHOR: &str = "ParachorSurfaceTension";
/// The class the liquid-liquid slot installs in model 0, and every slot's Firozabadi arm.
const FIROZABADI: &str = "FirozabadiRamleyInterfaceTension";

/// The interface tension of a phase pair by NeqSim's dispatch.
///
/// **The pair decides the slot**: `gas_oil` and `gas_aqueous` take Parachor in model 0,
/// and every other pair - including a reversed oil/gas one - takes Firozabadi-Ramley.
/// Models 3 and 5 are aliases of 0. The string selector's two calculators can be asked
/// for directly, and models 1, 2 and 4 are refused with the classes they need.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] for an unknown pair or selection, for the
///   three refused models, or for vectors that are not one entry per component.
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's range.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole, kelvins};
/// use azoth_eos::interface_properties::interface_properties;
///
/// // The capture's gas/oil pair: GAS first, OIL second, model 0's Parachor slot.
/// let r = interface_properties(
///     "gas_oil",
///     "model_0",
///     &[77.3, 191.7],
///     &[kelvins(190.56), kelvins(425.12)],
///     kilograms_per_cubic_meter(19.938328330315997),
///     kilograms_per_cubic_meter(549.3851668858096),
///     kilograms_per_mole(0.022959902730870334),
///     kilograms_per_mole(0.054107691623917306),
///     &[0.8356249351028912, 0.1643750648971088],
///     &[0.09542082642782022, 0.9045791735721797],
///     kelvins(300.0),
/// )?;
/// assert!((r.sigma.value - 0.009424889197268286).abs() < 1e-15);
/// assert_eq!(r.selected_model, "ParachorSurfaceTension");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T` is the symbol in the published equation
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn interface_properties(
    interface: &str,
    selection: &str,
    parachors: &[f64],
    tc: &[ThermodynamicTemperature],
    rho_1: MassDensity,
    rho_2: MassDensity,
    molar_mass_1: MolarMass,
    molar_mass_2: MolarMass,
    x_1: &[f64],
    x_2: &[f64],
    T: ThermodynamicTemperature,
) -> Result<InterfacePropertiesResult> {
    let spec = &model_gen::INTERFACE_PROPERTIES_SPEC;
    let mut warnings = Vec::new();

    let pair = match interface {
        "gas_oil" | "gas_aqueous" | "liquid_liquid" => interface,
        other => {
            return Err(AzothError::invalid_input(
                "interface",
                format!(
                    "`{other}` is not one of `gas_oil`, `gas_aqueous` or `liquid_liquid`; \
                     the pair decides which calculator slot runs"
                ),
            ));
        }
    };

    let use_firozabadi = match selection {
        "model_0" | "model_3" | "model_5" => pair == "liquid_liquid",
        "parachor" => false,
        "firozabadi_ramley" => true,
        "model_1" => {
            return Err(AzothError::invalid_input(
                "selection",
                "`model_1` installs `GTSurfaceTension` (`eos.gt_surface_tension`) in the \
                 gas-oil slot and `GTSurfaceTensionSimple` in the other two; all three are \
                 ported, and this dispatch does not carry the influence parameters the \
                 gas-oil calculator takes",
            ));
        }
        "model_2" => {
            return Err(AzothError::invalid_input(
                "selection",
                "`model_2` installs `LGTSurfaceTension` in all three slots. That calculator \
                 is `eos.lgt_surface_tension`, and this dispatch does not carry the influence \
                 parameters, reference component and interaction matrix it takes - so the arm \
                 is refused here rather than assumed",
            ));
        }
        "model_4" => {
            return Err(AzothError::invalid_input(
                "selection",
                "`model_4` installs `GTSurfaceTensionSimple` \
                 (`eos.gt_surface_tension_simple`), `LGTSurfaceTension` \
                 (`eos.lgt_surface_tension`) and Parachor in its three slots; all three are \
                 ported, and this dispatch does not carry the influence parameters, reference \
                 component and interaction matrix the two profile calculators take",
            ));
        }
        other => {
            return Err(AzothError::invalid_input(
                "selection",
                format!(
                    "`{other}` is not one of `model_0`, `model_3`, `model_5`, `parachor`, \
                     `firozabadi_ramley`, or the refused `model_1`, `model_2`, `model_4`"
                ),
            ));
        }
    };

    let count = x_1.len();
    if parachors.len() != count || tc.len() != count || x_2.len() != count {
        return Err(AzothError::invalid_input(
            "parachors",
            format!(
                "{} parachor(s), {} critical temperature(s), {} and {} mole fraction(s): \
                 the two phases share one component list, so every vector has one entry per \
                 component",
                parachors.len(),
                tc.len(),
                count,
                x_2.len()
            ),
        ));
    }

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "rho_1" => Some(rho_1.value),
            "rho_2" => Some(rho_2.value),
            "molar_mass_1" => Some(molar_mass_1.value),
            "molar_mass_2" => Some(molar_mass_2.value),
            "x_1" => x_1.iter().copied().reduce(f64::min),
            "x_2" => x_2.iter().copied().reduce(f64::min),
            _ => None,
        },
        &mut warnings,
    )?;

    let (sigma, selected_model) = if use_firozabadi {
        let pseudo_critical: f64 = x_1
            .iter()
            .zip(tc)
            .map(|(fraction, temperature)| fraction * temperature.value)
            .sum();
        let result =
            firozabadi_ramley_interface_tension(T, kelvins(pseudo_critical), rho_1, rho_2)?;
        warnings.extend(result.warnings);
        (result.sigma, FIROZABADI)
    } else {
        let result = parachor_mixture_surface_tension(
            parachors,
            rho_1,
            molar_mass_1,
            x_1,
            rho_2,
            molar_mass_2,
            x_2,
        )?;
        warnings.extend(result.warnings);
        (result.sigma, PARACHOR)
    };

    Ok(InterfacePropertiesResult {
        sigma,
        selected_model: selected_model.to_string(),
        warnings,
    })
}
