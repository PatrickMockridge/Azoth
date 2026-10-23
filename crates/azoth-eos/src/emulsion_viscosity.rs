//! `eos.emulsion_viscosity` - NeqSim's `EmulsionViscosityCalculator`.
//!
//! The Arirachakaran inversion point and one of six relative-viscosity models, with the
//! class's tightness and demulsifier corrections.

use azoth_core::units::{DynamicViscosity, MassDensity, pascal_seconds};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::{EmulsionType, EmulsionViscosityResult};
use crate::spec_gen;

#[derive(Clone, Copy, PartialEq, Eq)]
enum Model {
    Einstein,
    Taylor,
    Brinkman,
    PalRhodes,
    Woelflin,
    Richardson,
}

impl Model {
    fn parse(name: &str) -> Result<Self> {
        match name {
            "einstein" => Ok(Self::Einstein),
            "taylor" => Ok(Self::Taylor),
            "brinkman" => Ok(Self::Brinkman),
            "pal_rhodes" => Ok(Self::PalRhodes),
            "woelflin" => Ok(Self::Woelflin),
            "richardson" => Ok(Self::Richardson),
            other => Err(AzothError::invalid_input(
                "model",
                format!(
                    "`{other}` is not one of the models: einstein, taylor, brinkman, \
                     pal_rhodes, woelflin or richardson"
                ),
            )),
        }
    }
}

/// The emulsion viscosity, its inversion point and the phase that is continuous.
///
/// # Errors
/// * [`azoth_core::AzothError::OutOfRange`] if an input is outside the spec's hard range.
/// * [`azoth_core::AzothError::InvalidInput`] if the model name is unknown.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, pascal_seconds};
/// use azoth_eos::emulsion_viscosity::emulsion_viscosity;
///
/// let r = emulsion_viscosity(
///     pascal_seconds(0.010),
///     pascal_seconds(0.0005),
///     kilograms_per_cubic_meter(850.0),
///     kilograms_per_cubic_meter(1020.0),
///     0.30,
///     0.74,
///     "pal_rhodes",
///     0.5,
///     false,
///     0.0,
/// )?;
/// assert_eq!(r.emulsion_type.as_str(), "water_in_oil");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(clippy::too_many_arguments)] // The signature is the spec's declared inputs.
pub fn emulsion_viscosity(
    oil_viscosity: DynamicViscosity,
    water_viscosity: DynamicViscosity,
    oil_density: MassDensity,
    water_density: MassDensity,
    water_cut: f64,
    max_packing_fraction: f64,
    model: &str,
    tightness_factor: f64,
    demulsifier_present: bool,
    demulsifier_efficiency: f64,
) -> Result<EmulsionViscosityResult> {
    let spec = &spec_gen::EMULSION_VISCOSITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "oil_viscosity" => Some(oil_viscosity.value),
            "water_viscosity" => Some(water_viscosity.value),
            "oil_density" => Some(oil_density.value),
            "water_density" => Some(water_density.value),
            "water_cut" => Some(water_cut),
            "max_packing_fraction" => Some(max_packing_fraction),
            _ => None,
        },
        &mut warnings,
    )?;
    let model = Model::parse(model)?;

    let viscosity_ratio_phases = oil_viscosity.value / water_viscosity.value;
    let density_ratio = water_density.value / oil_density.value;
    let mut inversion = 1.0 / (1.0 + viscosity_ratio_phases.powf(0.25) * density_ratio.powf(0.5));
    // numerics-ok: the class clamps its own correlation to [0.20, 0.90], then shifts and
    // clamps again to [0.15, 0.95]; both are the class's transcription.
    inversion = inversion.clamp(0.20, 0.90);
    inversion = (inversion + (tightness_factor - 0.5) * 0.10).clamp(0.15, 0.95);

    let is_inverted = water_cut >= inversion;
    let (continuous, dispersed, phi) = if is_inverted {
        (water_viscosity.value, oil_viscosity.value, 1.0 - water_cut)
    } else {
        (oil_viscosity.value, water_viscosity.value, water_cut)
    };

    let relative_viscosity = if phi <= 0.0 {
        1.0
    } else if phi >= 1.0 {
        dispersed / continuous
    } else {
        match model {
            Model::Einstein => 1.0 + 2.5 * phi,
            Model::Taylor => {
                1.0 + 2.5 * phi * (dispersed + 0.4 * continuous) / (dispersed + continuous)
            }
            Model::Brinkman => (1.0 - phi).powf(-2.5),
            Model::PalRhodes => {
                if phi >= max_packing_fraction {
                    1.0e6
                } else {
                    (1.0 - phi / max_packing_fraction).powf(-2.5 * max_packing_fraction)
                }
            }
            Model::Woelflin => {
                let k = if tightness_factor < 0.33 {
                    2.5
                } else if tightness_factor < 0.67 {
                    4.0
                } else {
                    6.0
                };
                (k * phi).exp()
            }
            Model::Richardson => ((2.5 + 5.0 * tightness_factor) * phi).exp(),
        }
    };

    let mut effective = continuous * relative_viscosity;
    if tightness_factor > 0.5 {
        effective *= 1.0 + (tightness_factor - 0.5) * 2.0;
    }
    if demulsifier_present && demulsifier_efficiency > 0.0 {
        effective *= 1.0 - demulsifier_efficiency * 0.5;
    }

    Ok(EmulsionViscosityResult {
        effective_viscosity: pascal_seconds(effective),
        relative_viscosity,
        viscosity_ratio: effective / continuous,
        inversion_water_cut: inversion,
        emulsion_type: if is_inverted {
            EmulsionType::OilInWater
        } else {
            EmulsionType::WaterInOil
        },
        is_inverted,
        warnings,
    })
}
