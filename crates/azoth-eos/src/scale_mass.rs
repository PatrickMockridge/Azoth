//! `eos.scale_mass` - NeqSim's `ScaleMassCalculator`.
//!
//! The class's approximation: the limiting ion's concentration times `min(1 - 10^-SI, 0.99)`,
//! one arm per mineral and differing only in molar mass.

use azoth_core::units::{MolarConcentration, kilograms_per_cubic_meter, moles_per_cubic_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::results::ScaleMassResult;
use crate::spec_gen;

/// The five minerals' molar masses, in g/mol, as the class writes them.
const CACO3: f64 = 100.09;
const BASO4: f64 = 233.39;
const SRSO4: f64 = 183.68;
const CASO4: f64 = 136.14;
const FECO3: f64 = 115.85;

/// The mass a supersaturated water precipitates, by the class's approximation.
///
/// # Errors
/// * [`AzothError::InvalidInput`] for a name outside the five minerals.
/// * [`AzothError::OutOfRange`] if a concentration is negative.
///
/// # Example
/// ```
/// use azoth_core::units::{moles_per_cubic_meter};
/// use azoth_eos::scale_mass::scale_mass;
///
/// let r = scale_mass(
///     "caco3",
///     moles_per_cubic_meter(10.0),
///     moles_per_cubic_meter(20.0),
///     2.5,
/// )?;
/// assert!((r.precipitated_mass.value - 0.990891).abs() < 1e-12);
/// assert_eq!(r.fraction_excess, 0.99);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn scale_mass(
    mineral: &str,
    cation_concentration: MolarConcentration,
    anion_concentration: MolarConcentration,
    saturation_index: f64,
) -> Result<ScaleMassResult> {
    let spec = &spec_gen::SCALE_MASS_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "cation_concentration" => Some(cation_concentration.value),
            "anion_concentration" => Some(anion_concentration.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let molar_mass = match mineral {
        "caco3" => CACO3,
        "baso4" => BASO4,
        "srso4" => SRSO4,
        "caso4" => CASO4,
        "feco3" => FECO3,
        other => {
            return Err(AzothError::invalid_input(
                "mineral",
                format!("`{other}` is not one of `caco3`, `baso4`, `srso4`, `caso4` or `feco3`"),
            ));
        }
    };

    // The class's own units: mol/L and g/mol, with mg/L out.
    let cation = cation_concentration.value / 1000.0;
    let anion = anion_concentration.value / 1000.0;
    let limiting = cation.min(anion);
    let (fraction, moles) = if saturation_index <= 0.0 {
        (0.0, 0.0)
    } else {
        let fraction = (1.0 - 10.0f64.powf(-saturation_index)).min(0.99);
        (fraction, limiting * fraction)
    };
    let mass_mg_per_l = moles * molar_mass * 1000.0;

    Ok(ScaleMassResult {
        precipitated_mass: kilograms_per_cubic_meter(mass_mg_per_l * 1.0e-3),
        precipitated_moles: moles_per_cubic_meter(moles * 1000.0),
        fraction_excess: fraction,
        limiting_concentration: moles_per_cubic_meter(limiting * 1000.0),
        warnings,
    })
}
