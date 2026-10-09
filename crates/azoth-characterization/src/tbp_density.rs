//! `characterization.tbp_density` - a cut's specific gravity from its boiling point and molar
//! mass.
//!
//! Spec: `specs/models/characterization/tbp_density.toml`. Oracle:
//! `validation/neqsim/captures/tbp_closure_probe.tsv`.
//!
//! # One member of four
//!
//! `TbpClosure` is four closures and only the 1980 Riazi-Daubert pair can be inverted for
//! specific gravity. The other three are refused by name, and each refusal is a `[[unported]]`
//! row rather than a sentence, so `tools/check_unported.py` holds this half and the Python one to
//! the same three keys. The reasons are NeqSim's own and are worth reading - the 1987 pair turns
//! over in specific gravity, so two gravities reproduce one molar mass and the inverse has no
//! unique root at all.
//!
//! # The pair is not symmetric
//!
//! `calcMolarMass` divides by 1000 to reach kg/mol and this multiplies by 1000 to reach g/mol, so
//! the two directions cannot share a conversion. A port that did would be out by `10^6` on a
//! density, which is the same trap the terms in `tbp_cut_properties`' `assumptions` record.

use azoth_core::units::{MolarMass, ThermodynamicTemperature, kilograms_per_cubic_meter};
use azoth_core::{Result, apply_checks};

use crate::model_gen;
use crate::results::TbpDensityResult;
use crate::tbp_closure::TbpClosureKind;

/// The 1980 Riazi-Daubert pair, rearranged for specific gravity. `[g/cm3]`.
///
/// NeqSim writes the exponent as `-1.0 / -1.0164`; the two negatives cancel and it is written
/// here as the one positive power it is.
fn riazi_daubert_1980_molar_mass_to_density(boiling_point: f64, molar_mass: f64) -> f64 {
    let molar_mass_gmol = molar_mass * 1000.0;
    let rankine = boiling_point * 1.8;
    (4.5673e-5 * rankine.powf(2.1962) / molar_mass_gmol).powf(1.0 / 1.0164)
}

/// A cut's normal liquid density, from its normal boiling point and its molar mass.
///
/// # Errors
/// * [`azoth_core::AzothError::InvalidInput`] for the three closures that cannot be inverted,
///   each naming the class that would close it.
/// * [`azoth_core::AzothError::OutOfRange`] if `boiling_point` or `molar_mass` is not positive.
///
/// # Example
/// ```
/// use azoth_core::units::{kelvins, kilograms_per_mole};
/// use azoth_characterization::{TbpClosureKind, tbp_density};
///
/// let r = tbp_density(
///     TbpClosureKind::RiaziDaubert1980,
///     kelvins(300.0),
///     kilograms_per_mole(0.2),
/// )?;
/// assert!((r.density.value - 234.340430850704).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn tbp_density(
    closure: TbpClosureKind,
    boiling_point: ThermodynamicTemperature,
    molar_mass: MolarMass,
) -> Result<TbpDensityResult> {
    let spec = &model_gen::TBP_DENSITY_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "boiling_point" => Some(boiling_point.value),
            "molar_mass" => Some(molar_mass.value),
            _ => None,
        },
        &mut warnings,
    )?;

    // The three refusals, each the row the spec declares. A key is one string literal so the
    // checker reads it without following control flow, here and in Python alike.
    match closure {
        TbpClosureKind::RiaziDaubert1980 => {}
        TbpClosureKind::RiaziDaubert1987 => {
            return Err(crate::unported::refuse("closure=riazi_daubert_1987"));
        }
        TbpClosureKind::Soreide => return Err(crate::unported::refuse("closure=soreide")),
        TbpClosureKind::TbpModel => return Err(crate::unported::refuse("closure=tbp_model")),
    }

    // The boundary: the correlation answers in g/cm3 and this declares kg/m3.
    let density = riazi_daubert_1980_molar_mass_to_density(boiling_point.value, molar_mass.value);

    Ok(TbpDensityResult {
        density: kilograms_per_cubic_meter(density * 1000.0),
        warnings,
    })
}
