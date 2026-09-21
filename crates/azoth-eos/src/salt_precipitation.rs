//! `eos.salt_precipitation` - how much of one mineral precipitates from a brine.
//!
//! ```text
//! e  <-  the largest extent at which SR(e) = 1, up to min_i n_i / stoc_i
//! ```
//!
//! Spec: `specs/models/eos/salt_precipitation.toml`, which carries the bracket, the refusal
//! and the two captured states.
//!
//! NeqSim's `CalcSaltSatauration.precipitate()`, the per-mineral half of
//! `MultiSaltPrecipitation`'s complementarity loop: the loop holds one amount per mineral and
//! takes the largest `|SR - 1|` each step, so porting the loop is looping over this.
//!
//! # What the solve is
//!
//! The saturation ratio is wanted as a function of the extent the mineral has taken, which is
//! [`crate::scale_saturation_ratio`] with both ions' mole fractions moved down by `e` and the
//! whole renormalised - and with the activity coefficients [`crate::pitzer_phase`] gives for
//! the brine the extent *leaves*. That function is monotone in `e`, since taking ions out can
//! only lower the ion activity product, so the bracket is `[0, min_i n_i/stoc_i]`, the extent
//! at which an ion runs out, and the answer is where the ratio reaches one.
//!
//! **The brine is an electrolyte and not a cubic mixture**, which is why this takes names:
//! `mixture_of` refuses a cubic over an ion and rightly, so the coefficients come from the
//! Pitzer phase over the same composition rather than from a caller.
//!
//! **A mineral already under saturation takes nothing**, and its ratio is reported as it
//! stands rather than driven to one: that is the complementarity condition, and it is what
//! NeqSim's `if (initialSaturationRatio <= 1.0) return 0` says.

use azoth_core::units::{Pressure, ThermodynamicTemperature};
use azoth_core::{AzothError, Result, apply_checks};

use crate::databank;
use crate::model_gen;
use crate::pitzer_phase::pitzer_phase;
use crate::results::SaltPrecipitationResult;
use crate::scale_saturation_ratio;

/// The bisection's stopping rule on `|log10 SR|`, NeqSim's own.
const LOG10_TOLERANCE: f64 = 1.0e-8;

/// How far inside the maximum extent the bracket's upper end sits.
///
/// At the maximum an ion is exactly exhausted and the ratio is zero, and NeqSim brackets at
/// `maximum * (1 - 1e-12)` so that the upper end is a state the brine can still be in.
const BRACKET_INSET: f64 = 1.0e-12;

/// The solid one mineral takes from a brine, and the ratios that bound it.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the salt is not in the table, if the brine does not name
///   both of its ions and its water, or if `gammas` is not one per component.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
/// * [`AzothError::SolverNotConverged`] if the ratio is not under one at the maximum extent,
///   which is NeqSim's own `IllegalStateException` at the same place.
#[allow(clippy::too_many_arguments, non_snake_case)]
pub fn salt_precipitation(
    components: &[&str],
    salt: &str,
    T: ThermodynamicTemperature,
    P: Pressure,
    z: &[f64],
) -> Result<SaltPrecipitationResult> {
    let spec = &model_gen::SALT_PRECIPITATION_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T.value),
            "P" => Some(P.value),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = components.len();
    if z.len() != n {
        return Err(AzothError::invalid_input(
            "z",
            format!(
                "a brine of {n} components needs {n} mole fractions and got {}",
                z.len()
            ),
        ));
    }

    let record = databank::salt(salt).ok_or_else(|| {
        AzothError::invalid_input(
            "salt",
            format!("`{salt}` is not a row of `compsalt`, so it has no solubility product"),
        )
    })?;
    let index = |name: &str| -> Result<usize> {
        components
            .iter()
            .position(|candidate| candidate.trim().to_lowercase() == name.trim().to_lowercase())
            .ok_or_else(|| {
                AzothError::invalid_input(
                    "components",
                    format!(
                        "`{salt}` is built from `{name}`, and this brine does not carry it. A \
                     mineral the fluid cannot form is refused rather than reported as a zero, \
                     which would be an answer about a mineral that is not there"
                    ),
                )
            })
    };
    let (first, second, water) = (
        index(&record.cation)?,
        index(&record.anion)?,
        index("water")?,
    );

    let molar_mass_water = databank::entry("water", None)?.molar_mass.ok_or_else(|| {
        AzothError::invalid_input(
            "components",
            "water carries no molar mass, and a molality is per kilogram of it".to_string(),
        )
    })?;

    // The hydrogen ion, where the brine has one: `FeS`'s product is multiplied by its
    // molality and every other salt ignores it. Absent is the state NeqSim *skips* that salt
    // in, and the calc refuses it - so a brine without `H3O+` gets the refusal here too.
    let hydrogen = components
        .iter()
        .position(|candidate| candidate.trim().to_lowercase() == "h3o+")
        .map(|position| z[position] / (z[water] * molar_mass_water));

    // The ratio at one extent, which is the whole of the solve.
    let ratio = |extent: f64| -> Result<f64> {
        let mut moles = z.to_vec();
        moles[first] -= extent * record.cation_stoichiometry;
        moles[second] -= extent * record.anion_stoichiometry;
        let total: f64 = moles.iter().sum();
        let fraction = |position: usize| moles[position] / total;
        // **The coefficients are the brine's own**, re-solved at the composition the extent
        // leaves: `eos.pitzer_phase` over the same names and the same state.
        let x: Vec<f64> = (0..n).map(fraction).collect();
        let pitzer = pitzer_phase(components, T.value, &x)?;
        let result = scale_saturation_ratio(
            salt,
            fraction(first),
            fraction(second),
            fraction(water),
            pitzer.gamma[first],
            pitzer.gamma[second],
            fraction(water) * pitzer.gamma[water],
            hydrogen,
            T,
            P,
        )?;
        Ok(result.saturation_ratio)
    };

    let initial = ratio(0.0)?;
    // **Under saturation takes nothing**, and its ratio is the answer's.
    if initial <= 1.0 {
        return Ok(SaltPrecipitationResult {
            precipitated_moles: 0.0,
            initial_saturation_ratio: initial,
            final_saturation_ratio: initial,
            iterations: 0,
            extent_of_maximum: 0.0,
            warnings,
        });
    }

    // The extent at which an ion runs out, which is the bracket's upper end.
    let mut maximum = f64::INFINITY
        .min(z[first] / record.cation_stoichiometry)
        .min(z[second] / record.anion_stoichiometry);
    if record.water_stoichiometry > 0.0 {
        maximum = maximum.min(z[water] / record.water_stoichiometry);
    }
    if !maximum.is_finite() || maximum <= 0.0 {
        return Err(AzothError::out_of_range(
            "extent",
            maximum,
            format!("no finite positive precipitation extent is available for {salt}"),
        ));
    }

    let upper = maximum * (1.0 - BRACKET_INSET);
    let at_upper = ratio(upper)?;
    if at_upper.is_nan() || at_upper >= 1.0 {
        return Err(AzothError::SolverNotConverged {
            iterations: 1,
            residual: at_upper,
            tolerance: LOG10_TOLERANCE,
        });
    }

    let algorithm = crate::algorithm_of(spec)?;
    let mut low = 0.0;
    let mut high = upper;
    let mut iterations = 0;
    let mut extent = f64::NAN;
    for step in 1..=algorithm.max_iterations {
        iterations = step;
        extent = 0.5 * (low + high);
        let value = ratio(extent)?;
        if value.log10().abs() <= LOG10_TOLERANCE {
            break;
        }
        if value > 1.0 {
            low = extent;
        } else {
            high = extent;
        }
    }

    let final_ratio = ratio(extent)?;
    // **How much of the bracket the answer used.** One means an ion ran out - the brine had
    // no more of something - and less means the ratio crossed one first.
    let extent_of_maximum = extent / maximum;

    Ok(SaltPrecipitationResult {
        precipitated_moles: extent,
        initial_saturation_ratio: initial,
        final_saturation_ratio: final_ratio,
        iterations,
        extent_of_maximum,
        warnings,
    })
}
