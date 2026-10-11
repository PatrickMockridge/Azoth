//! `characterization.set_heavy_tbp_to_plus` - the TBP row a fluid's plus end begins at.
//!
//! Spec: `specs/models/characterization/set_heavy_tbp_to_plus.toml`. Oracle:
//! `validation/neqsim/captures/set_heavy_tbp_probe.tsv`.
//!
//! # The carbon number is in the name
//!
//! `PlusCharacterize.setHeavyTBPtoPlus` picks the row whose **name** carries the largest carbon
//! number, and it reads it out of the string with a positional slice:
//!
//! ```text
//! name.substring(3, 4).equals("_") ? name.substring(1, 3) : name.substring(1, 2)
//! ```
//!
//! So `C7_PC` reads `7`, `C10_PC` reads `10`, and a name whose fourth character is not the
// underscore reads a single character that will not parse - which the class catches along with
// everything else in the loop and **skips**. That is the one place this refuses where the class
// carries on: a plus end chosen from whatever happened to parse is a plausible answer and not the
// caller's.
//!
//! # The row below is the row before
//!
//! The class reads `component(compNumber - 1)` for `densLastTBP`, so it is the winner's neighbour
//! by **index** and not the row with the next carbon number. A slate whose names are out of order
//! takes its neighbour by position, and a winner at index zero reads index -1.

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::SetHeavyTbpToPlusResult;

/// Where the class slices for the underscore that says the carbon number is two digits.
const SUFFIX_MARKER: usize = 3;

/// The class's `0.03`, which it adds to a density in g/cm3. In the kg/m3 this declares, 30.
const INTERCEPT_STEP: f64 = 30.0;

/// `PlusCharacterize.setHeavyTBPtoPlus`'s parse: the carbon number a component's name carries.
///
/// `None` for a name the class's own slice would throw or fail to parse - which is what the
/// class's `catch` swallows, and what this refuses.
fn carbon_number(name: &str) -> Option<u32> {
    let chars: Vec<char> = name.chars().collect();
    // The class's first act is `name.substring(3, 4)`, which throws on a name of fewer than
    // four characters; the handler that catches a number that will not parse catches that too.
    let marker = *chars.get(SUFFIX_MARKER)?;
    let digits = if marker == '_' {
        chars.get(1..3)?
    } else {
        chars.get(1..2)?
    };
    digits.iter().collect::<String>().parse().ok()
}

/// The TBP row a fluid's plus end begins at, chosen by the carbon number in its name.
///
/// `names` are the fluid's TBP rows in the order it carries them, and the three vectors beside
/// them are those rows' data - each aligned with `names` entry for entry.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the vectors are not all the same length, if the slate is
///   empty, if a name carries no carbon number, or if the winner is the first row - the class
///   reads the row below it, and there is none.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::set_heavy_tbp_to_plus;
///
/// // The winner is the second row by name and not the last by position.
/// let result = set_heavy_tbp_to_plus(
///     &["C7_PC", "C20_PC", "C9_PC"],
///     &[
///         kilograms_per_mole(0.1),
///         kilograms_per_mole(0.35),
///         kilograms_per_mole(0.14),
///     ],
///     &[0.25, 0.25, 0.5],
///     &[
///         kilograms_per_cubic_meter(720.0),
///         kilograms_per_cubic_meter(900.0),
///         kilograms_per_cubic_meter(750.0),
///     ],
/// )?;
/// assert_eq!(result.plus_index, 1.0);
/// assert_eq!(result.plus_name, "C20_PC");
/// assert!((result.last_tbp_density.value - 720.0).abs() < 1e-12);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn set_heavy_tbp_to_plus(
    names: &[&str],
    molar_mass: &[MolarMass],
    mole_fraction: &[f64],
    normal_liquid_density: &[MassDensity],
) -> Result<SetHeavyTbpToPlusResult> {
    let spec = &model_gen::SET_HEAVY_TBP_TO_PLUS_SPEC;
    let mut warnings = Vec::new();

    apply_checks(
        spec.input_checks(),
        // The spec declares no range on these four columns: a bound on a vector is a bound on
        // each of its entries, and `apply_checks` takes one value per quantity. The length
        // agreement below is the check a slate of four aligned vectors actually needs.
        |_quantity| None,
        &mut warnings,
    )?;

    let count = names.len();
    if count == 0 {
        return Err(AzothError::invalid_input(
            "names",
            "a fluid of no TBP rows has no row for a plus end to begin at",
        ));
    }
    if molar_mass.len() != count
        || mole_fraction.len() != count
        || normal_liquid_density.len() != count
    {
        return Err(AzothError::invalid_input(
            "names",
            format!(
                "a fluid's TBP slate is four aligned vectors: {count} names against {} molar \
                 masses, {} mole fractions and {} normal liquid densities",
                molar_mass.len(),
                mole_fraction.len(),
                normal_liquid_density.len(),
            ),
        ));
    }

    // `plusCompNumber < firstPlusNumber`, strictly: a tie keeps the row that came first, which is
    // the class's own comparison against its running maximum.
    let mut winner: Option<(u32, usize)> = None;
    for (index, name) in names.iter().enumerate() {
        let Some(number) = carbon_number(name) else {
            return Err(AzothError::invalid_input(
                format!("names[{index}]"),
                format!("`{name}` carries no carbon number, so the walk has no place to put it"),
            ));
        };
        if winner.is_none_or(|(best, _)| best < number) {
            winner = Some((number, index));
        }
    }
    // The slate is not empty and every name parsed, so a winner exists.
    let (_, index) = winner.unwrap_or((0, 0));

    let Some(below) = index.checked_sub(1) else {
        return Err(AzothError::invalid_input(
            "names",
            format!(
                "the plus end is the first row, `{}`, and the row below it is what the class reads \
                 for `densLastTBP`",
                names[index],
            ),
        ));
    };

    Ok(SetHeavyTbpToPlusResult {
        plus_index: index as f64,
        plus_name: names[index].to_string(),
        plus_molar_mass: molar_mass[index],
        plus_mole_fraction: mole_fraction[index],
        plus_density: normal_liquid_density[index],
        last_tbp_density: normal_liquid_density[below],
        density_intercept: kilograms_per_cubic_meter(
            normal_liquid_density[below].value + INTERCEPT_STEP,
        ),
        warnings,
    })
}
