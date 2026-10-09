//! `characterization.characterise_plus_fraction` - a C7+ end characterised end to end.
//!
//! Spec: `specs/models/characterization/characterise_plus_fraction.toml`. Oracle:
//! `validation/neqsim/captures/plus_fraction_probe.tsv`.
//!
//! # The front door, and nothing of its own
//!
//! `Characterise.characterisePlusFraction` is four statements: initialise, find the plus fraction,
//! **replace the model if the plus fraction is heavier than the model can take**, and run the
//! lumping only if the split returned true. Everything between is
//! [`pedersen_plus_split`](crate::pedersen_plus_split),
//! [`whitson_gamma_split`](crate::whitson_gamma_split) and [`lumping`](crate::lumping); this
//! module is the wiring and the one rule that is the class's own.
//!
//! # The replacement is silent and this is the only thing that reports it
//!
//! `getMPlus() > getMaxPlusMolarMass()` is checked **before** the model is used, so a request for
//! Whitson Gamma above the threshold is swept away with the rest - and Whitson's own threshold is
//! `PedersenPlusModel`'s `0.605`, which its constructor inherits, not the 2.10 its sibling name
//! suggests. [`CharacterisePlusFractionResult::selected_model`] is what actually ran.
//!
//! # A false split is refused rather than skipped
//!
//! The class takes `characterizePlusFraction`'s boolean and, when it is false, leaves the fluid
//! with its single plus row and no lumps at all. An id that returns a table has no fluid to leave
//! unchanged, so the split's own refusal is propagated. That is a stated divergence.

use azoth_core::units::{MassDensity, MolarMass};
use azoth_core::{Result, apply_checks};

use crate::lumping::lumping;
use crate::model_gen;
use crate::pedersen_plus_split::pedersen_plus_split;
use crate::results::CharacterisePlusFractionResult;
use crate::whitson_gamma_split::{WhitsonDensityModel, whitson_gamma_split};

/// The first carbon number the gamma model's cut table begins at, which its override never moves.
const GAMMA_FIRST_CARBON_NUMBER: usize = 1;

/// Which plus model a caller asks for.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum PlusModel {
    /// `PedersenPlusModel`: 80 cuts from C20 and a maximum molar mass of 0.605 kg/mol.
    #[default]
    Pedersen,
    /// `PedersenHeavyOilPlusModel`: the same solve over 200 carbon numbers, and 2.10 kg/mol.
    PedersenHeavyOil,
    /// `WhitsonGammaModel`, which inherits the 0.605 threshold and the 80-cut range.
    WhitsonGamma,
}

impl PlusModel {
    /// `getMaxPlusMolarMass()`, which is the outer object's own field.
    #[must_use]
    pub const fn maximum_molar_mass(self) -> f64 {
        match self {
            Self::PedersenHeavyOil => 2.10,
            Self::Pedersen | Self::WhitsonGamma => 0.605,
        }
    }

    /// The last carbon number the model's split ends at, which its constructor sets.
    #[must_use]
    pub const fn last_carbon_number(self) -> usize {
        match self {
            Self::PedersenHeavyOil => 200,
            Self::Pedersen | Self::WhitsonGamma => 80,
        }
    }

    /// The spec's own spelling, which is also the spelling `selected_model` reports.
    #[must_use]
    pub fn as_str(self) -> &'static str {
        match self {
            Self::Pedersen => "pedersen",
            Self::PedersenHeavyOil => "pedersen_heavy_oil",
            Self::WhitsonGamma => "whitson_gamma",
        }
    }
}

impl std::str::FromStr for PlusModel {
    type Err = std::convert::Infallible;

    /// An unrecognised name is the class's own default, which `getModel("")` returns; the spec's
    /// vocabulary check is what should refuse.
    fn from_str(text: &str) -> std::result::Result<Self, Self::Err> {
        Ok(match text.trim() {
            "pedersen_heavy_oil" => Self::PedersenHeavyOil,
            "whitson_gamma" => Self::WhitsonGamma,
            _ => Self::Pedersen,
        })
    }
}

/// A C7+ end characterised end to end: model, split and lumps.
///
/// # Errors
/// * [`AzothError::InvalidInput`] or [`AzothError::OutOfRange`] from whichever of the three ids
///   refuses - including a plus fraction the split declines, which this raises rather than
///   skipping as the class does.
///
/// # Example
/// ```
/// use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
/// use azoth_characterization::{PlusModel, characterise_plus_fraction};
///
/// let characterised = characterise_plus_fraction(
///     kilograms_per_mole(0.4),
///     kilograms_per_cubic_meter(850.0),
///     0.1,
///     20,
///     Some(PlusModel::Pedersen),
///     None,
/// )?;
/// assert_eq!(characterised.fraction_of_heavy_end.len(), 7);
/// assert!((characterised.fraction_of_heavy_end.iter().sum::<f64>() - 1.0).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn characterise_plus_fraction(
    molar_mass: MolarMass,
    density: MassDensity,
    mole_fraction: f64,
    first_carbon_number: usize,
    plus_model: Option<PlusModel>,
    number_of_lumps: Option<usize>,
) -> Result<CharacterisePlusFractionResult> {
    let spec = &model_gen::CHARACTERISE_PLUS_FRACTION_SPEC;
    let mut warnings = Vec::new();
    apply_checks(spec.input_checks(), |_quantity| None, &mut warnings)?;

    let requested = plus_model.unwrap_or_default();
    // The class replaces the model before using it, so the request is discarded rather than
    // reported as an error - and `selected_model` is the only thing that says so.
    let selected = if molar_mass.value > requested.maximum_molar_mass() {
        PlusModel::PedersenHeavyOil
    } else {
        requested
    };
    let last_carbon_number = selected.last_carbon_number();

    // The two splits answer with one common shape: a table of cuts.
    let (cut_z, cut_molar_mass, cut_density) = match selected {
        PlusModel::WhitsonGamma => {
            // **The first carbon number is ignored here.** The gamma model's override never reads
            // it, so its cuts always begin at the outer object's own value of one.
            let split = whitson_gamma_split(
                molar_mass,
                density,
                mole_fraction,
                GAMMA_FIRST_CARBON_NUMBER,
                last_carbon_number,
                None,
                None,
                WhitsonDensityModel::Uop,
                false,
            )?;
            warnings.extend(split.warnings);
            (split.cut_z, split.cut_molar_mass, split.cut_density)
        }
        PlusModel::Pedersen | PlusModel::PedersenHeavyOil => {
            let split = pedersen_plus_split(
                molar_mass,
                density,
                mole_fraction,
                first_carbon_number,
                last_carbon_number,
            )?;
            warnings.extend(split.warnings);
            (split.cut_z, split.cut_molar_mass, split.cut_density)
        }
    };

    // **The lumping is handed the plus fraction's own two numbers**, not the cut table's sums: the
    // class reads them off the fluid and the difference moves a partition boundary on the heavy
    // rows.
    let grouped = lumping(
        molar_mass,
        mole_fraction,
        &cut_z,
        &cut_molar_mass,
        &cut_density,
        number_of_lumps,
    )?;
    warnings.extend(grouped.warnings);

    Ok(CharacterisePlusFractionResult {
        selected_model: selected,
        fraction_of_heavy_end: grouped.fraction_of_heavy_end,
        lump_molar_mass: grouped.lump_molar_mass,
        lump_density: grouped.lump_density,
        warnings,
    })
}
