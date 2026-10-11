//! `characterization.assay_mass_fractions` - an oil assay's declared fractions resolved to mass.
//!
//! Spec: `specs/models/characterization/assay_mass_fractions.toml`. Oracle:
//! `validation/neqsim/captures/oil_assay_probe.tsv`.
//!
//! # One basis, one closure, one answer
//!
//! An assay declares each cut's yield either as a mass fraction or as a liquid-volume fraction,
//! never as both and never as neither, and the whole assay shares one basis. The row is
//! normalised, and the closure - `|sum - 1| <= 1e-3` - is a **hard error**: a tenth of a percent
//! of closure is a broken assay and not a rounding artefact.
//!
//! A mass basis then resolves to itself. A volume basis multiplies each normalised volume
//! fraction by its cut's specific gravity and renormalises, which is what makes a light end lose
//! share to a heavy one.
//!
//! # The cut is a record
//!
//! NeqSim's `AssayCut` holds twelve nullable Doubles, and the one that matters here is the
//! density: it is optional *per cut*, which a parallel vector cannot say - a vector can be absent
//! per vector and not per entry. A volume basis needs every cut's, so the first that states none
//! is refused by position; a mass basis reads them for the bulk density alone, so a row may state
//! some and not others.
//!
//! # What is not here
//!
//! `OilAssayCharacterisation.apply` resolves a molar mass per cut, pulls standard components out
//! of the databank and adds a pseudo-component for each. The row this id takes is the two fields
//! it reads; the API gravity, the boiling point or Watson factor, the molar mass, the sulfur and
//! the nitrogen are still owed, and `apply` is the id that will want them.

use azoth_core::units::{MassDensity, kilograms_per_cubic_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::AssayMassFractionsResult;

/// `OilAssayCharacterisation.ASSAY_CLOSURE_TOLERANCE`, on the declared fractions' sum.
const ASSAY_CLOSURE_TOLERANCE: f64 = 1.0e-3;

/// One cut of an assay, as this id takes it.
///
/// **A record rather than a vector entry**, and the density is why: `Option` per *cut* has no
/// aligned-vector form. NeqSim's `AssayCut` holds ten more nullable fields beside these two; the
/// row grows as the ids that read them are ported.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct AssayCut {
    /// The cut's declared yield, on the assay's basis.
    pub declared_fraction: f64,
    /// The cut's own normal liquid density, absent when the assay does not state one.
    pub density: Option<MassDensity>,
}

/// Whether a cut's declared yield is a fraction of the assay's mass or of its liquid volume.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Default)]
pub enum AssayBasis {
    /// A mass fraction, which resolves to itself.
    #[default]
    Mass,
    /// A liquid-volume fraction, which needs a density per cut to become a mass fraction.
    Volume,
}

impl std::str::FromStr for AssayBasis {
    type Err = std::convert::Infallible;

    /// An unrecognised name is the mass basis; the spec's vocabulary check is what should refuse.
    fn from_str(text: &str) -> std::result::Result<Self, Self::Err> {
        Ok(match text.trim() {
            "volume" => Self::Volume,
            _ => Self::Mass,
        })
    }
}

/// An oil assay's declared yields, resolved to a mass basis.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the assay has no cuts, if a volume basis reaches a cut that
///   states no density, if the declared fractions do not sum to one within `1e-3`, or if the row
///   normalises to nothing.
/// * [`AzothError::OutOfRange`] if `total_declared_fraction` falls outside its declared range.
///
/// # Example
/// ```
/// use azoth_core::units::kilograms_per_cubic_meter;
/// use azoth_characterization::{AssayBasis, AssayCut, assay_mass_fractions};
///
/// let resolved = assay_mass_fractions(
///     AssayBasis::Volume,
///     &[
///         AssayCut {
///             declared_fraction: 0.5,
///             density: Some(kilograms_per_cubic_meter(700.0)),
///         },
///         AssayCut {
///             declared_fraction: 0.5,
///             density: Some(kilograms_per_cubic_meter(900.0)),
///         },
///     ],
/// )?;
/// assert!((resolved.mass_fraction[0] - 0.4375).abs() < 1e-15);
/// // `1 / (0.4375/700 + 0.5625/900)`, the mass-weighted harmonic mean in kg/m3.
/// assert!(
///     (resolved.bulk_density.expect("densities were given").value - 800.0).abs() < 1e-9
/// );
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn assay_mass_fractions(
    basis: AssayBasis,
    cuts: &[AssayCut],
) -> Result<AssayMassFractionsResult> {
    let spec = &model_gen::ASSAY_MASS_FRACTIONS_SPEC;
    let mut warnings = Vec::new();

    if cuts.is_empty() {
        return Err(AzothError::invalid_input(
            "cuts",
            "an assay of no cuts has no fractions to resolve",
        ));
    }

    let total: f64 = cuts.iter().map(|cut| cut.declared_fraction).sum();
    // The spec's one bound is on this sum, which is not an input - so it is a derived check and
    // the closure resolves it from what was computed rather than from what was handed over.
    apply_checks(
        spec.derived_checks(),
        |quantity| match quantity {
            "total_declared_fraction" => Some(total),
            _ => None,
        },
        &mut warnings,
    )?;
    if !total.is_finite() {
        return Err(AzothError::invalid_input(
            "cuts",
            format!("the declared fractions sum to {total}, which is no assay"),
        ));
    }
    if (total - 1.0).abs() > ASSAY_CLOSURE_TOLERANCE {
        return Err(AzothError::invalid_input(
            "cuts",
            format!(
                "assay fractions must sum to 1.0 within {ASSAY_CLOSURE_TOLERANCE}; supplied \
                 sum={total}"
            ),
        ));
    }

    let normalised: Vec<f64> = cuts.iter().map(|cut| cut.declared_fraction / total).collect();
    let mass_fraction = if basis == AssayBasis::Mass {
        normalised
    } else {
        // **Every cut needs a density, and the refusal names which one does not.** A volume
        // fraction has no mass to weight by without one, so a cut that states none is a hole in
        // the conversion rather than a cut without a bulk-density contribution.
        let mut relative = Vec::with_capacity(cuts.len());
        for (index, (fraction, cut)) in normalised.iter().zip(cuts).enumerate() {
            let Some(rho) = cut.density else {
                return Err(AzothError::invalid_input(
                    format!("cuts[{index}].density"),
                    "a volume-basis assay needs every cut's density to reach a mass fraction",
                ));
            };
            relative.push(fraction * rho.value);
        }
        let total_relative: f64 = relative.iter().sum();
        if !total_relative.is_finite() || total_relative <= 0.0 {
            return Err(AzothError::invalid_input(
                "cuts",
                "the volume-basis row has no mass to normalise: every cut's density is zero",
            ));
        }
        relative.iter().map(|mass| mass / total_relative).collect()
    };

    // The bulk density is the mass-weighted harmonic mean, and it is answerable only when **every**
    // cut carries a density. Collecting the options is that test: the `Option` is `Some` exactly
    // when nothing was absent, so a mass-basis assay that states a density for some cuts and not
    // others loses the bulk density rather than averaging over the ones it has.
    let densities: Option<Vec<f64>> = cuts
        .iter()
        .map(|cut| cut.density.map(|rho| rho.value))
        .collect();
    let bulk_density = densities.map(|values| {
        let reciprocal: f64 = mass_fraction
            .iter()
            .zip(&values)
            .map(|(mass, rho)| mass / rho)
            .sum();
        kilograms_per_cubic_meter(1.0 / reciprocal)
    });
    if let Some(bulk) = bulk_density
        && (!bulk.value.is_finite() || bulk.value <= 0.0)
    {
        return Err(AzothError::invalid_input(
            "bulk_density",
            format!(
                "the bulk density of these fractions and densities is {}",
                bulk.value
            ),
        ));
    }

    Ok(AssayMassFractionsResult {
        mass_fraction,
        total_declared_fraction: total,
        bulk_density,
        warnings,
    })
}
