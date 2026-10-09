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
//! # What is not here
//!
//! `OilAssayCharacterisation.apply` resolves a molar mass per cut, pulls standard components out
//! of the databank and adds a pseudo-component for each. `AssayCut` holds six optional fields per
//! cut, and a parallel vector cannot say "absent" per entry - only per vector - so the assay's
//! per-cut properties cannot be stated at this boundary at all. That is the `record_list` input
//! kind the plan costed as "ergonomics only" and deferred; for this id it is the difference
//! between the assay being expressible and not.

use azoth_core::units::{MassDensity, kilograms_per_cubic_meter};
use azoth_core::{AzothError, Result, apply_checks};

use crate::model_gen;
use crate::results::AssayMassFractionsResult;

/// `OilAssayCharacterisation.ASSAY_CLOSURE_TOLERANCE`, on the declared fractions' sum.
const ASSAY_CLOSURE_TOLERANCE: f64 = 1.0e-3;

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
/// * [`AzothError::InvalidInput`] if `density` is supplied and is not one entry per cut, if it is
///   needed and absent, if the declared fractions do not sum to one within `1e-3`, or if the row
///   normalises to nothing.
/// * [`AzothError::OutOfRange`] if `total_declared_fraction` falls outside its declared range.
///
/// # Example
/// ```
/// use azoth_core::units::kilograms_per_cubic_meter;
/// use azoth_characterization::{AssayBasis, assay_mass_fractions};
///
/// let resolved = assay_mass_fractions(
///     AssayBasis::Volume,
///     &[0.5, 0.5],
///     Some(&[
///         kilograms_per_cubic_meter(700.0),
///         kilograms_per_cubic_meter(900.0),
///     ]),
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
    declared_fraction: &[f64],
    density: Option<&[MassDensity]>,
) -> Result<AssayMassFractionsResult> {
    let spec = &model_gen::ASSAY_MASS_FRACTIONS_SPEC;
    let mut warnings = Vec::new();

    if declared_fraction.is_empty() {
        return Err(AzothError::invalid_input(
            "declared_fraction",
            "an assay of no cuts has no fractions to resolve",
        ));
    }
    if let Some(densities) = density
        && densities.len() != declared_fraction.len()
    {
        return Err(AzothError::invalid_input(
            "density",
            format!(
                "one density per cut is what the volume conversion needs: {} cuts against {} \
                 densities",
                declared_fraction.len(),
                densities.len()
            ),
        ));
    }
    if basis == AssayBasis::Volume && density.is_none() {
        return Err(AzothError::invalid_input(
            "density",
            "a volume-basis assay needs each cut's density to reach a mass fraction",
        ));
    }

    let total: f64 = declared_fraction.iter().sum();
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
            "declared_fraction",
            format!("the declared fractions sum to {total}, which is no assay"),
        ));
    }
    if (total - 1.0).abs() > ASSAY_CLOSURE_TOLERANCE {
        return Err(AzothError::invalid_input(
            "declared_fraction",
            format!(
                "assay fractions must sum to 1.0 within {ASSAY_CLOSURE_TOLERANCE}; supplied \
                 sum={total}"
            ),
        ));
    }

    let normalised: Vec<f64> = declared_fraction.iter().map(|f| f / total).collect();
    let mass_fraction = if basis == AssayBasis::Mass {
        normalised
    } else {
        // The densities are present: the guard above refused their absence.
        let densities = density.unwrap_or_default();
        let relative: Vec<f64> = normalised
            .iter()
            .zip(densities)
            .map(|(fraction, rho)| fraction * rho.value)
            .collect();
        let total_relative: f64 = relative.iter().sum();
        if !total_relative.is_finite() || total_relative <= 0.0 {
            return Err(AzothError::invalid_input(
                "density",
                "the volume-basis row has no mass to normalise: every cut's density is zero",
            ));
        }
        relative.iter().map(|mass| mass / total_relative).collect()
    };

    // The bulk density is the mass-weighted harmonic mean, and it is only answerable when every
    // cut carries a density - which a mass-basis assay need not state.
    let bulk_density = density.map(|densities| {
        let reciprocal: f64 = mass_fraction
            .iter()
            .zip(densities)
            .map(|(mass, rho)| mass / rho.value)
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
