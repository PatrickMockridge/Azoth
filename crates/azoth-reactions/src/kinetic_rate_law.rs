//! The kinetic rate law, from `ChemicalReaction.getRateFactor`.
//!
//! One method, two laws behind a selector, and **which one runs is decided by a field that the
//! class's own database-built reactions do not have**:
//!
//! ```text
//! LEGACY_TEMPERATURE_CORRELATION   2.576e9 * exp(-6024 / T) / 1000
//! REFERENCE_ARRHENIUS              rate * exp(-Ea / R * (1/T - 1/T_ref))
//! ```
//!
//! `getKineticRateLaw()` answers `LEGACY` when its field is **absent** - the class's comment is
//! explicit that an object serialised before the selector existed keeps the correlation - and
//! the measurement in the capture says every reaction `chemicalReactionInit` builds is in that
//! state. So the legacy law is what a fluid's reactions actually run, and the Arrhenius branch
//! is reached only by a caller who sets it.
//!
//! **The legacy law does not read the reaction's own parameters at all.** `2.576e9` and `6024`
//! are literals in the method, so every reaction of a fluid gets the same rate factor at a given
//! temperature - the capture's three reactions agree to the last digit at both temperatures - and
//! the stored `rateFactor` and `ACTENERGY` columns are never evaluated there. That is what makes
//! the Arrhenius branch the class documents a migration rather than a variant.
//!
//! The gas constant is [`crate::equilibrium_constant::GAS_CONSTANT`], which is
//! `ThermodynamicConstantsInterface.R` - the same `8.3144621` the equilibrium constant is built
//! with, and not the `8.314462` the RAND solver carries.

use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::equilibrium_constant::GAS_CONSTANT;

/// The legacy correlation's prefactor, from the method's own literal.
pub const LEGACY_PREFACTOR: f64 = 2.576e9;

/// The temperature it is written against, in K - `6024` as the method spells it.
pub const LEGACY_ACTIVATION_TEMPERATURE: f64 = 6024.0;

/// The divisor the result is scaled by, from the method's `/ 1000.0`.
pub const LEGACY_DIVISOR: f64 = 1000.0;

/// Which of the two laws is selected.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum KineticRateLaw {
    /// The historical correlation: one literal per temperature, the same for every reaction.
    LegacyTemperatureCorrelation,
    /// A reference rate, an activation energy and a reference temperature.
    ReferenceArrhenius,
}

impl KineticRateLaw {
    /// The law a selector names, with **the absent selector answering the legacy law** - which
    /// is the class's own fallback and the state its database-built reactions are in.
    ///
    /// # Errors
    /// [`AzothError::InvalidInput`] for a name that is neither.
    pub fn parse(name: &str) -> Result<Self> {
        match name.trim().to_lowercase().as_str() {
            "legacy" | "legacy_temperature_correlation" => Ok(Self::LegacyTemperatureCorrelation),
            "arrhenius" | "reference_arrhenius" => Ok(Self::ReferenceArrhenius),
            other => Err(AzothError::InvalidInput {
                field: "rate_law".to_string(),
                reason: format!("`{other}` is neither `legacy` nor `arrhenius`"),
            }),
        }
    }
}

/// A reaction's kinetic parameters: which law, and the three the Arrhenius one reads.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ReferenceKinetics {
    /// The selected law.
    pub law: KineticRateLaw,
    /// The rate at [`Self::reference_temperature`]. **Not read by the legacy law**, which the
    /// capture shows by returning the same number for three reactions whose stored rates agree.
    pub reference_rate: f64,
    /// The activation energy, in J/mol.
    pub activation_energy: f64,
    /// The temperature the reference rate is stated at, in K.
    pub reference_temperature: f64,
}

/// `getRateFactor`: the reaction's rate factor at a temperature.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a temperature that is not finite and positive - the class
/// throws there - and, on the Arrhenius branch, for parameters it refuses: a rate that is
/// negative or not finite, an activation energy that is not finite, or a reference temperature
/// that is not finite and positive.
pub fn rate_factor(kinetics: &ReferenceKinetics, temperature: f64) -> Result<f64> {
    if !temperature.is_finite() || temperature <= 0.0 {
        return Err(AzothError::InvalidInput {
            field: "temperature".to_string(),
            reason: format!("{temperature} is not a finite, positive temperature"),
        });
    }
    if kinetics.law == KineticRateLaw::LegacyTemperatureCorrelation {
        return Ok(
            LEGACY_PREFACTOR * (-LEGACY_ACTIVATION_TEMPERATURE / temperature).exp()
                / LEGACY_DIVISOR,
        );
    }
    if !kinetics.reference_rate.is_finite()
        || kinetics.reference_rate < 0.0
        || !kinetics.activation_energy.is_finite()
        || !kinetics.reference_temperature.is_finite()
        || kinetics.reference_temperature <= 0.0
    {
        return Err(AzothError::InvalidInput {
            field: "kinetics".to_string(),
            reason: format!(
                "the reference law needs a finite nonnegative rate, a finite energy and a \
                 positive temperature, and it has {}, {} and {}",
                kinetics.reference_rate, kinetics.activation_energy, kinetics.reference_temperature
            ),
        });
    }
    if kinetics.reference_rate == 0.0 {
        return Ok(0.0);
    }
    Ok(kinetics.reference_rate
        * (-kinetics.activation_energy / GAS_CONSTANT
            * (1.0 / temperature - 1.0 / kinetics.reference_temperature))
            .exp())
}

/// Result of `reactions.kinetic_rate_law`.
#[derive(Debug, Clone, PartialEq)]
pub struct KineticRateLawResult {
    /// The reaction's rate factor at `T`, by the selected law.
    pub rate_factor: f64,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for KineticRateLawResult {
    const CALC_ID: &'static str = "reactions.kinetic_rate_law";
    const FIELDS: &'static [&'static str] = &["rate_factor", "warnings"];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// `reactions.kinetic_rate_law`: the rate factor a reaction's own law answers at a state.
///
/// `law` is the selector's name, parsed by [`KineticRateLaw::parse`], and **the three
/// Arrhenius parameters are taken whether or not the legacy branch reads them** - NeqSim's
/// object carries them either way, and the legacy law ignores all three.
///
/// # Errors
/// [`AzothError::InvalidInput`] for a law this library does not carry, a temperature that is
/// not finite and positive, and on the reference branch an Arrhenius parameter the law
/// cannot use.
#[allow(non_snake_case)] // `T` is the symbol in the published equation
pub fn kinetic_rate_law(
    law: &str,
    T: f64,
    reference_rate: f64,
    activation_energy: f64,
    reference_temperature: f64,
) -> Result<KineticRateLawResult> {
    let spec = &crate::model_gen::KINETIC_RATE_LAW_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            "reference_temperature" => Some(reference_temperature),
            "reference_rate" => Some(reference_rate),
            _ => None,
        },
        &mut warnings,
    )?;

    let kinetics = ReferenceKinetics {
        law: KineticRateLaw::parse(law)?,
        reference_rate,
        activation_energy,
        reference_temperature,
    };
    Ok(KineticRateLawResult {
        rate_factor: rate_factor(&kinetics, T)?,
        warnings,
    })
}
