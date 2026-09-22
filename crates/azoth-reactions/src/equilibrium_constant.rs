//! `reactions.equilibrium_constant` - one reaction's `K` and its temperature derivative.
//!
//! ```text
//! ln K        = K1 + K2/T + K3*ln T + K4*T
//! d(ln K)/dT  = -K1/T**2 + K2/T + K3
//! dH          = d(ln K)/dT * R * T**2
//! ```
//!
//! Spec: `specs/calcs/reactions/equilibrium_constant.toml`, which carries the sources,
//! what the four coefficients are fitted to, and why the three reaction tables are not
//! interchangeable.
//!
//! **The four coefficients are a van't Hoff term, a heat-capacity term and a linear
//! term**, and they come from a named row rather than from the caller: `K1`, `K2` and
//! `K3` are dimensionless, K and 1/K respectively once divided by the powers of `T`
//! they are divided by, so they never cross the API and no coefficient unit is needed.
//! What crosses is the reaction's name and the temperature.
//!
//! The reference temperature in the row is **not read by the correlation**. It is the
//! temperature the fit is anchored at and the reference the rate law uses; `ln K` at
//! `T = Tref` is not any particular value.

use azoth_core::units::{MolarEnergy, ThermodynamicTemperature, joules_per_mole};
use azoth_core::{AzothError, CalcResult, Result, Warning, apply_checks};

use crate::databank::{ReactionDataSource, reaction};
use crate::spec_gen;

/// The gas constant `ThermodynamicConstantsInterface` fixes, in J/(mol*K).
///
/// **It is the SI value and the table's `ACTENERGY` is not in J/mol** - see the spec's
/// assumptions. Nothing in the correlation below reads the activation energy, so the two
/// never meet on this path.
pub const GAS_CONSTANT: f64 = 8.3144621;

/// Result of `reactions.equilibrium_constant`.
#[derive(Debug, Clone, PartialEq)]
pub struct EquilibriumConstantResult {
    /// `ln K` at the caller's temperature, from the four fitted coefficients.
    pub ln_k: f64,
    /// `K`, the exponential of `ln_k`. Dimensionless, because the correlation's
    /// standard state is what makes it so.
    pub k: f64,
    /// `d(ln K)/dT`, in 1/K. Van't Hoff's derivative, and the quantity the reaction heat
    /// is that times `R T**2`.
    pub ln_k_derivative: f64,
    /// `dH = d(ln K)/dT * R * T**2`, in J/mol.
    ///
    /// **A sign is a statement about the reaction as the table writes it**, which for
    /// `CO2water` means the forward direction is `CO2 + H2O -> H2CO3`. A negative value
    /// is exothermic in that direction and says nothing about the reverse.
    pub reaction_heat: MolarEnergy,
    /// The row's own citation, so which fit answered is visible in the output rather
    /// than only in the table.
    pub reference: String,
    /// Caveats.
    pub warnings: Vec<Warning>,
}

impl CalcResult for EquilibriumConstantResult {
    const CALC_ID: &'static str = "reactions.equilibrium_constant";
    const FIELDS: &'static [&'static str] = &[
        "ln_k",
        "k",
        "ln_k_derivative",
        "reaction_heat",
        "reference",
        "warnings",
    ];

    fn warnings(&self) -> &[Warning] {
        &self.warnings
    }
}

/// One reaction's `K`, its `d(ln K)/dT` and its heat of reaction at a temperature.
///
/// `source` selects which of the three reaction tables the row is read from, and the
/// three are on different standard states rather than being a fallback chain.
///
/// # Errors
/// [`AzothError::PropertyUnavailable`] where the source carries no row by that name -
/// a name the data could not answer, which is a different failure from a value the
/// caller supplied - and [`AzothError::OutOfRange`] for a temperature the spec bounds.
///
/// # Example
/// ```
/// use azoth_reactions::databank::ReactionDataSource;
/// use azoth_reactions::equilibrium_constant::equilibrium_constant;
/// use azoth_core::units::kelvins;
///
/// // The standard source's `CO2water` row at 298.15 K, from
/// // `validation/neqsim/captures/reaction_probe.tsv`.
/// let r = equilibrium_constant(
///     ReactionDataSource::Standard,
///     "CO2water",
///     kelvins(298.15),
/// )?;
/// assert!((r.ln_k + 14.633690407708798).abs() < 1e-12);
/// assert!((r.k - 4.412340363006617e-7).abs() < 1e-18);
/// assert!((r.reaction_heat.value - 9199.128566216192).abs() < 1e-6);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn equilibrium_constant(
    source: ReactionDataSource,
    reaction_name: &str,
    temperature: ThermodynamicTemperature,
) -> Result<EquilibriumConstantResult> {
    let spec = &spec_gen::EQUILIBRIUM_CONSTANT_SPEC;
    let mut warnings = Vec::new();

    let t = temperature.value;

    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(t),
            _ => None,
        },
        &mut warnings,
    )?;

    let row = reaction(source, reaction_name)?.ok_or_else(|| AzothError::PropertyUnavailable {
        fluid: reaction_name.to_string(),
        property: "reaction".to_string(),
        reason: format!(
            "the `{}` source carries no row by that name",
            source.identifier()
        ),
    })?;

    let [k1, k2, k3, k4] = row.coefficients;

    let ln_k = k1 + k2 / t + k3 * t.ln() + k4 * t;
    let k = ln_k.exp();
    let ln_k_derivative = -k2 / (t * t) + k3 / t + k4;
    let reaction_heat = joules_per_mole(ln_k_derivative * t * t * GAS_CONSTANT);

    apply_checks(
        spec.derived_checks(),
        |name| match name {
            "ln_k" => Some(ln_k),
            "k" => Some(k),
            _ => None,
        },
        &mut warnings,
    )?;

    Ok(EquilibriumConstantResult {
        ln_k,
        k,
        ln_k_derivative,
        reaction_heat,
        reference: row.reference,
        warnings,
    })
}

#[cfg(test)]
mod tests {
    use super::*;
    use azoth_core::units::kelvins;

    /// The van't Hoff relation itself: the heat is the derivative times `R T**2`, and
    /// this is that identity rather than a second implementation of it.
    #[test]
    fn the_heat_is_the_derivative_times_the_gas_constant_and_the_temperature() {
        let result =
            equilibrium_constant(ReactionDataSource::Standard, "CO2water", kelvins(298.15))
                .expect("the row exists");
        let from_identity = result.ln_k_derivative * 298.15 * 298.15 * GAS_CONSTANT;
        assert!((result.reaction_heat.value - from_identity).abs() < 1e-9);
    }

    /// A name the source does not carry is a data answer, not an input error.
    #[test]
    fn an_unknown_reaction_is_unavailable_rather_than_invalid() {
        let error = equilibrium_constant(ReactionDataSource::Standard, "no-such", kelvins(298.15))
            .expect_err("no row");
        assert!(
            matches!(error, AzothError::PropertyUnavailable { .. }),
            "{error:?}"
        );
    }

    /// Which source is read is part of the answer - but **not at the reference
    /// temperature**, where the two fits are anchored and cross within 0.01. The state
    /// that shows it is away from there, which is where `K3` and `K4` earn their place.
    #[test]
    fn the_sources_diverge_away_from_the_reference_temperature() {
        let at = |source, t| {
            equilibrium_constant(source, "CO2water", kelvins(t))
                .expect("both sources carry the row")
                .ln_k
        };

        let near =
            at(ReactionDataSource::Standard, 298.15) - at(ReactionDataSource::Pitzer, 298.15);
        assert!(
            near.abs() < 0.01,
            "the two fit the same reaction and nearly agree at the anchor: {near}"
        );

        let far = at(ReactionDataSource::Standard, 423.15) - at(ReactionDataSource::Pitzer, 423.15);
        assert!(
            far.abs() > 0.2,
            "125 K above the anchor they part company, and a port that shared one row \
             between the sources would be wrong by a factor of 1.24 in K: {far}"
        );
    }
}
