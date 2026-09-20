//! The molality-scale surface every electrolyte phase needs.
//!
//! An activity-coefficient phase works in mole fractions, and every electrolyte model
//! works in **molality**: `PhasePitzer`'s ionic strength is `1/2 sum m_i z_i^2`, its
//! osmotic coefficient is a sum over `m_i`, and `ComponentGePitzer.getMolality` divides
//! the mole count by the solvent's mass. This module is that change of scale, and it is
//! here rather than in each model because the arithmetic is the same for all of them.
//!
//! # What is *not* here
//!
//! **The osmotic coefficient itself.** `Pitzer`'s, `DesmukhMather`'s and
//! `KentEisenberg`'s are each a different formula over the same molalities, so only the
//! relation from a coefficient to a water activity is shared - [`ln_water_activity`].
//!
//! # The molality NeqSim means
//!
//! `Component.getMolality` is `getMolarity / (density / 1e3)`, which is moles per
//! kilogram of *solution*, and its own comment says `// return mol/kg`. **The two
//! electrolyte models that read a molality override it**: `ComponentGePitzer` and
//! `ComponentDesmukhMather` both return `n_i / getSolventWeight()`, which is moles per
//! kilogram of *water* and is the textbook quantity. This module implements the
//! override, because it is the one the models that need it use.

use azoth_core::{AzothError, Result};

/// The substance a phase's molality is measured against.
///
/// NeqSim matches the component's *name*, case-sensitively, in
/// `PhasePitzer.getSolventWeight`. Matched without regard to case here, which is the
/// rule the rest of this library resolves names by.
const SOLVENT: &str = "water";

/// One mole of mixture's worth of the molality-scale composition.
///
/// Per mole of mixture rather than per phase, because that is what the mole fractions
/// determine: the total mole count cancels out of `n_i / m_solvent`, so a phase's
/// molalities are a function of its composition and not of its size.
#[derive(Debug, Clone, PartialEq)]
pub struct ElectrolyteComposition {
    /// Each component's molality `n_i / m_solvent`, in mol/kg.
    ///
    /// **The solvent's own entry is its reciprocal molar mass**, about `55.5 mol/kg` for
    /// water: `n_w / m_w` where `m_w = n_w M_w`. NeqSim computes the same, and a model
    /// that summed this column as "the solute" would include the solvent.
    pub molality: Vec<f64>,
    /// The mass of solvent per mole of mixture, in kg.
    pub solvent_mass: f64,
    /// `I = 1/2 sum m_i z_i^2`, in mol/kg.
    pub ionic_strength: f64,
}

/// The mass of solvent per mole of mixture, in kg: `PhasePitzer.getSolventWeight`.
///
/// Summed over every component named [`SOLVENT`], which is NeqSim's rule - a mixture with
/// two water-named components weights both.
#[must_use]
pub fn solvent_mass(names: &[&str], x: &[f64], molar_mass: &[f64]) -> f64 {
    names
        .iter()
        .zip(x)
        .zip(molar_mass)
        .filter(|((name, _), _)| name.eq_ignore_ascii_case(SOLVENT))
        .map(|((_, &x), &mass)| x * mass)
        .sum()
}

/// `I = 1/2 sum m_i z_i^2`, in mol/kg: `PhasePitzer.getIonicStrength`.
///
/// Every component contributes, and a neutral one contributes zero because its charge is
/// zero - so this needs no test of which components are ions.
#[must_use]
pub fn ionic_strength(molality: &[f64], charge: &[f64]) -> f64 {
    0.5 * molality
        .iter()
        .zip(charge)
        .map(|(&m, &z)| m * z * z)
        .sum::<f64>()
}

/// The molality-scale composition of a phase, from its mole fractions.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the vectors disagree in length, or if the mixture
///   carries no solvent.
///
/// **The refusal is a divergence from NeqSim, and a deliberate one.**
/// `ComponentGePitzer.getMolality` returns `0.0` when the phase has no solvent weight,
/// so every molality is zero, the ionic strength is zero, and Pitzer evaluates as a
/// solution of nothing - a finite answer with no symptom. A phase with no water is not a
/// phase this surface describes, so it is refused.
pub fn composition(
    names: &[&str],
    x: &[f64],
    molar_mass: &[f64],
    charge: &[f64],
) -> Result<ElectrolyteComposition> {
    let n = names.len();
    if x.len() != n || molar_mass.len() != n || charge.len() != n {
        return Err(AzothError::invalid_input(
            "components",
            format!(
                "an electrolyte composition needs one entry per component: {} names, {} mole \
                 fractions, {} molar masses and {} charges",
                n,
                x.len(),
                molar_mass.len(),
                charge.len()
            ),
        ));
    }

    let mass = solvent_mass(names, x, molar_mass);
    if !mass.is_finite() || mass <= 0.0 {
        return Err(AzothError::invalid_input(
            "components",
            "the mixture carries no solvent, so every molality would be zero and the ionic \
             strength with it. NeqSim returns 0.0 here and evaluates the model as a solution \
             of nothing; a phase without water is not one this surface describes"
                .to_string(),
        ));
    }

    let molality: Vec<f64> = x.iter().map(|&xi| xi / mass).collect();
    Ok(ElectrolyteComposition {
        ionic_strength: ionic_strength(&molality, charge),
        molality,
        solvent_mass: mass,
    })
}

/// `ln(a_w) = -phi M_w sum_i m_i`, the relation from an osmotic coefficient to a water
/// activity, with `M_w` in kg/mol.
///
/// NeqSim's own comment (`ComponentGePitzer`, above its `lnaw`): `// Water activity:
/// ln(a_w) = -phi * M_w * sumM / 1000`, where the `1000` converts its `18.015 g/mol` to
/// kg/mol. The two agree - the databank's water molar mass is `0.018015 kg/mol` exactly -
/// so this takes the molar mass rather than restating the constant, and a table whose
/// water row moved would move this with it.
///
/// `sum_molalities` is the sum the model's own `phi` was built over, which is **not** the
/// solvent: `PhasePitzer` sums the charged components and the neutral solutes it has
/// interactions for, and returns `phi = 1` when that sum is below `1e-12`.
#[must_use]
pub fn ln_water_activity(
    osmotic_coefficient: f64,
    sum_molalities: f64,
    water_molar_mass: f64,
) -> f64 {
    -osmotic_coefficient * water_molar_mass * sum_molalities
}

/// `gamma_w = a_w / x_w`: NeqSim's conversion from the water activity it computes to the
/// mole-fraction activity coefficient the fugacity kernel takes.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if the water mole fraction is not positive, where the
///   ratio is not a number. NeqSim guards the same case and returns `gamma = 1`; refusing
///   is the house rule for a state the expression is not defined at.
pub fn water_activity_coefficient(ln_a_w: f64, x_water: f64) -> Result<f64> {
    if !x_water.is_finite() || x_water <= 0.0 {
        return Err(AzothError::invalid_input(
            "components",
            "the water mole fraction is not positive, so `a_w / x_w` is not a number. \
             NeqSim returns `gamma = 1` here, which is the ideal-solution value rather than \
             the one this state has"
                .to_string(),
        ));
    }
    Ok(ln_a_w.exp() / x_water)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Names, mole fractions, molar masses in kg/mol, charges - a sodium chloride solution.
    fn brine() -> (&'static [&'static str], [f64; 3], [f64; 3], [f64; 3]) {
        (
            &["water", "na+", "cl-"],
            [0.9, 0.05, 0.05],
            [0.018015, 0.02299, 0.03545],
            [0.0, 1.0, -1.0],
        )
    }

    /// **The molality is `n_i / m_solvent`, and the solvent's own entry is not zero.**
    ///
    /// The solvent's molality is its reciprocal molar mass, which is the part of this a
    /// reader is most likely to assume away - a sum over "the molalities" that included
    /// it would carry `55.5 mol/kg` of water into an ionic strength it does not belong to.
    #[test]
    fn the_composition_is_per_kilogram_of_water() {
        let (names, x, mass, charge) = brine();
        let c = composition(names, &x, &mass, &charge).expect("a brine");

        // m_solvent = 0.9 * 0.018015 = 0.0162135 kg per mole of mixture.
        assert!((c.solvent_mass - 0.9 * 0.018015).abs() < 1.0e-15);
        // The solutes: 0.05 / 0.0162135 = 3.0839 mol/kg each.
        assert!((c.molality[1] - 0.05 / (0.9 * 0.018015)).abs() < 1.0e-12);
        assert!((c.molality[1] - c.molality[2]).abs() < 1.0e-12, "1:1");
        // The solvent: 1 / 0.018015 = 55.51 mol/kg.
        assert!((c.molality[0] - 1.0 / 0.018015).abs() < 1.0e-9);

        // I = 1/2 (m_na * 1 + m_cl * 1) = the solute molality for a 1:1 salt.
        assert!((c.ionic_strength - c.molality[1]).abs() < 1.0e-12);
    }

    /// The ionic strength is a valence-weighted sum, so a 2:1 salt is not its molality.
    #[test]
    fn the_ionic_strength_weights_by_charge_squared() {
        let names = ["water", "ca++", "cl-", "cl-"];
        let x = [0.9, 0.033, 0.033, 0.033];
        let mass = [0.018015, 0.04008, 0.03545, 0.03545];
        let charge = [0.0, 2.0, -1.0, -1.0];
        let c = composition(&names, &x, &mass, &charge).expect("a brine");

        let m_calcium = 0.033 / (0.9 * 0.018015);
        let m_chloride = m_calcium;
        let expected = 0.5 * (m_calcium * 4.0 + m_chloride * 1.0 + m_chloride * 1.0);
        assert!((c.ionic_strength - expected).abs() < 1.0e-12);
        // And it is not simply the total molality, which is what makes the weighting
        // worth a test.
        assert!(c.ionic_strength > 1.5 * m_calcium);
    }

    /// A neutral component contributes nothing, without being excluded by a test.
    #[test]
    fn a_neutral_solute_does_not_reach_the_ionic_strength() {
        let names = ["water", "co2"];
        let x = [0.99, 0.01];
        let mass = [0.018015, 0.04401];
        let charge = [0.0, 0.0];
        let c = composition(&names, &x, &mass, &charge).expect("a solution");
        assert!(c.molality[1] > 0.0, "the solute has a molality");
        assert_eq!(c.ionic_strength, 0.0, "and no ionic strength");
    }

    /// A phase without water is refused rather than evaluated as a solution of nothing.
    #[test]
    fn a_phase_without_a_solvent_is_refused() {
        let names = ["methane", "co2"];
        let x = [0.6, 0.4];
        let mass = [0.016043, 0.04401];
        let charge = [0.0, 0.0];
        let error = composition(&names, &x, &mass, &charge).expect_err("no solvent");
        assert!(error.to_string().contains("no solvent"), "{error}");

        // And a mixture with the vectors out of step is refused for the other reason.
        assert!(composition(&names, &x, &mass[..1], &charge).is_err());
    }

    /// The water activity relation, against NeqSim's own constants.
    ///
    /// `phi = 1`, `sum m = 55.51` (pure water) gives `ln a_w = -1 * 0.018015 * 55.51 =
    /// -1.0000`, so `a_w = 0.3679`. That is the *pure water* limit through this
    /// expression, and the interesting part is that the expression is applied to the sum
    /// the model built - so this is a test of the arithmetic rather than of a physical
    /// value.
    #[test]
    fn the_water_activity_is_the_osmotic_relation() {
        let ln_a_w = ln_water_activity(1.0, 1.0 / 0.018015, 0.018015);
        assert!((ln_a_w + 1.0).abs() < 1.0e-12, "ln a_w = {ln_a_w}");
        assert!((ln_a_w.exp() - (-1.0f64).exp()).abs() < 1.0e-12);

        // And `a_w / x_w`, which is the coefficient the kernel takes.
        let gamma = water_activity_coefficient(ln_a_w, 0.9).expect("a positive mole fraction");
        assert!((gamma - ln_a_w.exp() / 0.9).abs() < 1.0e-15);
    }

    /// A zero water mole fraction is refused, where NeqSim returns the ideal value.
    #[test]
    fn a_vanished_water_mole_fraction_is_refused() {
        let error = water_activity_coefficient(-1.0, 0.0).expect_err("not a number");
        assert!(error.to_string().contains("mole fraction"), "{error}");
    }
}
