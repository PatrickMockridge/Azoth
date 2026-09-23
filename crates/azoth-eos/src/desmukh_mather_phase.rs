//! `eos.desmukh_mather_phase` - the activity coefficients of a Desmukh-Mather phase.
//!
//! Spec: `specs/models/eos/desmukh_mather_phase.toml`. A *direct* model: no iteration, so
//! no `algorithm` block.
//!
//! NeqSim's `ComponentDesmukhMather.getGamma` is a Debye-Huckel term plus a pairwise sum,
//! followed by a change of scale:
//!
//! ```text
//! ln gamma*_i = -A z_i^2 sqrt(I) / (1 + B d_i 1e-10 sqrt(I)) + 2 sum_j beta_ij m_j
//! gamma_i     = m_i M_solvent exp(ln gamma*_i) / x_i
//! ```
//!
//! with `A = 1.174`, `B = 3.32384e9`, `d_i` the component's `DeshMatIonicDiameter` in
//! angstrom, and `beta_ij = aij + bij T`. The last line is the mole-fraction conversion
//! `ComponentDesmukhMather` does *inside* the activity coefficient rather than beside it,
//! and `getMolality` is the same `n_i / m_solvent` the other electrolyte models use.
//!
//! # What the pair sum carries
//!
//! **Four rows of the interaction table**, all of them one chemistry: `MDEA+`/`CO2`,
//! `HCO3-`/`MDEA`, `CO3--`/`MDEA` and `HCO3-`/`MDEA+`. `bij` is zero on every row, so
//! `beta` does not vary with temperature at all, and every other pair evaluates to zero -
//! either from a zero row or from no row, which NeqSim cannot tell apart and this can.
//!
//! # The two things it does that the other electrolyte models do not
//!
//! **Its solvent is selected by `REFERENCESTATETYPE`, not by the name `water`.**
//! `PhaseDesmukhMather.getSolventWeight` sums `n_i M_i` over every component whose
//! reference state is `solvent` - which for an amine brine includes the amine, so the
//! molar mass it scales by is a *mixture* mean. `PhasePitzer`'s matches the name.
//!
//! **And its fugacity coefficient is the `(gamma / gamma^infinity) H / P` branch**, with
//! `gamma^infinity` a second evaluation of the same expression at a two-component
//! reference state. Measured, that ratio is one to eleven digits on this data - so the
//! branch *looks* different from Pitzer's `gamma H (m/x) / P` and is the same number, and
//! the two agree only because the reference state's pair terms are absent.

use azoth_core::{AzothError, Result, apply_checks};

use crate::results::DesmukhMatherPhaseResult;

/// The Debye-Huckel `A` of the expression, dimensionless.
pub const A: f64 = 1.174;

/// The gas constant the Poynting correction uses, J/(mol K). NeqSim's `Component.R`.
pub const R: f64 = 8.314_462_618;

/// The `B` that scales the ionic diameter, in reciprocal metres.
///
/// NeqSim multiplies the diameter - which the database carries in angstrom - by `1e-10` at
/// the point of use, so `B d 1e-10` is the coefficient of `sqrt(I)`.
pub const B: f64 = 3.32384e9;

/// The `aij` of the `MDEA`/`MDEA` diagonal, which is **not in the interaction table**.
///
/// `PhaseDesmukhMather.getParameters` assigns it when both components are named `MDEA` and
/// reads nothing for that pair, so the one value the model's own author fitted is in the
/// code and not in the data.
pub const MDEA_DIAGONAL: f64 = -0.0828487;

/// The fugacity coefficient an ion gets, dimensionless: `ComponentDesmukhMather.fugcoef`.
///
/// The third of the tranche's three insoluble-ion constants and the smallest: NeqSim gives
/// `1e12` in `ComponentGePitzer` and `1e8` in `ComponentKentEisenberg`. A `1e-15` says the
/// ion is *absent* from the vapour rather than sparingly present, which is a different
/// claim from the other two.
pub const INSOLUBLE_ION: f64 = 1.0e-15;

/// The activity coefficients and fugacity coefficients of a Desmukh-Mather phase.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if a component is not in the databank, if `x` is not a
///   composition, or if the mixture carries no solvent.
/// * [`AzothError::PropertyUnavailable`] if a `solvent` component carries no Antoine
///   correlation or a neutral solute one carries no Henry row.
/// * [`AzothError::OutOfRange`] if `T` or `P` is not positive.
///
/// # Example
/// ```
/// use azoth_eos::desmukh_mather_phase::desmukh_mather_phase;
///
/// let r = desmukh_mather_phase(
///     &["water", "MDEA+", "Cl-", "CO2"],
///     313.15,
///     500_000.0,
///     &[0.89, 0.04, 0.04, 0.03],
/// )?;
/// assert!((r.ionic_strength - 2.494_799_901_46).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
#[allow(non_snake_case)] // `T`, `P` and `x` are the symbols in the chemistry
pub fn desmukh_mather_phase(
    names: &[&str],
    T: f64,
    P: f64,
    x: &[f64],
) -> Result<DesmukhMatherPhaseResult> {
    let spec = &crate::model_gen::DESMUKH_MATHER_PHASE_SPEC;
    let mut warnings = Vec::new();
    apply_checks(
        spec.input_checks(),
        |quantity| match quantity {
            "T" => Some(T),
            "P" => Some(P),
            _ => None,
        },
        &mut warnings,
    )?;

    let n = names.len();
    if x.len() != n {
        return Err(AzothError::invalid_input(
            "x",
            format!("a mixture of {n} components has {} mole fractions", x.len()),
        ));
    }
    if let Some(bad) = x.iter().position(|&value| value < 0.0) {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "x[{bad}] is {} but a mole fraction cannot be negative",
                x[bad]
            ),
        ));
    }
    let sum: f64 = x.iter().sum();
    if (sum - 1.0).abs() > 1.0e-09 {
        return Err(AzothError::invalid_input(
            "x",
            format!(
                "the mole fractions sum to {sum}, not to one. Renormalising them here \
                 would make a composition error invisible in every number downstream, so \
                 it is refused instead"
            ),
        ));
    }

    let entries: Vec<crate::databank::Entry> = names
        .iter()
        .map(|name| crate::databank::entry(name, None))
        .collect::<Result<Vec<_>>>()?;
    let charge: Vec<f64> = entries.iter().map(|e| e.ionic_charge).collect();
    let diameter: Vec<f64> = entries.iter().map(|e| e.deshmukh_mather_diameter).collect();

    // **The solvent is whatever `REFERENCESTATETYPE` says it is**, which is NeqSim's rule
    // here and not `PhasePitzer`'s. An amine tagged `solvent` joins water in the sum.
    let solvent_weight: f64 = entries
        .iter()
        .zip(x)
        .filter(|(entry, _)| entry.reference_state == crate::databank::SOLVENT)
        .map(|(entry, &xi)| xi * entry.molar_mass.unwrap_or_default())
        .sum();
    let solvent_moles: f64 = entries
        .iter()
        .zip(x)
        .filter(|(entry, _)| entry.reference_state == crate::databank::SOLVENT)
        .map(|(_, &xi)| xi)
        .sum();
    if solvent_weight <= 0.0 || solvent_moles <= 0.0 {
        return Err(AzothError::invalid_input(
            "components",
            "the mixture carries no `solvent`-reference component, so its molalities would \
             all be zero. NeqSim divides by a zero solvent weight and returns a phase of \
             nothing"
                .to_string(),
        ));
    }
    let solvent_molar_mass = solvent_weight / solvent_moles;

    let molality: Vec<f64> = x.iter().map(|&xi| xi / solvent_weight).collect();
    let ionic_strength = 0.5
        * molality
            .iter()
            .zip(&charge)
            .map(|(&m, &z)| m * z * z)
            .sum::<f64>();

    let ln_gamma: Vec<f64> = (0..n)
        .map(|i| {
            let star = ln_gamma_star(names, &molality, &charge, &diameter, T, i);
            let value = molality[i].max(0.0) * solvent_molar_mass * star.exp() / x[i].max(1.0e-30);
            // **NeqSim's own fallback**: `gamma` is set to `exp(ln gamma)` when the
            // conversion is not finite or not positive, which never happens for a present
            // component. Reproduced rather than improved on.
            if value > 0.0 && value.is_finite() {
                value.ln()
            } else {
                star
            }
        })
        .collect();

    let mut ln_phi = Vec::with_capacity(n);
    for (i, entry) in entries.iter().enumerate() {
        // **Four branches, and the first is a name test.** `ComponentDesmukhMather.fugcoef`
        // asks for `water` *by name* before it asks what the reference state is, so a
        // solvent-reference component that is not water takes the second branch and gets
        // no Poynting correction - which is a third of a percent at 313 K and 5 bar.
        let solvent = entry.reference_state == crate::databank::SOLVENT;
        let neutral = entry.ionic_charge == 0.0;
        let coefficient = if entry.name.eq_ignore_ascii_case("water") {
            let p0 = crate::antoine_vapor_pressure::saturation_pressure(entry, T, &mut warnings)?;
            // `exp(v (P - P0) / R T)`, with `v` the liquid molar volume. NeqSim converts
            // the pressure difference from bar with a `1e5` because its pressures are in
            // bar; everything here is SI already, so the factor is absent.
            let molar_volume = 1.0e-3 * entry.molar_mass.unwrap_or_default();
            let poynting = (molar_volume / (R * T) * (P - p0)).exp();
            ln_gamma[i].exp() * p0 / P * poynting
        } else if neutral && solvent {
            let p0 = crate::antoine_vapor_pressure::saturation_pressure(entry, T, &mut warnings)?;
            ln_gamma[i].exp() * p0 / P
        } else if neutral && entry.reference_state == crate::databank::SOLUTE {
            let infinite =
                ln_gamma_infinite_dilution(&entry.name, entry.ionic_charge, diameter[i], T);
            crate::henry::effective_coefficient(entry, T)? * 1.0e5 / P
                * (ln_gamma[i] - infinite).exp()
        } else {
            INSOLUBLE_ION
        };
        if coefficient <= 0.0 || !coefficient.is_finite() {
            return Err(AzothError::property_unavailable(
                entry.name.clone(),
                "fugacity coefficient".to_string(),
                format!(
                    "the branch its reference state selects gives {coefficient}, which is \
                     not a coefficient. NeqSim returns it without comment"
                ),
            ));
        }
        ln_phi.push(coefficient.ln());
    }

    Ok(DesmukhMatherPhaseResult {
        gamma: ln_gamma.iter().map(|value| value.exp()).collect(),
        ln_gamma,
        molality,
        ionic_strength,
        solvent_molar_mass,
        ln_phi,
        warnings,
    })
}

/// `ln gamma*`, the expression before the mole-fraction conversion.
///
/// The ionic strength is computed here rather than passed: it is `1/2 sum m_i z_i^2` over
/// the same three slices, so a caller computing it separately could disagree with this.
fn ln_gamma_star(
    names: &[&str],
    molality: &[f64],
    charge: &[f64],
    diameter: &[f64],
    temperature: f64,
    component: usize,
) -> f64 {
    let ionic_strength = 0.5
        * molality
            .iter()
            .zip(charge)
            .map(|(&m, &z)| m * z * z)
            .sum::<f64>();
    let sqrt_i = ionic_strength.sqrt();
    let debye_huckel = -A * charge[component] * charge[component] * sqrt_i
        / (1.0 + B * diameter[component] * 1.0e-10 * sqrt_i);

    let mut pairs = 0.0;
    for (j, &m) in molality.iter().enumerate() {
        // **NeqSim's loop skips neither the diagonal nor water.** It sums over every
        // component and relies on the diagonal's `beta` being zero - true except for
        // `MDEA`/`MDEA`, which `beta` handles - and `ComponentDesmukhMather` tests for
        // water by name before adding its term, which is the one thing the loop does skip.
        if j == component || names[j].eq_ignore_ascii_case("water") {
            continue;
        }
        pairs += 2.0 * beta(names[component], names[j], temperature) * m;
    }
    debye_huckel + pairs
}

/// `beta_ij = aij + bij T`, either order round.
///
/// The `MDEA`/`MDEA` diagonal comes from [`MDEA_DIAGONAL`], because the interaction table
/// carries no row for it and `PhaseDesmukhMather` states the number in code.
fn beta(first: &str, second: &str, temperature: f64) -> f64 {
    if first.eq_ignore_ascii_case("MDEA") && second.eq_ignore_ascii_case("MDEA") {
        return MDEA_DIAGONAL;
    }
    match crate::databank::desmukh_mather_pair(first, second) {
        Some((aij, bij)) => aij + bij * temperature,
        // NeqSim's own reader leaves its arrays at zero when the query finds nothing, so
        // an absent pair and a fitted zero evaluate the same there. They are different
        // here and the model states which it uses.
        None => 0.0,
    }
}

/// `gamma^infinity` for a solute: the same expression at the two-component reference state
/// `Phase.initRefPhases` builds.
///
/// That phase is the solute at `1e-10` mol and **water at 10 mol**, so its solvent weight
/// is `10 M_water` and its only pair is solute-with-water. Measured, the ratio this gives
/// is one to eleven digits on the brines the tranche tested, so the solute branch of
/// `fugcoef` and Pitzer's `gamma H (m/x) / P` are the same number here - and they agree
/// because the reference state's pair terms are absent rather than because they must.
fn ln_gamma_infinite_dilution(name: &str, charge: f64, diameter: f64, temperature: f64) -> f64 {
    let water = crate::pitzer_phase::WATER_MOLAR_MASS;
    let weight = 10.0 * water;
    let names = [name, "water"];
    let molality = [1.0e-10 / weight, 10.0 / weight];
    let charge = [charge, 0.0];
    let diameter = [diameter, 0.0];
    let star = ln_gamma_star(&names, &molality, &charge, &diameter, temperature, 0);
    // **The mole fraction, not the mole count.** The reference phase holds `1e-10` mol of
    // solute beside 10 mol of water, so its fraction is `1e-11` and the first draft of
    // this read `1e-10` - a factor of ten on every solute's `gamma^infinity` and so on
    // every solute's fugacity coefficient.
    let reference_fraction = 1.0e-10 / (10.0 + 1.0e-10);
    let value = molality[0] * water * star.exp() / reference_fraction;
    if value > 0.0 && value.is_finite() {
        value.ln()
    } else {
        star
    }
}
