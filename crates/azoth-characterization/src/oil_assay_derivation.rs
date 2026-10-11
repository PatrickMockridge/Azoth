//! `characterization.oil_assay_derivation` - an assay's cuts derived to molar masses and moles.
//!
//! Spec: `specs/models/characterization/oil_assay_derivation.toml`. Oracle:
//! `validation/neqsim/captures/oil_assay_derivation_probe.tsv`.
//!
//! # The arithmetic, and not the mutation
//!
//! `OilAssayCharacterisation.apply` resolves the cuts' mass fractions and then derives one component
//! per cut with a positive yield. In NeqSim that component goes into a `SystemInterface`; here the
//! answer is the table - each component's name, molar mass and moles - and the fluid is the
//! caller's step, which is what every id in this namespace does with the thing a class mutates.
//!
//! # Three routes to one molar mass, in this order
//!
//! A **standard component** is a substance the database knows, and its molar mass comes out of
//! `azoth_eos::databank` by name. A **stated molar mass** is used as it is. Otherwise the closed
//! form `5.805e-5 * Tb^2.3776 / SG^0.9371` is evaluated, in g/mol and g/cm3 - the class's own
//! scale - on a boiling point the cut states or one resolved from its Watson factor, and a gravity
//! it states or one resolved from its API gravity.
//!
//! # The gate is the class's own
//!
//! `|sum(moles * M) - total| / max(total, 1e-30)` must be within `1e-10`, and the value is
//! published. It closes by construction - each cut's moles are its mass over its molar mass - so
//! the number says what the arithmetic did rather than whether it converged.

use std::collections::HashSet;

use azoth_core::units::{MassDensity, MolarMass, kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, Result, apply_checks};

use crate::assay_mass_fractions::{AssayBasis, AssayCut, assay_mass_fractions};
use crate::model_gen;
use crate::results::OilAssayDerivationResult;

/// `apply`'s own gate, on the mass it rebuilt against the mass it started from.
const MASS_BALANCE_TOLERANCE: f64 = 1.0e-10;

/// The divisor's floor, the class's `Math.max(totalAssayMass, 1.0e-30)`.
const MASS_BALANCE_FLOOR: f64 = 1.0e-30;

/// The `141.5` and `131.5` of `AssayCut.resolveDensity`'s API relation.
const API_NUMERATOR: f64 = 141.5;
const API_OFFSET: f64 = 131.5;

/// The `5.805e-5` and the two exponents of `AssayCut.resolveMolarMass`.
const MOLAR_MASS_COEFFICIENT: f64 = 5.805e-5;
const BOILING_POINT_EXPONENT: f64 = 2.3776;
const DENSITY_EXPONENT: f64 = 0.9371;

/// The `1.8` of `AssayCut.resolveAverageBoilingPoint`'s Watson relation.
const WATSON_DIVISOR: f64 = 1.8;

/// One cut of an assay, as the derivation takes it.
///
/// The derivation reads the yields and the properties a cut's molar mass is derived from, and every
/// route to one is a field here. NeqSim's `AssayCut` holds the same ones among its twelve nullable
/// Doubles, and this row grows as the ids that read more of them are ported.
#[derive(Debug, Clone, PartialEq)]
pub struct AssayCutDerivation {
    /// The cut's name, which is what its pseudo-component is named after.
    pub name: String,
    /// The cut's declared yield, on the assay's basis.
    pub declared_fraction: f64,
    /// The cut's normal liquid density, where it states one.
    pub density: Option<MassDensity>,
    /// The cut's API gravity, where it states one and no density.
    pub api_gravity: Option<f64>,
    /// The cut's molar mass, where it is known rather than correlated.
    pub molar_mass: Option<MolarMass>,
    /// The cut's average normal boiling point, where it states one.
    pub average_boiling_point: Option<f64>,
    /// The cut's Watson characterization factor, where it states one and no boiling point.
    pub watson_characterization_factor: Option<f64>,
    /// The name of the database substance this cut *is*, where it is one.
    pub standard_component_name: Option<String>,
}

/// `AssayCut.resolveDensity`: the stated gravity, or the one its API gravity implies.
///
/// g/cm3, which is the scale every correlation in the class is in and the scale the API relation
/// is written in. `None` where the cut states neither, which the derivation refuses.
fn specific_gravity(cut: &AssayCutDerivation) -> Option<f64> {
    if let Some(density) = cut.density {
        return Some(density.value / 1000.0);
    }
    cut.api_gravity
        .map(|api| API_NUMERATOR / (api + API_OFFSET))
}

/// `AssayCut.resolveAverageBoilingPoint`: the stated boiling point, or the one its Watson factor
/// implies.
fn average_boiling_point(cut: &AssayCutDerivation, gravity: f64) -> Option<f64> {
    if let Some(boiling_point) = cut.average_boiling_point {
        return Some(boiling_point);
    }
    cut.watson_characterization_factor
        .map(|watson| (watson * gravity).powi(3) / WATSON_DIVISOR)
}

/// `AssayCut.resolveMolarMass`'s closed form, in g/mol from a boiling point in K and a gravity in
/// g/cm3 - the inverse of `tbp_cut_properties`' `TBPBaseModel.calcTB`.
fn correlated_molar_mass(boiling_point: f64, gravity: f64) -> f64 {
    MOLAR_MASS_COEFFICIENT * boiling_point.powf(BOILING_POINT_EXPONENT)
        / gravity.powf(DENSITY_EXPONENT)
        / 1000.0
}

/// An assay's cuts derived to one component each, by `OilAssayCharacterisation.apply`'s arithmetic.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if two cuts share a name, if two standard components name one
///   substance, if a standard component is on a volume basis or also states pseudo-component
///   properties, if a cut reaches the correlation with neither a boiling point nor an API gravity,
///   if a standard component is not in the database, or if the mass it rebuilt misses the total by
///   more than `1e-10` relative.
///
/// # Example
/// ```
/// use azoth_core::units::kilograms_per_cubic_meter;
/// use azoth_characterization::{AssayBasis, AssayCutDerivation, oil_assay_derivation};
///
/// let derived = oil_assay_derivation(
///     AssayBasis::Mass,
///     &[AssayCutDerivation {
///         name: "cut_a".to_string(),
///         declared_fraction: 1.0,
///         density: Some(kilograms_per_cubic_meter(780.0)),
///         api_gravity: None,
///         molar_mass: None,
///         average_boiling_point: Some(400.0),
///         watson_characterization_factor: None,
///         standard_component_name: None,
///     }],
///     1000.0,
/// )?;
/// assert_eq!(derived.component_name, vec!["cut_a_PC"]);
/// assert!((derived.molar_mass[0].value - 0.112_610_323_807_509).abs() < 1e-15);
/// assert!((derived.total_moles - 1000.0 / 0.112_610_323_807_509).abs() < 1e-9);
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn oil_assay_derivation(
    basis: AssayBasis,
    cuts: &[AssayCutDerivation],
    total_assay_mass: f64,
) -> Result<OilAssayDerivationResult> {
    let spec = &model_gen::OIL_ASSAY_DERIVATION_SPEC;
    let mut warnings = Vec::new();

    apply_checks(spec.input_checks(), |_quantity| None, &mut warnings)?;

    // `validateUniqueCutNames`: by name, case-sensitively.
    let mut seen: HashSet<&str> = HashSet::new();
    for cut in cuts {
        if !seen.insert(cut.name.as_str()) {
            return Err(AzothError::invalid_input(
                "cuts",
                format!("two cuts carry the name `{}`, and a component cannot be named twice", cut.name),
            ));
        }
    }
    // `validateUniqueStandardComponentNames`: by lower-cased name, because a substance is one
    // substance whatever case its name was written in.
    let mut substances: HashSet<String> = HashSet::new();
    for cut in cuts {
        if let Some(name) = &cut.standard_component_name
            && !substances.insert(name.to_lowercase())
        {
            return Err(AzothError::invalid_input(
                "cuts",
                format!("two cuts name the standard component `{name}`, which is one substance"),
            ));
        }
    }

    // The resolution, over the gravity each cut actually has rather than over the field it was
    // stated in: `apply` reaches `resolveDensity` for every cut it derives.
    let mut resolved: Vec<Option<f64>> = Vec::with_capacity(cuts.len());
    let mut rows: Vec<AssayCut> = Vec::with_capacity(cuts.len());
    for cut in cuts {
        let gravity = specific_gravity(cut);
        resolved.push(gravity);
        rows.push(AssayCut {
            declared_fraction: cut.declared_fraction,
            // The assay kernel's own resolution is over SI densities; a cut that states no gravity
            // passes none, which its mass branch never reads and its volume branch refuses.
            density: gravity.map(|value| kilograms_per_cubic_meter(value * 1000.0)),
        });
    }
    let assay = assay_mass_fractions(basis, &rows)?;
    let mass_fractions = &assay.mass_fraction;

    let mut component_name = Vec::new();
    let mut derived_fraction = Vec::new();
    let mut molar_mass = Vec::new();
    let mut moles = Vec::new();
    let mut reconstructed = 0.0;
    for (index, cut) in cuts.iter().enumerate() {
        let mass_fraction = mass_fractions[index];
        // The class's own `!(massFraction > 0.0)`, written the way these crates write a positive
        // test: a NaN and a non-positive both skip, and an infinity does not.
        if mass_fraction.is_nan() || mass_fraction <= 0.0 {
            continue;
        }
        let (name, cut_molar_mass) = match &cut.standard_component_name {
            Some(substance) => {
                if basis == AssayBasis::Volume {
                    return Err(AzothError::invalid_input(
                        "cuts",
                        format!(
                            "`{}` is the standard component `{substance}`, which needs a mass-basis \
                             assay fraction",
                            cut.name
                        ),
                    ));
                }
                if cut.density.is_some()
                    || cut.api_gravity.is_some()
                    || cut.molar_mass.is_some()
                    || cut.average_boiling_point.is_some()
                    || cut.watson_characterization_factor.is_some()
                {
                    return Err(AzothError::invalid_input(
                        "cuts",
                        format!(
                            "`{}` is the standard component `{substance}` and also states \
                             pseudo-component properties, which is one cut wanting two definitions",
                            cut.name
                        ),
                    ));
                }
                let entry = azoth_eos::databank::entry(substance, None).map_err(|_| {
                    AzothError::invalid_input(
                        "cuts",
                        format!("`{substance}` is not a substance the databank carries"),
                    )
                })?;
                let mass = entry.molar_mass.filter(|value| value.is_finite() && *value > 0.0);
                let Some(mass) = mass else {
                    return Err(AzothError::invalid_input(
                        "cuts",
                        format!("the databank carries no molar mass for `{substance}`"),
                    ));
                };
                (substance.clone(), mass)
            }
            None => {
                let Some(gravity) = resolved[index] else {
                    return Err(AzothError::invalid_input(
                        format!("cuts[{index}].density"),
                        format!(
                            "`{}` states neither a density nor an API gravity, and the mass the \
                             moles are taken over needs one",
                            cut.name
                        ),
                    ));
                };
                let mass = match cut.molar_mass {
                    Some(stated) => stated.value,
                    None => {
                        let Some(boiling_point) = average_boiling_point(cut, gravity) else {
                            return Err(AzothError::invalid_input(
                                format!("cuts[{index}].average_boiling_point"),
                                format!(
                                    "`{}` states neither an average boiling point nor a Watson \
                                     factor, and one of them is what the correlation is evaluated on",
                                    cut.name
                                ),
                            ));
                        };
                        correlated_molar_mass(boiling_point, gravity)
                    }
                };
                (format!("{}_PC", cut.name), mass)
            }
        };
        let cut_moles = total_assay_mass * mass_fraction / cut_molar_mass;
        reconstructed += cut_moles * cut_molar_mass;
        component_name.push(name);
        derived_fraction.push(mass_fraction);
        molar_mass.push(kilograms_per_mole(cut_molar_mass));
        moles.push(cut_moles);
    }

    let mass_balance_error = (reconstructed - total_assay_mass).abs()
        / total_assay_mass.abs().max(MASS_BALANCE_FLOOR);
    if !mass_balance_error.is_finite() || mass_balance_error > MASS_BALANCE_TOLERANCE {
        return Err(AzothError::invalid_input(
            "total_assay_mass",
            format!(
                "the derived components come to {reconstructed} against a total of \
                 {total_assay_mass}, which is {mass_balance_error} relative and past the class's \
                 {MASS_BALANCE_TOLERANCE}"
            ),
        ));
    }

    let total_moles = moles.iter().sum();
    apply_checks(spec.derived_checks(), |_quantity| None, &mut warnings)?;

    Ok(OilAssayDerivationResult {
        component_name,
        mass_fraction: derived_fraction,
        molar_mass,
        moles,
        total_moles,
        mass_balance_error,
        warnings,
    })
}
