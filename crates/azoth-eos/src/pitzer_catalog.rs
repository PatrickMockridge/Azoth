//! The PHREEQC Pitzer catalogue, and the rule that chooses between it and the legacy CSV.
//!
//! # Two datasets, one phase
//!
//! `PhasePitzer.loadParametersFromDatabase` tries
//! `PitzerParameterDatasets.tryApplyCompletePhreeqcPitzerCatalog` **and returns if it
//! succeeds**, reading the `pitzerparameters` table - `PitzerParameters.csv`, the 30 rows
//! this library has always had - only when it fails. Measured, the choice is per-mixture:
//! `water + Na+ + Cl-` takes the catalogue and `water + Na+ + HCO3-` takes the CSV,
//! because the catalogue has no `C0` row for that pair.
//!
//! **And the two disagree on pairs they share.** Na+/Cl- is `0.07534 / 0.2769 / 0.00148`
//! under PHREEQC against the CSV's `0.0765 / 0.2664 / 0.00127` - Cphi by 16%. So which
//! dataset applies is not a detail: it is the difference between two parameter sets for
//! the same fluid.
//!
//! # The coverage rule, exactly
//!
//! The catalogue applies only when it covers the phase's whole topology, and the
//! requirements are pairwise:
//!
//! | topology | required |
//! |---|---|
//! | an opposite-sign ion pair | `B0`, `B1`, `C0` |
//! | a same-sign ion pair | `THETA` |
//! | three ions with one or two positive | `PSI` |
//! | a neutral with another neutral, itself included, and with each ion | `LAMBDA` |
//! | a neutral with a cation and an anion | `ZETA` |
//!
//! `B2` is the one that is *optional*: `applyCatalogIonRows` looks it up rather than
//! requiring it, and applies it only where the catalogue carries one.
//!
//! **One requirement can never be met and is not a data question.** A phase with no ions
//! returns `false` before the catalogue is consulted at all, so an ion-free mixture always
//! takes the CSV - which is why `water + CO2` does despite the catalogue carrying
//! `CO2|CO2` and `CO2|Na+`.
//!
//! # The species names
//!
//! PHREEQC's, canonicalised: a numeric charge suffix becomes NeqSim's repeated-sign form
//! (`Ba+2` -> `Ba++`, `SO4-2` -> `SO4--`) and nothing else changes, so `B(OH)4-` is left
//! alone. It is a rewrite and not a table, which is what makes it safe - a species the
//! catalogue carries and the databank does not still compiles, and the model refuses it
//! later rather than here.

use std::collections::HashMap;
use std::sync::OnceLock;

/// A Pitzer parameter family, as `PhreeqcPitzerParameterCatalog.Family` names it.
///
/// `Mu`, `Eta` and `Alphas` are declared and **empty in the shipped catalogue** - the
/// enum is wider than the file - so they are listed to be looked up and found absent,
/// which is a different answer from a family that does not exist.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Family {
    /// The binary constant term of an opposite-sign ion pair.
    B0,
    /// The binary `1/T`-like term of the same pair.
    B1,
    /// The binary term a 2:2 pair carries and a 1:1 pair does not. Optional.
    B2,
    /// The third binary term.
    C0,
    /// The same-sign ion mixing parameter.
    Theta,
    /// The ternary mixing parameter.
    Psi,
    /// A neutral's interaction with another species.
    Lambda,
    /// A neutral's interaction with a cation and an anion.
    Zeta,
    /// Declared upstream, absent from the catalogue.
    Mu,
    /// Declared upstream, absent from the catalogue.
    Eta,
    /// Declared upstream, absent from the catalogue.
    Alphas,
}

impl Family {
    /// Every family, in the enum's own order.
    pub const ALL: [Family; 11] = [
        Family::B0,
        Family::B1,
        Family::B2,
        Family::C0,
        Family::Theta,
        Family::Psi,
        Family::Lambda,
        Family::Zeta,
        Family::Mu,
        Family::Eta,
        Family::Alphas,
    ];

    /// How many species a row of this family names.
    #[must_use]
    pub fn species_count(self) -> usize {
        match self {
            Family::Psi | Family::Zeta | Family::Mu | Family::Eta => 3,
            _ => 2,
        }
    }

    /// The name the catalogue and the compiled table write.
    #[must_use]
    pub fn name(self) -> &'static str {
        match self {
            Family::B0 => "B0",
            Family::B1 => "B1",
            Family::B2 => "B2",
            Family::C0 => "C0",
            Family::Theta => "THETA",
            Family::Psi => "PSI",
            Family::Lambda => "LAMBDA",
            Family::Zeta => "ZETA",
            Family::Mu => "MU",
            Family::Eta => "ETA",
            Family::Alphas => "ALPHAS",
        }
    }
}

/// One row of the catalogue: its family, its species and the six temperature coefficients.
#[derive(Debug, Clone, PartialEq)]
pub struct CatalogRow {
    /// Which parameter this row is.
    pub family: Family,
    /// The species, canonicalised and **in the order the catalogue wrote them**.
    pub species: Vec<String>,
    /// `a0`..`a5` of `PitzerTemperatureFunction`, whose units are the family's own.
    pub a: [f64; 6],
}

/// One species name in NeqSim's spelling: `PhreeqcPitzerParameterCatalog.canonicalSpeciesName`.
///
/// A numeric charge suffix becomes the repeated-sign form and nothing else changes. The
/// order of the tests matters: `+3` before `+2`, or `Fe+3` would come back `Fe+` with a
/// stray `3`.
#[must_use]
pub fn canonical_species(name: &str) -> String {
    for (suffix, sign) in [("+3", "+++"), ("-3", "---"), ("+2", "++"), ("-2", "--")] {
        if let Some(stem) = name.strip_suffix(suffix) {
            return format!("{stem}{sign}");
        }
    }
    name.to_string()
}

/// The catalogue's lookup key: the canonical species **sorted** and joined with `|`.
///
/// Sorted, so a row is found whichever order its species are written in - the catalogue
/// writes `Ba+2 Cl-` and `Cl- H+`, and both must be found by a caller holding the names
/// in either order.
#[must_use]
pub fn species_key(species: &[&str]) -> String {
    let mut canonical: Vec<String> = species.iter().map(|name| canonical_species(name)).collect();
    canonical.sort();
    canonical.join("|")
}

/// The compiled catalogue, by `(family, species_key)`.
fn rows() -> &'static HashMap<(Family, String), [f64; 6]> {
    fn parse() -> HashMap<(Family, String), [f64; 6]> {
        let text = crate::databank::PITZER_PHREEQC_CSV;
        // By name rather than by position, off the file's own header: a column inserted in
        // the middle should resolve to the wrong *name* rather than shift every value.
        let index: HashMap<&str, usize> = text
            .lines()
            .next()
            .unwrap_or("")
            .split(',')
            .enumerate()
            .map(|(position, name)| (name.trim(), position))
            .collect();
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(text.as_bytes());
        let mut out = HashMap::new();
        for record in reader.records().flatten() {
            let field = |name: &str| -> &str {
                index
                    .get(name)
                    .and_then(|position| record.get(*position))
                    .unwrap_or("")
                    .trim()
            };
            let Some(family) = Family::ALL
                .into_iter()
                .find(|candidate| candidate.name() == field("family"))
            else {
                continue;
            };
            let key = field("species_key").to_string();
            if key.is_empty() {
                continue;
            }
            let mut a = [0.0; 6];
            for (index, slot) in a.iter_mut().enumerate() {
                *slot = field(&format!("a{index}")).parse().unwrap_or_default();
            }
            out.insert((family, key), a);
        }
        out
    }
    static ROWS: OnceLock<HashMap<(Family, String), [f64; 6]>> = OnceLock::new();
    ROWS.get_or_init(parse)
}

/// Every row of one family, in the compiled table's order.
#[must_use]
pub fn family_rows(family: Family) -> Vec<CatalogRow> {
    let mut out: Vec<CatalogRow> = rows()
        .iter()
        .filter(|((row_family, _), _)| *row_family == family)
        .map(|((_, key), a)| CatalogRow {
            family,
            species: key.split('|').map(str::to_string).collect(),
            a: *a,
        })
        .collect();
    // The map has no order, so sort by key to make this a function of the data rather
    // than of the hasher.
    out.sort_by(|left, right| left.species.cmp(&right.species));
    out
}

/// One row's six coefficients, or `None` where the catalogue carries no such row.
///
/// `find` and not `require`: the caller decides whether an absent row is a fallback, an
/// optional term or an error, which is exactly the distinction `B2` and the coverage rule
/// are built on.
#[must_use]
pub fn find(family: Family, species: &[&str]) -> Option<[f64; 6]> {
    rows().get(&(family, species_key(species))).copied()
}

/// One component as the selection rule sees it.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Species<'a> {
    /// The databank name, which is the catalogue's species name after canonicalisation.
    pub name: &'a str,
    /// Moles in the phase. A component below `1e-20` is not in the topology.
    pub moles: f64,
    /// The ionic charge, in units of the elementary charge.
    pub charge: f64,
    /// The component's formula, for the hydrocarbon test.
    pub formula: &'a str,
    /// NeqSim's `componentType == "HC"` clause of `isHydrocarbon`.
    pub hydrocarbon: bool,
}

/// Which dataset covers a phase's topology.
#[derive(Debug, Clone, PartialEq)]
pub enum Selection {
    /// The PHREEQC catalogue covers every required pair, and is in force.
    Phreeqc,
    /// It does not, and `PitzerParameters.csv` is. The reason is carried because the two
    /// datasets give different numbers for the same pair, so *why* is part of the answer.
    Legacy(LegacyReason),
}

/// Why the catalogue was not applied.
#[derive(Debug, Clone, PartialEq)]
pub enum LegacyReason {
    /// The phase has no ion above the activity threshold, so the catalogue is never
    /// consulted: `tryApplyCompletePhreeqcPitzerCatalog` returns `false` before it looks.
    NoIons,
    /// A required row is absent. The string is the family and the species it wanted, in
    /// the shape `require`'s own message has.
    Uncovered(String),
}

/// Whether a neutral component is a hydrocarbon the automatic catalogue excludes.
///
/// `PitzerParameterDatasets.isHydrocarbonForAutomaticCatalog`: NeqSim's own flag, or **a
/// formula of carbon, hydrogen and digits with both present**. The formula clause is what
/// the method exists for - a database component such as methane keeps the type `normal`
/// in a GE phase, so the type alone would let it into a topology the catalogue has no
/// rows for.
///
/// **Two of NeqSim's three flag clauses have no counterpart here.** `isIsTBPfraction` and
/// `isPlusFraction` are set through the API rather than read from the databank, and the
/// components they mark are typed `HC` in the table anyway, so the type clause covers
/// them. Recording that rather than inventing a flag.
#[must_use]
pub fn is_hydrocarbon(species: &Species<'_>) -> bool {
    if species.hydrocarbon {
        return true;
    }
    if species.formula.is_empty() {
        return false;
    }
    let mut carbon = false;
    let mut hydrogen = false;
    for character in species.formula.chars() {
        match character {
            'C' => carbon = true,
            'H' => hydrogen = true,
            other if other.is_ascii_digit() => {}
            _ => return false,
        }
    }
    carbon && hydrogen
}

/// Moles below which a component is not in the topology.
const ACTIVE_MOLES: f64 = 1.0e-20;

/// Charge below which a component is a neutral rather than an ion.
const ACTIVE_CHARGE: f64 = 0.5;

/// Which dataset covers a phase's topology, and why not when the answer is the CSV.
///
/// A pure function of the composition, which is the point: NeqSim decides this inside
/// `loadParametersFromDatabase` against a live `PhasePitzer`, and the rule is about the
/// *topology* rather than the phase's state, so it can be answered - and tested - without
/// one.
#[must_use]
pub fn select_dataset(species: &[Species<'_>]) -> Selection {
    let active: Vec<&Species<'_>> = species.iter().filter(|s| s.moles > ACTIVE_MOLES).collect();
    let ions: Vec<&&Species<'_>> = active
        .iter()
        .filter(|s| s.charge.abs() >= ACTIVE_CHARGE)
        .collect();

    // The first thing `tryApplyCompletePhreeqcPitzerCatalog` does, and the reason an
    // ion-free brine takes the CSV however well the catalogue covers it.
    if ions.is_empty() {
        return Selection::Legacy(LegacyReason::NoIons);
    }

    let neutrals: Vec<&&Species<'_>> = active
        .iter()
        .filter(|s| {
            s.charge.abs() < ACTIVE_CHARGE
                && !s.name.eq_ignore_ascii_case("water")
                && !is_hydrocarbon(s)
        })
        .collect();

    let missing = |family: Family, names: &[&str]| -> Option<LegacyReason> {
        find(family, names)
            .is_none()
            .then(|| LegacyReason::Uncovered(format!("{} for {}", family.name(), names.join(", "))))
    };

    for (first, one) in ions.iter().enumerate() {
        for other in ions.iter().skip(first + 1) {
            let names = [one.name, other.name];
            if one.charge * other.charge < 0.0 {
                for family in [Family::B0, Family::B1, Family::C0] {
                    if let Some(reason) = missing(family, &names) {
                        return Selection::Legacy(reason);
                    }
                }
            } else if let Some(reason) = missing(Family::Theta, &names) {
                return Selection::Legacy(reason);
            }
        }
    }

    // A triple with one or two positives, which is the shape PSI describes: two of one
    // sign and one of the other.
    for (first, one) in ions.iter().enumerate() {
        for (second, other) in ions.iter().enumerate().skip(first + 1) {
            for third in ions.iter().skip(second + 1) {
                let positives = [one, other, third]
                    .iter()
                    .filter(|s| s.charge > 0.0)
                    .count();
                if positives == 1 || positives == 2 {
                    let names = [one.name, other.name, third.name];
                    if let Some(reason) = missing(Family::Psi, &names) {
                        return Selection::Legacy(reason);
                    }
                }
            }
        }
    }

    for (position, neutral) in neutrals.iter().enumerate() {
        for other in neutrals.iter().skip(position) {
            let names = [neutral.name, other.name];
            if let Some(reason) = missing(Family::Lambda, &names) {
                return Selection::Legacy(reason);
            }
        }
        for ion in &ions {
            let names = [neutral.name, ion.name];
            if let Some(reason) = missing(Family::Lambda, &names) {
                return Selection::Legacy(reason);
            }
        }
        for cation in ions.iter().filter(|s| s.charge > 0.0) {
            for anion in ions.iter().filter(|s| s.charge < 0.0) {
                let names = [neutral.name, cation.name, anion.name];
                if let Some(reason) = missing(Family::Zeta, &names) {
                    return Selection::Legacy(reason);
                }
            }
        }
    }

    Selection::Phreeqc
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::databank;

    fn species<'a>(
        ions: &'a [(&'a str, f64)],
        neutrals: &'a [(&'a str, &'a str)],
    ) -> Vec<Species<'a>> {
        let mut out = vec![Species {
            name: "water",
            moles: 0.9,
            charge: 0.0,
            formula: "H2O",
            hydrocarbon: false,
        }];
        out.extend(neutrals.iter().map(|(name, formula)| Species {
            name,
            moles: 0.01,
            charge: 0.0,
            formula,
            hydrocarbon: false,
        }));
        out.extend(ions.iter().map(|(name, charge)| Species {
            name,
            moles: 0.05,
            charge: *charge,
            formula: "",
            hydrocarbon: false,
        }));
        out
    }

    /// The catalogue compiles, and its rows are found by the key NeqSim looks up by.
    #[test]
    fn the_catalogue_is_read() {
        // The probe's own reading: `b0(Na+,Cl-) = 0.07534` under PHREEQC.
        let sodium_chloride = find(Family::B0, &["Na+", "Cl-"]).expect("the commonest pair");
        assert!(
            (sodium_chloride[0] - 0.07534).abs() < 1.0e-12,
            "{sodium_chloride:?}"
        );

        // **Either order round**, because the key is sorted - the catalogue writes
        // `Cl- Na+` and a caller holds the cation first.
        assert_eq!(find(Family::B0, &["Cl-", "Na+"]), Some(sodium_chloride));

        // **And PHREEQC's spelling of a charge is canonicalised**, so a caller holding
        // NeqSim's `Ba++` finds the row the catalogue wrote as `Ba+2`.
        assert!(find(Family::B0, &["Ba++", "Cl-"]).is_some());
        assert_eq!(canonical_species("Ba+2"), "Ba++");
        assert_eq!(canonical_species("SO4-2"), "SO4--");
        assert_eq!(canonical_species("Fe+3"), "Fe+++");
        // A species with no numeric suffix is untouched, brackets and all.
        assert_eq!(canonical_species("B(OH)4-"), "B(OH)4-");

        // **The families the enum declares and the file does not carry** are empty, which
        // is a different answer from a family that does not exist.
        assert!(family_rows(Family::Mu).is_empty());
        assert!(family_rows(Family::Eta).is_empty());
        assert!(family_rows(Family::Alphas).is_empty());
        assert_eq!(family_rows(Family::B0).len(), 54);
        assert_eq!(
            Family::ALL
                .iter()
                .map(|f| family_rows(*f).len())
                .sum::<usize>(),
            268,
            "the file's whole row count"
        );
    }

    /// **An ion-free phase never reaches the catalogue**, however well it covers it.
    ///
    /// This is the first thing `tryApplyCompletePhreeqcPitzerCatalog` does, and it is why
    /// `water + CO2` takes the CSV although the catalogue carries `CO2|CO2` and
    /// `CO2|Na+`. Measured on the real system; reproduced here without one.
    #[test]
    fn an_ion_free_phase_falls_back_before_the_catalogue_is_consulted() {
        let mixture = species(&[], &[("CO2", "CO2")]);
        assert!(
            find(Family::Lambda, &["CO2", "CO2"]).is_some(),
            "the catalogue does cover it, which is what makes the refusal about the ions"
        );
        assert_eq!(
            select_dataset(&mixture),
            Selection::Legacy(LegacyReason::NoIons)
        );
    }

    /// A brine the catalogue covers takes it; one it does not falls back, and says which
    /// row was missing.
    #[test]
    fn the_coverage_rule_decides_the_dataset() {
        assert_eq!(
            select_dataset(&species(&[("Na+", 1.0), ("Cl-", -1.0)], &[])),
            Selection::Phreeqc
        );

        // **Hydrogen carbonate is the measured case**: the catalogue has `B0` and `B1` for
        // Na+/HCO3- and no `C0`, so the whole dataset is abandoned rather than completed
        // from the other one.
        assert!(find(Family::B0, &["Na+", "HCO3-"]).is_some());
        assert!(find(Family::B1, &["Na+", "HCO3-"]).is_some());
        assert!(find(Family::C0, &["Na+", "HCO3-"]).is_none());
        assert_eq!(
            select_dataset(&species(&[("Na+", 1.0), ("HCO3-", -1.0)], &[])),
            Selection::Legacy(LegacyReason::Uncovered("C0 for Na+, HCO3-".to_string()))
        );
    }

    /// A 2:2 pair takes the catalogue's optional `B2` without it being required.
    #[test]
    fn an_optional_family_does_not_gate_the_selection() {
        // `Mg++/SO4--` has a `B2`; the selection does not ask for one either way.
        assert!(find(Family::B2, &["Mg++", "SO4--"]).is_some());
        let selected = select_dataset(&species(&[("Mg++", 2.0), ("SO4--", -2.0)], &[]));
        assert_eq!(selected, Selection::Phreeqc);
    }

    /// The hydrocarbon exclusion is the formula test, which is what it exists for.
    #[test]
    fn a_hydrocarbon_is_excluded_from_the_neutral_topology() {
        let methane = Species {
            name: "methane",
            moles: 0.1,
            charge: 0.0,
            formula: "CH4",
            hydrocarbon: false,
        };
        assert!(
            is_hydrocarbon(&methane),
            "a C-H formula is a hydrocarbon even when the type flag is unset, which is \
             the case the method exists for"
        );
        // `CO2` has an oxygen, so it is an active neutral - and the catalogue covers it.
        let carbon_dioxide = Species {
            name: "CO2",
            moles: 0.1,
            charge: 0.0,
            formula: "CO2",
            hydrocarbon: false,
        };
        assert!(!is_hydrocarbon(&carbon_dioxide));
        // And the flag is the other clause.
        let typed = Species {
            name: "default",
            moles: 0.1,
            charge: 0.0,
            formula: "",
            hydrocarbon: true,
        };
        assert!(is_hydrocarbon(&typed));
    }

    /// A component below the activity threshold is not in the topology.
    #[test]
    fn a_trace_component_does_not_enter_the_topology() {
        let mut mixture = species(&[("Na+", 1.0), ("Cl-", -1.0)], &[]);
        mixture.push(Species {
            name: "Li+",
            moles: 1.0e-25,
            charge: 1.0,
            formula: "",
            hydrocarbon: false,
        });
        // Lithium has no catalogue rows at all, so if the trace entry counted this would
        // fall back.
        assert_eq!(select_dataset(&mixture), Selection::Phreeqc);
    }

    /// The databank's own rows resolve to the catalogue's species, which is what makes
    /// the two tables addressable from each other.
    #[test]
    fn the_databanks_ions_are_the_catalogues_species() {
        for name in ["na+", "cl-", "ca++", "so4--", "hco3-", "mg++", "k+"] {
            let entry = databank::entry(name, None).expect("the databank has it");
            assert_eq!(
                canonical_species(&entry.name),
                entry.name,
                "{name} is already in the canonical spelling"
            );
        }
    }
}
