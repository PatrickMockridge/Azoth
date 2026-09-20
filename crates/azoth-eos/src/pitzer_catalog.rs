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
//! # What the chosen dataset covers
//!
//! [`coverage`] is the second question, and a different one: not which dataset *would*
//! apply but whether the one that loaded defines every interaction the topology needs.
//! **The legacy CSV carries no same-sign rows at all**, so its theta and psi sets are empty
//! and every mixed topology it is asked to cover is incomplete.
//!
//! Measured, and the reason the audit is not merely a diagnostic:
//! `water + Na+ + K+ + Cl- + HCO3-` falls back to the CSV and then **throws from `init(1)`**,
//! because `validateParameterCoverageOncePerState` refuses on the first
//! `getExcessGibbsEnergy` after a level-zero init. A single cation and a single anion are
//! exempt from that refusal, so `water + NH4+ + Cl-` initializes and evaluates an absent
//! pair at **zero** - which is what [`require_complete`] exists to make a caller able to
//! refuse instead.
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

use azoth_core::{AzothError, Result};

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

/// The catalogue's lookup key: the canonical species **sorted**, **folded to lower case**,
/// and joined with `|`.
///
/// Sorted, so a row is found whichever order its species are written in - the catalogue
/// writes `Ba+2 Cl-` and `Cl- H+`, and both must be found by a caller holding the names in
/// either order.
///
/// **Folded, because the catalogue's namespace is NeqSim's and this library's is not.**
/// The catalogue writes `Na+` and `HCO3-`, which are `getComponentName()`'s spellings,
/// while the component table writes `na+` and `hco3-`. Case is the only thing that differs
/// between the two, and resolving names without regard to case is this library's rule, so
/// folding here is what makes a databank name addressable against the catalogue at all.
#[must_use]
pub fn species_key(species: &[&str]) -> String {
    let mut canonical: Vec<String> = species
        .iter()
        .map(|name| canonical_species(name).to_lowercase())
        .collect();
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
///
/// `PhasePitzer`'s own test, written inline in each place that needs it; named here because
/// the model classifies its components by the same rule the selection does.
pub const ACTIVE_CHARGE: f64 = 0.5;

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

/// The identity NeqSim gives the legacy dataset: `PhasePitzer.DEFAULT_PARAMETER_DATASET_ID`.
pub const LEGACY_DATASET_ID: &str = "neqsim-legacy-pitzer-parameters-v1";

/// The identity NeqSim gives the catalogue:
/// `PitzerParameterDatasets.PHREEQC_PITZER_CATALOG_ID`.
///
/// It carries the commit the vendored file came from, so it is a literal rather than a
/// composed string - the hash is part of the identity, not a path.
pub const PHREEQC_DATASET_ID: &str =
    "usgs-phreeqc-pitzer-b0b3be767158ccc3322d2c816625cf470045e67e-catalog-v1";

impl Selection {
    /// The identity of the dataset this selection loads.
    #[must_use]
    pub fn dataset_id(&self) -> &'static str {
        match self {
            Selection::Phreeqc => PHREEQC_DATASET_ID,
            Selection::Legacy(_) => LEGACY_DATASET_ID,
        }
    }
}

/// Molality above which an ion is in the audited topology: `ACTIVE_ION_MOLALITY`.
///
/// **A different threshold from the selection rule's `1e-20` moles**, and deliberately:
/// `tryApplyCompletePhreeqcPitzerCatalog` tests an absolute mole count while
/// `activeIonIndexes` tests a *molality*, so the audit's topology is the narrower one. A
/// trace ion can therefore be absent from the audit and still have forced the fallback.
pub const ACTIVE_ION_MOLALITY: f64 = 1.0e-8;

/// The four species the primary-salt audit excludes, by
/// `PhasePitzer.isPrimarySaltCoverageSpecies`.
///
/// They are the acid-base species the reaction solver creates, whose trial molalities are
/// not a stable input-brine topology. [`reaction_coverage`] is the variant that keeps
/// them, and `SystemPitzer` exposes it separately for that reason.
const REACTION_SPECIES: [&str; 4] = ["H3O+", "OH-", "HCO3-", "CO3--"];

/// Whether a component is covered by the primary-salt audit rather than the reaction one.
#[must_use]
pub fn is_primary_salt_species(name: &str) -> bool {
    !REACTION_SPECIES
        .iter()
        .any(|species| species.eq_ignore_ascii_case(name))
}

/// What the loaded dataset fails to cover, in the shape `PitzerParameterCoverage` reports.
#[derive(Debug, Clone, PartialEq)]
pub struct Coverage {
    /// The identity of the dataset the audit ran against.
    pub dataset_id: &'static str,
    /// Active cations, sorted.
    pub active_cations: Vec<String>,
    /// Active anions, sorted.
    pub active_anions: Vec<String>,
    /// Absent cation-anion pairs.
    pub missing_binary: Vec<String>,
    /// Absent same-sign pairs.
    pub missing_theta: Vec<String>,
    /// Absent cation-cation-anion and anion-anion-cation triples.
    pub missing_psi: Vec<String>,
}

impl Coverage {
    /// Whether every interaction the active topology needs is defined.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.missing_binary.is_empty()
            && self.missing_theta.is_empty()
            && self.missing_psi.is_empty()
    }

    /// The diagnostic, byte for byte as `formatDiagnostic` writes it.
    ///
    /// Java's `List.toString`, which is the form a reader of NeqSim's exceptions has seen.
    #[must_use]
    pub fn diagnostic(&self) -> String {
        format!(
            "Pitzer parameter coverage incomplete for dataset '{}': activeCations={}, \
             activeAnions={}, missingBinary={}, missingTheta={}, missingPsi={}",
            self.dataset_id,
            java_list(&self.active_cations),
            java_list(&self.active_anions),
            java_list(&self.missing_binary),
            java_list(&self.missing_theta),
            java_list(&self.missing_psi)
        )
    }
}

/// `requireCompletePitzerParameterCoverage`: the audit's refusal.
///
/// # Errors
/// [`AzothError::InvalidInput`] with [`Coverage::diagnostic`] as the reason, where an
/// interaction the topology needs is absent from the loaded dataset. **A refusal rather
/// than a zero**: an absent same-sign or ternary parameter is not a fitted ideal solution,
/// so evaluating it as one would return a number with no symptom.
///
/// **NeqSim does not always enforce this.** `validateParameterCoverageOncePerState` calls
/// it only when [`has_mixed_primary_salt_topology`] holds, so a single-salt brine with an
/// absent binary pair initializes and evaluates the pair at zero. The audit still reports
/// it incomplete, and a caller wanting NeqSim's behaviour asks for the audit rather than
/// the refusal.
pub fn require_complete(coverage: &Coverage) -> Result<()> {
    if coverage.is_complete() {
        return Ok(());
    }
    Err(AzothError::invalid_input(
        "composition",
        coverage.diagnostic(),
    ))
}

/// Whether more than one active primary-salt cation or anion is present.
///
/// The condition under which NeqSim enforces the audit at all: a single cation and a
/// single anion keep the "established binary behavior", and a mixed topology does not.
#[must_use]
pub fn has_mixed_primary_salt_topology(species: &[Species<'_>], solvent_mass: f64) -> bool {
    let active_moles = ACTIVE_ION_MOLALITY * solvent_mass;
    let (mut cations, mut anions) = (0_usize, 0_usize);
    for one in species {
        let charge = one.charge;
        if charge == 0.0 || !is_primary_salt_species(one.name) || one.moles <= active_moles {
            continue;
        }
        if charge > 0.0 {
            cations += 1;
        } else {
            anions += 1;
        }
        if cations > 1 || anions > 1 {
            return true;
        }
    }
    false
}

/// The coverage audit of a phase's primary-salt ions: `getPitzerParameterCoverage`.
#[must_use]
pub fn coverage(species: &[Species<'_>], selection: &Selection, solvent_mass: f64) -> Coverage {
    audit(species, selection, solvent_mass, false)
}

/// The same audit including the reaction solver's acid-base species:
/// `getPitzerReactionParameterCoverage`.
#[must_use]
pub fn reaction_coverage(
    species: &[Species<'_>],
    selection: &Selection,
    solvent_mass: f64,
) -> Coverage {
    audit(species, selection, solvent_mass, true)
}

fn audit(
    species: &[Species<'_>],
    selection: &Selection,
    solvent_mass: f64,
    include_reaction_species: bool,
) -> Coverage {
    let active_moles = ACTIVE_ION_MOLALITY * solvent_mass;
    let active = |one: &&Species<'_>, positive: bool| {
        let charge = one.charge;
        let wanted = if positive { charge > 0.0 } else { charge < 0.0 };
        wanted
            && (include_reaction_species || is_primary_salt_species(one.name))
            && one.moles > active_moles
    };
    let cations: Vec<&Species<'_>> = species.iter().filter(|s| active(s, true)).collect();
    let anions: Vec<&Species<'_>> = species.iter().filter(|s| active(s, false)).collect();

    let mut missing_binary = Vec::new();
    for cation in &cations {
        for anion in &anions {
            if !binary_is_defined(selection, cation.name, anion.name) {
                missing_binary.push(report_pair_key(cation.name, anion.name));
            }
        }
    }

    let mut missing_theta = Vec::new();
    let mut missing_psi = Vec::new();
    mixed_interactions(
        &cations,
        &anions,
        selection,
        &mut missing_theta,
        &mut missing_psi,
    );
    mixed_interactions(
        &anions,
        &cations,
        selection,
        &mut missing_theta,
        &mut missing_psi,
    );

    Coverage {
        dataset_id: selection.dataset_id(),
        active_cations: sorted(cations.iter().map(|s| s.name.to_string()).collect()),
        active_anions: sorted(anions.iter().map(|s| s.name.to_string()).collect()),
        missing_binary: sorted(missing_binary),
        missing_theta: sorted(missing_theta),
        missing_psi: sorted(missing_psi),
    }
}

/// The absent theta/psi definitions of every same-sign pair against every opposite ion.
fn mixed_interactions(
    same_sign: &[&Species<'_>],
    opposite: &[&Species<'_>],
    selection: &Selection,
    missing_theta: &mut Vec<String>,
    missing_psi: &mut Vec<String>,
) {
    for (first, one) in same_sign.iter().enumerate() {
        for other in same_sign.iter().skip(first + 1) {
            let names = [one.name, other.name];
            if !mixed_is_defined(selection, Family::Theta, &names) {
                missing_theta.push(report_pair_key(one.name, other.name));
            }
            for third in opposite {
                let triple = [one.name, other.name, third.name];
                if !mixed_is_defined(selection, Family::Psi, &triple) {
                    missing_psi.push(report_psi_key(one.name, other.name, third.name));
                }
            }
        }
    }
}

/// Whether the loaded dataset defines a cation-anion pair.
///
/// **Per pair and not per family**, which is NeqSim's own model: the loader calls
/// `setBinaryParameters(i, j, b0, b1, c)` once for a pair and records one key, so the
/// three families are defined together or not at all.
fn binary_is_defined(selection: &Selection, first: &str, second: &str) -> bool {
    match selection {
        Selection::Phreeqc => [Family::B0, Family::B1, Family::C0]
            .iter()
            .all(|family| find(*family, &[first, second]).is_some()),
        Selection::Legacy(_) => crate::databank::pitzer_pair(first, second).is_some(),
    }
}

/// Whether the loaded dataset defines a same-sign or ternary interaction.
///
/// **False for the legacy dataset however the row reads.** Its loader calls no theta or
/// psi setter at all, so `definedThetaPairs` and `definedPsiTuples` stay empty and every
/// mixed topology is incomplete - which is what makes a mixed legacy brine uninitializable
/// rather than merely inaccurate.
fn mixed_is_defined(selection: &Selection, family: Family, names: &[&str]) -> bool {
    matches!(selection, Selection::Phreeqc) && find(family, names).is_some()
}

fn sorted(mut values: Vec<String>) -> Vec<String> {
    values.sort();
    values
}

/// `PitzerParameterCoverage.immutableSortedCopy`'s rendering: Java's `List.toString`.
fn java_list(values: &[String]) -> String {
    format!("[{}]", values.join(", "))
}

/// `PhasePitzer.pairKey`: the two names sorted and joined with `|`.
///
/// **Not [`species_key`]**, although the two agree on two names. This one is a *report*
/// key and never a lookup - the catalogue's own key sorts all three of a psi row's species
/// and this one must not, or the diagnostic would name a triple NeqSim never writes.
fn report_pair_key(first: &str, second: &str) -> String {
    if first <= second {
        format!("{first}|{second}")
    } else {
        format!("{second}|{first}")
    }
}

/// `PhasePitzer.psiKey`: the sorted same-sign pair, then the opposite ion unsorted.
fn report_psi_key(first: &str, second: &str, opposite: &str) -> String {
    format!("{}|{opposite}", report_pair_key(first, second))
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

    /// **The databank's own names address the catalogue's rows**, which is what the fold is
    /// for and was not true before it.
    ///
    /// The catalogue spells a species the way NeqSim's `getComponentName()` does and the
    /// component table spells it lower case, so a case-sensitive key finds nothing here -
    /// and every brine silently falls back to the CSV, where the numbers differ. Asserting
    /// the two spellings resolve to the *same row* is the property; asserting either one
    /// alone is not.
    #[test]
    fn the_databanks_ions_address_the_catalogue() {
        let named = |key: &str| {
            databank::entry(key, None)
                .expect("the databank has it")
                .name
        };
        let sodium = named("na+");
        let chloride = named("cl-");
        let hydrogen_carbonate = named("hco3-");
        assert_eq!(sodium, "na+", "the table's own spelling is lower case");

        let from_the_table = find(Family::B0, &[&sodium, &chloride]);
        assert!(from_the_table.is_some());
        assert_eq!(
            from_the_table,
            find(Family::B0, &["Na+", "Cl-"]),
            "two spellings of one species are one row"
        );
        // And the PHREEQC spelling the file itself carries, for a species whose charge
        // suffix is numeric there.
        assert_eq!(
            find(Family::B0, &["Ba+2", "Cl-"]),
            find(Family::B0, &["ba++", "cl-"])
        );
        // Hydrogen carbonate is the measured fallback: `B0` and `B1` are there and `C0` is
        // not, which is what takes `water + Na+ + HCO3-` to the CSV.
        assert!(find(Family::B0, &[&sodium, &hydrogen_carbonate]).is_some());
        assert!(find(Family::C0, &[&sodium, &hydrogen_carbonate]).is_none());
    }

    /// One mole of mixture's worth, with the solvent mass the audit's threshold scales by.
    fn brine<'a>(entries: &'a [(&'a str, f64, f64)]) -> (Vec<Species<'a>>, f64) {
        let water = 0.9;
        let mut out = vec![Species {
            name: "water",
            moles: water,
            charge: 0.0,
            formula: "H2O",
            hydrocarbon: false,
        }];
        out.extend(entries.iter().map(|(name, charge, moles)| Species {
            name,
            moles: *moles,
            charge: *charge,
            formula: "",
            hydrocarbon: false,
        }));
        (out, water * 0.018_015)
    }

    /// **The catalogue covers a mixed brine**, so the audit accepts it and nothing is
    /// missing in any of the three families. Measured on the real system.
    #[test]
    fn a_mixed_catalogue_brine_is_complete() {
        let (mixture, mass) = brine(&[
            ("Na+", 1.0, 0.01),
            ("K+", 1.0, 0.01),
            ("Cl-", -1.0, 0.01),
            ("SO4--", -2.0, 0.01),
        ]);
        let selection = select_dataset(&mixture);
        assert_eq!(selection, Selection::Phreeqc);
        let audit = coverage(&mixture, &selection, mass);
        assert!(audit.is_complete(), "{audit:?}");
        assert_eq!(audit.active_cations, ["K+", "Na+"]);
        assert_eq!(audit.active_anions, ["Cl-", "SO4--"]);
        assert_eq!(audit.dataset_id, PHREEQC_DATASET_ID);
        assert!(require_complete(&audit).is_ok());
    }

    /// **A pair with no parameters is refused, not evaluated at zero.** This is the
    /// sabotage check: `NH4+`/`Cl-` appears in neither dataset, and the audit must say so
    /// rather than let the model compute a substance as an ideal solution.
    #[test]
    fn a_pair_with_no_parameters_is_refused() {
        let (mixture, mass) = brine(&[("NH4+", 1.0, 0.01), ("Cl-", -1.0, 0.01)]);
        let selection = select_dataset(&mixture);
        assert_eq!(
            selection,
            Selection::Legacy(LegacyReason::Uncovered("B0 for NH4+, Cl-".to_string())),
            "the catalogue has no ammonium rows at all"
        );
        assert!(databank::pitzer_pair("NH4+", "Cl-").is_none());
        assert!(find(Family::B0, &["NH4+", "Cl-"]).is_none());

        let audit = coverage(&mixture, &selection, mass);
        assert!(!audit.is_complete());
        assert_eq!(audit.missing_binary, ["Cl-|NH4+"]);
        assert!(audit.missing_theta.is_empty());
        assert!(audit.missing_psi.is_empty());

        // **The refusal, and its exact wording** - the string is NeqSim's, because a
        // reader comparing a port's failure against NeqSim's exception must not have to
        // tell two formats apart.
        let refused = require_complete(&audit).expect_err("a missing pair is a refusal");
        assert_eq!(
            refused.to_string(),
            "invalid input `composition`: Pitzer parameter coverage incomplete for dataset \
             'neqsim-legacy-pitzer-parameters-v1': activeCations=[NH4+], activeAnions=[Cl-], \
             missingBinary=[Cl-|NH4+], missingTheta=[], missingPsi=[]"
        );

        // **And NeqSim would still initialize it.** A single cation and a single anion is
        // not a mixed topology, and `validateParameterCoverageOncePerState` only enforces
        // the audit when it is - so this brine evaluates the pair at zero. The audit
        // reports it; the phase does not refuse it.
        assert!(!has_mixed_primary_salt_topology(&mixture, mass));
    }

    /// A mixed brine on the legacy dataset cannot pass, because that dataset carries no
    /// same-sign rows at all. Measured: NeqSim throws from `init(1)` here.
    #[test]
    fn a_mixed_legacy_brine_is_incomplete() {
        let (mixture, mass) = brine(&[
            ("Na+", 1.0, 0.01),
            ("K+", 1.0, 0.01),
            ("Cl-", -1.0, 0.01),
            ("HCO3-", -1.0, 0.01),
        ]);
        let selection = select_dataset(&mixture);
        assert!(
            matches!(selection, Selection::Legacy(_)),
            "hydrogen carbonate is in no catalogue family"
        );
        let audit = coverage(&mixture, &selection, mass);

        // **`HCO3-` is absent from `activeAnions`** although it is the larger of the two
        // anions, because the primary-salt audit excludes the reaction species. The
        // topology it audits is `Na+`/`K+`/`Cl-`.
        assert_eq!(audit.active_cations, ["K+", "Na+"]);
        assert_eq!(audit.active_anions, ["Cl-"]);
        assert!(
            audit.missing_binary.is_empty(),
            "the CSV carries both chlorides"
        );
        assert_eq!(audit.missing_theta, ["K+|Na+"]);
        assert_eq!(audit.missing_psi, ["K+|Na+|Cl-"]);
        assert!(has_mixed_primary_salt_topology(&mixture, mass));
        assert_eq!(
            require_complete(&audit).unwrap_err().to_string(),
            "invalid input `composition`: Pitzer parameter coverage incomplete for dataset \
             'neqsim-legacy-pitzer-parameters-v1': activeCations=[K+, Na+], activeAnions=[Cl-], \
             missingBinary=[], missingTheta=[K+|Na+], missingPsi=[K+|Na+|Cl-]"
        );
    }

    /// **The reaction variant is a different observable.** The same brine is complete
    /// under the primary-salt audit and incomplete under the reaction one, which is what
    /// `SystemPitzer` exposes both for.
    #[test]
    fn the_reaction_audit_keeps_the_species_the_primary_audit_drops() {
        let (mixture, mass) = brine(&[
            ("Na+", 1.0, 0.01),
            ("HCO3-", -1.0, 0.01),
            ("CO3--", -2.0, 0.01),
        ]);
        let selection = select_dataset(&mixture);

        let primary = coverage(&mixture, &selection, mass);
        assert!(primary.active_anions.is_empty());
        assert!(primary.is_complete());

        let reaction = reaction_coverage(&mixture, &selection, mass);
        assert_eq!(reaction.active_anions, ["CO3--", "HCO3-"]);
        assert!(!reaction.is_complete());
        assert_eq!(reaction.missing_theta, ["CO3--|HCO3-"]);
        assert_eq!(reaction.missing_psi, ["CO3--|HCO3-|Na+"]);
    }

    /// **The audit's threshold is a molality and the selection's is a mole count**, so a
    /// trace ion can have chosen the dataset and still not be audited. The psi and theta
    /// loops see only what the binary loop saw.
    #[test]
    fn a_trace_ion_leaves_the_audited_topology() {
        let (mixture, mass) = brine(&[
            ("Na+", 1.0, 0.01),
            ("Cl-", -1.0, 0.01),
            ("K+", 1.0, 1.0e-11),
        ]);
        let molality = 1.0e-11 / mass;
        assert!(
            molality < ACTIVE_ION_MOLALITY && molality > 1.0e-20,
            "k+ is under the audit's threshold and over the selection's: {molality:e}"
        );

        // The selection sees it, and the catalogue covers it, so the answer is still the
        // catalogue - the two thresholds agreeing here is a coincidence of the data.
        let selection = select_dataset(&mixture);
        assert_eq!(selection, Selection::Phreeqc);
        let audit = coverage(&mixture, &selection, mass);
        assert_eq!(audit.active_cations, ["Na+"]);
        assert!(!has_mixed_primary_salt_topology(&mixture, mass));
    }

    /// A same-sign pair on the catalogue is covered, and the psi tuple is reported with
    /// **only its first two species sorted** - `psiKey`'s own order, not the catalogue's
    /// all-sorted, folded lookup key.
    #[test]
    fn the_reported_keys_are_neqsims() {
        assert_eq!(report_pair_key("Na+", "K+"), "K+|Na+");
        assert_eq!(report_pair_key("Cl-", "Na+"), "Cl-|Na+");
        assert_eq!(report_psi_key("Na+", "K+", "Cl-"), "K+|Na+|Cl-");
        // The lookup key sorts all three, which is a triple NeqSim's diagnostic never
        // contains - and it folds case, which the reported key must not because the
        // diagnostic has to read the way NeqSim's exception does.
        assert_eq!(species_key(&["Na+", "K+", "Cl-"]), "cl-|k+|na+");
        assert_ne!(
            report_psi_key("Na+", "K+", "Cl-"),
            species_key(&["Na+", "K+", "Cl-"])
        );
    }
}
