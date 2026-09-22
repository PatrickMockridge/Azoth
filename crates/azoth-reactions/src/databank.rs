//! The reaction data: elements, stoichiometry, and the reaction rows themselves.
//!
//! Five compiled tables under `data/reactions/`, written by
//! `tools/gen_reaction_data.py` from the vendored copies of NeqSim's own `element`,
//! `STOCCOEFDATA` and `REACTIONDATA` tables, and embedded with `include_str!` so a
//! wheel carries what it reads.
//!
//! # The element table is what makes a component reactive
//!
//! `Component.getElements()` builds an `Element` from a row of this table, keyed by
//! component name, and the element matrix is built from those rows. **It is not
//! parsed from the name**: a substance the table does not carry has no composition,
//! so it cannot enter the balance and takes no part in a reaction. That is a
//! property of the data rather than of the caller, and it is why the table's 80
//! components bound what a fluid can react at all - MEG, DEG and TEG have no row.
//!
//! # Three sources, and they are not interchangeable
//!
//! `REACTIONDATA`, `REACTIONDATAPITZER` and `REACTIONDATAKENTEISENBERG` hold the same
//! reaction names on different standard states, and the Pitzer source has an extra
//! `ValidationStatus` column the other two do not. `ChemicalReactionDataSource`'s
//! javadoc says the parameters of one source are not valid for another's model, so
//! the source is an input rather than a default that silently differs.

use std::collections::HashMap;
use std::sync::OnceLock;

use azoth_core::AzothError;
use csv::ReaderBuilder;

// **One column is read from the component databank rather than from `data/reactions/`.**
// The electroneutrality row of the element matrix holds each substance's ionic charge, and
// NeqSim reads that from the component's own `COMP.csv` row rather than from any reaction
// table. Vendoring it again under `data/reactions/` would make two copies of one upstream
// column, so the compiled component table is read here instead - one field, by name.
const COMPONENTS_CSV: &str = include_str!("../../../data/components/components.csv");
const ELEMENTS_CSV: &str = include_str!("../../../data/reactions/elements.csv");
const STOICHIOMETRY_CSV: &str = include_str!("../../../data/reactions/stoichiometry.csv");
const STANDARD_CSV: &str = include_str!("../../../data/reactions/REACTIONDATA.csv");
const PITZER_CSV: &str = include_str!("../../../data/reactions/REACTIONDATAPITZER.csv");
const KENT_EISENBERG_CSV: &str =
    include_str!("../../../data/reactions/REACTIONDATAKENTEISENBERG.csv");

/// The compiled paths, for the Python side's agreement test.
pub const ELEMENTS_PATH: &str = "data/reactions/elements.csv";
/// The compiled stoichiometry table's path.
pub const STOICHIOMETRY_PATH: &str = "data/reactions/stoichiometry.csv";
/// The compiled component databank's path, read for the ionic charge alone.
pub const COMPONENTS_PATH: &str = "data/components/components.csv";

/// Which table a reaction's `K` is read from.
///
/// The three members hold the same reactions on different standard states, so this
/// selects a convention and not a fallback chain. `Pitzer` is the only one whose rows
/// carry a validation status, and the molality basis is the convention that goes with it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ReactionDataSource {
    /// `reactiondata`, the general set an electrolyte EOS system uses.
    Standard,
    /// `reactiondatapitzer`, on the molality standard state, with validated rows.
    Pitzer,
    /// `reactiondatakenteisenberg`, apparent equilibrium constants.
    KentEisenberg,
}

impl ReactionDataSource {
    /// The table's own identifier, as `ChemicalReactionDataSource` spells it.
    #[must_use]
    pub fn identifier(self) -> &'static str {
        match self {
            Self::Standard => "standard",
            Self::Pitzer => "pitzer",
            Self::KentEisenberg => "kent-eisenberg",
        }
    }

    /// The compiled text this source reads.
    #[must_use]
    pub fn text(self) -> &'static str {
        match self {
            Self::Standard => STANDARD_CSV,
            Self::Pitzer => PITZER_CSV,
            Self::KentEisenberg => KENT_EISENBERG_CSV,
        }
    }
}

impl std::str::FromStr for ReactionDataSource {
    type Err = AzothError;

    /// **A source this library does not carry is refused**, rather than falling back to
    /// the general set: the three are different standard states, so a fallback would
    /// answer a question the caller did not ask.
    fn from_str(value: &str) -> Result<Self, Self::Err> {
        match value.trim().to_ascii_lowercase().as_str() {
            "standard" => Ok(Self::Standard),
            "pitzer" => Ok(Self::Pitzer),
            "kent-eisenberg" | "kent_eisenberg" => Ok(Self::KentEisenberg),
            other => Err(AzothError::InvalidInput {
                field: "reaction_data_source".to_string(),
                reason: format!("`{other}` is not one of `standard`, `pitzer`, `kent-eisenberg`"),
            }),
        }
    }
}

/// One element's count in one component.
#[derive(Debug, Clone, PartialEq)]
pub struct ElementRow {
    /// The substance, as NeqSim spells it.
    pub component: String,
    /// The element symbol.
    pub element: String,
    /// How many atoms of it the formula carries, which is a count and not a rate.
    pub count: f64,
}

/// One component's stoichiometric coefficient in one reaction.
#[derive(Debug, Clone, PartialEq)]
pub struct StoichiometryRow {
    /// The reaction's name.
    pub reaction: String,
    /// The substance.
    pub component: String,
    /// The coefficient, negative for a reactant and positive for a product.
    pub coefficient: f64,
}

/// One reaction's fitted constants, in the units the table holds them.
#[derive(Debug, Clone, PartialEq)]
pub struct ReactionRow {
    /// The reaction's name, which is what selects it.
    pub name: String,
    /// The four coefficients of `ln K = K1 + K2/T + K3 ln T + K4 T`.
    pub coefficients: [f64; 4],
    /// The reference temperature the fit is anchored at, in K. **The `ln K` expression
    /// does not read it**: it is the temperature the coefficients were fitted at, and
    /// the rate law's own reference is the same column.
    pub reference_temperature: f64,
    /// The pre-exponential rate factor, read by the reference-Arrhenius rate law.
    pub rate_factor: f64,
    /// The activation energy, in the table's own units, read by the same law.
    pub activation_energy: f64,
    /// The literature citation for the fit, which is a per-row provenance string.
    pub reference: String,
    /// Whether NeqSim loads the row at all. A zero here is why the one combustion row
    /// in the table is dormant.
    pub use_reaction: bool,
    /// The Pitzer source's evidence column. `None` for the other two sources, which do
    /// not carry it rather than carrying an empty one.
    pub validation_status: Option<String>,
}

/// Every row of the element table.
///
/// # Errors
/// Fails if the compiled table is malformed, which is a build-time invariant.
pub fn elements() -> Result<&'static [ElementRow], AzothError> {
    static TABLE: OnceLock<Result<Vec<ElementRow>, AzothError>> = OnceLock::new();
    match TABLE.get_or_init(|| {
        let mut out = Vec::new();
        for row in reader(ELEMENTS_CSV)?.records() {
            let row = row.map_err(malformed)?;
            out.push(ElementRow {
                component: column(&row, 1)?.to_string(),
                element: column(&row, 2)?.to_string(),
                count: number(&row, 3)?,
            });
        }
        Ok(out)
    }) {
        Ok(table) => Ok(table),
        Err(error) => Err(error.clone()),
    }
}

/// The elements one component's formula is made of, in the table's own order.
///
/// **`None` where the table has no row for the name**, which is the state that keeps a
/// component out of the element matrix. An empty vector is a different answer and does
/// not occur: every row states a count.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn element_composition(component: &str) -> Result<Option<Vec<(String, f64)>>, AzothError> {
    let mut found: Vec<(String, f64)> = Vec::new();
    for row in elements()? {
        if row.component == component {
            found.push((row.element.clone(), row.count));
        }
    }
    Ok(if found.is_empty() { None } else { Some(found) })
}

/// The reactions a fluid can run: every reaction whose name the table carries and whose
/// `use_reaction` flag is set.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn reactions(source: ReactionDataSource) -> Result<Vec<ReactionRow>, AzothError> {
    let mut out = Vec::new();
    for row in reader(source.text())?.records() {
        let row = row.map_err(malformed)?;
        let use_reaction = number(&row, 11)? != 0.0;
        out.push(ReactionRow {
            name: column(&row, 2)?.to_string(),
            coefficients: [
                number(&row, 3)?,
                number(&row, 4)?,
                number(&row, 5)?,
                number(&row, 6)?,
            ],
            reference_temperature: number(&row, 7)?,
            rate_factor: number(&row, 8)?,
            activation_energy: number(&row, 9)?,
            reference: column(&row, 10)?.to_string(),
            use_reaction,
            // The Pitzer table's thirteenth column; the other two tables stop at twelve.
            validation_status: row.get(12).map(str::to_string),
        });
    }
    Ok(out)
}

/// One reaction by name, or `None` where the source does not carry it.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn reaction(source: ReactionDataSource, name: &str) -> Result<Option<ReactionRow>, AzothError> {
    Ok(reactions(source)?.into_iter().find(|row| row.name == name))
}

/// Every reaction name a source carries, loaded or not.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn reaction_names(source: ReactionDataSource) -> Result<Vec<String>, AzothError> {
    Ok(reactions(source)?.into_iter().map(|row| row.name).collect())
}

/// The stoichiometric coefficients of one reaction, keyed by component.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn stoichiometry(reaction: &str) -> Result<Vec<(String, f64)>, AzothError> {
    let mut out = Vec::new();
    for row in stoichiometry_rows()? {
        if row.reaction == reaction {
            out.push((row.component.clone(), row.coefficient));
        }
    }
    Ok(out)
}

/// Every row of the stoichiometry table.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn stoichiometry_rows() -> Result<&'static [StoichiometryRow], AzothError> {
    static TABLE: OnceLock<Result<Vec<StoichiometryRow>, AzothError>> = OnceLock::new();
    match TABLE.get_or_init(|| {
        let mut out = Vec::new();
        for row in reader(STOICHIOMETRY_CSV)?.records() {
            let row = row.map_err(malformed)?;
            out.push(StoichiometryRow {
                reaction: column(&row, 2)?.to_string(),
                component: column(&row, 3)?.to_string(),
                coefficient: number(&row, 4)?,
            });
        }
        Ok(out)
    }) {
        Ok(table) => Ok(table),
        Err(error) => Err(error.clone()),
    }
}

/// Every component the element table carries, which is the set that can react.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn reactive_components() -> Result<Vec<String>, AzothError> {
    let mut names: Vec<String> = Vec::new();
    for row in elements()? {
        if !names.contains(&row.component) {
            names.push(row.component.clone());
        }
    }
    Ok(names)
}

/// One component, by name, as the element table keys it.
///
/// # Errors
/// Fails if the compiled table is malformed.
pub fn composition_table() -> Result<HashMap<String, Vec<(String, f64)>>, AzothError> {
    let mut out: HashMap<String, Vec<(String, f64)>> = HashMap::new();
    for row in elements()? {
        out.entry(row.component.clone())
            .or_default()
            .push((row.element.clone(), row.count));
    }
    Ok(out)
}

/// One component's ionic charge, in elementary charges.
///
/// **`None` where the component databank has no row for the name.** That is the same
/// kind of absence as a missing element row and is refused by the caller rather than
/// defaulted to zero: a substance treated as neutral when it is charged is a charge row
/// that does not balance, and the answer would still look like a number.
///
/// # Errors
/// Fails if the compiled component table is malformed.
pub fn ionic_charge(component: &str) -> Result<Option<f64>, AzothError> {
    Ok(ionic_charges()?
        .get(&component.trim().to_lowercase())
        .copied())
}

/// Every component's ionic charge, keyed by lowercased name.
///
/// **The two tables spell their names differently and this is where that is settled.**
/// `element.csv` and the reaction tables keep NeqSim's spelling, so the element matrix is
/// built over `CO2` and `H3O+`; the component databank is lowercased by its generator,
/// which is also how `azoth-eos` looks a component up. The column is found by header
/// rather than by position, because the component table carries ~140 of them.
fn ionic_charges() -> Result<&'static HashMap<String, f64>, AzothError> {
    static TABLE: OnceLock<Result<HashMap<String, f64>, AzothError>> = OnceLock::new();
    match TABLE.get_or_init(|| {
        let mut parsed = reader(COMPONENTS_CSV)?;
        let header = parsed.headers().map_err(malformed)?.clone();
        let name_column = header
            .iter()
            .position(|column| column == "name")
            .ok_or_else(|| absent_column("name"))?;
        let charge_column = header
            .iter()
            .position(|column| column == "ioniccharge")
            .ok_or_else(|| absent_column("ioniccharge"))?;

        let mut out = HashMap::new();
        for row in parsed.records() {
            let row = row.map_err(malformed)?;
            out.insert(
                column(&row, name_column)?.trim().to_lowercase(),
                number(&row, charge_column)?,
            );
        }
        Ok(out)
    }) {
        Ok(table) => Ok(table),
        Err(error) => Err(error.clone()),
    }
}

/// One component's standard-state formation properties, in the units the table holds them.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct FormationProperties {
    /// The Gibbs energy of formation, in J/mol.
    pub gibbs_energy_of_formation: f64,
    /// The ideal-gas enthalpy of formation, in J/mol. Read by the reactive flash's
    /// reference potential, and **not** by the ideal-gas enthalpy: NeqSim's `getHID`
    /// multiplies this column by zero (`Component.java:1631`).
    pub enthalpy_of_formation: f64,
    /// The ideal-gas absolute entropy, in J/(mol*K).
    pub absolute_entropy: f64,
}

/// One component's formation properties, or `None` where the databank has no row for it.
///
/// **A zero is returned as a zero rather than read as absent.** The table has no blanks,
/// and a substance in its own standard state has a zero Gibbs energy of formation and a
/// zero enthalpy of formation because that is the definition; `oxygen`, `nitrogen`,
/// `hydrogen` and `argon` are the ones here. `H+` is zero in all three columns, which is
/// the aqueous standard state's own convention. The same columns also carry zeros that
/// *are* missing data, and
/// a value alone cannot tell the two apart: `formic acid`'s Gibbs energy of formation is
/// `0.0` against a fitted `-378700` enthalpy, and `ethylene`'s absolute entropy is `0.0`
/// against a Gibbs energy of formation that is right. NeqSim reads all of them as
/// numbers and so does this, which is why the ambiguity is recorded here instead of being
/// resolved by a rule the data does not support.
///
/// Names are matched as the component databank spells them - lowercased and trimmed -
/// which is how [`ionic_charge`] settles the same difference between the two tables.
///
/// # Errors
/// Fails if the compiled component table is malformed.
pub fn formation_properties(component: &str) -> Result<Option<FormationProperties>, AzothError> {
    Ok(formation_table()?
        .get(&component.trim().to_lowercase())
        .copied())
}

/// Every component's formation properties, keyed by lowercased name.
fn formation_table() -> Result<&'static HashMap<String, FormationProperties>, AzothError> {
    static TABLE: OnceLock<Result<HashMap<String, FormationProperties>, AzothError>> =
        OnceLock::new();
    match TABLE.get_or_init(|| {
        let mut parsed = reader(COMPONENTS_CSV)?;
        let header = parsed.headers().map_err(malformed)?.clone();
        let indexed = |name: &str| -> Result<usize, AzothError> {
            header
                .iter()
                .position(|column| column == name)
                .ok_or_else(|| absent_column(name))
        };
        let name_column = indexed("name")?;
        let gibbs_column = indexed("gibbsenergyofformation")?;
        let enthalpy_column = indexed("enthalpyofformation")?;
        let entropy_column = indexed("absoluteentropy")?;

        let mut out = HashMap::new();
        for row in parsed.records() {
            let row = row.map_err(malformed)?;
            out.insert(
                column(&row, name_column)?.trim().to_lowercase(),
                FormationProperties {
                    gibbs_energy_of_formation: number(&row, gibbs_column)?,
                    enthalpy_of_formation: number(&row, enthalpy_column)?,
                    absolute_entropy: number(&row, entropy_column)?,
                },
            );
        }
        Ok(out)
    }) {
        Ok(table) => Ok(table),
        Err(error) => Err(error.clone()),
    }
}

fn absent_column(name: &str) -> AzothError {
    AzothError::InvalidInput {
        field: "component_data".to_string(),
        reason: format!("the compiled component table has no `{name}` column"),
    }
}

fn reader(text: &'static str) -> Result<csv::Reader<&'static [u8]>, AzothError> {
    Ok(ReaderBuilder::new()
        .has_headers(true)
        .from_reader(text.as_bytes()))
}

fn column(record: &csv::StringRecord, index: usize) -> Result<&str, AzothError> {
    record.get(index).ok_or_else(|| AzothError::InvalidInput {
        field: "reaction_data".to_string(),
        reason: format!(
            "the compiled table has {} column(s), so column {index} is absent",
            record.len()
        ),
    })
}

fn number(record: &csv::StringRecord, index: usize) -> Result<f64, AzothError> {
    let raw = column(record, index)?;
    raw.trim()
        .parse::<f64>()
        .map_err(|_| AzothError::InvalidInput {
            field: "reaction_data".to_string(),
            reason: format!("{raw:?} is not a number"),
        })
}

fn malformed(error: csv::Error) -> AzothError {
    AzothError::InvalidInput {
        field: "reaction_data".to_string(),
        reason: error.to_string(),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_element_table_loads_and_bounds_what_can_react() {
        let names = reactive_components().expect("the table parses");
        assert_eq!(names.len(), 80, "the element table's component count");
        // **The glycols are absent**, which is a property of the data and not of the
        // caller: a P9 inhibitor fluid is not reactive.
        for absent in ["MEG", "DEG", "TEG"] {
            assert!(
                element_composition(absent).expect("parses").is_none(),
                "{absent} should have no element row"
            );
        }
        assert_eq!(
            element_composition("CO2").expect("parses"),
            Some(vec![("C".to_string(), 1.0), ("O".to_string(), 2.0)])
        );
    }

    #[test]
    fn the_formation_properties_are_read_by_name_in_the_databanks_own_spelling() {
        // The element table spells it `CO2` and the reaction tables spell it `CO2`; the
        // component databank spells it `co2`, and the lookup settles that.
        let co2 = formation_properties("CO2")
            .expect("parses")
            .expect("the databank carries it");
        assert_eq!(
            co2,
            FormationProperties {
                gibbs_energy_of_formation: -394_359.0,
                enthalpy_of_formation: -393_509.0,
                absolute_entropy: 213.8,
            }
        );
        assert!(
            formation_properties("no-such-substance")
                .expect("parses")
                .is_none()
        );
    }

    #[test]
    fn a_zero_formation_property_is_returned_as_a_value() {
        // Oxygen's zeros are the standard state's definition, and `H+` is zero in all
        // three columns by the aqueous convention. Both are answers.
        let oxygen = formation_properties("oxygen")
            .expect("parses")
            .expect("the databank carries it");
        assert_eq!(oxygen.gibbs_energy_of_formation, 0.0);
        assert_eq!(oxygen.enthalpy_of_formation, 0.0);
        assert_eq!(oxygen.absolute_entropy, 205.1);

        // And the same column carries a zero that is *not* an answer: formic acid's
        // fitted enthalpy of formation is a number while its Gibbs energy of formation is
        // zero, and nothing in the value says which kind of zero it is.
        let formic_acid = formation_properties("formic acid")
            .expect("parses")
            .expect("the databank carries it");
        assert_eq!(formic_acid.gibbs_energy_of_formation, 0.0);
        assert_eq!(formic_acid.enthalpy_of_formation, -378_700.0);
    }

    #[test]
    fn the_three_sources_carry_the_same_reactions_on_different_constants() {
        let standard = reaction(ReactionDataSource::Standard, "CO2water")
            .expect("parses")
            .expect("the standard source carries it");
        let pitzer = reaction(ReactionDataSource::Pitzer, "CO2water")
            .expect("parses")
            .expect("the Pitzer source carries it");

        // The same reaction, and the constants differ: that is the standard state, and a
        // port that shared one set between the two would be silently wrong.
        assert_ne!(standard.coefficients[0], pitzer.coefficients[0]);
        assert_eq!(standard.validation_status, None, "only Pitzer carries it");
        assert_eq!(pitzer.validation_status.as_deref(), Some("VALIDATED"));
    }

    #[test]
    fn a_source_this_library_does_not_carry_is_refused() {
        assert!(
            "phreeqc".parse::<ReactionDataSource>().is_err(),
            "a fourth source is not one of the three"
        );
    }

    #[test]
    fn the_stoichiometry_is_keyed_by_reaction_name() {
        let rows = stoichiometry("CO2water").expect("parses");
        assert!(!rows.is_empty(), "CO2water has coefficients");
        assert!(
            rows.iter().any(|(component, _)| component == "CO2"),
            "the reaction names CO2: {rows:?}"
        );
    }

    #[test]
    fn a_reaction_name_the_source_does_not_carry_is_absent() {
        assert!(
            reaction(ReactionDataSource::Standard, "not-a-reaction")
                .expect("parses")
                .is_none()
        );
    }
}
