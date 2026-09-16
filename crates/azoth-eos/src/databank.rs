//! The component databank: a mixture built from substance names.
//!
//! `data/components/components.csv` and `kij.csv` are generated from NeqSim's `COMP.csv`
//! and `INTER.csv` by `tools/gen_databank.py`, embedded with `include_str!` so a wheel
//! cannot find a different file or none. The Python side opens the same file from the
//! repository and `python/tests/test_data_agreement.py` compares the embedded bytes
//! against it.
//!
//! A keycard overrides the embedded tables by name, parameter by parameter: a card naming
//! only `omega` keeps the shipped `Tc` and `Pc`, and a name the databank does not have is
//! added - with no heat-capacity coefficients, because a card supplies the parameters a
//! cubic needs. An overlay is a value a caller passes and not a store; [`crate::card`]
//! reads a file and produces one.

use std::collections::HashMap;
use std::sync::OnceLock;

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, Result};

use crate::mixture::{Component, Mixture};
use crate::molar_enthalpy_entropy::IdealGasModel;

/// The compiled component table, generated from NeqSim's `COMP.csv`.
const COMPONENTS_CSV: &str = include_str!("../../../data/components/components.csv");

/// The compiled interaction table, generated from NeqSim's `INTER.csv`.
const KIJ_CSV: &str = include_str!("../../../data/components/kij.csv");

/// The compiled UNIFAC tables, generated from NeqSim's `UNIFACcomp.csv`,
/// `UNIFACGroupParam.csv` and `UNIFACInterParam.csv`.
const UNIFAC_COMP_CSV: &str = include_str!("../../../data/components/UNIFACcomp.csv");
const UNIFAC_GROUP_CSV: &str = include_str!("../../../data/components/UNIFACGroupParam.csv");
const UNIFAC_INTER_CSV: &str = include_str!("../../../data/components/UNIFACInterParam.csv");

/// Repo-relative path of the component table, which is how Python addresses the same
/// file. A constant rather than a string restated at the call site for the reason the
/// whole databank exists: two copies of a path can disagree.
pub const COMPONENTS_PATH: &str = "data/components/components.csv";

/// Repo-relative path of the interaction table.
pub const KIJ_PATH: &str = "data/components/kij.csv";

/// Repo-relative paths of the three UNIFAC tables.
pub const UNIFAC_COMP_PATH: &str = "data/components/UNIFACcomp.csv";
pub const UNIFAC_GROUP_PATH: &str = "data/components/UNIFACGroupParam.csv";
pub const UNIFAC_INTER_PATH: &str = "data/components/UNIFACInterParam.csv";

/// The exact bytes this build embedded for the component table.
///
/// Exposed so the Python side can compare bytes rather than parsed values - parsed
/// values agree across two *different* files, which is what a stale copy bundled into
/// a wheel looks like. `python/tests/test_data_agreement.py` is the caller.
#[must_use]
pub fn embedded_components() -> &'static str {
    COMPONENTS_CSV
}

/// The exact bytes this build embedded for the interaction table.
#[must_use]
pub fn embedded_kij() -> &'static str {
    KIJ_CSV
}

/// The exact bytes this build embedded for each UNIFAC table, so the Python side can
/// compare bytes rather than parsed values.
#[must_use]
pub fn embedded_unifac_comp() -> &'static str {
    UNIFAC_COMP_CSV
}

#[must_use]
pub fn embedded_unifac_group() -> &'static str {
    UNIFAC_GROUP_CSV
}

#[must_use]
pub fn embedded_unifac_inter() -> &'static str {
    UNIFAC_INTER_CSV
}

/// One substance's constants, in the units the compiled table holds them in.
#[derive(Debug, Clone, PartialEq)]
pub struct Entry {
    /// The name it is looked up by, lower case.
    pub name: String,
    /// Critical temperature, in K.
    pub tc: f64,
    /// Critical pressure, in Pa.
    pub pc: f64,
    /// Acentric factor, dimensionless.
    pub omega: f64,
    /// The `Cp` polynomial's five coefficients, in J/(mol*K**n). `None` for a substance
    /// an overlay added: a card supplies the parameters a cubic needs. `mixture_of`
    /// refuses such a name rather than defaulting to zeros.
    pub cp: Option<[f64; 5]>,
    /// Molar mass, in kg/mol.
    pub molar_mass: Option<f64>,
    /// Critical molar volume, in m³/mol.
    pub critical_volume: Option<f64>,
    /// Dipole moment, in debye.
    pub dipole: Option<f64>,
}

impl Entry {
    /// The cubic's record for this substance.
    pub fn component(&self) -> Result<Component> {
        Ok(Component::new(kelvins(self.tc), pascals(self.pc), self.omega)?
            .with_molar_mass(self.molar_mass))
    }
}

/// One substance as an overlay states it: each parameter named, or none.
///
/// Every field is optional, and that is the rule a keycard follows rather than a
/// convenience: a card naming only `omega` keeps the shipped `Tc` and `Pc`. A record
/// that replaced the shipped one whole would make a user correcting one value restate
/// the others, and lose them silently if they did not.
///
/// **Three fields, and that is a closed list.** It widens when the component model
/// does, and `specs/schema/component.schema.json` is where that is declared. Until
/// then this is the second place after `keycard.COMPONENT_PARAMETERS` that says which
/// parameters a *cubic* reads, and the two are held together by a test.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct ComponentOverride {
    /// Critical temperature, in K.
    pub tc: Option<f64>,
    /// Critical pressure, in Pa.
    pub pc: Option<f64>,
    /// Acentric factor, dimensionless.
    pub omega: Option<f64>,
    /// The five ideal-gas heat-capacity coefficients, all or none.
    pub cp: Option<[f64; 5]>,
}

impl ComponentOverride {
    /// Whether this override names every parameter a cubic needs. Only asked of a
    /// substance the table does not have; one it does have is completed from the table.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.tc.is_some() && self.pc.is_some() && self.omega.is_some()
    }
}

/// A keycard's data, as a value a caller passes. [`crate::card::Card::overlay`] is where
/// a card file becomes one; a caller with the values in hand builds one with
/// [`Overlay::new`], so a carded lookup needs no file. Nothing holds one.
#[derive(Debug, Clone, Default, PartialEq)]
pub struct Overlay {
    /// Overrides by lower-cased name, and additions the table does not have.
    components: HashMap<String, ComponentOverride>,
    /// Overrides by lower-cased pair, stored both ways round.
    kij: HashMap<(String, String), f64>,
}

impl Overlay {
    /// An overlay that changes nothing.
    #[must_use]
    pub fn new() -> Self {
        Self::default()
    }

    /// Override one substance's parameters, or add a substance by name.
    pub fn set_component(&mut self, name: &str, parameters: ComponentOverride) -> &mut Self {
        self.components
            .insert(name.trim().to_lowercase(), parameters);
        self
    }

    /// Override one pair's interaction parameter.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if both names are the same substance. `Mixture::new`
    ///   refuses a non-zero diagonal too, but only for a pair whose *both* names are in
    ///   the mixture being built - so a self-pair on a name nothing builds would be stored
    ///   and read by nothing, which is the failure this crate refuses everywhere.
    pub fn set_kij(&mut self, first: &str, second: &str, value: f64) -> Result<&mut Self> {
        let a = first.trim().to_lowercase();
        let b = second.trim().to_lowercase();
        if a == b {
            return Err(AzothError::invalid_input(
                "kij",
                format!("`{a}` does not interact with itself"),
            ));
        }
        // Stored both ways, exactly as the table is, so a caller need not know which
        // name came first.
        self.kij.insert((a.clone(), b.clone()), value);
        self.kij.insert((b, a), value);
        Ok(self)
    }

    /// This overlay's statement about a substance, or `None` if it makes none.
    #[must_use]
    pub fn component(&self, name: &str) -> Option<&ComponentOverride> {
        self.components.get(&name.trim().to_lowercase())
    }

    /// This overlay's interaction parameter for a pair, or `None` if it states none.
    ///
    /// `None` and `Some(0.0)` are different answers, and the difference is a caller's
    /// deliberate reset to ideal mixing. Anything that filters a zero here silently
    /// undoes it.
    #[must_use]
    pub fn kij(&self, first: &str, second: &str) -> Option<f64> {
        self.kij
            .get(&(first.trim().to_lowercase(), second.trim().to_lowercase()))
            .copied()
    }

    /// Every pair this overlay states, each once, with the lower name first.
    ///
    /// The de-duplication is the opposite of the *storage*, which keeps both orderings
    /// so a caller need not know which name came first. This is the display direction:
    /// one row per pair.
    #[must_use]
    pub fn kij_pairs(&self) -> Vec<(String, String)> {
        let mut out: Vec<(String, String)> = self
            .kij
            .keys()
            .filter(|(first, second)| first < second)
            .cloned()
            .collect();
        out.sort();
        out
    }

    /// The names this overlay adds to the table, in no particular order.
    #[must_use]
    pub fn component_names(&self) -> Vec<&str> {
        self.components.keys().map(String::as_str).collect()
    }

    /// Whether this overlay changes nothing.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.components.is_empty() && self.kij.is_empty()
    }
}

/// One interaction-pair record from `INTER.csv`, keyed by the ordered pair of names
/// exactly as they appear in the file.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Interaction {
    /// The cubic binary interaction parameter, symmetric.
    kij: f64,
    /// The NRTL non-randomness parameter, symmetric.
    alpha: f64,
    /// The NRTL energy parameter for the ordered pair `(first, second)`: the `g_ij` in
    /// `tau_ij = g_ij / T`, in Kelvin. Directional - `g_ij` differs from `g_ji`.
    gij: f64,
}

/// The two parsed tables: substances by name, and interaction parameters by pair.
type Tables = (
    HashMap<String, Entry>,
    HashMap<(String, String), Interaction>,
);

/// The parsed tables, read once.
fn tables() -> &'static Tables {
    static TABLES: OnceLock<Tables> = OnceLock::new();
    TABLES.get_or_init(|| {
        (
            parse_components().expect("the embedded component table should parse"),
            parse_kij().expect("the embedded interaction table should parse"),
        )
    })
}

/// The column index of `name` in a header row, or a named failure.
fn column(header: &csv::StringRecord, name: &str) -> Result<usize> {
    header
        .iter()
        .position(|field| field == name)
        .ok_or_else(|| AzothError::InvalidInput {
            field: "databank".to_string(),
            reason: format!(
                "the compiled table has no `{name}` column; its header is {}",
                header.iter().collect::<Vec<_>>().join(", ")
            ),
        })
}

/// A field of one record, parsed as a number.
fn number(record: &csv::StringRecord, index: usize, column: &str, row: usize) -> Result<f64> {
    let raw = record.get(index).unwrap_or("").trim();
    raw.parse::<f64>().map_err(|_| AzothError::InvalidInput {
        field: "databank".to_string(),
        reason: format!("row {row}: `{column}` is {raw:?}, which is not a number"),
    })
}

/// A `csv` failure as this crate's error.
///
/// Unreachable in practice - the table is embedded, so a malformed one is a build
/// defect rather than a caller condition - and it is a `Result` anyway, because this
/// crate does not panic on data and because the table is a build artefact whose parse
/// is fallible in principle. It is not here for a keycard: a card is read by
/// [`crate::card`], which reports what is wrong with the document, and an overlay
/// arrives here already resolved.
fn csv_failure(error: csv::Error) -> AzothError {
    AzothError::InvalidInput {
        field: "databank".to_string(),
        reason: format!("the embedded table could not be read: {error}"),
    }
}

fn parse_components() -> Result<HashMap<String, Entry>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(COMPONENTS_CSV.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let mut index = HashMap::new();
    for name in [
        "name",
        "tc_k",
        "pc_pa",
        "acentric_factor",
        "cpa",
        "cpb",
        "cpc",
        "cpd",
        "cpe",
        "molar_mass_kg_per_mol",
        "critical_volume_m3_per_mol",
        "dipole_moment_debye",
    ] {
        index.insert(name, column(&header, name)?);
    }

    let mut out = HashMap::new();
    for (offset, record) in reader.records().enumerate() {
        let record = record.map_err(csv_failure)?;
        let name = record
            .get(index["name"])
            .unwrap_or("")
            .trim()
            .to_lowercase();
        if name.is_empty() {
            continue;
        }
        let row = offset + 2; // one for the header, one to count from one
        let mut cp = [0.0; 5];
        for (slot, key) in cp.iter_mut().zip(["cpa", "cpb", "cpc", "cpd", "cpe"]) {
            *slot = number(&record, index[key], key, row)?;
        }
        out.insert(
            name.clone(),
            Entry {
                name,
                tc: number(&record, index["tc_k"], "tc_k", row)?,
                pc: number(&record, index["pc_pa"], "pc_pa", row)?,
                omega: number(&record, index["acentric_factor"], "acentric_factor", row)?,
                // `Some` for everything the table carries: it holds the polynomial for
                // every row it has, and `mixture_of` refuses a name without one.
                cp: Some(cp),
                molar_mass: Some(number(
                    &record,
                    index["molar_mass_kg_per_mol"],
                    "molar_mass_kg_per_mol",
                    row,
                )?),
                critical_volume: Some(number(
                    &record,
                    index["critical_volume_m3_per_mol"],
                    "critical_volume_m3_per_mol",
                    row,
                )?),
                dipole: Some(number(
                    &record,
                    index["dipole_moment_debye"],
                    "dipole_moment_debye",
                    row,
                )?),
            },
        );
    }
    Ok(out)
}

fn parse_kij() -> Result<HashMap<(String, String), Interaction>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(KIJ_CSV.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let mut index = HashMap::new();
    for name in [
        "component_a",
        "component_b",
        "kij_pr",
        "nrtlalpha",
        "nrtlgij",
        "nrtlgji",
    ] {
        index.insert(name, column(&header, name)?);
    }

    let mut out = HashMap::new();
    for (offset, record) in reader.records().enumerate() {
        let record = record.map_err(csv_failure)?;
        let a = record
            .get(index["component_a"])
            .unwrap_or("")
            .trim()
            .to_lowercase();
        let b = record
            .get(index["component_b"])
            .unwrap_or("")
            .trim()
            .to_lowercase();
        if a.is_empty() || b.is_empty() {
            continue;
        }
        let raw = record.get(index["kij_pr"]).unwrap_or("").trim();
        if raw.is_empty() {
            continue;
        }
        let value = raw.parse::<f64>().map_err(|_| AzothError::InvalidInput {
            field: "databank".to_string(),
            reason: format!("row {}: `kij_pr` is {raw:?}", offset + 2),
        })?;
        let row = offset + 2;
        let alpha = number(&record, index["nrtlalpha"], "nrtlalpha", row)?;
        let gij = number(&record, index["nrtlgij"], "nrtlgij", row)?;
        let gji = number(&record, index["nrtlgji"], "nrtlgji", row)?;
        // `kij` and `alpha` are symmetric, stored both ways round so a caller need not
        // know which name came first. `gij` is directional, so the reversed key carries
        // the reversed energy.
        out.insert(
            (a.clone(), b.clone()),
            Interaction {
                kij: value,
                alpha,
                gij,
            },
        );
        out.insert(
            (b, a),
            Interaction {
                kij: value,
                alpha,
                gij: gji,
            },
        );
    }
    Ok(out)
}

/// One substance's constants, or a failure naming it.
///
/// `overlay` is the card this call reads; `None` means the data this crate ships. An
/// override is applied here, so every caller resolves a name the same way. Owned rather
/// than `&'static`, because a substance an overlay adds is in no table.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if neither the table nor the overlay carries
///   the substance, or if the overlay adds one without every parameter a cubic reads.
pub fn entry(name: &str, overlay: Option<&Overlay>) -> Result<Entry> {
    let key = name.trim().to_lowercase();
    let base = tables().0.get(&key).cloned();
    let over = overlay.and_then(|o| o.component(&key));

    match (base, over) {
        (None, None) => Err(AzothError::property_unavailable(
            key,
            "critical constants".to_string(),
            "neither the databank nor the keycard has it. Names come from NeqSim's COMP.csv, \
             carried in data/components/components.csv; a keycard adds one by name"
                .to_string(),
        )),
        (Some(base), None) => Ok(base),
        (None, Some(over)) => {
            if !over.is_complete() {
                let mut missing = Vec::new();
                if over.tc.is_none() {
                    missing.push("Tc");
                }
                if over.pc.is_none() {
                    missing.push("Pc");
                }
                if over.omega.is_none() {
                    missing.push("omega");
                }
                return Err(AzothError::property_unavailable(
                    key,
                    "critical constants".to_string(),
                    format!(
                        "the keycard adds it but is missing {missing:?}. A substance the \
                         databank does not have needs every parameter a cubic reads, because \
                         completing it from a similar one would be inventing data"
                    ),
                ));
            }
            Ok(Entry {
                name: key,
                // Checked complete above; the defaults are unreachable rather than
                // meaningful, and a zero here would be a critical constant.
                tc: over.tc.unwrap_or_default(),
                pc: over.pc.unwrap_or_default(),
                omega: over.omega.unwrap_or_default(),
                // The card may also supply the polynomial, in which case the substance
                // has an enthalpy path; without it, it is a cubic only.
                cp: over.cp,
                // A card states the parameters a cubic reads; it carries no molar mass,
                // critical volume or dipole, so a card-added substance has none.
                molar_mass: None,
                critical_volume: None,
                dipole: None,
            })
        }
        (Some(base), Some(over)) => Ok(Entry {
            // Parameter by parameter: what the overlay names, else what ships. A card
            // overriding one value does not restate, and does not lose, the others.
            tc: over.tc.unwrap_or(base.tc),
            pc: over.pc.unwrap_or(base.pc),
            omega: over.omega.unwrap_or(base.omega),
            cp: over.cp.or(base.cp),
            molar_mass: base.molar_mass,
            critical_volume: base.critical_volume,
            dipole: base.dipole,
            name: base.name,
        }),
    }
}

/// The binary interaction parameter for a pair, or zero.
///
/// Zero rather than a failure: an absent pair is the ideal-mixture default, which is
/// what NeqSim's own reader substitutes. An overlay's zero wins over a fitted value,
/// because overriding a pair back to ideal mixing is a caller stating something.
#[must_use]
pub fn kij(first: &str, second: &str, overlay: Option<&Overlay>) -> f64 {
    if let Some(value) = overlay.and_then(|o| o.kij(first, second)) {
        return value;
    }
    tables()
        .1
        .get(&(first.trim().to_lowercase(), second.trim().to_lowercase()))
        .map_or(0.0, |interaction| interaction.kij)
}

/// The Wilke-Chang association parameter for a solvent, by name.
///
/// A name resolves through a lowercased lookup against NeqSim's table; a solvent
/// that is not listed - a non-associated one, most hydrocarbons - is 1.0. NeqSim's
/// table mixes cases: the six uppercase keys (`MEG`, `DEG`, `TEG`, `MDEA`, `MEA`,
/// `DEA`) are never matched by the lowercased lookup, so they fall through to 1.0
/// exactly as NeqSim's own `getAssociationParameter` does.
#[must_use]
pub fn wilke_chang_phi(name: &str) -> f64 {
    match name.trim().to_lowercase().as_str() {
        "water" | "h2o" | "d2o" => 2.26,
        "methanol" => 1.9,
        "ethanol" => 1.5,
        "1-propanol" | "2-propanol" => 1.2,
        "1-butanol" | "n-butanol" => 1.0,
        "acetic acid" => 1.3,
        "formic acid" => 1.6,
        _ => 1.0,
    }
}

/// Every substance name available, sorted: the table plus whatever an overlay adds.
#[must_use]
pub fn names(overlay: Option<&Overlay>) -> Vec<String> {
    let mut out: Vec<String> = tables().0.keys().cloned().collect();
    if let Some(overlay) = overlay {
        out.extend(overlay.components.keys().cloned());
    }
    out.sort();
    out.dedup();
    out
}

/// Every substance, ordered by name.
///
/// Ordered rather than in file order because the two implementations must be able to
/// compare them one for one, and a sort is the one ordering both can reproduce without
/// agreeing on how the file happens to be laid out.
#[must_use]
pub fn all_entries() -> Vec<&'static Entry> {
    let mut out: Vec<&Entry> = tables().0.values().collect();
    out.sort_by(|left, right| left.name.cmp(&right.name));
    out
}

/// Every interaction pair the databank carries, ordered, each pair once.
///
/// The table is stored both ways round so a caller need not know which name came first;
/// this un-does that by keeping only the ordering where the first name sorts lower.
#[must_use]
pub fn all_kij() -> Vec<(String, String, f64)> {
    let mut out: Vec<(String, String, f64)> = tables()
        .1
        .iter()
        .filter(|((a, b), _)| a < b)
        .map(|((a, b), interaction)| (a.clone(), b.clone(), interaction.kij))
        .collect();
    out.sort_by(|left, right| (&left.0, &left.1).cmp(&(&right.0, &right.1)));
    out
}

/// The NRTL non-randomness matrix for a list of names, flattened row-major.
///
/// `alpha[i][j] = alpha_ij`, symmetric; the diagonal and any pair the table does not
/// carry are `0.0`, the ideal-mixture default.
#[must_use]
pub fn nrtl_alpha(names: &[&str]) -> Vec<f64> {
    let table = &tables().1;
    let n = names.len();
    let mut out = vec![0.0; n * n];
    for (i, a) in names.iter().enumerate() {
        for (j, b) in names.iter().enumerate() {
            if i != j {
                out[i * n + j] = table
                    .get(&(a.trim().to_lowercase(), b.trim().to_lowercase()))
                    .map_or(0.0, |interaction| interaction.alpha);
            }
        }
    }
    out
}

/// The NRTL energy matrix `Dij` for a list of names, flattened row-major.
///
/// `dij[i][j] = g_ij`, the Kelvin energy in `tau_ij = g_ij / T`; directional, so
/// `dij[i][j]` and `dij[j][i]` differ in general. The diagonal and any absent pair are
/// `0.0`.
#[must_use]
pub fn nrtl_dij(names: &[&str]) -> Vec<f64> {
    let table = &tables().1;
    let n = names.len();
    let mut out = vec![0.0; n * n];
    for (i, a) in names.iter().enumerate() {
        for (j, b) in names.iter().enumerate() {
            if i != j {
                out[i * n + j] = table
                    .get(&(a.trim().to_lowercase(), b.trim().to_lowercase()))
                    .map_or(0.0, |interaction| interaction.gij);
            }
        }
    }
    out
}

/// The resolved UNIFAC inputs for a mixture, each matrix flattened row-major.
///
/// `groups` is `N x G` (one row per component, one column per group), `group_r` and
/// `group_q` are length `G`, and `aij` is `G x G` (`a_{mn}` in Kelvin). `G` is the
/// union of the named components' subgroups, sorted by subgroup number, with absent
/// groups counted zero.
#[derive(Debug, Clone, PartialEq)]
pub struct UnifacParameters {
    /// Per-component group counts, `N x G` row-major.
    pub groups: Vec<f64>,
    /// The volume `R` of each group, length `G`.
    pub group_r: Vec<f64>,
    /// The surface area `Q` of each group, length `G`.
    pub group_q: Vec<f64>,
    /// The main-group interaction matrix, `G x G` row-major, in Kelvin.
    pub aij: Vec<f64>,
}

/// The parsed UNIFAC tables: group constants by subgroup, main-group interactions,
/// and per-component group memberships.
struct UnifacTables {
    /// Subgroup number -> `(R, Q, main group)`.
    group: HashMap<i64, (f64, f64, i64)>,
    /// `(main group, main group)` -> `a_mn`.
    aij: HashMap<(i64, i64), f64>,
    /// Lower-cased name -> `[(subgroup, count)]`.
    members: HashMap<String, Vec<(i64, i64)>>,
}

fn unifac_tables() -> &'static UnifacTables {
    static TABLES: OnceLock<UnifacTables> = OnceLock::new();
    TABLES.get_or_init(|| parse_unifac().expect("the embedded UNIFAC tables should parse"))
}

/// One integer field of a record, parsed strictly.
fn integer(record: &csv::StringRecord, index: usize, column: &str, row: usize) -> Result<i64> {
    let raw = record.get(index).unwrap_or("").trim();
    raw.parse::<i64>().map_err(|_| AzothError::InvalidInput {
        field: "databank".to_string(),
        reason: format!("row {row}: `{column}` is {raw:?}, which is not an integer"),
    })
}

fn parse_unifac() -> Result<UnifacTables> {
    let mut group = HashMap::new();
    {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(UNIFAC_GROUP_CSV.as_bytes());
        let header = reader.headers().map_err(csv_failure)?.clone();
        let mut idx = HashMap::new();
        for name in ["secondary", "volumer", "surfareaq", "main"] {
            idx.insert(name, column(&header, name)?);
        }
        for (offset, record) in reader.records().enumerate() {
            let record = record.map_err(csv_failure)?;
            let row = offset + 2;
            group.insert(
                integer(&record, idx["secondary"], "secondary", row)?,
                (
                    number(&record, idx["volumer"], "volumer", row)?,
                    number(&record, idx["surfareaq"], "surfareaq", row)?,
                    integer(&record, idx["main"], "main", row)?,
                ),
            );
        }
    }

    let mut aij = HashMap::new();
    {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(UNIFAC_INTER_CSV.as_bytes());
        let header = reader.headers().map_err(csv_failure)?.clone();
        let main_idx = column(&header, "maingroup")?;
        let mut n_idx = HashMap::new();
        for n in 1..=64 {
            n_idx.insert(n, column(&header, &format!("n{n}"))?);
        }
        for (offset, record) in reader.records().enumerate() {
            let record = record.map_err(csv_failure)?;
            let row = offset + 2;
            let main = integer(&record, main_idx, "maingroup", row)?;
            for (n, &idx) in &n_idx {
                let raw = record.get(idx).unwrap_or("").trim();
                let value = if raw.is_empty() {
                    0.0
                } else {
                    raw.parse::<f64>().map_err(|_| AzothError::InvalidInput {
                        field: "databank".to_string(),
                        reason: format!("row {row}: `n{n}` is {raw:?}, which is not a number"),
                    })?
                };
                aij.insert((main, *n), value);
            }
        }
    }

    let mut members = HashMap::new();
    {
        let mut reader = csv::ReaderBuilder::new()
            .has_headers(true)
            .from_reader(UNIFAC_COMP_CSV.as_bytes());
        let header = reader.headers().map_err(csv_failure)?.clone();
        let name_idx = column(&header, "name")?;
        let mut sub_idx = HashMap::new();
        for s in 1..=140 {
            sub_idx.insert(s, column(&header, &format!("sub{s}"))?);
        }
        for (offset, record) in reader.records().enumerate() {
            let record = record.map_err(csv_failure)?;
            let name = record.get(name_idx).unwrap_or("").trim().to_lowercase();
            if name.is_empty() {
                continue;
            }
            let mut subs = Vec::new();
            for (s, &idx) in &sub_idx {
                let count = integer(&record, idx, &format!("sub{s}"), offset + 2)?;
                if count > 0 {
                    subs.push((*s, count));
                }
            }
            members.insert(name, subs);
        }
    }

    Ok(UnifacTables {
        group,
        aij,
        members,
    })
}

/// The UNIFAC inputs for a list of names, resolved from the vendored group tables.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty.
/// * [`AzothError::PropertyUnavailable`] if a name has no UNIFAC group assignment.
#[allow(clippy::missing_panics_doc)]
pub fn unifac_parameters(names: &[&str]) -> Result<UnifacParameters> {
    let tables = unifac_tables();
    let n = names.len();
    if n == 0 {
        return Err(AzothError::invalid_input(
            "components",
            "a mixture needs at least one component",
        ));
    }
    let mut union: Vec<i64> = Vec::new();
    let mut resolved: Vec<&Vec<(i64, i64)>> = Vec::with_capacity(n);
    for name in names {
        let key = name.trim().to_lowercase();
        let subs = tables.members.get(&key).ok_or_else(|| {
            AzothError::property_unavailable(
                key,
                "UNIFAC group assignment".to_string(),
                "not in UNIFACcomp.csv; a UNIFAC activity coefficient needs a group \
                 decomposition for every component"
                    .to_string(),
            )
        })?;
        for &(s, _) in subs {
            if !union.contains(&s) {
                union.push(s);
            }
        }
        resolved.push(subs);
    }
    union.sort_unstable();
    let g = union.len();

    let mut group_r = vec![0.0; g];
    let mut group_q = vec![0.0; g];
    let mut aij = vec![0.0; g * g];
    for (k, &s) in union.iter().enumerate() {
        let (r, q, main) = *tables.group.get(&s).expect("a subgroup with a member row");
        group_r[k] = r;
        group_q[k] = q;
        for (m, &t) in union.iter().enumerate() {
            let (_, _, main_t) = *tables.group.get(&t).expect("a subgroup with a member row");
            aij[k * g + m] = tables.aij.get(&(main, main_t)).copied().unwrap_or(0.0);
        }
    }

    let mut groups = vec![0.0; n * g];
    for (i, subs) in resolved.iter().enumerate() {
        for &(s, count) in subs.iter() {
            let k = union
                .iter()
                .position(|&x| x == s)
                .expect("a member subgroup in the union");
            groups[i * g + k] = count as f64;
        }
    }

    Ok(UnifacParameters {
        groups,
        group_r,
        group_q,
        aij,
    })
}

/// A mixture and its ideal-gas model, built from substance names, which come back
/// together because a mixture without heat-capacity coefficients cannot produce an
/// enthalpy.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty.
/// * [`AzothError::PropertyUnavailable`] if a name is in neither source, or if one has
///   no heat-capacity coefficients - which is what an overlay-added substance is.
/// * Propagates [`Component::new`]'s range checks.
pub fn mixture_of(names: &[&str], overlay: Option<&Overlay>) -> Result<(Mixture, IdealGasModel)> {
    if names.is_empty() {
        return Err(AzothError::invalid_input(
            "components",
            "a mixture needs at least one component",
        ));
    }
    let entries: Vec<Entry> = names
        .iter()
        .map(|name| entry(name, overlay))
        .collect::<Result<_>>()?;

    let missing: Vec<&str> = entries
        .iter()
        .filter(|e| e.cp.is_none())
        .map(|e| e.name.as_str())
        .collect();
    if !missing.is_empty() {
        return Err(AzothError::property_unavailable(
            missing.join(", "),
            "heat-capacity coefficients".to_string(),
            "the databank carries them for every substance it ships; one a keycard adds \
             needs its own, because a cubic needs `Tc`, `Pc` and `omega` and an enthalpy \
             needs the polynomial as well"
                .to_string(),
        ));
    }

    let components = entries
        .iter()
        .map(Entry::component)
        .collect::<Result<Vec<_>>>()?;

    // Flattened row-major, symmetric with a zero diagonal - the shape `Mixture::new`
    // validates. The diagonal is zero because a component does not interact with
    // itself, and neither source has a self-pair to look up.
    let n = entries.len();
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let value = kij(&entries[i].name, &entries[j].name, overlay);
            matrix[i * n + j] = value;
            matrix[j * n + i] = value;
        }
    }

    // Every entry has a polynomial: the filter above refuses the mixture otherwise.
    let coefficient = |index: usize| -> Vec<f64> {
        entries
            .iter()
            .map(|e| e.cp.map_or(0.0, |cp| cp[index]))
            .collect()
    };
    let ideal_gas = IdealGasModel {
        cp_a: coefficient(0),
        cp_b: coefficient(1),
        cp_c: coefficient(2),
        cp_d: coefficient(3),
        cp_e: coefficient(4),
    };

    Ok((Mixture::new(components, matrix)?, ideal_gas))
}
