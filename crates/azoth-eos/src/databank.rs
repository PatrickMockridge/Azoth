//! The component databank: a mixture built from substance names.
//!
//! Until this existed, every calculation that needed a mixture took its critical
//! constants, acentric factors, heat-capacity coefficients and interaction parameters
//! as **parallel vectors supplied by the caller** - nine of them, retyped into every
//! spec file. A number in a spec file is a number nobody can check against anything,
//! and the whole point of vendoring NeqSim's tables is that it does not have to be.
//!
//! **Where the data comes from.** `data/components/components.csv` and `kij.csv` are
//! generated from NeqSim's `COMP.csv` and `INTER.csv` by `tools/gen_databank.py`, which
//! carries every column it can and records the unit of each in `databank/manifest.yaml`.
//! They are embedded with `include_str!` rather than read at runtime, which is the
//! arrangement `azoth-hydraulics` already uses for the fittings registry and for the
//! same reason: a wheel that had to find a path would find a different file or none.
//! The Python side opens the same file from the repository, and
//! `python/tests/test_data_agreement.py` compares the embedded bytes against it - so
//! the two implementations read one file rather than two copies that are supposed to
//! match.
//!
//! **What a keycard can change.** The embedded tables are the baseline. A user's
//! keycard overrides them by name, parameter by parameter: a card naming only `omega`
//! keeps the shipped `Tc` and `Pc`, and a name the databank does not have is added -
//! with no heat-capacity coefficients, because a card supplies the parameters a cubic
//! needs and a polynomial is not one of them.
//!
//! An overlay is a value a caller passes, and not a store: nothing here holds one, so
//! two cards in one process are two calls and neither answer depends on what was
//! passed before it. **This crate never parses a card.** A keycard is YAML, the
//! workspace takes no YAML dependency for the reason `tools/gen_registry.py` gives
//! about the spec tree, and `python/src/azoth/keycard.py` reads and validates the file -
//! an [`Overlay`] is built from the values it resolved to.

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

/// Repo-relative path of the component table, which is how Python addresses the same
/// file. A constant rather than a string restated at the call site for the reason the
/// whole databank exists: two copies of a path can disagree.
pub const COMPONENTS_PATH: &str = "data/components/components.csv";

/// Repo-relative path of the interaction table.
pub const KIJ_PATH: &str = "data/components/kij.csv";

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
    /// The `Cp` polynomial's five coefficients, in J/(mol*K**n).
    ///
    /// `None` for a substance an overlay *added*, and that is the whole of why this is
    /// optional: a keycard supplies the parameters a **cubic** needs, and a heat-capacity
    /// polynomial is not one of them. A mixture needs one, so `mixture_of` refuses such
    /// a name rather than defaulting to zeros - which would be a zero heat capacity
    /// wearing the shape of a polynomial.
    pub cp: Option<[f64; 5]>,
}

impl Entry {
    /// The cubic's record for this substance.
    pub fn component(&self) -> Result<Component> {
        Component::new(kelvins(self.tc), pascals(self.pc), self.omega)
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
}

impl ComponentOverride {
    /// Whether this override names every parameter a cubic needs.
    ///
    /// Only asked of a substance the table does not have: one it *does* have is
    /// completed from the table rather than required to be complete.
    #[must_use]
    pub fn is_complete(&self) -> bool {
        self.tc.is_some() && self.pc.is_some() && self.omega.is_some()
    }
}

/// A keycard's data, as a value a caller passes.
///
/// **Built from values, never from a file.** A keycard is YAML, this workspace takes no
/// YAML dependency, and this crate does not read one: `python/src/azoth/keycard.py`
/// reads and validates the file, and an overlay is built from what it resolved to. A
/// Rust caller builds one directly.
///
/// Nothing holds one. Two overlays in one process are two calls, and neither answer
/// depends on what was passed before it - which is the property a module-level card
/// cannot have, whoever set it.
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

/// The two parsed tables: substances by name, and interaction parameters by pair.
type Tables = (HashMap<String, Entry>, HashMap<(String, String), f64>);

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
/// is fallible in principle. It is not here for a keycard: a keycard is YAML, this
/// crate never reads one, and an overlay is built from values rather than parsed.
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
            },
        );
    }
    Ok(out)
}

fn parse_kij() -> Result<HashMap<(String, String), f64>> {
    let mut reader = csv::ReaderBuilder::new()
        .has_headers(true)
        .from_reader(KIJ_CSV.as_bytes());
    let header = reader.headers().map_err(csv_failure)?.clone();
    let mut index = HashMap::new();
    for name in ["component_a", "component_b", "kij_pr"] {
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
        // Stored both ways, so a caller need not know which name came first.
        out.insert((a.clone(), b.clone()), value);
        out.insert((b, a), value);
    }
    Ok(out)
}

/// One substance's constants, or a failure naming it.
///
/// `overlay` is the card this call reads, and `None` means the data this crate ships.
/// The two are one path rather than two: an override is applied *here*, so every caller
/// resolves a name the same way whether a card is in play or not.
///
/// Returns an owned [`Entry`] rather than a `&'static` one, and that is forced rather
/// than chosen: a substance an overlay *adds* is in no table, so there is nothing static
/// to borrow. The clone is a `String` and six numbers, once per component per mixture.
///
/// # Errors
/// * [`AzothError::PropertyUnavailable`] if neither the table nor the overlay carries
///   the substance, or if the overlay adds one without every parameter a cubic reads.
///   Both are refused rather than approximated: a mixture silently missing a component,
///   or holding one completed from a similar substance, is a wrong answer with every
///   symptom of a right one.
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
                // A card supplies the parameters a cubic needs, and a polynomial is not
                // one of them.
                cp: None,
            })
        }
        (Some(base), Some(over)) => Ok(Entry {
            // Parameter by parameter: what the overlay names, else what ships. A card
            // overriding one value does not restate, and does not lose, the others.
            tc: over.tc.unwrap_or(base.tc),
            pc: over.pc.unwrap_or(base.pc),
            omega: over.omega.unwrap_or(base.omega),
            cp: base.cp,
            name: base.name,
        }),
    }
}

/// The binary interaction parameter for a pair, or zero.
///
/// Zero rather than a failure: an absent pair is the ideal-mixture default, which is
/// what NeqSim's own reader substitutes. A pair either source *does* carry is never
/// silently ignored.
///
/// **An overlay's zero wins over a fitted value**, which is the opposite of how an
/// absent pair reads, and deliberately so: overriding a fitted pair back to ideal
/// mixing is a caller stating something, and treating that zero as "no opinion" would
/// undo it with no symptom.
#[must_use]
pub fn kij(first: &str, second: &str, overlay: Option<&Overlay>) -> f64 {
    if let Some(value) = overlay.and_then(|o| o.kij(first, second)) {
        return value;
    }
    tables()
        .1
        .get(&(first.trim().to_lowercase(), second.trim().to_lowercase()))
        .copied()
        .unwrap_or(0.0)
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
        .map(|((a, b), value)| (a.clone(), b.clone(), *value))
        .collect();
    out.sort_by(|left, right| (&left.0, &left.1).cmp(&(&right.0, &right.1)));
    out
}

/// A mixture and its ideal-gas model, built from substance names.
///
/// The two come back together because they are one object in practice: a mixture
/// without the heat-capacity coefficients cannot produce an enthalpy, and building them
/// from two separate lookups invites a call that names different components in each.
///
/// `overlay` is the card this call reads, and `None` means the data this crate ships.
///
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty.
/// * [`AzothError::PropertyUnavailable`] if a name is in neither source, if an overlay
///   adds one without every parameter a cubic reads, or if one has no heat-capacity
///   coefficients - which is what an overlay-added substance is, and why an enthalpy
///   cannot be produced from one.
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
