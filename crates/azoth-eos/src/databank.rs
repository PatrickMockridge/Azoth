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
    pub cp: [f64; 5],
}

impl Entry {
    /// The cubic's record for this substance.
    pub fn component(&self) -> Result<Component> {
        Component::new(kelvins(self.tc), pascals(self.pc), self.omega)
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
                cp,
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
/// # Errors
/// * [`AzothError::InvalidInput`] if the databank carries no such substance. A name it
///   does not have is refused rather than approximated: a mixture silently missing a
///   component is a wrong answer with every symptom of a right one.
pub fn entry(name: &str) -> Result<&'static Entry> {
    let key = name.trim().to_lowercase();
    tables()
        .0
        .get(&key)
        .ok_or_else(|| AzothError::InvalidInput {
            field: "components".to_string(),
            reason: format!(
                "the databank has no component `{key}`. Names come from NeqSim's COMP.csv, \
                 carried in data/components/components.csv; a keycard adds one by name"
            ),
        })
}

/// The binary interaction parameter for a pair, or zero.
///
/// Zero rather than a failure: an absent pair is the ideal-mixture default, which is
/// what NeqSim's own reader substitutes. A pair the databank *does* carry is never
/// silently ignored.
pub fn kij(first: &str, second: &str) -> f64 {
    tables()
        .1
        .get(&(first.trim().to_lowercase(), second.trim().to_lowercase()))
        .copied()
        .unwrap_or(0.0)
}

/// Every substance name in the databank, in file order.
pub fn names() -> Vec<String> {
    let mut out: Vec<String> = tables().0.keys().cloned().collect();
    out.sort();
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
/// # Errors
/// * [`AzothError::InvalidInput`] if `names` is empty or one of them is not in the
///   databank.
/// * Propagates [`Component::new`]'s range checks.
pub fn mixture_of(names: &[&str]) -> Result<(Mixture, IdealGasModel)> {
    if names.is_empty() {
        return Err(AzothError::invalid_input(
            "components",
            "a mixture needs at least one component",
        ));
    }
    let entries: Vec<&Entry> = names
        .iter()
        .map(|name| entry(name))
        .collect::<Result<_>>()?;

    let components = entries
        .iter()
        .map(|e| e.component())
        .collect::<Result<Vec<_>>>()?;

    // Flattened row-major, symmetric with a zero diagonal - the shape `Mixture::new`
    // validates. The diagonal is zero because a component does not interact with
    // itself, and the databank has no self-pair to look up.
    let n = entries.len();
    let mut matrix = vec![0.0; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let value = kij(&entries[i].name, &entries[j].name);
            matrix[i * n + j] = value;
            matrix[j * n + i] = value;
        }
    }

    let ideal_gas = IdealGasModel {
        cp_a: entries.iter().map(|e| e.cp[0]).collect(),
        cp_b: entries.iter().map(|e| e.cp[1]).collect(),
        cp_c: entries.iter().map(|e| e.cp[2]).collect(),
        cp_d: entries.iter().map(|e| e.cp[3]).collect(),
        cp_e: entries.iter().map(|e| e.cp[4]).collect(),
    };

    Ok((Mixture::new(components, matrix)?, ideal_gas))
}
