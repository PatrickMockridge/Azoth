//! The fittings registry: equivalent-length coefficients read from
//! `data/fittings/crane_k_factors.csv`.
//!
//! The CSV is embedded with `include_str!` rather than read at runtime, for two
//! reasons. It removes any question of which file was loaded, and it means the
//! Rust and Python implementations necessarily read byte-identical data - the
//! Python side opens the same file, and a test compares the parsed coefficient
//! for every row across the two languages.
//!
//! **Every coefficient in that file is currently an estimated dummy value, not
//! engineering data.** See the file's header. [`Fitting::is_estimated`] reports
//! it, and `crane_k_factors` turns it into a warning on every affected result so
//! a computed pressure drop cannot be mistaken for a design-grade one.

use azoth_core::{AzothError, Result};

use crate::provenance::VerifyStatus;

/// The embedded registry. Path is relative to this source file.
const FITTINGS_CSV: &str = include_str!("../../../data/fittings/crane_k_factors.csv");

/// One row of the registry.
#[derive(Debug, Clone, PartialEq)]
pub struct Fitting {
    /// Stable identifier used by callers and specs.
    pub id: String,
    /// Coarse category: `bend`, `valve`, and so on.
    pub family: String,
    /// Human-readable name.
    pub name: String,
    /// Equivalent length ratio, `L_eq / D`, for fully turbulent flow.
    pub n_ld: f64,
    /// What the coefficient is multiplied by, conventionally `f_t`.
    pub f_t_basis: String,
    /// Where the value came from, or a statement that it came from nowhere.
    pub citation: String,
    /// How far the value can be trusted.
    pub status: VerifyStatus,
    /// The document the value was read from, in a fetchable form.
    ///
    /// `arweave:<txid>` is preferred: an Arweave transaction ID is the hash of
    /// its content, so the document is immutable, independently timestamped, and
    /// fetchable byte-for-byte by anyone. That makes a single number's
    /// provenance auditable rather than a matter of trusting whoever typed it.
    pub source_ref: Option<String>,
    /// Where inside that document to look, e.g. "Table 2, 90 deg elbow".
    pub source_locator: Option<String>,
}

impl Fitting {
    /// True when this coefficient is a placeholder rather than a measurement.
    ///
    /// The same concept the fluid tables carry, under the name that reads right
    /// here: for a coefficient, "estimated" is what a placeholder means.
    #[must_use]
    pub fn is_estimated(&self) -> bool {
        self.status.is_placeholder()
    }
}

/// Parse the embedded registry. Comments and blank lines are skipped; the
/// header block is documentation and provenance, not decoration.
fn parse() -> Result<Vec<Fitting>> {
    let body: String = FITTINGS_CSV
        .lines()
        .filter(|line| !line.trim_start().starts_with('#') && !line.trim().is_empty())
        .collect::<Vec<_>>()
        .join("\n");

    let mut reader = csv::Reader::from_reader(body.as_bytes());
    let mut out = Vec::new();
    for record in reader.records() {
        let record = record.map_err(|e| {
            AzothError::invalid_input("fittings", format!("malformed CSV row: {e}"))
        })?;
        let field = |name: &str| -> Result<&str> {
            record.get(record_field_index(name)).ok_or_else(|| {
                AzothError::invalid_input("fittings", format!("row is missing column `{name}`"))
            })
        };
        let n_ld: f64 = field("n_ld")?
            .trim()
            .parse()
            .map_err(|e| AzothError::invalid_input("n_ld", format!("not a number: {e}")))?;
        out.push(Fitting {
            id: field("fitting_id")?.to_string(),
            family: field("family")?.to_string(),
            name: field("name")?.to_string(),
            n_ld,
            f_t_basis: field("f_t_basis")?.to_string(),
            citation: field("citation")?.to_string(),
            status: VerifyStatus::parse(field("verify_status")?)?,
            source_ref: optional(field("source_ref")?),
            source_locator: optional(field("source_locator")?),
        });
    }
    Ok(out)
}

/// Column order as declared in the CSV header.
const COLUMNS: [&str; 9] = [
    "fitting_id",
    "family",
    "name",
    "n_ld",
    "f_t_basis",
    "citation",
    "verify_status",
    "source_ref",
    "source_locator",
];

/// An empty CSV field means absent, not an empty string.
fn optional(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_string())
}

fn record_field_index(name: &str) -> usize {
    COLUMNS.iter().position(|c| *c == name).unwrap_or(0)
}

/// The registry, parsed once.
///
/// # Errors
/// Returns an error if the embedded CSV is malformed. That is a build-time
/// invariant rather than a user error, and a test asserts it holds, but it is
/// reported rather than panicked on because library code does not panic.
pub fn registry() -> Result<&'static [Fitting]> {
    use std::sync::OnceLock;
    static REGISTRY: OnceLock<Result<Vec<Fitting>>> = OnceLock::new();
    match REGISTRY.get_or_init(parse) {
        Ok(rows) => Ok(rows.as_slice()),
        Err(e) => Err(e.clone()),
    }
}

/// Look up one fitting by id.
///
/// # Errors
/// Returns [`AzothError::UnknownFitting`] if the id is not in the registry.
/// An unknown id is an error rather than a skip: silently treating it as zero
/// loss would under-report pressure drop, which is the dangerous direction to be
/// wrong in.
pub fn find(id: &str) -> Result<&'static Fitting> {
    registry()?
        .iter()
        .find(|f| f.id == id)
        .ok_or_else(|| AzothError::UnknownFitting { id: id.to_string() })
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_embedded_registry_parses() {
        let rows = registry().expect("embedded registry must parse");
        assert!(!rows.is_empty());
    }

    #[test]
    fn every_row_is_well_formed() {
        for row in registry().unwrap() {
            assert!(!row.id.is_empty(), "empty id");
            assert!(row.n_ld > 0.0, "{} has non-positive n_ld", row.id);
            assert!(!row.name.is_empty(), "{} has no name", row.id);
            assert!(!row.citation.is_empty(), "{} has no citation", row.id);
        }
    }

    #[test]
    fn ids_are_unique() {
        let rows = registry().unwrap();
        let unique: std::collections::HashSet<_> = rows.iter().map(|r| &r.id).collect();
        assert_eq!(unique.len(), rows.len(), "duplicate fitting id");
    }

    #[test]
    fn an_estimated_row_says_so_in_its_citation() {
        // The marker and the prose must agree, or a reader skimming the citation
        // would not realise the number is a placeholder. spec_lint enforces this
        // too; this catches it at test time as well.
        for row in registry().unwrap() {
            if row.is_estimated() {
                assert!(
                    row.citation.to_uppercase().contains("DUMMY"),
                    "{} is marked estimated but its citation does not say so",
                    row.id
                );
            }
        }
    }

    #[test]
    fn the_fittings_the_slice_needs_are_present() {
        // The CLI example and the crane spec worked example both rely on these.
        for id in ["90_elbow", "gate_valve_open"] {
            find(id).unwrap_or_else(|e| panic!("{id} missing from registry: {e}"));
        }
    }

    #[test]
    fn unknown_id_is_an_error_not_a_silent_zero() {
        let err = find("no_such_fitting").unwrap_err();
        assert!(matches!(err, AzothError::UnknownFitting { .. }));
        assert!(err.to_string().contains("no_such_fitting"));
    }
}
