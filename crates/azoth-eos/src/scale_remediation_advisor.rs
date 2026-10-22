//! `eos.scale_remediation_advisor` - NeqSim's `ScaleRemediationAdvisor`.
//!
//! A lookup over a vendored pipe-delimited knowledge base: a deposit name is normalised
//! through the class's alias table, its rows are selected (or the `UNKNOWN` fallback is),
//! and they are ordered by descending effectiveness.

use azoth_core::Result;
use azoth_core::units::kelvins;

use crate::results::{RemediationOptionRow, ScaleRemediationAdvisorResult};

/// The knowledge base's repository path, which the Python side reads.
pub const REMEDIATION_PATH: &str = "data/scale/scale_remediation.csv";

/// The embedded knowledge base. Path is relative to this source file.
const REMEDIATION_CSV: &str = include_str!("../../../data/scale/scale_remediation.csv");

/// One row of the knowledge base, in its own units: degrees Celsius.
struct RemediationOption {
    scale_type: String,
    category: String,
    dissolver: String,
    concentration: String,
    method: String,
    effectiveness: String,
    min_temperature_c: f64,
    max_temperature_c: f64,
    cautions: String,
    standard: String,
}

impl RemediationOption {
    fn effectiveness_rank(&self) -> i32 {
        if self.effectiveness.eq_ignore_ascii_case("HIGH") {
            3
        } else if self.effectiveness.eq_ignore_ascii_case("MEDIUM") {
            2
        } else if self.effectiveness.eq_ignore_ascii_case("LOW") {
            1
        } else {
            0
        }
    }
}

/// The class's alias table, lower-case alias to canonical key.
const ALIASES: &[(&str, &str)] = &[
    ("calcite", "CaCO3"),
    ("aragonite", "CaCO3"),
    ("calcium carbonate", "CaCO3"),
    ("siderite", "FeCO3"),
    ("iron carbonate", "FeCO3"),
    ("mackinawite", "FeS"),
    ("pyrrhotite", "FeS"),
    ("troilite", "FeS"),
    ("iron sulfide", "FeS"),
    ("iron sulphide", "FeS"),
    ("barite", "BaSO4"),
    ("barytes", "BaSO4"),
    ("barium sulfate", "BaSO4"),
    ("barium sulphate", "BaSO4"),
    ("celestite", "SrSO4"),
    ("strontium sulfate", "SrSO4"),
    ("strontium sulphate", "SrSO4"),
    ("gypsum", "CaSO4"),
    ("anhydrite", "CaSO4"),
    ("calcium sulfate", "CaSO4"),
    ("calcium sulphate", "CaSO4"),
    ("halite", "NaCl"),
    ("salt", "NaCl"),
    ("dithiazine", "Dithiazine"),
    ("triazine solids", "Dithiazine"),
    ("scavenger solids", "Dithiazine"),
    ("asphaltene", "Asphaltene"),
    ("asphaltenes", "Asphaltene"),
    ("wax", "Wax"),
    ("paraffin", "Wax"),
    ("naphthenate", "Naphthenate"),
    ("calcium naphthenate", "Naphthenate"),
];

/// The candidate remedies for a deposit.
///
/// # Errors
/// This calculation has no range-checked input; the signature returns [`Result`] for the
/// registry's uniform contract.
///
/// # Example
/// ```
/// use azoth_eos::scale_remediation_advisor::scale_remediation_advisor;
///
/// let r = scale_remediation_advisor("calcite")?;
/// assert_eq!(r.canonical_key, "CaCO3");
/// assert_eq!(r.options[0].effectiveness, "HIGH");
/// # Ok::<(), azoth_core::AzothError>(())
/// ```
pub fn scale_remediation_advisor(deposit: &str) -> Result<ScaleRemediationAdvisorResult> {
    let canonical_key = canonical_key(deposit);
    let options = recommend_for(&canonical_key);

    Ok(ScaleRemediationAdvisorResult {
        canonical_key,
        options: options
            .iter()
            .map(|option| RemediationOptionRow {
                scale_type: option.scale_type.clone(),
                category: option.category.clone(),
                dissolver: option.dissolver.clone(),
                concentration: option.concentration.clone(),
                method: option.method.clone(),
                effectiveness: option.effectiveness.clone(),
                min_temperature: kelvins(option.min_temperature_c + 273.15),
                max_temperature: kelvins(option.max_temperature_c + 273.15),
                cautions: option.cautions.clone(),
                standard: option.standard.clone(),
            })
            .collect(),
        warnings: Vec::new(),
    })
}

/// The class's `canonicalKey`: trim, map the empty string to `UNKNOWN`, consult the alias
/// table case-insensitively, and otherwise keep the trimmed input as spelled.
#[must_use]
pub fn canonical_key(deposit: &str) -> String {
    let trimmed = deposit.trim();
    if trimmed.is_empty() {
        return "UNKNOWN".to_owned();
    }
    let lowered = trimmed.to_lowercase();
    ALIASES
        .iter()
        .find(|(alias, _)| *alias == lowered)
        .map_or_else(|| trimmed.to_owned(), |(_, key)| (*key).to_owned())
}

fn recommend_for(key: &str) -> Vec<&'static RemediationOption> {
    let registry = registry();
    let mut matches: Vec<&RemediationOption> = registry
        .iter()
        .filter(|option| option.scale_type.eq_ignore_ascii_case(key))
        .collect();
    if matches.is_empty() {
        matches = registry
            .iter()
            .filter(|option| option.scale_type.eq_ignore_ascii_case("UNKNOWN"))
            .collect();
    }
    // A stable sort, as `List.sort` is: equal ranks keep the CSV's order.
    matches.sort_by_key(|option| std::cmp::Reverse(option.effectiveness_rank()));
    matches
}

fn registry() -> &'static [RemediationOption] {
    use std::sync::OnceLock;
    static REGISTRY: OnceLock<Vec<RemediationOption>> = OnceLock::new();
    REGISTRY.get_or_init(parse)
}

fn parse() -> Vec<RemediationOption> {
    let mut options = Vec::new();
    for line in REMEDIATION_CSV.lines() {
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let parts: Vec<&str> = line.split('|').collect();
        if parts.len() < 10 {
            continue;
        }
        options.push(RemediationOption {
            scale_type: parts[0].trim().to_owned(),
            category: parts[1].trim().to_owned(),
            dissolver: parts[2].trim().to_owned(),
            concentration: parts[3].trim().to_owned(),
            method: parts[4].trim().to_owned(),
            effectiveness: parts[5].trim().to_owned(),
            min_temperature_c: parse_double(parts[6].trim(), 0.0),
            max_temperature_c: parse_double(parts[7].trim(), 200.0),
            cautions: parts[8].trim().to_owned(),
            standard: parts[9].trim().to_owned(),
        });
    }
    options
}

fn parse_double(value: &str, fallback: f64) -> f64 {
    value.parse().unwrap_or(fallback)
}

/// The path the Python side opens, so the data-agreement test compares the two sides.
#[must_use]
pub fn embedded_path() -> &'static str {
    REMEDIATION_PATH
}

/// The knowledge base, so the data-agreement test can compare both sides row by row.
#[must_use]
pub fn embedded_csv() -> &'static str {
    REMEDIATION_CSV
}

/// The parsed registry, as `(scale_type, category, effectiveness, min C, max C)`.
#[must_use]
pub fn registry_rows() -> Vec<(String, String, String, f64, f64)> {
    registry()
        .iter()
        .map(|option| {
            (
                option.scale_type.clone(),
                option.category.clone(),
                option.effectiveness.clone(),
                option.min_temperature_c,
                option.max_temperature_c,
            )
        })
        .collect()
}
