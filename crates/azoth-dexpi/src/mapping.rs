//! The azoth palette <-> DEXPI Process mapping: types and the generated table's readers.
//!
//! `specs/dexpi/process_units.toml` is the one hand-written statement of the mapping;
//! `tools/gen_dexpi_mapping.py` compiles it into [`process_mapping_gen`], and the tables below
//! are how the rest of the crate reads them. Adding a palette entry without a row fails the
//! generator, so a new unit operation cannot leave the exchange silently incomplete.

use crate::process_mapping_gen::{PROCESS_TYPE_CANDIDATES, PROCESS_UNITS};

/// What exchanging a unit as a DEXPI Process step costs.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Disposition {
    /// The type carries what the unit does.
    Faithful,
    /// The type is the nearest reviewed counterpart, and the row's `note` says what is lost.
    Lossy,
    /// No reviewed type exists; the export fails naming the unit.
    Refused,
}

/// One palette entry's mapping.
#[derive(Debug, Clone, Copy)]
pub struct UnitMapping {
    pub azoth_id: &'static str,
    pub disposition: Disposition,
    /// The DEXPI Process types this unit may be exchanged as. Empty for a refusal.
    pub dexpi_types: &'static [&'static str],
    /// What a lossy exchange loses, or why a refusal is a refusal.
    pub note: &'static str,
}

/// Every mapping row, in palette-id order.
#[must_use]
pub fn units() -> &'static [UnitMapping] {
    PROCESS_UNITS
}

/// The row for one azoth palette id.
#[must_use]
pub fn for_azoth(id: &str) -> Option<&'static UnitMapping> {
    PROCESS_UNITS.iter().find(|row| row.azoth_id == id)
}

/// The mapping's counts: faithful, lossy and refused.
///
/// The documentation's claims are measured from the same table, so these are the numbers a
/// page may quote rather than count by hand.
#[must_use]
pub fn counts() -> (usize, usize, usize) {
    (
        crate::process_mapping_gen::FAITHFUL_COUNT,
        crate::process_mapping_gen::LOSSY_COUNT,
        crate::process_mapping_gen::REFUSED_COUNT,
    )
}

/// The azoth units a DEXPI Process type may stand for, in palette-id order.
#[must_use]
pub fn candidates_for(dexpi_type: &str) -> &'static [&'static str] {
    PROCESS_TYPE_CANDIDATES
        .iter()
        .find_map(|(ty, ids)| (*ty == dexpi_type).then_some(*ids))
        .unwrap_or(&[])
}
