//! The mapping table covers the shipped palette exactly once.

mod common;

use std::collections::BTreeSet;

use azoth_dexpi::mapping;

#[test]
fn every_palette_entry_has_exactly_one_row() {
    let palette = common::palette();
    let ids: BTreeSet<&str> = palette.iter().map(|spec| spec.id.as_str()).collect();
    let rows: BTreeSet<&str> = mapping::units().iter().map(|row| row.azoth_id).collect();
    assert_eq!(
        rows.len(),
        mapping::units().len(),
        "no row is declared twice"
    );
    assert_eq!(ids, rows, "the table and the palette name the same entries");
}

#[test]
fn every_mapped_type_is_in_the_reviewed_profile() {
    for row in mapping::units() {
        for dexpi_type in row.dexpi_types {
            assert!(
                azoth_dexpi::profile::supports(dexpi_type),
                "{} names an unreviewed type {dexpi_type}",
                row.azoth_id
            );
        }
    }
}

#[test]
fn every_refused_row_has_a_reason_and_every_lossy_row_has_a_note() {
    for row in mapping::units() {
        match row.disposition {
            mapping::Disposition::Refused => {
                assert!(row.dexpi_types.is_empty(), "{}", row.azoth_id);
                assert!(!row.note.is_empty(), "{} carries no reason", row.azoth_id);
            }
            mapping::Disposition::Lossy => {
                assert!(!row.dexpi_types.is_empty(), "{}", row.azoth_id);
                assert!(!row.note.is_empty(), "{} carries no note", row.azoth_id);
            }
            mapping::Disposition::Faithful => {
                assert!(!row.dexpi_types.is_empty(), "{}", row.azoth_id);
            }
        }
    }
}

#[test]
fn the_counts_agree_with_the_table() {
    let (faithful, lossy, refused) = mapping::counts();
    let total = mapping::units().len();
    assert_eq!(faithful + lossy + refused, total);
    assert_eq!(
        mapping::units()
            .iter()
            .filter(|row| row.disposition == mapping::Disposition::Faithful)
            .count(),
        faithful
    );
    assert_eq!(
        mapping::units()
            .iter()
            .filter(|row| row.disposition == mapping::Disposition::Lossy)
            .count(),
        lossy
    );
    assert_eq!(
        mapping::units()
            .iter()
            .filter(|row| row.disposition == mapping::Disposition::Refused)
            .count(),
        refused
    );
}

#[test]
fn the_reverse_index_lists_only_mappable_units() {
    let mut types = BTreeSet::new();
    for row in mapping::units() {
        types.extend(row.dexpi_types.iter().copied());
    }
    for dexpi_type in types {
        let candidates = mapping::candidates_for(dexpi_type);
        let unique: BTreeSet<&str> = candidates.iter().copied().collect();
        assert_eq!(
            unique.len(),
            candidates.len(),
            "{dexpi_type} lists a unit twice"
        );
        for row in mapping::units() {
            if row.dexpi_types.contains(&dexpi_type) {
                assert!(
                    candidates.contains(&row.azoth_id),
                    "{dexpi_type} is missing {}",
                    row.azoth_id
                );
            }
        }
    }
}
