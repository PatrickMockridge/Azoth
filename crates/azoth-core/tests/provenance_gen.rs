//! The generated provenance table, held to what it claims about itself.
//!
//! Every registered id's block is emitted by `tools/gen_registry.py`, so what is worth
//! testing is not the values - those are hashes of committed bytes, and the Python suite
//! recomputes every one against the file it names - but the table's *shape*: that nothing
//! was truncated, that an id appears once, and that each block's own spelling of its id and
//! its status agree with the entry it is filed under. A table that had lost an entry would
//! still compile, still be ordered, and still look right.

use azoth_core::provenance::VerificationStatus;
use azoth_core::provenance_gen::{ALL_PROVENANCE, provenance};

#[test]
fn the_table_is_not_empty_and_every_id_is_found_by_its_own_name() {
    assert!(!ALL_PROVENANCE.is_empty(), "the table carries nothing");

    for entry in ALL_PROVENANCE {
        let found = provenance(entry.calc_id).expect("an entry is not found by its own id");
        assert_eq!(found.calc_id, entry.calc_id);
        assert!(!entry.name.is_empty(), "{} has no name", entry.calc_id);
        assert!(
            !entry.source.is_empty(),
            "{} names no source",
            entry.calc_id
        );
    }
}

#[test]
fn the_table_is_sorted_and_has_no_duplicate_ids() {
    // Sorted because the file says it is, and a reader checking it against the registry by
    // eye depends on that. Duplicates the lookup would hide: `find` returns the first, so a
    // second entry for one id would be unreachable rather than reported.
    let mut previous: Option<&str> = None;
    for entry in ALL_PROVENANCE {
        if let Some(last) = previous {
            assert!(
                last < entry.calc_id,
                "{:?} does not follow {last:?}: the table is out of order or repeats an id",
                entry.calc_id
            );
        }
        previous = Some(entry.calc_id);
    }
}

#[test]
fn every_status_is_one_the_type_parses() {
    // The table and the vocabulary are generated and hand-written respectively, so this is
    // the seam where they could disagree - and a status that did not parse would be a block
    // whose verification claim nothing could read.
    for entry in ALL_PROVENANCE {
        assert_eq!(
            VerificationStatus::parse(entry.verification.as_str()),
            Some(entry.verification),
            "{} carries a status the vocabulary does not know",
            entry.calc_id
        );
    }
}

#[test]
fn every_block_names_a_spec_and_both_implementations() {
    for entry in ALL_PROVENANCE {
        let named = [
            entry.spec_path,
            entry.spec_sha256,
            entry.python_path,
            entry.python_sha256,
            entry.rust_path,
            entry.rust_sha256,
        ];
        for value in named {
            assert!(
                !value.is_empty(),
                "{} names a file or hash it does not carry",
                entry.calc_id
            );
        }
    }
}

#[test]
fn an_id_this_build_does_not_ship_has_no_block() {
    assert!(provenance("hydraulics.not_a_calculation").is_none());
    assert!(provenance("").is_none());
}
