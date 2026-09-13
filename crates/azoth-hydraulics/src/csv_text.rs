//! Reading the shared data files, where the format is CSV with a comment banner.
//!
//! Written once because it was about to be written twice: the fittings registry and
//! the fluid tables are both `#`-bannered CSV and both needed the banner stripped.
//!
//! It carried an `optional` helper for the provenance columns, which parsed an
//! empty field as an absent value. Those columns are gone - the databank's
//! provenance is institutional and a per-row URL was fiction - so the helper went
//! with them.

/// The leading `#` comment banner and blank lines, removed.
///
/// The banner is not decoration: it carries the copyright statement and the
/// meaning of `verify_status`, so the parser has to tolerate it rather than the
/// file having to be machine-only.
#[must_use]
pub(crate) fn body(raw: &str) -> String {
    raw.lines()
        .filter(|line| !line.trim_start().starts_with('#') && !line.trim().is_empty())
        .collect::<Vec<_>>()
        .join("\n")
}
