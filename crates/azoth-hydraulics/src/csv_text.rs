//! Reading the shared data files, where the format is text and empty means absent.
//!
//! The fittings registry and the fluid tables are both CSV with `#` comment
//! banners, and both carry provenance columns that a placeholder row leaves empty.
//! That makes "empty field" and "absent value" the same thing in both, and it is a
//! distinction worth having once rather than twice: `Option<String>` says a value
//! is not there, where `String::new()` would say it is there and says nothing -
//! which is exactly the difference between a row that has no source and a row
//! whose source is an empty string.

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

/// An empty CSV field means absent, not an empty string.
#[must_use]
pub(crate) fn optional(raw: &str) -> Option<String> {
    let trimmed = raw.trim();
    (!trimmed.is_empty()).then(|| trimmed.to_string())
}
