//! What this namespace does not carry, as a lookup rather than a sentence.
//!
//! **A refusal is a fact about the declaration, so it is declared once and read here.** The
//! values this namespace's models do not carry are rows in the spec's `[[unported]]` array,
//! generated into [`crate::unported_gen::UNPORTED`]; a call site names its row by key and gets
//! back the NeqSim symbol that would close the gap.
//!
//! **The shape is what makes the two implementations comparable**, and it is the same argument
//! `azoth-process::unported` makes at length: a refusal written as a sentence is a sentence in
//! Rust and a different sentence in Python, and the pair drifts unobserved. Here it matters for
//! three rows that a wrong reading would turn into three plausible numbers - `TbpClosure`'s
//! members that cannot be inverted for density.
//!
//! A miss is an error as well, and it names the key that is missing: a refusal whose row was
//! never declared has to be a message somebody reads rather than a silent pass.

use azoth_core::AzothError;

use crate::unported_gen::UNPORTED;

/// One value a model declares as not carried.
#[derive(Debug, Clone, Copy)]
pub struct Unported {
    /// How a call site names the row: `parameter=value`, or `parameter@other=member` with the
    /// conditions sorted and joined by `&`.
    pub key: &'static str,
    /// The model the row belongs to.
    pub model: &'static str,
    /// The NeqSim class or member that would close the gap.
    pub class: &'static str,
    /// The capture that measures the refusal, or empty where none does.
    pub capture: &'static str,
}

/// The row for a key, or `None` where the declaration does not carry one.
#[must_use]
pub fn row(key: &str) -> Option<&'static Unported> {
    UNPORTED.iter().find(|row| row.key == key)
}

/// The input a key names - the part before the `=` or the `@`.
fn parameter_of(key: &str) -> &str {
    key.split(['=', '@']).next().unwrap_or(key)
}

/// Refuse a value the declaration marks as not carried.
///
/// `key` is exactly the row's own key, and it is a literal at every call site on purpose - a key
/// the compiler does not have to compute is a key the checker can read without following control
/// flow, in this language and in the other.
#[must_use]
pub fn refuse(key: &str) -> AzothError {
    let parameter = parameter_of(key);
    let Some(row) = row(key) else {
        return AzothError::invalid_input(
            parameter,
            format!(
                "`{key}` is refused, and no `[[unported]]` row declares it. A refusal is a claim \
                 about the declaration, so it belongs in the `[[unported]]` array of the spec that \
                 carries this input, and `tools/check_unported.py` is what holds the two together."
            ),
        );
    };
    let close = format!("`{}` is the class that would close it", row.class);
    let measured = if row.capture.is_empty() {
        String::new()
    } else {
        format!(" `{}` measures what it does.", row.capture)
    };
    AzothError::invalid_input(
        parameter,
        format!("`{key}` is not ported: {close}.{measured}"),
    )
}
