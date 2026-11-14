//! Identities: the XML ids DEXPI requires, and the azoth names they must not lose.
//!
//! DEXPI's `xsd:ID` is `[A-Za-z_][A-Za-z_0-9]*`, and azoth's instance ids are free-form -
//! they may carry a hyphen (`10-VA-001`), which is a legal azoth id and an illegal XML id.
//! The two are not made to agree by mangling one into the other. **Structural XML ids are
//! counters** (`ProcessStep1`, `Stream1`, `Port1`), so they are unique by construction; the
//! azoth name travels in the standard `Identifier`/`Label` fields, which are free strings and
//! preserve it exactly.

/// Whether a string is a legal DEXPI `xsd:ID`.
#[must_use]
pub fn is_xml_id(text: &str) -> bool {
    let mut chars = text.chars();
    let Some(first) = chars.next() else {
        return false;
    };
    (first.is_ascii_alphabetic() || first == '_')
        && chars.all(|value| value.is_ascii_alphanumeric() || value == '_')
}

/// A string that can be an azoth instance id and survive an endpoint's grammar.
///
/// **The dot and the brackets are the endpoint's**, not XML's: `instance.port[i]` is how a
/// connection addresses a stream, so an instance id carrying either would make the endpoint
/// ambiguous. **A hyphen is not**: `V-101` is a legal azoth id and a tag a P&ID actually
/// carries, so it survives; whitespace and anything outside `[A-Za-z0-9_-]` does not.
#[must_use]
pub fn sanitize_instance_id(text: &str) -> String {
    let mut sanitized = String::with_capacity(text.len() + 1);
    for (index, value) in text.chars().enumerate() {
        let accepted = if index == 0 {
            value.is_ascii_alphabetic() || value == '_' || value == '-'
        } else {
            value.is_ascii_alphanumeric() || value == '_' || value == '-'
        };
        sanitized.push(if accepted { value } else { '_' });
    }
    if sanitized.is_empty() {
        sanitized.push('_');
    }
    sanitized
}

/// A DEXPI `name` from an azoth flowsheet id: `[A-Za-z_][A-Za-z_0-9]*`, no hyphen.
///
/// **Stricter than an instance id**, because `Model/@name` is an XSD `name` - a dot and a
/// hyphen are both illegal there, while an azoth id may carry either.
#[must_use]
pub fn model_name(id: &str) -> String {
    let mut sanitized = String::with_capacity(id.len() + 1);
    for (index, value) in id.chars().enumerate() {
        let accepted = if index == 0 {
            value.is_ascii_alphabetic() || value == '_'
        } else {
            value.is_ascii_alphanumeric() || value == '_'
        };
        sanitized.push(if accepted { value } else { '_' });
    }
    if sanitized.is_empty() || sanitized == "_" {
        "AzothFlowsheet".to_string()
    } else {
        sanitized
    }
}

/// A counter that mints structural ids, with the collision-free guarantee a counter gives.
#[derive(Debug, Default)]
pub struct Mint {
    next: usize,
}

impl Mint {
    #[must_use]
    pub fn new() -> Self {
        Self { next: 0 }
    }

    /// The next id: `prefix` followed by a one-based count.
    pub fn next(&mut self, prefix: &str) -> String {
        self.next += 1;
        format!("{prefix}{}", self.next)
    }
}

/// Make `candidate` unique against `taken`, suffixing `_2`, `_3`, ... as needed.
pub fn uniquify(candidate: &str, taken: &mut std::collections::BTreeSet<String>) -> String {
    if taken.insert(candidate.to_string()) {
        return candidate.to_string();
    }
    let mut suffix = 2usize;
    loop {
        let alternative = format!("{candidate}_{suffix}");
        if taken.insert(alternative.clone()) {
            return alternative;
        }
        suffix += 1;
    }
}
