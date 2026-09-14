//! A keycard, read: the file a user supplies to extend or override what azoth ships.
//!
//! `keycard.example.toml` is the template, `specs/schema/keycard.schema.json` the
//! contract, and `python/src/azoth/keycard.py` the other implementation of it.
//!
//! Every section is marked `runtime` or `compiled` in the schema, and that annotation is
//! what this reader is organised around. `components` and `kij` become an [`Overlay`],
//! which is what the databank resolves a name against. `coefficients`, `models`,
//! `keyholder`, `fittings` and `fluids` are carried as the data they are; nothing in this
//! crate interprets them.
//!
//! An unknown section, parameter or unit is refused rather than skipped: a value nothing
//! reads is data that looks in use and is not.

use std::collections::BTreeMap;
use std::path::Path;

use azoth_core::unit_vocab_gen::{dimension, si_factor};
use azoth_core::{AzothError, Result};
use serde::Deserialize;

use crate::databank::{ComponentOverride, Overlay};

/// The keycard format version this reader understands. A card declaring anything else
/// is refused rather than guessed at.
pub const SCHEMA_VERSION: i64 = 2;

/// The component parameters a cubic reads, and the canonical unit each is stated in.
///
/// A card is held to the unit's **dimension** rather than its name: any unit in the
/// vocabulary of the same dimension is converted, and one of another dimension is
/// refused. That stops a temperature being read as a pressure, which is a plausible
/// wrong answer with no symptom.
pub const COMPONENT_PARAMETERS: &[(&str, &str)] =
    &[("Tc", "K"), ("Pc", "Pa"), ("omega", "dimensionless")];

/// The model vocabularies this build implements, each the only member it admits, and the
/// same lists as the schema's enums. A shape is listed only when the code runs it.
pub const MODEL_KINDS: &[&str] = &["cubic_eos"];
/// See [`MODEL_KINDS`].
pub const MODEL_SHAPES: &[&str] = &["peng_robinson"];
/// See [`MODEL_KINDS`].
pub const MODEL_ALPHAS: &[&str] = &["peng_robinson"];
/// See [`MODEL_KINDS`].
pub const MODEL_MIXING_RULES: &[&str] = &["classical_kij"];

/// Who is asserting the right to use a card's values.
///
/// Nothing reads this - no result carries it and no check requires it - and it is here
/// so that the file says whose it is. `licence` is a fetchable reference to the terms a
/// value was used under, optional because a card of public data has no licence to name.
#[derive(Debug, Clone, PartialEq, Eq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Keyholder {
    /// The person or organisation, as they want to be named.
    pub name: String,
    /// A fetchable reference to the licence itself, if the holder named one.
    pub licence: Option<String>,
}

/// One parameter of one substance, as the card states it.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Parameter {
    /// The number, in the unit beside it.
    pub value: f64,
    /// The unit, as a name from `specs/vocabulary/vocabulary.toml`.
    pub unit: String,
    /// Where the value came from, if the holder said. Never validated: provenance is the
    /// holder's, and a field a tool could check would be a field people fill in.
    pub citation: Option<String>,
}

/// One binary interaction parameter.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct KijRow {
    /// One substance of the pair.
    pub component_a: String,
    /// The other. The two are unordered: `a`/`b` and `b`/`a` name one pair.
    pub component_b: String,
    /// The parameter, dimensionless.
    pub value: f64,
    /// The temperature it was fitted at, where the source says.
    ///
    /// Carried and not read, which is how the format states it: nothing here
    /// interpolates an interaction parameter with temperature, and recording one is
    /// disclosure rather than a requirement.
    pub temperature_k: Option<f64>,
    /// Where the value came from, if the holder said.
    pub citation: Option<String>,
}

/// One coefficient a card supplies for a calculation's argument.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Coefficient {
    /// The number, in the unit beside it.
    pub value: f64,
    /// The unit, as a name from `specs/vocabulary/vocabulary.toml`.
    pub unit: String,
    /// Which convention the value is stated in, where the quantity is not a true ratio.
    pub convention: Option<String>,
    /// Where the value came from, if the holder said.
    pub citation: Option<String>,
}

impl Coefficient {
    /// The value as an SI base magnitude, the conversion the Python side's
    /// `azoth.core.units.to_si` makes, so a coefficient declared in `mm` means the same
    /// thing in both languages. The unit is checked when the card is read.
    #[must_use]
    pub fn si_value(&self) -> f64 {
        si_factor(&self.unit).map_or(self.value, |factor| self.value * factor)
    }
}

/// One declarative model definition: named choices, and nothing to execute.
#[derive(Debug, Clone, PartialEq, Deserialize)]
#[serde(deny_unknown_fields)]
pub struct Model {
    /// The model family.
    pub kind: String,
    /// The cubic's two constants, named as a pair rather than as loose numbers.
    pub shape: String,
    /// The alpha function.
    pub alpha: String,
    /// The mixing rule.
    pub mixing_rule: String,
    /// The substances, resolved against this card's `components` section.
    pub components: Vec<String>,
    /// The extra parameters an alpha function needs, which the one this build
    /// implements does not. Present and empty in the format, and refused when non-empty
    /// rather than stored: a fitted alpha is a different equation.
    #[serde(default)]
    pub alpha_parameters: BTreeMap<String, f64>,
}

/// A keycard document, as it is written. Private: what a caller gets is the [`Card`],
/// which is the same thing with every parameter resolved and every name checked.
#[derive(Debug, Deserialize)]
#[serde(deny_unknown_fields)]
struct Document {
    schema_version: i64,
    keyholder: Option<Keyholder>,
    components: Option<BTreeMap<String, BTreeMap<String, Parameter>>>,
    kij: Option<Vec<KijRow>>,
    coefficients: Option<BTreeMap<String, BTreeMap<String, Coefficient>>>,
    models: Option<BTreeMap<String, Model>>,
    /// The `compiled` stage, as the document wrote it. Carried rather than dropped so
    /// that opening a card loses nothing, and interpreted by nothing here.
    fittings: Option<Vec<toml::Value>>,
    /// See `fittings`.
    fluids: Option<toml::Table>,
}

/// A card, read and checked.
///
/// The type is the answer to "what did the file say", and it holds no file: a card is a
/// value, and reading one twice gives two values rather than one that something else can
/// change underneath a caller.
#[derive(Debug, Clone, PartialEq)]
pub struct Card {
    /// Whose card this is, if the file says.
    pub keyholder: Option<Keyholder>,
    /// The interaction parameters as written, each pair once, in file order.
    pub kij: Vec<KijRow>,
    /// The coefficients this card supplies, by calculation id and then input name.
    pub coefficients: BTreeMap<String, BTreeMap<String, Coefficient>>,
    /// The model definitions, by name.
    pub models: BTreeMap<String, Model>,
    /// The `compiled` fitting rows, exactly as the document wrote them.
    pub fittings: Vec<toml::Value>,
    /// The `compiled` fluid tables, exactly as the document wrote them.
    pub fluids: toml::Table,
    /// The component and interaction overrides, built while the document is checked, so
    /// building one cannot fail: every refusal has already happened by the time a `Card`
    /// exists.
    overlay: Overlay,
}

impl Card {
    /// Read a card from its text.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if the document does not parse as TOML, declares a
    ///   format version this build does not read, carries a section, a parameter, a unit
    ///   or a model choice this build does not know, or pairs a substance with itself.
    pub fn from_toml(text: &str) -> Result<Self> {
        let document: Document = toml::from_str(text).map_err(|error| {
            AzothError::invalid_input("<document>", format!("does not parse as a card: {error}"))
        })?;

        if document.schema_version != SCHEMA_VERSION {
            return Err(AzothError::invalid_input(
                "schema_version",
                format!(
                    "is {}; this build reads {SCHEMA_VERSION}. Refusing rather than \
                     guessing: version 1 had no `components`, `coefficients` or `models`, \
                     so one of those in a version-1 file would be silently dropped.",
                    document.schema_version
                ),
            ));
        }

        let mut overlay = Overlay::new();
        resolve_components(&mut overlay, document.components.as_ref())?;
        resolve_kij(&mut overlay, document.kij.as_ref())?;
        check_coefficients(document.coefficients.as_ref())?;
        check_models(document.models.as_ref())?;

        Ok(Self {
            keyholder: document.keyholder,
            kij: document.kij.unwrap_or_default(),
            coefficients: document.coefficients.unwrap_or_default(),
            models: document.models.unwrap_or_default(),
            fittings: document.fittings.unwrap_or_default(),
            fluids: document.fluids.unwrap_or_default(),
            overlay,
        })
    }

    /// Read a card from a file.
    ///
    /// # Errors
    /// * [`AzothError::InvalidInput`] if the file cannot be read, and as
    ///   [`Card::from_toml`] otherwise.
    pub fn from_path(path: &Path) -> Result<Self> {
        let text = std::fs::read_to_string(path).map_err(|error| {
            AzothError::invalid_input(
                path.display().to_string(),
                format!("cannot be read: {error}"),
            )
        })?;
        Self::from_toml(&text)
    }

    /// The overrides a Rust-native caller passes to `databank::entry(name, Some(...))`.
    #[must_use]
    pub fn overlay(&self) -> &Overlay {
        &self.overlay
    }

    /// This card's statement about one substance, or `None` if it makes none.
    #[must_use]
    pub fn component(&self, name: &str) -> Option<&ComponentOverride> {
        self.overlay.component(name)
    }

    /// The interaction parameter this card states for a pair, in either order.
    #[must_use]
    pub fn kij_for(&self, first: &str, second: &str) -> Option<f64> {
        let a = first.trim().to_lowercase();
        let b = second.trim().to_lowercase();
        if a == b {
            return None;
        }
        self.kij
            .iter()
            .find(|row| {
                let x = row.component_a.trim().to_lowercase();
                let y = row.component_b.trim().to_lowercase();
                (x == a && y == b) || (x == b && y == a)
            })
            .map(|row| row.value)
    }

    /// The coefficient this card supplies for one calculation's argument, or `None`.
    #[must_use]
    pub fn coefficient(&self, calc_id: &str, name: &str) -> Option<&Coefficient> {
        self.coefficients
            .get(calc_id)
            .and_then(|args| args.get(name))
    }

    /// A model definition by name, or `None`.
    #[must_use]
    pub fn model(&self, name: &str) -> Option<&Model> {
        self.models.get(name)
    }
}

/// Add every component the document names to the overlay, converting each parameter.
fn resolve_components(
    overlay: &mut Overlay,
    components: Option<&BTreeMap<String, BTreeMap<String, Parameter>>>,
) -> Result<()> {
    let Some(components) = components else {
        return Ok(());
    };

    for (name, parameters) in components {
        let mut override_ = ComponentOverride::default();
        for (parameter, body) in parameters {
            let field = format!("components.{name}.{parameter}");
            let Some((_, canonical)) = COMPONENT_PARAMETERS
                .iter()
                .find(|(accepted, _)| accepted == parameter)
            else {
                return Err(AzothError::invalid_input(
                    field,
                    format!(
                        "is not a parameter this build reads. Accepted: {:?}. Refused \
                         rather than stored: a value nothing reads is data that looks in \
                         use and is not.",
                        COMPONENT_PARAMETERS
                            .iter()
                            .map(|(name, _)| *name)
                            .collect::<Vec<_>>()
                    ),
                ));
            };
            let value = converted(body.value, &body.unit, canonical, &field)?;
            match parameter.as_str() {
                "Tc" => override_.tc = Some(value),
                "Pc" => override_.pc = Some(value),
                "omega" => override_.omega = Some(value),
                // Unreachable: the lookup above admits these three and no others.
                other => {
                    return Err(AzothError::invalid_input(
                        field,
                        format!("`{other}` has no reader in this build"),
                    ));
                }
            }
        }
        overlay.set_component(name, override_);
    }
    Ok(())
}

/// Add every interaction pair the document states to the overlay, which refuses a pair
/// naming one substance twice.
fn resolve_kij(overlay: &mut Overlay, kij: Option<&Vec<KijRow>>) -> Result<()> {
    let Some(kij) = kij else {
        return Ok(());
    };

    for (index, row) in kij.iter().enumerate() {
        let field = format!("kij[{index}]");
        if row.component_a.trim().is_empty() || row.component_b.trim().is_empty() {
            return Err(AzothError::invalid_input(
                field,
                "needs a non-empty `component_a` and `component_b`",
            ));
        }
        overlay.set_kij(&row.component_a, &row.component_b, row.value)?;
    }
    Ok(())
}

/// Check the units on every coefficient. The values are carried as written.
fn check_coefficients(
    coefficients: Option<&BTreeMap<String, BTreeMap<String, Coefficient>>>,
) -> Result<()> {
    let Some(coefficients) = coefficients else {
        return Ok(());
    };

    for (calc_id, arguments) in coefficients {
        for (name, body) in arguments {
            let field = format!("coefficients.{calc_id}.{name}");
            if dimension(&body.unit).is_none() {
                return Err(AzothError::invalid_input(
                    field,
                    format!(
                        "declares the unit {:?}, which is not in the vocabulary this \
                         format uses. The vocabulary spells kelvin `K`.",
                        body.unit
                    ),
                ));
            }
            // The dimension is the *spec's* to declare for an input - a coefficient is
            // whatever the calculation it belongs to takes - so the check here is that
            // the unit is one this build can convert at all, and the conversion is
            // `Coefficient::si_value`.
        }
    }
    Ok(())
}

/// Check every choice a model definition makes.
fn check_models(models: Option<&BTreeMap<String, Model>>) -> Result<()> {
    let Some(models) = models else {
        return Ok(());
    };

    for (name, model) in models {
        let field = format!("models.{name}");
        for (key, chosen, accepted) in [
            ("kind", &model.kind, MODEL_KINDS),
            ("shape", &model.shape, MODEL_SHAPES),
            ("alpha", &model.alpha, MODEL_ALPHAS),
            ("mixing_rule", &model.mixing_rule, MODEL_MIXING_RULES),
        ] {
            if !accepted.contains(&chosen.as_str()) {
                return Err(AzothError::invalid_input(
                    format!("{field}.{key}"),
                    format!(
                        "is {chosen:?}, which this build does not implement. Accepted: \
                         {accepted:?}. Refused rather than accepted-and-ignored: a model \
                         declaring a variant and evaluated as another is a wrong answer \
                         with no symptom."
                    ),
                ));
            }
        }

        if model.components.is_empty() {
            return Err(AzothError::invalid_input(
                field,
                "needs a non-empty `components` list",
            ));
        }

        if !model.alpha_parameters.is_empty() {
            return Err(AzothError::invalid_input(
                format!("{field}.alpha_parameters"),
                format!(
                    "gives {:?}, but the alpha function this build implements \
                     (`peng_robinson`) takes none. A fitted alpha is a different \
                     equation, not a parameterisation of this one.",
                    model.alpha_parameters.keys().collect::<Vec<_>>()
                ),
            ));
        }
    }
    Ok(())
}

/// One declared value in whichever unit the card used, as an SI base magnitude.
fn converted(value: f64, unit: &str, canonical: &str, field: &str) -> Result<f64> {
    let Some(declared) = dimension(unit) else {
        return Err(AzothError::invalid_input(
            field,
            format!(
                "declares the unit {unit:?}, which is not in the vocabulary this format \
                 uses. Accepted: {canonical} and anything of the same dimension."
            ),
        ));
    };
    let expected = dimension(canonical).ok_or_else(|| {
        // A defect in this crate rather than in the card: the canonical name is written
        // here and the vocabulary does not have it, which only a test can catch.
        AzothError::invalid_input(
            "<vocabulary>",
            format!("`{canonical}` is not a unit this build knows"),
        )
    })?;
    if declared != expected {
        return Err(AzothError::invalid_input(
            field,
            format!(
                "is declared in {unit:?}, which is not a {canonical}. A `Tc` meant in \
                 another dimension is a factor with no symptom, so the dimension is \
                 checked rather than the number."
            ),
        ));
    }
    // The unit is in the vocabulary and `dimension` and `si_factor` are built from the
    // same table, so this factor exists. A missing one is a defect in the generated
    // table, and it is reported rather than passed through as the declared number.
    si_factor(unit)
        .ok_or_else(|| {
            AzothError::invalid_input(field, format!("`{unit}` has no conversion in this build"))
        })
        .map(|factor| value * factor)
}

#[cfg(test)]
mod tests {
    use super::COMPONENT_PARAMETERS;

    /// The parameter set `specs/schema/component.schema.json` declares, name and unit.
    ///
    /// Rust does not read the schema, so the expected list is written here and the
    /// Python-side test holds both to the schema. A parameter is added in all three
    /// places or none.
    const DECLARED: &[(&str, &str)] = &[("Tc", "K"), ("Pc", "Pa"), ("omega", "dimensionless")];

    #[test]
    fn component_parameters_match_the_declaration() {
        assert_eq!(COMPONENT_PARAMETERS, DECLARED);
    }
}
