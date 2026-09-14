//! The keycard reader: what a card says, and what it refuses.
//!
//! The file these tests read is the shipped template, from the repository root, by
//! path. Reading the file the documentation tells a user to copy is the point: a
//! template that this reader cannot take teaches the wrong shape, exactly as one that
//! `tools/check_user_data.py` rejects does.

use std::path::PathBuf;

use azoth_eos::card::{COMPONENT_PARAMETERS, Card, MODEL_KINDS, SCHEMA_VERSION};
use azoth_eos::databank;

/// The template, at the path the README tells a user to copy it from.
fn template_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .join("../../keycard.example.toml")
        .canonicalize()
        .expect("keycard.example.toml is at the repository root")
}

fn template() -> Card {
    Card::from_path(&template_path()).expect("the template is a card this build reads")
}

/// A card carrying one section, as text. Written here rather than through a writer:
/// there is no TOML writer in this repository, and a test that needed one would be
/// testing something other than the reader.
fn a_card(body: &str) -> String {
    format!("schema_version = {SCHEMA_VERSION}\n{body}")
}

#[test]
fn the_shipped_template_loads() {
    let card = template();

    let holder = card
        .keyholder
        .as_ref()
        .expect("the template names a keyholder");
    assert_eq!(holder.name, "Example Engineering Ltd");
    assert!(holder.licence.is_none(), "the template names no licence");

    // The one parameter the template states, converted from `K` to the SI base
    // magnitude the databank works in. The value is the template's own, and it is
    // deliberately not methane's real critical temperature.
    let methane = card.component("methane").expect("the card names methane");
    assert_eq!(methane.tc, Some(190.0));
    assert_eq!(
        methane.pc, None,
        "a card naming one parameter keeps the rest"
    );
    assert_eq!(methane.omega, None);

    // A pair, in either order, because a card may write it either way.
    assert_eq!(card.kij_for("methane", "n-butane"), Some(0.0135));
    assert_eq!(card.kij_for("n-butane", "methane"), Some(0.0135));
    assert_eq!(card.kij_for("methane", "methane"), None);

    // A coefficient, carried with the unit it was declared in and converted beside it.
    let cd = card
        .coefficient("hydraulics.orifice_flow", "Cd")
        .expect("the template supplies a discharge coefficient");
    assert_eq!(cd.value, 0.61);
    assert_eq!(cd.unit, "dimensionless");
    assert_eq!(cd.si_value(), 0.61);
    assert_eq!(cd.convention.as_deref(), Some("iso_5167_corner"));
    assert!(
        card.coefficient("hydraulics.orifice_flow", "nope")
            .is_none()
    );

    // The compiled sections, carried as written and interpreted by nothing.
    assert_eq!(card.fittings.len(), 1);
    assert!(card.fluids.contains_key("example_fluid"));
    assert!(
        card.models.is_empty(),
        "the template's models are a comment"
    );
}

#[test]
fn a_card_reaches_the_databank() {
    let card = template();
    let shipped = databank::entry("methane", None).expect("methane ships");

    let carded = databank::entry("methane", Some(card.overlay())).expect("methane");
    assert_eq!(carded.tc, 190.0, "the card's Tc, not the table's");
    assert_eq!(
        carded.pc, shipped.pc,
        "a parameter the card does not name is kept"
    );
    assert_eq!(carded.omega, shipped.omega);

    let pair = databank::kij("methane", "n-butane", Some(card.overlay()));
    assert_eq!(pair, 0.0135);
}

#[test]
fn a_coefficient_is_carried_as_written_and_converted_beside_it() {
    // A coefficient is whatever the calculation it belongs to takes, so the card's unit
    // is not checked against a dimension the way a component parameter's is - it is
    // checked against the vocabulary, and the conversion is offered beside the value.
    // `mm` is in the vocabulary and is a thousandth of a metre, which is the whole
    // content of the test.
    let card = Card::from_toml(&a_card(
        "[coefficients.\"hydraulics.choked_flow_area\".d]\nvalue = 50.0\nunit = \"mm\"\n",
    ))
    .expect("a valid card");
    let d = card
        .coefficient("hydraulics.choked_flow_area", "d")
        .expect("the card states it");

    assert_eq!(d.value, 50.0, "what the card wrote is disclosure");
    assert_eq!(d.unit, "mm");
    assert_eq!(
        d.si_value(),
        0.05,
        "what a calculation takes is an SI magnitude"
    );
}

#[test]
fn a_unit_of_the_wrong_dimension_is_refused() {
    // A `Tc` in `Pa` is a number that would be read as kelvin. The dimension is checked
    // rather than the number, because the wrong one has no symptom: every result shifts
    // and none of them looks wrong.
    let error = Card::from_toml(&a_card(
        "[components.ethane.Tc]\nvalue = 1.0\nunit = \"Pa\"\n",
    ))
    .expect_err("`Pa` is not a temperature");
    let message = error.to_string();

    assert!(message.contains("components.ethane.Tc"), "{message}");
    assert!(message.contains("not a K"), "{message}");
}

#[test]
fn a_unit_outside_the_vocabulary_is_refused() {
    // `kelvin` is the spelling people reach for, and the vocabulary spells it `K`.
    let error = Card::from_toml(&a_card(
        "[components.ethane.Tc]\nvalue = 1.0\nunit = \"kelvin\"\n",
    ))
    .expect_err("`kelvin` is not in the vocabulary");

    assert!(error.to_string().contains("kelvin"), "{error}");
}

#[test]
fn a_parameter_nothing_reads_is_refused() {
    let error = Card::from_toml(&a_card(
        "[components.ethane.critical_pressure]\nvalue = 1.0\nunit = \"Pa\"\n",
    ))
    .expect_err("nothing reads `critical_pressure`");
    let message = error.to_string();

    assert!(
        message.contains("components.ethane.critical_pressure"),
        "{message}"
    );
    for (parameter, _) in COMPONENT_PARAMETERS {
        assert!(
            message.contains(parameter),
            "the message names {parameter}: {message}"
        );
    }
}

#[test]
fn a_substance_paired_with_itself_is_refused() {
    // A self-pair would silently rescale that substance's attraction.
    let error = Card::from_toml(&a_card(
        "[[kij]]\ncomponent_a = \"methane\"\ncomponent_b = \"methane\"\nvalue = 0.1\n",
    ))
    .expect_err("a substance does not interact with itself");

    assert!(
        error.to_string().contains("does not interact with itself"),
        "{error}"
    );
}

#[test]
fn a_version_this_build_does_not_read_is_refused() {
    let error =
        Card::from_toml("schema_version = 1\n[components.methane.Tc]\nvalue = 1.0\nunit = \"K\"\n")
            .expect_err("version 1 is not this format");

    assert!(error.to_string().contains("schema_version"), "{error}");
}

#[test]
fn a_section_this_build_does_not_know_is_refused() {
    // Refused rather than skipped: a misspelled `component` would otherwise be a whole
    // section of a card that looks in use and is read by nothing.
    let error = Card::from_toml(&a_card(
        "[component.methane.Tc]\nvalue = 1.0\nunit = \"K\"\n",
    ))
    .expect_err("`component` is not a section");

    assert!(error.to_string().contains("component"), "{error}");
}

#[test]
fn a_model_declaring_a_variant_this_build_lacks_is_refused() {
    // Accepted-and-ignored would be the worst outcome: a card naming Soave-Redlich-Kwong
    // and evaluated as Peng-Robinson gives a wrong answer with no symptom.
    let error = Card::from_toml(&a_card(
        "[models.vendor_gas]\nkind = \"cubic_eos\"\nshape = \"soave_redlich_kwong\"\n\
         alpha = \"peng_robinson\"\nmixing_rule = \"classical_kij\"\ncomponents = [\"methane\"]\n",
    ))
    .expect_err("this build implements one shape");
    let message = error.to_string();

    assert!(message.contains("models.vendor_gas.shape"), "{message}");
    for shape in azoth_eos::card::MODEL_SHAPES {
        assert!(message.contains(shape), "{message}");
    }
    for kind in MODEL_KINDS {
        assert!(!message.is_empty(), "{kind}");
    }
}

#[test]
fn a_model_with_an_alpha_parameter_is_refused() {
    let error = Card::from_toml(&a_card(
        "[models.vendor_gas]\nkind = \"cubic_eos\"\nshape = \"peng_robinson\"\n\
         alpha = \"peng_robinson\"\nmixing_rule = \"classical_kij\"\ncomponents = [\"methane\"]\n\
         [models.vendor_gas.alpha_parameters]\nc1 = 0.5\n",
    ))
    .expect_err("the alpha function this build implements takes none");

    assert!(error.to_string().contains("alpha_parameters"), "{error}");
}

#[test]
fn a_document_that_is_not_toml_is_refused_as_such() {
    let error = Card::from_toml("schema_version: 2\n").expect_err("that is not TOML");

    assert!(error.to_string().contains("does not parse"), "{error}");
}
