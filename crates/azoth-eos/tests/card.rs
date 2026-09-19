//! The keycard reader: what a card says, and what it refuses.
//!
//! The file these tests read is the shipped template, from the repository root, by
//! path. Reading the file the documentation tells a user to copy is the point: a
//! template that this reader cannot take teaches the wrong shape, exactly as one that
//! `tools/check_user_data.py` rejects does.

use std::path::PathBuf;

use azoth_eos::Cubic;
use azoth_eos::association::SiteScheme;
use azoth_eos::card::{
    CoefficientValue,
    {
        ASSOCIATION_PARAMETERS, ASSOCIATION_SCHEMES, COMPONENT_PARAMETERS, Card, MODEL_KINDS,
        SCHEMA_VERSION,
    },
};
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
    assert_eq!(cd.value, CoefficientValue::Scalar(0.61));
    assert_eq!(cd.unit, "dimensionless");
    assert_eq!(cd.si_value().scalar(), Some(0.61));
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

    let pair = databank::kij("methane", "n-butane", Cubic::Pr, Some(card.overlay()));
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

    assert_eq!(
        d.value,
        CoefficientValue::Scalar(50.0),
        "what the card wrote is disclosure"
    );
    assert_eq!(d.unit, "mm");
    assert_eq!(
        d.si_value().scalar(),
        Some(0.05),
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

/// A coefficient may be a matrix, and reading it keeps its shape and its unit.
///
/// `eos.uniquac_activity_coefficients`'s `aij` is the one that needs this: NeqSim carries
/// no UNIQUAC interaction table, so a card is the only place the parameter can come from
/// that is not the caller's own argument.
#[test]
fn a_matrix_coefficient_keeps_its_shape_and_its_unit() {
    let card = Card::from_toml(&a_card(
        "[coefficients.\"eos.uniquac_activity_coefficients\".aij]\n\
         value = [[0.0, -71.0], [209.0, 0.0]]\nunit = \"K\"\n",
    ))
    .expect("a valid card");
    let aij = card
        .coefficient("eos.uniquac_activity_coefficients", "aij")
        .expect("the card states it");

    assert_eq!(aij.unit, "K");
    assert_eq!(aij.si_value().scalar(), None, "a matrix is not a scalar");
    assert_eq!(
        aij.si_value(),
        CoefficientValue::List(vec![
            CoefficientValue::List(vec![
                CoefficientValue::Scalar(0.0),
                CoefficientValue::Scalar(-71.0),
            ]),
            CoefficientValue::List(vec![
                CoefficientValue::Scalar(209.0),
                CoefficientValue::Scalar(0.0),
            ]),
        ])
    );
}

/// A matrix in millimetres converts every entry, and a vector is a value too.
#[test]
fn a_vector_is_a_value_and_a_unit_scales_every_entry() {
    let card = Card::from_toml(&a_card(
        "[coefficients.\"eos.uniquac_activity_coefficients\".aij]\n\
         value = [[0.0, 2.0], [3.0, 0.0]]\nunit = \"mm\"\n",
    ))
    .expect("a valid card");
    let aij = card
        .coefficient("eos.uniquac_activity_coefficients", "aij")
        .expect("the card states it");
    assert_eq!(
        aij.si_value(),
        CoefficientValue::List(vec![
            CoefficientValue::List(vec![
                CoefficientValue::Scalar(0.0),
                CoefficientValue::Scalar(0.002),
            ]),
            CoefficientValue::List(vec![
                CoefficientValue::Scalar(0.003),
                CoefficientValue::Scalar(0.0),
            ]),
        ])
    );
}

/// A ragged list of lists is not a matrix, and is refused rather than read row by row
/// against a matrix the calculation expects.
#[test]
fn a_ragged_coefficient_is_refused() {
    let err = Card::from_toml(&a_card(
        "[coefficients.\"eos.uniquac_activity_coefficients\".aij]\n\
         value = [[0.0, -71.0], [209.0]]\nunit = \"K\"\n",
    ))
    .expect_err("a ragged matrix is refused");
    let message = format!("{err}");
    assert!(
        message.contains("rectangular"),
        "the error should say what is wrong with the shape: {message}"
    );
}

/// An empty value is refused: a coefficient with no numbers is not a default.
#[test]
fn an_empty_coefficient_is_refused() {
    let err = Card::from_toml(&a_card(
        "[coefficients.\"eos.uniquac_activity_coefficients\".aij]\nvalue = []\nunit = \"K\"\n",
    ))
    .expect_err("an empty value is refused");
    assert!(
        format!("{err}").contains("rectangular"),
        "an empty value is not a matrix: {err}"
    );
}

/// One substance's association, stated in SI, as the overlay and the databank resolve it.
///
/// The card's numbers are the published CPA set for water - `a0 = 0.12277`,
/// `b = 1.4515e-5`, `eps = 16655`, `kappa_AB = 0.0692` - and the databank's are NeqSim's
/// internal scale, which is the same numbers a hundred thousand times larger. Reading the
/// databank's back is what says the crossing happens once and in the right direction.
#[test]
fn an_association_stated_in_si_reaches_the_databank_in_its_own_scale() {
    let card = Card::from_toml(&a_card(
        "[associations.water]\nscheme = \"4C\"\n\
         [associations.water.energy]\nvalue = 16655.0\nunit = \"J/mol\"\n\
         [associations.water.volume_srk]\nvalue = 0.0692\nunit = \"dimensionless\"\n\
         [associations.water.a_srk]\nvalue = 0.12277\nunit = \"Pa*m**6/mol**2\"\n\
         [associations.water.b_srk]\nvalue = 1.4515e-5\nunit = \"m**3/mol\"\n",
    ))
    .expect("a valid card");

    let stated = card
        .association_for("water")
        .expect("the card states water's association");
    assert_eq!(stated.scheme, SiteScheme::FourC);
    assert_eq!(stated.energy, Some(16655.0), "J/mol is the canonical unit");
    assert_eq!(stated.a_srk, Some(0.12277));

    let shipped = databank::entry("water", None).expect("water ships");
    let carded = databank::entry("water", Some(card.overlay())).expect("water");
    let base = shipped
        .association
        .as_ref()
        .expect("water ships associating");
    let merged = carded.association.as_ref().expect("the card's association");

    assert_eq!(merged.scheme, SiteScheme::FourC);
    assert_eq!(merged.energy, 16655.0, "J/mol is the table's unit too");
    assert_eq!(
        merged.volume_srk, 0.0692,
        "kappa_AB is a ratio and crosses no scale"
    );
    assert_eq!(merged.a_srk, base.a_srk, "0.12277 SI is the table's 12277");
    assert_eq!(merged.b_srk, base.b_srk, "and 1.4515e-5 is its 1.4515");
    assert_eq!(
        merged.m_srk, base.m_srk,
        "a parameter the card does not name is the table's"
    );
    assert_eq!(
        merged.a_pr, base.a_pr,
        "and the other family's set is untouched"
    );
    assert_eq!(
        merged.sites, base.sites,
        "water's scheme and its stated count agree"
    );
}

/// A parameter the card names is the card's, and one it does not is the table's - the rule
/// a card follows for the cubic's parameters too, applied to the association's.
#[test]
fn an_association_is_overridden_parameter_by_parameter() {
    let card = Card::from_toml(&a_card(
        "[associations.methanol]\nscheme = \"2B\"\n\
         [associations.methanol.energy]\nvalue = 20000.0\nunit = \"J/mol\"\n",
    ))
    .expect("a valid card");

    let shipped = databank::entry("methanol", None).expect("methanol ships");
    let base = shipped
        .association
        .as_ref()
        .expect("methanol ships associating");
    let carded = databank::entry("methanol", Some(card.overlay())).expect("methanol");
    let record = carded.association.as_ref().expect("the card states one");

    assert_eq!(record.scheme, SiteScheme::TwoB, "the card's");
    assert_eq!(record.energy, 20000.0, "the card's, against 24591 shipped");
    assert_eq!(record.a_srk, base.a_srk, "the table's");
    assert_eq!(
        record.volume_pr, base.volume_pr,
        "and the other family's too"
    );
}

/// A new substance can carry an association, which is the point of the section being
/// separate from `components`: a holder's own fluid has no shipped association to inherit.
#[test]
fn a_card_added_substance_can_carry_an_association() {
    let card = Card::from_toml(&a_card(
        "[components.lab_solvent.Tc]\nvalue = 500.0\nunit = \"K\"\n\
         [components.lab_solvent.Pc]\nvalue = 4.0e6\nunit = \"Pa\"\n\
         [components.lab_solvent.omega]\nvalue = 0.5\nunit = \"dimensionless\"\n\
         [associations.lab_solvent]\nscheme = \"2B\"\n\
         [associations.lab_solvent.energy]\nvalue = 20000.0\nunit = \"J/mol\"\n",
    ))
    .expect("a valid card");

    let added = databank::entry("lab_solvent", Some(card.overlay())).expect("the card adds it");
    let record = added.association.as_ref().expect("the card states one");
    assert_eq!(record.scheme, SiteScheme::TwoB);
    assert_eq!(record.energy, 20000.0);
    assert_eq!(
        record.volume_srk, 0.0,
        "no table row to inherit a kappa_AB from, so the card's is the whole of it"
    );
}

/// A scheme this build does not implement, and a parameter nothing reads. Both are
/// refused rather than stored.
#[test]
fn an_association_this_build_cannot_use_is_refused() {
    let unknown_scheme = Card::from_toml(&a_card("[associations.water]\nscheme = \"3B\"\n"))
        .expect_err("`3B` is not a scheme");
    let message = format!("{unknown_scheme}");
    assert!(message.contains("associations.water.scheme"), "{message}");
    for scheme in ASSOCIATION_SCHEMES {
        assert!(
            message.contains(scheme),
            "the message names {scheme}: {message}"
        );
    }

    let unknown_parameter = Card::from_toml(&a_card(
        "[associations.water]\nscheme = \"4C\"\n\
         [associations.water.epsilon]\nvalue = 1.0\nunit = \"J/mol\"\n",
    ))
    .expect_err("`epsilon` is not the parameter's name");
    let message = format!("{unknown_parameter}");
    assert!(
        message.contains("associations.water.epsilon"),
        "the refusal names the key that is wrong: {message}"
    );
    for (parameter, _) in ASSOCIATION_PARAMETERS {
        assert!(
            message.contains(parameter),
            "the message names {parameter}: {message}"
        );
    }
}

/// A unit of the wrong dimension, refused for the same reason a component's is.
#[test]
fn an_association_parameter_in_the_wrong_unit_is_refused() {
    let error = Card::from_toml(&a_card(
        "[associations.water]\nscheme = \"4C\"\n\
         [associations.water.a_srk]\nvalue = 0.12277\nunit = \"dimensionless\"\n",
    ))
    .expect_err("an attraction is not a ratio");

    assert!(
        format!("{error}").contains("Pa*m**6/mol**2"),
        "the refusal names the unit it wanted: {error}"
    );
}

/// The shipped baseline card resolves to the table it was generated from.
///
/// **Nothing read this file before the association's fitted values went into it**, and a
/// generated artefact nothing reads is where a factor of a hundred thousand lives without
/// a symptom: the components' parameters are copied through unchanged, but `acpa_*` and
/// `bcpa_*` cross from NeqSim's internal scale to SI in `tools/gen_keycard.py` and back in
/// [`AssociationOverride::applied_to`]. Only a round trip through both says the two
/// crossings agree - and neither the generator's own output nor the reader's own input can
/// say it alone.
#[test]
fn the_baseline_card_resolves_to_the_table_it_was_generated_from() {
    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../databank/keycard.toml");
    let card = Card::from_path(&path).expect("the baseline card is one this build reads");

    let mut names = card.overlay().component_names();
    names.sort_unstable();
    assert_eq!(names.len(), 348, "the table's own count");

    let mut compared = 0;
    for name in names {
        let shipped = databank::entry(name, None).expect("the card names shipped substances");
        let carded = databank::entry(name, Some(card.overlay())).expect(name);

        // Every parameter the card copies through is the table's exactly: the unit it is
        // declared in is the unit it is stored in, so the conversion is an identity and a
        // difference would be a defect rather than a rounding.
        assert_eq!(carded.tc, shipped.tc, "{name}: Tc");
        assert_eq!(carded.pc, shipped.pc, "{name}: Pc");
        assert_eq!(carded.omega, shipped.omega, "{name}: omega");
        assert_eq!(carded.cp, shipped.cp, "{name}: Cp");
        assert_eq!(carded.antoine, shipped.antoine, "{name}: Antoine");
        assert_eq!(carded.reference_state, shipped.reference_state, "{name}");
        // **A card never changes a substance's class.** It is the class that decides
        // whether a cubic may be built at all, and the baseline card carries NeqSim's
        // ion rows - filler `Tc`, `Pc` and `omega` and all - so an overlay that reset it
        // would hand a cubic the filler the refusal exists to keep out of one.
        assert_eq!(carded.class, shipped.class, "{name}: class");

        let (Some(record), Some(base)) = (&carded.association, &shipped.association) else {
            assert_eq!(
                carded.association, shipped.association,
                "{name}: a shipped association survives the card"
            );
            continue;
        };
        compared += 1;
        assert_eq!(record.scheme, base.scheme, "{name}");
        assert_eq!(record.sites, base.sites, "{name}");
        assert_eq!(record.energy, base.energy, "{name}");
        assert_eq!(record.volume_srk, base.volume_srk, "{name}");
        assert_eq!(record.m_srk, base.m_srk, "{name}");
        assert_eq!(record.volume_pr, base.volume_pr, "{name}");
        assert_eq!(record.m_pr, base.m_pr, "{name}");
        assert_eq!(record.racket_z, base.racket_z, "{name}");
        assert_eq!(record.volume_correction, base.volume_correction, "{name}");

        // The four that cross a scale, and the only tolerance in this test. A decimal
        // shift here is a factor of a hundred thousand rather than a rounding, so the
        // tolerance is there to admit the last bit of a double and nothing else.
        for (crossed, table) in [
            (record.a_srk, base.a_srk),
            (record.b_srk, base.b_srk),
            (record.a_pr, base.a_pr),
            (record.b_pr, base.b_pr),
        ] {
            assert!(
                (crossed - table).abs() <= 1e-9 * table.abs().max(1.0),
                "{name}: {crossed} is not {table}"
            );
        }
    }

    assert_eq!(
        compared, 167,
        "the associating substances the table carries"
    );
}

/// A `kij` row's associating columns reach the column an associating mixture reads.
///
/// The row states three parameters of one pair, and the two columns are separate fits:
/// water/methanol is `-0.153` in `cpakij_SRK`, so a card overriding only the SRK one must
/// leave the PR one at the table's value rather than carrying the override across.
#[test]
fn a_kij_row_states_the_associating_columns_beside_the_classical_one() {
    use azoth_eos::association::AssociationCubic;

    let card = Card::from_toml(&a_card(
        "[[kij]]\ncomponent_a = \"water\"\ncomponent_b = \"methanol\"\nvalue = -0.0789\n\
         cpa_value_srk = -0.08\ncpa_value_pr = -0.31\n",
    ))
    .expect("a valid card");
    let overlay = card.overlay();
    let names = ["water", "methanol"];

    assert_eq!(
        databank::cpa_kij(&names, AssociationCubic::Srk, Some(overlay))[1],
        -0.08
    );
    assert_eq!(
        databank::cpa_kij(&names, AssociationCubic::Pr, Some(overlay))[1],
        -0.31
    );
    // And the classical column is the card's too, not either of the other two.
    let (classical, _) = databank::mixture_of(&names, Cubic::Pr, Some(overlay)).expect("a mixture");
    assert!((classical.kij(0, 1) - -0.0789).abs() < 1.0e-12);
}

/// A `kij` row still needs its classical value, and a card that omits one is refused
/// rather than read as stating the associating column alone.
#[test]
fn a_kij_row_without_the_classical_value_is_refused() {
    let error = Card::from_toml(&a_card(
        "[[kij]]\ncomponent_a = \"water\"\ncomponent_b = \"methanol\"\ncpa_value_srk = -0.08\n",
    ))
    .expect_err("`value` is required");

    assert!(format!("{error}").contains("value"), "{error}");
}

/// The shipped baseline card's interaction rows are the table's, in all three columns.
///
/// The companion of the components test, for the half of the card that is keyed by a pair
/// rather than a name. A row the generator emits with a `cpa_value_*` is asserting a value
/// for a column, so a crossing that was applied in the wrong direction - or a row written
/// against the wrong column - would show here and nowhere else.
#[test]
fn the_baseline_cards_pairs_are_the_tables() {
    use azoth_eos::association::AssociationCubic;

    let path = PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../databank/keycard.toml");
    let card = Card::from_path(&path).expect("the baseline card is one this build reads");

    let mut pairs = card.overlay().kij_pairs();
    pairs.sort();
    assert!(pairs.len() > 400, "the shipped pairs are the table's");

    let mut associating = 0;
    for (first, second) in pairs {
        let names = [first.as_str(), second.as_str()];
        assert_eq!(
            databank::kij(&first, &second, Cubic::Pr, Some(card.overlay())),
            databank::kij(&first, &second, Cubic::Pr, None),
            "{first}/{second}: the classical column"
        );
        for family in [AssociationCubic::Srk, AssociationCubic::Pr] {
            let carded = databank::cpa_kij(&names, family, Some(card.overlay()))[1];
            let table = databank::cpa_kij(&names, family, None)[1];
            assert_eq!(carded, table, "{first}/{second}: the {family:?} column");
            if table != 0.0 {
                associating += 1;
            }
        }
    }

    // Counted over the shipped pairs at both families, of which the generator writes only
    // the non-zero ones. A generator that dropped the column entirely would leave this at
    // zero and the loop above comparing nothing but zeros.
    assert_eq!(
        associating, 711,
        "the non-zero associating cells the table has"
    );
}
