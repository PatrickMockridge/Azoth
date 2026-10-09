//! Spec-driven tests for the `characterization.characterise_plus_fraction` model.
//!
//! Every expected value is a sweep row of `validation/neqsim/captures/plus_fraction_probe.tsv`:
//! the lumping model's `fractionOfHeavyEnd`, the added components' molar masses and gravities, and
//! the model `Characterise` actually ran.

use azoth_characterization::{
    CharacterisePlusFractionResult, PlusModel, characterise_plus_fraction, model_gen,
};
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_core::{AzothError, spec::TestCase};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.characterise_plus_fraction";

fn call(case: &TestCase) -> CharacterisePlusFractionResult {
    characterise_plus_fraction(
        kilograms_per_mole(common::input(case, "molar_mass")),
        kilograms_per_cubic_meter(common::input(case, "density")),
        common::input(case, "mole_fraction"),
        common::input(case, "first_carbon_number") as usize,
        case.string("plus_model")
            .map(|name| name.parse().expect("the case names a model")),
        case.input("number_of_lumps").map(|n| n as usize),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        if let Some(expected) = case.expected_string("selected_model") {
            assert_eq!(result.selected_model.as_str(), expected, "{context}::selected_model");
        }
        for (name, actual, expected) in [
            (
                "fraction_of_heavy_end",
                result.fraction_of_heavy_end.clone(),
                case.expected_vector("fraction_of_heavy_end"),
            ),
            (
                "lump_molar_mass",
                result.lump_molar_mass.iter().map(|q| q.value).collect(),
                case.expected_vector("lump_molar_mass"),
            ),
            (
                "lump_density",
                result.lump_density.iter().map(|q| q.value).collect(),
                case.expected_vector("lump_density"),
            ),
        ] {
            let Some(expected) = expected else { continue };
            assert_eq!(actual.len(), expected.len(), "{context}::{name}: the lump count");
            for (index, (got, want)) in actual.iter().zip(expected).enumerate() {
                common::assert_close(
                    *got,
                    *want,
                    case.tolerance,
                    &format!("{context}::{name}[{index}]"),
                );
            }
        }
        common::assert_consistent(&result, &context);
    }
}

/// The three models' own thresholds, which are the class's fields rather than a name's promise.
#[test]
fn the_thresholds_are_the_classes_own_fields() {
    assert!((PlusModel::Pedersen.maximum_molar_mass() - 0.605).abs() < 1e-15);
    assert!((PlusModel::PedersenHeavyOil.maximum_molar_mass() - 2.10).abs() < 1e-15);
    // **`WhitsonGammaModel` extends `PedersenPlusModel`**, so its constructor sets `0.605` - and
    // not the 2.10 the heavy-oil name next to it suggests.
    assert!((PlusModel::WhitsonGamma.maximum_molar_mass() - 0.605).abs() < 1e-15);
    assert_eq!(PlusModel::Pedersen.last_carbon_number(), 80);
    assert_eq!(PlusModel::PedersenHeavyOil.last_carbon_number(), 200);
}

#[test]
fn the_replacement_is_silent_and_the_heavy_model_is_the_same_run() {
    // 0.65 kg/mol is above Pedersen's 0.605 and below the heavy model's 2.10, so the request for
    // `pedersen` is replaced - and asking for the heavy model outright is then the same run.
    let run = |model: PlusModel| {
        characterise_plus_fraction(
            kilograms_per_mole(0.65),
            kilograms_per_cubic_meter(900.0),
            0.1,
            20,
            Some(model),
            None,
        )
        .expect("computes")
    };
    let replaced = run(PlusModel::Pedersen);
    let asked = run(PlusModel::PedersenHeavyOil);
    assert_eq!(replaced.selected_model, PlusModel::PedersenHeavyOil);
    assert_eq!(asked.selected_model, PlusModel::PedersenHeavyOil);
    assert_eq!(replaced.fraction_of_heavy_end, asked.fraction_of_heavy_end);
    assert_eq!(replaced.lump_molar_mass, asked.lump_molar_mass);
}

#[test]
fn a_gamma_request_above_the_threshold_is_swept_away_too() {
    // The replacement fires before the model is used, so it is not a fallback *within* the gamma
    // path - and `selected_model` is the only thing that says so.
    let swept = characterise_plus_fraction(
        kilograms_per_mole(0.65),
        kilograms_per_cubic_meter(900.0),
        0.1,
        20,
        Some(PlusModel::WhitsonGamma),
        None,
    )
    .expect("computes");
    assert_eq!(swept.selected_model, PlusModel::PedersenHeavyOil);
    // And below the threshold it is the gamma model that runs, with its own range.
    let gamma = characterise_plus_fraction(
        kilograms_per_mole(0.4),
        kilograms_per_cubic_meter(850.0),
        0.1,
        20,
        Some(PlusModel::WhitsonGamma),
        None,
    )
    .expect("computes");
    assert_eq!(gamma.selected_model, PlusModel::WhitsonGamma);
}

#[test]
fn a_split_that_declines_is_refused_rather_than_skipped() {
    // `PedersenPlusModel.characterizePlusFraction` compares the table's first cut against the
    // plus fraction's molar mass and returns false when the plus fraction is lighter. The class
    // then leaves the fluid with its single plus row and no lumps; an id that returns a table has
    // no fluid to leave unchanged, so the split's own message comes out.
    for model in [PlusModel::Pedersen, PlusModel::PedersenHeavyOil] {
        let refused = characterise_plus_fraction(
            kilograms_per_mole(0.2),
            kilograms_per_cubic_meter(780.0),
            0.1,
            20,
            Some(model),
            None,
        );
        let error = refused.expect_err("a plus fraction under C20's 275 g/mol is refused");
        assert!(
            error.to_string().contains("275"),
            "the split's own refusal is propagated: {error}"
        );
        assert!(matches!(error, AzothError::InvalidInput { .. }));
    }
}
