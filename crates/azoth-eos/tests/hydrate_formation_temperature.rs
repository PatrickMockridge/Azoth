//! Spec-driven tests for the `eos.hydrate_formation_temperature` model.

use azoth_core::AzothError;
use azoth_core::units::pascals;
use azoth_eos::hydrate::HydrateModel;
use azoth_eos::{Cubic, hydrate, hydrate_formation_temperature, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hydrate_formation_temperature";

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> azoth_eos::mixture::Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    hydrate::hydrate_mixture_of(
        names,
        Cubic::Srk,
        None,
        case.string("hydrate_model")
            .and_then(|text| text.parse().ok())
            .unwrap_or(HydrateModel::Pvtsim),
    )
    .expect("the case's components resolve for a hydrate")
    .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let result = hydrate_formation_temperature(
            &mixture,
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.temperature.value,
            common::expected(case, "temperature"),
            case.tolerance,
            &format!("{context} (T)"),
        );
        // The structure is an output of the same comparison the temperature is, so a case
        // that pinned only the temperature would pass with the wrong cages behind it.
        assert_eq!(
            Some(result.structure.as_str()),
            case.expected_string("structure"),
            "{context}: the stable structure"
        );
        assert!(
            result.residual.abs() < 1e-6,
            "{context}: the residual at the answer is {}",
            result.residual
        );
    }
}

/// A fluid with no water, or with nothing that occupies a cage, is refused.
#[test]
fn a_fluid_that_cannot_form_a_hydrate_is_refused() {
    let no_water = hydrate::hydrate_mixture_of(
        &["methane", "ethane"],
        Cubic::Srk,
        None,
        HydrateModel::Pvtsim,
    )
    .expect_err("no water, no hydrate");
    assert!(
        matches!(no_water, AzothError::InvalidInput { .. }),
        "{no_water:?}"
    );

    let no_guest = hydrate::hydrate_mixture_of(
        &["water", "methanol"],
        Cubic::Srk,
        None,
        HydrateModel::Pvtsim,
    )
    .expect_err("nothing in the mixture occupies a cage");
    assert!(
        matches!(no_guest, AzothError::InvalidInput { .. }),
        "{no_guest:?}"
    );
}

/// A mixture resolved without the hydrate's tables is refused rather than given a default.
#[test]
fn a_plain_mixture_is_refused() {
    let (plain, _) = azoth_eos::databank::mixture_of(&["methane", "water"], Cubic::Srk, None)
        .expect("the pair resolves");
    let error = hydrate_formation_temperature(&plain, pascals(1.0e7), &[0.9, 0.1])
        .expect_err("the mixture carries no guests");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}
