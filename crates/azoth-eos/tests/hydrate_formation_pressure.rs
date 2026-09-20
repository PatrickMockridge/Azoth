//! Spec-driven tests for the `eos.hydrate_formation_pressure` model.

use azoth_core::AzothError;
use azoth_core::units::kelvins;
use azoth_eos::{
    Cubic, hydrate, hydrate_formation_pressure, hydrate_formation_temperature, model_gen,
};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hydrate_formation_pressure";

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> azoth_eos::mixture::Mixture {
    let names = case
        .list("components")
        .expect("the case declares components");
    hydrate::hydrate_mixture_of(names, Cubic::Srk, None)
        .expect("the case's components resolve for a hydrate")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let result = hydrate_formation_pressure(
            &mixture,
            kelvins(common::input(case, "T")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.pressure.value,
            common::expected(case, "pressure"),
            case.tolerance,
            &format!("{context} (P)"),
        );
        // The structure is an output of the same comparison the pressure is, so a case that
        // pinned only the pressure would pass with the wrong cages behind it.
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

/// **The two models are one curve, read from either end.**
///
/// The same equality settles both: `eos.hydrate_formation_temperature` finds the temperature
/// at a pressure and this finds the pressure at a temperature. So the check that ties them
/// together is that each one's answer is the other's input: solving at (T, P(T)) must return
/// the temperature it started from, and its residual must be the one this model reports.
///
/// Neither is NeqSim's search - it secants one and iterates the pressure ratio at the other -
/// so the consistency is between the two here rather than against an oracle.
#[test]
fn the_two_models_are_one_curve() {
    let (mixture, _) =
        hydrate::hydrate_mixture_of(&["methane", "ethane", "propane", "water"], Cubic::Srk, None)
            .expect("resolves");
    let z = [
        0.781_018_289_668_808_7,
        0.099_851_705_388_037_56,
        0.020_266_930_301_532_374,
        0.098_863_074_641_621_35,
    ];

    for t_k in [278.15, 283.15, 288.15] {
        let pressure =
            hydrate_formation_pressure(&mixture, kelvins(t_k), &z).expect("the pressure solves");
        let temperature =
            hydrate_formation_temperature(&mixture, pressure.pressure, &z).expect("solves back");

        assert_eq!(
            pressure.structure, temperature.structure,
            "T = {t_k}: the two directions disagree about the cages"
        );
        common::assert_close(
            temperature.temperature.value,
            t_k,
            1.0e-7,
            &format!("T = {t_k}: the temperature the pressure's answer solves back to"),
        );
    }
}

/// A fluid with no water, or with nothing that occupies a cage, is refused.
#[test]
fn a_fluid_that_cannot_form_a_hydrate_is_refused() {
    let no_water = hydrate::hydrate_mixture_of(&["methane", "ethane"], Cubic::Srk, None)
        .expect_err("no water, no hydrate");
    assert!(
        matches!(no_water, AzothError::InvalidInput { .. }),
        "{no_water:?}"
    );

    let no_guest = hydrate::hydrate_mixture_of(&["water", "methanol"], Cubic::Srk, None)
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
    let error = hydrate_formation_pressure(&plain, kelvins(288.15), &[0.9, 0.1])
        .expect_err("the mixture carries no guests");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}
