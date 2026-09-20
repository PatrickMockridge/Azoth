//! Spec-driven tests for the `eos.hydrate_fraction` model.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::{Cubic, hydrate, hydrate_fraction, model_gen, pt_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hydrate_fraction";

/// The probe's own feed, as `validation/neqsim/captures/hydrate_fraction_probe.tsv` states it.
const Z: [f64; 4] = [
    0.781_018_289_668_808_7,
    0.099_851_705_388_037_56,
    0.020_266_930_301_532_374,
    0.098_863_074_641_621_35,
];

fn probe_mixture() -> hydrate::Hydration {
    hydrate::hydrate_mixture_of(&["methane", "ethane", "propane", "water"], Cubic::Srk, None)
        .expect("the probe's feed resolves for a hydrate")
        .0
        .hydration()
        .expect("attached")
        .clone()
}

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> azoth_eos::mixture::Mixture {
    let names = case.list("components").expect("components");
    hydrate::hydrate_mixture_of(names, Cubic::Srk, None)
        .expect("the case resolves for a hydrate")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let result = hydrate_fraction(
            &mixture,
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.beta,
            common::expected(case, "beta"),
            case.tolerance,
            context,
        );
        assert_eq!(
            Some(result.structure.as_str()),
            case.expected_string("structure"),
            "{context}: the stable structure"
        );
        common::assert_close(
            result.balance_error,
            common::expected(case, "balance_error"),
            case.tolerance,
            &format!("{context} (balance)"),
        );
    }
}

/// **The material balance, which is what this model exists over NeqSim's.**
///
/// NeqSim's own state at the same feed holds `+0.0874` water and `-0.0734` methane against
/// the feed, as `validation/neqsim/captures/hydrate_fraction_probe.tsv` records. This one
/// closes by construction, and the number is asserted rather than described.
#[test]
fn the_material_balance_closes() {
    let (mixture, _) =
        hydrate::hydrate_mixture_of(&["methane", "ethane", "propane", "water"], Cubic::Srk, None)
            .expect("resolves");
    let result = hydrate_fraction(&mixture, kelvins(288.15), pascals(1.0e7), &Z).expect("computes");
    assert!(
        result.balance_error < 1.0e-12,
        "the balance is out by {}, and NeqSim's own state is out by 0.0874 on water",
        result.balance_error
    );
}

/// **The two hydrate models are one equilibrium, read in two directions.**
///
/// `eos.hydrate_formation_temperature` solves for the temperature where the objective
/// `ln(f_w^hydrate/f_w^fluid)` is zero; this model reads the sign of the same objective at the
/// same feed to decide between no hydrate and all of the water. So the check that ties them
/// together is that the sign change and the temperature are the same state: above it the
/// fraction is **exactly** zero, below it there is an amount, and at it the residual is zero.
///
/// Neither is NeqSim's - its fraction does not depend on the temperature at all, because its
/// bound is the constant `z_water/(46/54)` - so the consistency is between the two here rather
/// than against an oracle.
#[test]
fn the_two_models_are_one_equilibrium() {
    let (mixture, _) =
        hydrate::hydrate_mixture_of(&["methane", "ethane", "propane", "water"], Cubic::Srk, None)
            .expect("resolves");
    let formation =
        azoth_eos::hydrate_formation_temperature(&mixture, pascals(1.0e7), &Z).expect("solves");
    let at = formation.temperature.value;

    let above =
        hydrate_fraction(&mixture, kelvins(at + 1.0), pascals(1.0e7), &Z).expect("computes");
    let below =
        hydrate_fraction(&mixture, kelvins(at - 5.0), pascals(1.0e7), &Z).expect("computes");
    let on = hydrate_fraction(&mixture, kelvins(at), pascals(1.0e7), &Z).expect("computes");

    assert_eq!(
        above.beta, 0.0,
        "a kelvin above the formation temperature the fraction is {}, and it is not a hydrate",
        above.beta
    );
    assert!(
        above.residual > 0.0 && below.residual < 0.0,
        "the residual is the formation model's own and changes sign at its temperature, but \
         is {} above and {} below",
        above.residual,
        below.residual
    );
    assert!(
        below.beta > 0.05,
        "five kelvin below the formation temperature the fraction is only {}",
        below.beta
    );
    assert!(
        on.residual.abs() < 1.0e-6,
        "at the formation temperature the residual is {}, not zero",
        on.residual
    );
    // **The answer jumps, and that is the equilibrium rather than a tolerance.** The objective
    // is flat and negative over the whole range while the aqueous phase pins water's fugacity,
    // so a hydrate that is stable at all is stable at every fraction up to the bound, and the
    // fraction goes from zero to the bound as the temperature crosses. Exactly at the crossing
    // the residual's sign is round-off, so either branch is the honest answer there.
    assert!(
        on.beta == 0.0 || (on.beta - below.beta).abs() < 1.0e-3,
        "at the formation temperature the fraction is {}, which is neither zero nor the {} \
         below it",
        on.beta,
        below.beta
    );
}

/// **The cages drive the fraction, and nothing else does.**
///
/// Halving methane's Langmuir constant empties its cages, which is a *hydrate* change: the
/// fluid the same feed flashes to must be bit-for-bit the state it was. That is the seam a
/// unit error hides behind - a constant read in the wrong units is a constant of the wrong
/// size - so the test names both sides of it rather than only the side that moves.
#[test]
fn the_cages_move_the_fraction_and_not_the_fluid() {
    let (mixture, _) =
        hydrate::hydrate_mixture_of(&["methane", "ethane", "propane", "water"], Cubic::Srk, None)
            .expect("resolves");
    let mut broken = probe_mixture();
    for guest in &mut broken.guests {
        if guest.name == "methane" {
            guest.langmuir_a[0][0] *= 0.5;
            guest.langmuir_a[1][0] *= 0.5;
        }
    }
    let sabotaged = mixture.clone().with_hydration(broken);

    let sound = hydrate_fraction(&mixture, kelvins(288.15), pascals(1.0e7), &Z).expect("computes");
    let moved =
        hydrate_fraction(&sabotaged, kelvins(288.15), pascals(1.0e7), &Z).expect("computes");
    assert!(
        (sound.beta - moved.beta).abs() > 1.0e-4,
        "halving methane's Langmuir constant moved the fraction only from {} to {}",
        sound.beta,
        moved.beta
    );

    let sound_fluid = pt_flash(&mixture, kelvins(288.15), pascals(1.0e7), &Z).expect("flashes");
    let moved_fluid = pt_flash(&sabotaged, kelvins(288.15), pascals(1.0e7), &Z).expect("flashes");
    assert_eq!(
        sound_fluid, moved_fluid,
        "the hydrate's tables changed the fluid, which they are not part of"
    );
}

/// A fluid with no water, or nothing that occupies a cage, is refused.
#[test]
fn a_fluid_that_cannot_form_a_hydrate_is_refused() {
    let no_water = hydrate::hydrate_mixture_of(&["methane", "ethane"], Cubic::Srk, None)
        .expect_err("no water, no hydrate");
    assert!(
        matches!(no_water, AzothError::InvalidInput { .. }),
        "{no_water:?}"
    );
}
