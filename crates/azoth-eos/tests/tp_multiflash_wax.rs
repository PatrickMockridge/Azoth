//! Spec-driven tests for the `eos.tp_multiflash_wax` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::mixture::{Component, Mixture};
use azoth_eos::{Cubic, databank, model_gen, tp_multiflash_wax};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.tp_multiflash_wax";

fn mixture_from_case(case: &azoth_core::spec::TestCase) -> Mixture {
    let names = case.list("components").expect("components");
    databank::mixture_of(names, Cubic::Srk, None)
        .expect("the case's components resolve")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let mixture = mixture_from_case(case);
        let result = tp_multiflash_wax(
            &mixture,
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.wax_fraction,
            common::expected(case, "wax_fraction"),
            case.tolerance,
            context,
        );
        assert_eq!(
            result.phase_count,
            common::expected(case, "phase_count") as u32,
            "{context}: the phase count"
        );
    }
}

/// **A material balance, which is what the amount has to be.**
///
/// `sum_p beta_p x_ip` must be the feed for every component, and the wax phase is the one
/// that a model can get wrong while still landing on a plausible fraction - the hydrate
/// family's instrument found exactly that, at `+0.0874` on water. The fractions are the ones
/// the solve reports, so this tests the state it assembled and not the iteration's algebra.
#[test]
fn the_state_balances() {
    let names = ["methane", "n-heptane", "nc14", "nc20"];
    let (mixture, _) = databank::mixture_of(&names, Cubic::Srk, None).expect("resolves");
    let z = [0.7, 0.1, 0.1, 0.1];

    for t_k in [290.0, 275.0, 265.0, 255.0, 245.0] {
        let result = tp_multiflash_wax(&mixture, kelvins(t_k), pascals(5.0e5), &z).expect("runs");
        let mut worst: f64 = 0.0;
        for (i, &total) in z.iter().enumerate() {
            let mut counted = 0.0;
            for (phase, &beta) in result.beta.iter().enumerate() {
                counted += beta * result.x[phase][i];
            }
            worst = worst.max((counted - total).abs());
        }
        assert!(
            worst < 1.0e-08,
            "T = {t_k}: the balance is out by {worst} over the {} phase(s)",
            result.phase_count
        );
        // And the fractions are a split of the feed, not a set of unrelated numbers.
        let total: f64 = result.beta.iter().sum();
        assert!(
            (total - 1.0).abs() < 1.0e-08,
            "T = {t_k}: beta sums to {total}"
        );
    }
}

/// **The wax amount rises as the temperature falls, and only the oil feeds it.**
///
/// The captured states give the shape: at 265, 255 and 245 K the wax fraction is `0.0707`,
/// `0.1320` and `0.1542` while the gas holds at `0.6914`, `0.6931` and `0.6932` and the oil
/// falls from `0.2379` to `0.1526`. A model that precipitated out of the gas would move the
/// wrong phase, and one that had the sign of the fusion term backwards would rise with the
/// temperature instead.
#[test]
fn the_wax_takes_from_the_oil_and_rises_as_it_cools() {
    let names = ["methane", "n-heptane", "nc14", "nc20"];
    let (mixture, _) = databank::mixture_of(&names, Cubic::Srk, None).expect("resolves");
    let z = [0.7, 0.1, 0.1, 0.1];

    let mut previous = -1.0;
    for (t_k, expected) in [
        (265.0, 0.070_673_392_486_494_9),
        (255.0, 0.132_007_556_310_333),
        (245.0, 0.154_193_075_421_103),
    ] {
        let result = tp_multiflash_wax(&mixture, kelvins(t_k), pascals(5.0e5), &z).expect("runs");
        assert!(
            result.wax_fraction > previous,
            "T = {t_k}: the fraction is {} and it was {previous} at the warmer state",
            result.wax_fraction
        );
        previous = result.wax_fraction;
        // `1e-7`: NeqSim's own fraction solve stops at `ans.norm2() > 1e-6`, so this is finer
        // than the criterion the capture was taken at.
        // **`1e-6`, which is NeqSim's own criterion and therefore the finest this can be
        // held to**: its fraction solve stops at `ans.norm2() > 1e-6`, so the three states
        // agree to `2.95e-8`, `2.12e-7` and `1.38e-7` relative - at or below the line its own
        // iteration stopped at. A tighter tolerance would pin where NeqSim stopped rather
        // than the state it stopped at.
        assert!((result.wax_fraction / expected - 1.0).abs() < 1.0e-06);
        assert_eq!(result.phase_count, 3);
    }

    // Above the appearance the answer is the two-phase flash and `converged` says so.
    for t_k in [290.0, 275.0] {
        let result = tp_multiflash_wax(&mixture, kelvins(t_k), pascals(5.0e5), &z).expect("runs");
        assert_eq!(result.wax_fraction, 0.0, "T = {t_k}");
        assert_eq!(result.phase_count, 2, "T = {t_k}");
        assert!(
            !result.converged,
            "T = {t_k}: the flag must say which was solved"
        );
    }
}

/// A fluid with no wax former is answered by the two-phase flash rather than refused.
#[test]
fn a_feed_with_no_wax_former_is_two_phase() {
    let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], Cubic::Srk, None)
        .expect("the pair resolves");
    let result =
        tp_multiflash_wax(&mixture, kelvins(280.0), pascals(5.0e5), &[0.7, 0.3]).expect("runs");
    assert_eq!(result.wax_fraction, 0.0);
    assert_eq!(result.phase_count, 2);
    assert!(
        result.converged,
        "nothing was in the set that could not be solved"
    );
}

/// A wax former with no melt data is refused rather than given a zero.
#[test]
fn a_wax_former_with_no_melt_data_is_refused() {
    // A caller-supplied component that says it is a wax former and states nothing else.
    let (base, _) = databank::mixture_of(&["methane", "nc14"], Cubic::Srk, None).expect("resolves");
    // **By index, because a `Component` carries no name** - it is critical constants and
    // nothing else, and the name it was resolved from is the databank's key and not a field.
    let components: Vec<Component> = base
        .components()
        .iter()
        .enumerate()
        .map(|(index, component)| component.clone().with_wax_data(index == 1, 0.0, 0.0))
        .collect();
    let mixture = Mixture::new(components, vec![0.0; 4])
        .expect("two components")
        .with_cubic(Cubic::Srk);
    let error = tp_multiflash_wax(&mixture, kelvins(250.0), pascals(5.0e5), &[0.7, 0.3])
        .expect_err("a wax former with no heat of fusion has no fusion term");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}
