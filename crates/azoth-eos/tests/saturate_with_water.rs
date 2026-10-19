//! `eos.saturate_with_water`: the source's check value, the phase state it reports, and the
//! two refusals.
//!
//! The case file carries the scalar assertion the registry can express. This carries what it
//! cannot: the boolean the generator refuses as an expectation, the phase state, and the
//! refusals - which is what the generator's own message asks for.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::saturate_with_water;

fn names(list: &[&str]) -> Vec<String> {
    list.iter().map(|name| (*name).to_string()).collect()
}

/// `SaturateWithWaterTest.testRun4`, the source's own published value.
#[test]
fn the_pr_check_value_reproduces() {
    let result = saturate_with_water(
        &names(&["methane", "n-heptane"]),
        &[2.0, 75.0],
        kelvins(293.15),
        pascals(150.0e5),
        "pr",
        false,
    )
    .expect("the fluid is in the databank");

    assert!(
        (result.water_mole_fraction - 2.4301370485671443e-4).abs() < 1e-9,
        "{}",
        result.water_mole_fraction
    );
    // **The boolean the case cannot carry.** Flooding must find a free aqueous phase, or the
    // answer is the feed returned unchanged and the number above is meaningless.
    assert!(
        result.aqueous_present,
        "the flood forms a free aqueous phase"
    );
    assert_eq!(result.phase_count, 1, "the saturated oil is one phase");
    assert!(
        (result.phase_fractions.iter().sum::<f64>() - 1.0).abs() < 1e-12,
        "the reported fractions are renormalised over what remains"
    );
    assert!(result.warnings.is_empty(), "{:?}", result.warnings);
}

/// A water-free feed is the source's own case: the class adds the trace, and so does this.
#[test]
fn the_trace_of_water_is_added_when_the_feed_has_none() {
    let without = saturate_with_water(
        &names(&["methane", "n-heptane"]),
        &[2.0, 75.0],
        kelvins(293.15),
        pascals(150.0e5),
        "pr",
        false,
    )
    .expect("the fluid is in the databank");
    let with = saturate_with_water(
        &names(&["methane", "n-heptane", "water"]),
        &[2.0, 75.0, 0.0],
        kelvins(293.15),
        pascals(150.0e5),
        "pr",
        false,
    )
    .expect("the fluid is in the databank");

    // The two differ by the trace's share of the total, which is 1% of a feed: the saturated
    // composition they produce is the same amount of dissolved water over a slightly
    // different total, so the fractions differ by about that much and by nothing else.
    let relative = (without.water_mole_fraction - with.water_mole_fraction).abs()
        / without.water_mole_fraction;
    assert!(relative < 2.0e-2, "{without:?} against {with:?}");
}

/// The CPA state, which the source's own `SystemSrkCPAstatoil` asserts.
///
/// The fluid is the generic CPA mixture plus `PhaseSrkCPAs.calc_g`'s Chapman distribution;
/// the Carnahan-Starling one moves this answer by 15%, and the case file carries the same
/// bound on the source's own number.
#[test]
fn an_associating_fluid_is_the_statoil_cpa() {
    let result = saturate_with_water(
        &names(&["methane", "n-heptane"]),
        &[2.0, 75.0],
        kelvins(293.15),
        pascals(150.0e5),
        "srk",
        true,
    )
    .expect("the Statoil CPA is built");
    let expected = 2.846_502_497_491_981_6e-4;
    let relative = (result.water_mole_fraction - expected).abs() / expected;
    assert!(relative < 5.0e-3, "{relative}");
}

#[test]
fn a_cubic_this_library_does_not_carry_is_refused() {
    let error = saturate_with_water(
        &names(&["methane", "n-heptane"]),
        &[2.0, 75.0],
        kelvins(293.15),
        pascals(150.0e5),
        "not-a-cubic",
        false,
    )
    .expect_err("an unknown cubic is a refusal");
    assert!(error.to_string().contains("not-a-cubic"), "{error}");
}
