//! Spec-driven tests for the `eos.hydrate_equilibrium_line` model.
//!
//! The oracle is `validation/neqsim/captures/hydrate_equilibrium_probe.tsv`: NeqSim's
//! `HydrateEquilibriumLine` over the same fluid and the same bounds. The capture also carries
//! the same ten pressures solved **independently**, which is what makes the seeding's absence
//! checkable rather than asserted - the two vectors agree to every printed digit.
//!
//! **The agreement across the grid is not uniform**: measured against NeqSim it runs from
//! `7.6e-10` at 111.6 bara to `1.6e-7` at 177.9, so the bar here is the case's own `1e-6`
//! rather than the tightest the best point manages.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::units::pascals;
use azoth_eos::hydrate::{self, HydrateModel};
use azoth_eos::{Cubic, hydrate_equilibrium_line, hydrate_formation_temperature, model_gen};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.hydrate_equilibrium_line";

/// The probe's bounds, in bara: `validation/neqsim/HydrateEquilibriumProbe.java`.
const P_MIN_BARA: f64 = 1.0;

const P_MAX_BARA: f64 = 200.0;

/// The capture's ten temperatures, in grid order, at the probe's feed.
const TEMPERATURES: [f64; 10] = [
    258.099582916960,
    282.316876654403,
    287.762981872390,
    290.688987372908,
    292.552200127104,
    293.861737270018,
    294.863816418626,
    295.690077826984,
    296.410250673584,
    297.061002920514,
];

/// The probe's own feed: methane/ethane/propane/water at `79.0 10.10 2.050 10.0` moles.
const Z: [f64; 4] = [
    0.781_018_289_668_808_7,
    0.099_851_705_388_037_56,
    0.020_266_930_301_532_374,
    0.098_863_074_641_621_35,
];

fn mixture(model: HydrateModel) -> azoth_eos::mixture::Mixture {
    hydrate::hydrate_mixture_of(
        &["methane", "ethane", "propane", "water"],
        Cubic::Srk,
        None,
        model,
    )
    .expect("resolves")
    .0
}

#[test]
fn the_line_reproduces_the_capture() {
    let fluid = mixture(HydrateModel::Pvtsim);
    let line = hydrate_equilibrium_line(
        &fluid,
        pascals(P_MIN_BARA * 1.0e5),
        pascals(P_MAX_BARA * 1.0e5),
        &Z,
    )
    .expect("the grid solves");

    assert_eq!(line.temperature.len(), 10, "the grid is ten points");
    assert_eq!(line.pressure.len(), 10);

    // **The grid**, which is NeqSim's own stepping: `min + dp i` over ten points, landing on
    // the maximum at the last one. A port that used `max - dp i` would walk the other way and
    // still hit both ends.
    for (point, pressure) in line.pressure.iter().enumerate() {
        let want = (P_MIN_BARA + (P_MAX_BARA - P_MIN_BARA) / 9.0 * point as f64) * 1.0e5;
        common::assert_close(
            *pressure,
            want,
            1.0e-6,
            &format!("grid point {point}'s pressure"),
        );
        common::assert_close(
            line.temperature[point],
            TEMPERATURES[point],
            1.0e-6,
            &format!("grid point {point}'s temperature"),
        );
    }
}

/// **The curve rises with the pressure**, which is the shape the model is used for and the
/// premise NeqSim's seeding rests on.
#[test]
fn the_temperature_rises_with_the_pressure() {
    let fluid = mixture(HydrateModel::Pvtsim);
    let line = hydrate_equilibrium_line(
        &fluid,
        pascals(P_MIN_BARA * 1.0e5),
        pascals(P_MAX_BARA * 1.0e5),
        &Z,
    )
    .expect("the grid solves");
    for pair in line.temperature.windows(2) {
        assert!(
            pair[1] > pair[0],
            "the line is not monotone: {} then {}",
            pair[0],
            pair[1]
        );
    }
}

/// A bound that does not make a grid is refused rather than solved around.
#[test]
fn a_grid_that_is_not_a_grid_is_refused() {
    let fluid = mixture(HydrateModel::Pvtsim);
    for (low, high) in [(200.0, 1.0), (1.0, 1.0), (0.0, 200.0), (-1.0, 200.0)] {
        let error =
            hydrate_equilibrium_line(&fluid, pascals(low * 1.0e5), pascals(high * 1.0e5), &Z)
                .expect_err("not a grid");
        // A non-positive bound is refused by the spec's own range check and a reversed pair
        // by this model's, so the two kinds are both refusals and the test says so.
        assert!(
            matches!(
                error,
                azoth_core::AzothError::InvalidInput { .. }
                    | azoth_core::AzothError::OutOfRange { .. }
            ),
            "{low} to {high}: {error:?}"
        );
    }
}

/// Every case in the spec runs, and the first and last points are what the source states.
#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());
    for case in spec.cases {
        let names = case
            .list("components")
            .expect("the case declares components");
        let fluid = hydrate::hydrate_mixture_of(names, Cubic::Srk, None, HydrateModel::Pvtsim)
            .expect("the case's components resolve for a hydrate")
            .0;
        let line = hydrate_equilibrium_line(
            &fluid,
            pascals(common::input(case, "P_min")),
            pascals(common::input(case, "P_max")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));
        assert_eq!(line.temperature.len(), 10);
        common::assert_close(
            line.temperature[0],
            258.099582916960,
            1.0e-6,
            &format!("{}::{} first point", spec.id, case.id),
        );
        common::assert_close(
            line.temperature[9],
            297.061002920514,
            1.0e-6,
            &format!("{}::{} last point", spec.id, case.id),
        );
    }
}

/// **The `guo_finch` route is not this grid over another model, and this records what it is.**
///
/// Its reference is the empty lattice's vapour pressure with no reference-water-fugacity term,
/// and NeqSim only ever reaches it through a `CPA-SRK-EOS`-named system - so on a fluid this
/// library can build it solves a *different* equilibrium rather than a better one. Measured
/// against the same feed: at 100 bara `pvtsim` gives `293.2277` K and `guo_finch` `259.3619`,
/// which is the electrolyte fluid's absence showing through and not a correction.
///
/// It also fails to bracket at the ends of this grid, at 1 and 200 bara. Pinned here rather
/// than left in prose, so that the assumptions and the behaviour cannot drift apart.
#[test]
fn the_guo_finch_route_is_a_different_equilibrium() {
    let pvtsim = hydrate_equilibrium_line(
        &mixture(HydrateModel::Pvtsim),
        pascals(P_MIN_BARA * 1.0e5),
        pascals(P_MAX_BARA * 1.0e5),
        &Z,
    )
    .expect("the grid solves");

    let guo_finch = hydrate_equilibrium_line(
        &mixture(HydrateModel::GuoFinch),
        pascals(P_MIN_BARA * 1.0e5),
        pascals(P_MAX_BARA * 1.0e5),
        &Z,
    );
    assert!(
        matches!(
            guo_finch,
            Err(azoth_core::AzothError::SolverNotConverged { .. })
        ),
        "the Guo-Finch grid used to refuse at 1 bara; it now gives {guo_finch:?}"
    );

    // At the middle of the grid it does bracket, and lower - which is the difference being
    // recorded rather than a tolerance being met.
    let at_100_bara =
        hydrate_formation_temperature(&mixture(HydrateModel::GuoFinch), pascals(100.0 * 1.0e5), &Z)
            .expect("the middle of the grid brackets");
    assert!(
        at_100_bara.temperature.value < pvtsim.temperature[5] - 20.0,
        "the Guo-Finch temperature at 100 bara is {}, against PVTsim's {}",
        at_100_bara.temperature.value,
        pvtsim.temperature[5]
    );
}
