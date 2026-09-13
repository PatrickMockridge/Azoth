//! Spec-driven tests for `eos.pr_departure`.

use azoth_core::AzothError;
use azoth_eos::pr_departure;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.pr_departure";

/// The three outputs the spec declares, in order.
const FIELDS: [&str; 3] = ["ln_phi", "h_dep_rt", "s_dep_r"];

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PrDepartureResult {
    pr_departure(
        common::input(case, "a_reduced"),
        common::input(case, "b_reduced"),
        common::input(case, "z"),
        common::input(case, "kappa"),
        common::input(case, "Tr"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
}

fn field(result: &azoth_eos::PrDepartureResult, name: &str) -> f64 {
    match name {
        "ln_phi" => result.ln_phi,
        "h_dep_rt" => result.h_dep_rt,
        _ => result.s_dep_r,
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        match case.kind {
            "worked_example" | "reference" => {
                let result = call(case);
                for name in FIELDS {
                    common::assert_close(
                        field(&result, name),
                        common::expected(case, name),
                        case.tolerance,
                        &format!("{}::{} ({name})", spec.id, case.id),
                    );
                }
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "a_reduced" => Some(common::input(case, "a_reduced")),
                        "b_reduced" => Some(common::input(case, "b_reduced")),
                        "z" => Some(common::input(case, "z")),
                        "kappa" => Some(common::input(case, "kappa")),
                        "Tr" => Some(common::input(case, "Tr")),
                        "z_minus_b_reduced" => {
                            Some(common::input(case, "z") - common::input(case, "b_reduced"))
                        }
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
            }
            "property" => match case.property {
                Some("consistency_with") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 4,
        "expected several active cases, ran {executed}"
    );
}

/// `h_dep_rt - s_dep_r` equals `ln_phi`, swept rather than spot-checked.
///
/// The identity is exact in real arithmetic - the two `psi` terms cancel, and the
/// spec's `notes` show the three lines - so any disagreement is rounding.
/// It is the only check that catches a sign error in either departure function:
/// flipping one leaves `ln_phi` untouched and moves the other to a value that is
/// still plausible as an enthalpy or entropy departure.
///
/// Asserted at 1e-12, looser than the identity deserves algebraically and for a
/// measured reason: on the liquid root both departures are near -6.6 and cancel to
/// -0.32, losing about two digits. 1.3e-15 is what that case achieves.
#[test]
fn the_gibbs_identity_holds() {
    for a in [0.05, 0.20206500174625697, 0.4572355289213822] {
        for b in [0.005, 0.02431127309496514, 0.07779607390388846] {
            for z in [0.05, 0.3, 0.79, 0.95] {
                if z <= b {
                    continue;
                }
                let r = pr_departure(a, b, z, 0.60282728832, 0.8).unwrap();
                let identity = r.h_dep_rt - r.s_dep_r;
                assert!(
                    (identity - r.ln_phi).abs() < 1e-12,
                    "A = {a}, B = {b}, z = {z}: h - s = {identity} against ln_phi = {}",
                    r.ln_phi
                );
            }
        }
    }
}

/// As `A` and `B` go to zero at `z = 1`, all three outputs go to zero - linearly.
///
/// The one check that constrains all three expressions using no reference values at
/// all. It also pins `ln(z - B)` approaching `ln(1) = 0` rather than diverging,
/// which a stray additive constant in either departure function would break.
#[test]
fn the_low_pressure_limit_is_the_ideal_gas() {
    let mut previous: Option<[f64; 3]> = None;
    for exponent in [3.0, 5.0, 7.0, 9.0] {
        let b = 10.0_f64.powf(-exponent);
        let a = 8.0 * b;
        let r = pr_departure(a, b, 1.0, 0.6, 0.8).unwrap();
        let magnitudes = [r.ln_phi.abs(), r.h_dep_rt.abs(), r.s_dep_r.abs()];

        if let Some(previous) = previous {
            // Two decades of B should shrink each output by about two decades. The
            // band is generous - a factor of 50 to 200 rather than exactly 100 -
            // because what is claimed is *linearity in B*, not a particular
            // coefficient. Measured, the factor is 99.99 per decade.
            for (name, (before, after)) in FIELDS.iter().zip(
                previous
                    .iter()
                    .copied()
                    .zip(magnitudes)
                    .collect::<Vec<(f64, f64)>>(),
            ) {
                let shrink = before / after;
                assert!(
                    (50.0..=200.0).contains(&shrink),
                    "{name} should shrink roughly 100-fold per decade of B, shrank {shrink}"
                );
            }
        }
        previous = Some(magnitudes);
    }

    // And at the smallest B the outputs are genuinely negligible. Threshold taken
    // from the measurement - 1.2e-8 at B = 1e-9 - rather than picked to look small.
    let b = 1e-9_f64;
    let r = pr_departure(8.0 * b, b, 1.0, 0.6, 0.8).unwrap();
    for (name, magnitude) in FIELDS
        .iter()
        .zip([r.ln_phi.abs(), r.h_dep_rt.abs(), r.s_dep_r.abs()])
    {
        assert!(
            magnitude < 1e-7,
            "{name} should be negligible at B = 1e-9, but is {magnitude:e}"
        );
    }
}

#[test]
fn a_non_positive_b_reduced_is_an_error() {
    for b in [0.0, -0.01] {
        let err = pr_departure(0.2, b, 0.79, 0.6, 0.8).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("b_reduced"), "for B = {b}");
    }
}

#[test]
fn a_z_at_or_below_b_reduced_is_an_error_naming_the_difference() {
    // The bound is on `z - B`, not on `z`, because that is the quantity
    // `ln(z - B)` needs to be positive. `eos.pr_z_factor` already discards roots at
    // or below B, so this is for the caller who supplied `z` some other way - and
    // it names the failure rather than letting a domain error do it.
    for z in [0.02431127309496514, 0.01, 0.0, -1.0] {
        let err = pr_departure(0.2, 0.02431127309496514, z, 0.6, 0.8).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }));
        assert_eq!(err.field(), Some("z_minus_b_reduced"), "for z = {z}");
    }
}

#[test]
fn kappa_and_tr_do_not_reach_the_fugacity_coefficient() {
    // The quiet failure the spec warns about, asserted so it stays visible: `kappa`
    // and `Tr` enter only through `psi`, which the `ln_phi` expression does not use.
    // A caller who mixes states therefore gets an *unchanged* fugacity coefficient
    // and shifted departures, rather than an obviously wrong set.
    let base = pr_departure(
        0.20206500174625697,
        0.02431127309496514,
        0.79,
        0.60282728832,
        0.8,
    )
    .unwrap();
    for (kappa, tr) in [(0.2, 0.8), (0.60282728832, 1.4), (0.9, 0.5)] {
        let other =
            pr_departure(0.20206500174625697, 0.02431127309496514, 0.79, kappa, tr).unwrap();
        assert_eq!(
            other.ln_phi.to_bits(),
            base.ln_phi.to_bits(),
            "ln_phi should not depend on kappa or Tr, but ({kappa}, {tr}) moved it"
        );
        assert!(
            (other.h_dep_rt - base.h_dep_rt).abs() > 1e-9,
            "h_dep_rt should move with ({kappa}, {tr})"
        );
    }
}
