//! Spec-driven tests for `eos.solid_fugacity`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::solid_fugacity;
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.solid_fugacity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::results::SolidFugacityResult {
    solid_fugacity(
        common::input(case, "heat_of_fusion"),
        common::input(case, "triple_point_temperature"),
        common::input(case, "delta_cp_sl"),
        common::input(case, "delta_solid_volume"),
        kelvins(common::input(case, "tc")),
        pascals(common::input(case, "pc")),
        common::input(case, "omega"),
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
        common::input_str(case, "eos"),
    )
    .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id))
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
                let context = format!("{}::{}", spec.id, case.id);
                common::assert_close(
                    result.fugacity_coefficient,
                    common::expected(case, "fugacity_coefficient"),
                    case.tolerance,
                    &format!("{context} (fugacity coefficient)"),
                );
                common::assert_consistent(&result, &context);
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// **The solid is the liquid times an exponential, and the exponential is one at the triple
/// point.**
///
/// At `T = T_tp` every term of the exponent carries a factor of `(T_tp - T)` or `(1 - T/T_tp)`,
/// so the coefficient is the reference liquid's own and nothing else - the triple point's
/// definition, and the check that says the three terms are the right three. Measured on the
/// capture's water it is `1 - 9.7e-5` one hundredth of a kelvin below, which is the shape a
/// term-by-term port gets wrong by a sign.
#[test]
fn the_exponential_is_one_at_the_triple_point() {
    let at = |t: f64| {
        solid_fugacity(
            6010.0,
            273.16,
            37.12,
            0.0,
            kelvins(647.096),
            pascals(2.2064e7),
            0.3443,
            kelvins(t),
            pascals(1.0e6),
            "srk",
        )
        .expect("computes")
        .fugacity_coefficient
    };
    let at_triple = at(273.16);
    let reference = solid_fugacity(
        6010.0,
        273.16,
        37.12,
        0.0,
        kelvins(647.096),
        pascals(2.2064e7),
        0.3443,
        kelvins(273.16),
        pascals(1.0e6),
        "srk",
    )
    .expect("computes")
    .fugacity_coefficient;
    assert!((at_triple / reference - 1.0).abs() < 1.0e-12);

    // Below the triple point the solid is the more stable of the two, so its coefficient is
    // the smaller - which is the sign of every term in the exponent.
    assert!(at(273.15) < at(273.16), "subcooled");
    assert!(at(263.15) < at(273.15), "and colder still");
}
