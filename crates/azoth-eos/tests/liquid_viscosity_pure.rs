//! Spec-driven tests for `eos.liquid_viscosity_pure`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::liquid_viscosity_pure::{LiquidViscosityLadder, liquid_viscosity_pure};
use azoth_eos::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "eos.liquid_viscosity_pure";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::LiquidViscosityPureResult {
    let form: LiquidViscosityLadder = common::input_str(case, "form").parse().unwrap();
    liquid_viscosity_pure(
        form,
        common::input(case, "model") as u32,
        common::input(case, "l1"),
        common::input(case, "l2"),
        common::input(case, "l3"),
        common::input(case, "l4"),
        kelvins(common::input(case, "Tc")),
        pascals(common::input(case, "Pc")),
        common::input(case, "omega"),
        kelvins(common::input(case, "T")),
        pascals(common::input(case, "P")),
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
                common::assert_close(
                    result.mu.value,
                    common::expected(case, "mu"),
                    case.tolerance,
                    &format!("{}::{} (mu)", spec.id, case.id),
                );
                common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
                common::assert_warnings_agree_with_spec(
                    spec,
                    &result.warnings,
                    |quantity| match quantity {
                        "Tc" => Some(common::input(case, "Tc")),
                        "Pc" => Some(common::input(case, "Pc")),
                        "T" => Some(common::input(case, "T")),
                        "P" => Some(common::input(case, "P")),
                        "mu" => Some(result.mu.value),
                        _ => None,
                    },
                    &format!("{}::{}", spec.id, case.id),
                );
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

/// **The two ladders differ in exactly one branch, and this is it.**
///
/// Every other model answers the same under both; model 2 is zero under the common-phase class
/// and `exp(l1 + l2/T)` under the liquid one. The test moves only the `form` and holds everything
/// else, so a port that had merged the ladders would fail here and nowhere else.
#[test]
fn the_ladders_agree_everywhere_but_model_two() {
    let run = |form: LiquidViscosityLadder, model: u32| {
        liquid_viscosity_pure(
            form,
            model,
            -7.811,
            3140.0,
            0.0,
            0.0,
            kelvins(425.12),
            pascals(3796000.0),
            0.2002,
            kelvins(300.0),
            pascals(2000000.0),
        )
        .expect("the ladder runs")
        .mu
        .value
    };
    // Model 2: the branch they disagree on.
    assert_eq!(run(LiquidViscosityLadder::CommonPhase, 2), 0.0);
    assert!(run(LiquidViscosityLadder::Liquid, 2) > 1.0e-3);
    // Every other model agrees, including the sentinels.
    for model in [1, 3, 4, 0] {
        assert_eq!(
            run(LiquidViscosityLadder::CommonPhase, model),
            run(LiquidViscosityLadder::Liquid, model),
            "model {model} is the same under both ladders"
        );
    }
}
