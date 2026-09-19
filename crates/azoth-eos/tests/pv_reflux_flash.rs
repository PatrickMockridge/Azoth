//! Spec-driven tests for the `eos.pv_reflux_flash` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::pv_reflux_flash::{RefluxPhase, pv_reflux_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pv_reflux_flash";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PvRefluxFlashResult {
    let (mixture, _) = databank::mixture_of(
        case.list("components").expect("components"),
        Cubic::Pr,
        None,
    )
    .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    let phase = match common::input_str(case, "phase") {
        "vapour" => RefluxPhase::Vapour,
        _ => RefluxPhase::Liquid,
    };
    pv_reflux_flash(
        &mixture,
        pascals(common::input(case, "P")),
        common::input(case, "reflux"),
        phase,
        kelvins(common::input(case, "temperature")),
        case.vector("z").expect("z"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model");
    assert!(!spec.cases.is_empty(), "the model should have cases");
    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.t.value,
            common::expected(case, "T"),
            case.tolerance,
            &format!("{context} (T)"),
        );
        assert!(
            result.residual <= spec.algorithm.expect("a procedure").tolerance,
            "{context}: the ratio residual is {:e}, beyond the declared tolerance",
            result.residual
        );
        common::assert_consistent(&result, context);
    }
}

/// The phase named decides the answer, and the two ratios are reciprocals.
///
/// A model that ignored `phase` would answer the two halves of this with the same
/// temperature. At 330 K and 25 bar the vapour ratio is 0.187 and the liquid's 5.337, so
/// asking for each must give two different temperatures - and the one at which the
/// liquid ratio is the *reciprocal* of the vapour's is the same temperature.
#[test]
fn the_phase_named_decides_the_answer() {
    let (mixture, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("the pair");
    let z = [0.6, 0.4];
    let p = pascals(2_500_000.0);

    let vapour = pv_reflux_flash(
        &mixture,
        p,
        0.187_358_600_133_836_1,
        RefluxPhase::Vapour,
        kelvins(330.0),
        &z,
    )
    .expect("the vapour ratio has a temperature");
    let liquid = pv_reflux_flash(
        &mixture,
        p,
        5.337_358_409_412_052,
        RefluxPhase::Liquid,
        kelvins(330.0),
        &z,
    )
    .expect("the liquid ratio has a temperature");
    assert!(
        (vapour.t.value - liquid.t.value).abs() < 1.0e-4,
        "reciprocal ratios should be the same state: {} against {}",
        vapour.t.value,
        liquid.t.value
    );

    // A vapour ratio of one is the same request as a liquid ratio of one - equal phase
    // amounts - so the two must agree where neither is a round trip.
    let half_vapour = pv_reflux_flash(&mixture, p, 1.0, RefluxPhase::Vapour, kelvins(330.0), &z)
        .expect("a ratio of one");
    let half_liquid = pv_reflux_flash(&mixture, p, 1.0, RefluxPhase::Liquid, kelvins(330.0), &z)
        .expect("a ratio of one");
    // Loose because the two are separate searches that stop on their own residual:
    // each lands within the declared tolerance of the same ratio, and the two
    // tolerances buy a few microkelvin of difference in where that is.
    assert!(
        (half_vapour.t.value - half_liquid.t.value).abs() < 1.0e-4,
        "a ratio of one is equal phase amounts, whichever phase it is of: {} against {}",
        half_vapour.t.value,
        half_liquid.t.value
    );
    assert!((half_vapour.beta.expect("a split") - 0.5).abs() < 1.0e-6);
}
