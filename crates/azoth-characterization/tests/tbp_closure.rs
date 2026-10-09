//! Spec-driven tests for the `characterization.tbp_closure` model.
//!
//! Every expected value is a row of `validation/neqsim/captures/tbp_closure_probe.tsv`, which
//! `TbpClosureProbe` printed from the pinned jar.

use azoth_characterization::{TbpClosureKind, TbpClosureResult, TbpModel, model_gen, tbp_closure};
use azoth_core::units::{kelvins, kilograms_per_cubic_meter};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.tbp_closure";

fn call(case: &azoth_core::spec::TestCase) -> TbpClosureResult {
    let closure: TbpClosureKind = common::input_str(case, "closure")
        .parse()
        .expect("the case names a closure");
    let model = case.string("model").map(|name| {
        name.parse::<TbpModel>()
            .expect("the case names a TBP model")
    });
    tbp_closure(
        closure,
        kelvins(common::input(case, "boiling_point")),
        kilograms_per_cubic_meter(common::input(case, "density")),
        model,
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        common::assert_close(
            result.molar_mass.value,
            common::expected(case, "molar_mass"),
            case.tolerance,
            &format!("{}::{} (molar_mass)", spec.id, case.id),
        );
        common::assert_consistent(&result, &format!("{}::{}", spec.id, case.id));
    }
}

#[test]
fn the_1987_pair_is_not_the_1980_one() {
    // The two closed forms are different correlations over the same pair, so a port that
    // wired one into the other's branch would pass every case that named only one of them.
    let b1980 = tbp_closure(
        TbpClosureKind::RiaziDaubert1980,
        kelvins(500.0),
        kilograms_per_cubic_meter(850.0),
        None,
    )
    .expect("in range");
    let b1987 = tbp_closure(
        TbpClosureKind::RiaziDaubert1987,
        kelvins(500.0),
        kilograms_per_cubic_meter(850.0),
        None,
    )
    .expect("in range");
    assert!(
        (b1980.molar_mass.value - b1987.molar_mass.value).abs() > 1.0e-4,
        "the two pairs must differ: {} against {}",
        b1980.molar_mass.value,
        b1987.molar_mass.value
    );
}

#[test]
fn the_tbp_model_member_needs_a_model() {
    let refused = tbp_closure(
        TbpClosureKind::TbpModel,
        kelvins(500.0),
        kilograms_per_cubic_meter(850.0),
        None,
    );
    assert!(refused.is_err(), "a null model must be refused");
}

#[test]
fn an_unattainable_boiling_point_is_refused_rather_than_narrowed() {
    // 5000 K is far above anything the bracket reaches, so the endpoints share a sign and the
    // search must refuse instead of halving toward a bound the caller asked for.
    let refused = tbp_closure(
        TbpClosureKind::Soreide,
        kelvins(5000.0),
        kilograms_per_cubic_meter(850.0),
        None,
    );
    assert!(
        refused.is_err(),
        "5000 K is not attainable and must be refused"
    );
}
