//! Spec-driven tests for the `eos.pure_saturation` model.

use azoth_core::units::{kelvins, pascals};
use azoth_core::{AzothError, CalcResult};
use azoth_eos::{model_gen, pure_saturation};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pure_saturation";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::PureSaturationResult {
    pure_saturation(
        kelvins(common::input(case, "Tc")),
        pascals(common::input(case, "Pc")),
        common::input(case, "omega"),
        kelvins(common::input(case, "T")),
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
            result.p_sat.value,
            common::expected(case, "p_sat"),
            case.tolerance,
            &format!("{}::{} (p_sat)", spec.id, case.id),
        );
        assert!(
            result.iterations > 0,
            "{}::{}: the search should have taken at least one step",
            spec.id,
            case.id
        );
        assert!(
            result.residual <= spec.algorithm.expect("a procedure").tolerance,
            "{}::{}: residual {:e} exceeds the declared tolerance {:e}",
            spec.id,
            case.id,
            result.residual,
            spec.algorithm.expect("a procedure").tolerance
        );
    }
}

/// The residual vanishes at the returned pressure, which is what the search is for.
///
/// The check that does not depend on knowing the answer, and the one that matters
/// for a model: `p_sat` is only meaningful if the two fugacities actually agree
/// there, and a bisection that took a wrong branch or stopped early would return a
/// plausible pressure with a residual nowhere near zero.
///
/// Recomputed here from the kernels rather than read from the result, so a result
/// that reported a residual it had not achieved would fail.
#[test]
fn the_fugacities_agree_at_the_returned_pressure() {
    for (tc, pc, w, t) in [
        (369.83, 4_248_000.0, 0.1523, 300.0),
        (304.13, 7_377_000.0, 0.2239, 280.0),
        (425.12, 3_796_000.0, 0.2002, 350.0),
        (190.56, 4_599_200.0, 0.01142, 150.0),
    ] {
        let r = pure_saturation(kelvins(tc), pascals(pc), w, kelvins(t)).unwrap();

        // Rebuild the state at the answer, through the same kernels the model uses.
        let tr = t / tc;
        let pr = r.p_sat.value / pc;
        let kappa = azoth_eos::pr_kappa(w).unwrap().kappa;
        let ab = azoth_eos::pr_alpha_ab(kappa, tr, pr).unwrap();
        let z = azoth_eos::pr_z_factor(ab.a_reduced, ab.b_reduced).unwrap();
        let liquid =
            azoth_eos::pr_departure(ab.a_reduced, ab.b_reduced, z.z_min, kappa, tr).unwrap();
        let vapour =
            azoth_eos::pr_departure(ab.a_reduced, ab.b_reduced, z.z_max, kappa, tr).unwrap();

        assert!(
            (liquid.ln_phi - vapour.ln_phi).abs() < 1e-9,
            "Tc={tc}, T={t}: the fugacities differ by {:e} at the returned pressure",
            liquid.ln_phi - vapour.ln_phi
        );
        // And the reported `ln_phi` is the one the kernels produce there.
        assert!((r.ln_phi - liquid.ln_phi).abs() < 1e-12);
    }
}

/// The answer lies strictly inside the bracket the scan found.
///
/// The bracket's upper end is the spinodal, and this asserts the search never
/// returns a pressure at or above it - where there is no liquid branch and the
/// fugacity equality has no meaning. A scan whose upper bound overshot, or a
/// bisection that walked past the root, would return a number here that passed every
/// point-value case.
#[test]
fn the_answer_is_below_the_spinodal() {
    let (tc, pc, w, t) = (369.83, 4_248_000.0, 0.1523, 300.0);
    let r = pure_saturation(kelvins(tc), pascals(pc), w, kelvins(t)).unwrap();
    let tr = t / tc;
    let kappa = azoth_eos::pr_kappa(w).unwrap().kappa;

    // Walk up until the cubic stops having a liquid branch, which is the spinodal.
    let bracket = model_gen::PURE_SATURATION_SPEC
        .algorithm
        .expect("a procedure")
        .bracket
        .expect("this model's scheme brackets");
    let mut spinodal = bracket.lower;
    for step in 0..bracket.steps {
        let fraction = f64::from(step) / f64::from(bracket.steps - 1);
        let pr = bracket.lower + (bracket.upper - bracket.lower) * fraction;
        let ab = azoth_eos::pr_alpha_ab(kappa, tr, pr).unwrap();
        if azoth_eos::pr_z_factor(ab.a_reduced, ab.b_reduced)
            .unwrap()
            .root_structure
            == azoth_eos::RootStructure::Three
        {
            spinodal = pr;
        }
    }
    assert!(
        r.p_sat.value < spinodal * pc,
        "the answer {} should be below the spinodal {}",
        r.p_sat.value,
        spinodal * pc
    );
    assert!(r.p_sat.value > 0.0);
}

#[test]
fn a_temperature_at_or_above_the_critical_is_refused() {
    // The bound is on the ratio. At exactly `Tc` the two roots have merged and the
    // residual is zero everywhere, so the search would converge on nothing - hence
    // exclusive, and hence this asserting the boundary as well as past it.
    for t in [369.83, 370.0, 500.0] {
        let err =
            pure_saturation(kelvins(369.83), pascals(4_248_000.0), 0.1523, kelvins(t)).unwrap_err();
        assert!(matches!(err, AzothError::OutOfRange { .. }), "for T = {t}");
        assert_eq!(err.field(), Some("t_over_tc"));
    }
}

#[test]
fn a_non_positive_input_is_an_error() {
    for (tc, pc, t) in [
        (0.0, 4_248_000.0, 300.0),
        (369.83, 0.0, 300.0),
        (369.83, 4_248_000.0, 0.0),
        (369.83, 4_248_000.0, -1.0),
    ] {
        assert!(
            pure_saturation(kelvins(tc), pascals(pc), 0.1523, kelvins(t)).is_err(),
            "Tc={tc}, Pc={pc}, T={t}"
        );
    }
}

#[test]
fn the_result_is_clean_at_an_ordinary_state() {
    // No warning bounds are declared for this model, so a valid state is silent.
    // Asserted so a bound added without a rationale shows up here.
    let r = pure_saturation(
        kelvins(369.83),
        pascals(4_248_000.0),
        0.1523,
        kelvins(300.0),
    )
    .unwrap();
    assert!(r.is_clean(), "unexpected warnings: {:?}", r.warnings);
}
