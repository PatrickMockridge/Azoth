//! Spec-driven tests for the `eos.freezing_point` model.

// The oracle values are NeqSim's ground truth, so their full digits are the point.
#![allow(clippy::excessive_precision)]

use azoth_core::AzothError;
use azoth_core::units::pascals;
use azoth_eos::{freezing_point, model_gen, parahydrogen_solid};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.freezing_point";

fn names(case: &azoth_core::spec::TestCase) -> Vec<String> {
    case.list("components")
        .expect("the case declares components")
        .iter()
        .map(|name| (*name).to_string())
        .collect()
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let result = freezing_point(&names(case), pascals(common::input(case, "P")))
            .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.temperature.value,
            common::expected(case, "temperature"),
            case.tolerance,
            &format!("{context} (T)"),
        );
        // The residual is the equation being solved, so at the answer it is the solver's own
        // tolerance rather than a property of the answer - which is why the cases do not pin
        // it. Checked here so a solve that returned its starting temperature by accident
        // would not pass.
        assert!(
            result.residual.abs() < 1e-9,
            "{context}: the residual at the answer is {}",
            result.residual
        );
    }
}

/// The calibration is the triple point's own definition, and the raw solid does not satisfy it.
///
/// NeqSim computes the two constants at first use by comparing the raw solid against the
/// para-Leachman **liquid** at the triple point, and this is the check that they are what that
/// comparison gives: with them, the solid's Gibbs energy at the triple point is the liquid's.
#[test]
fn the_calibration_meets_the_liquid_at_the_triple_point() {
    let calibration = freezing_point::calibration();
    let t_tp = parahydrogen_solid::TRIPLE_POINT_TEMPERATURE;
    let p_tp = parahydrogen_solid::TRIPLE_POINT_PRESSURE;

    let fluid =
        freezing_point::residual(t_tp, p_tp, freezing_point::FluidRoot::Liquid, &calibration)
            .expect("the triple point has a dense root");
    assert!(
        fluid.abs() < 1e-9,
        "the calibrated residual at the triple point is {fluid}"
    );

    // And the raw solid misses it by the amount measured from the capture: `-3.12` J/mol of
    // Gibbs energy, which is `gibbs_shift` to four figures.
    let raw = parahydrogen_solid::properties(t_tp, p_tp);
    let uncalibrated = freezing_point::residual(
        t_tp,
        p_tp,
        freezing_point::FluidRoot::Liquid,
        &freezing_point::Calibration {
            gibbs_shift: 0.0,
            entropy_shift: 0.0,
        },
    )
    .expect("the triple point has a dense root");
    assert!(
        (uncalibrated - fluid).abs() > 1e-3,
        "the raw and calibrated residuals are the same to 1e-3, so the shift does nothing: \
         {uncalibrated} against {fluid}"
    );
    assert!(
        (calibration.gibbs_shift + 3.122_683).abs() < 1e-5,
        "the Gibbs shift is {}, and the capture's offset is -3.122683",
        calibration.gibbs_shift
    );
    let _ = raw;
}

/// A substance with no solid equation here is refused rather than approximated.
#[test]
fn a_substance_without_a_solid_equation_is_refused() {
    let error = freezing_point(&["methane".to_string()], pascals(1.0e5))
        .expect_err("methane has no solid equation in this crate");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );

    let two = freezing_point(
        &["para-hydrogen".to_string(), "argon".to_string()],
        pascals(1.0e5),
    )
    .expect_err("a freezing point is a pure substance's");
    assert!(matches!(two, AzothError::InvalidInput { .. }), "{two:?}");
}

/// The root rule is NeqSim's: a gas below the triple-point pressure, a liquid at and above it.
#[test]
fn the_fluid_root_follows_the_triple_point_pressure() {
    let calibration = freezing_point::calibration();
    let t_tp = parahydrogen_solid::TRIPLE_POINT_TEMPERATURE;
    let below = parahydrogen_solid::TRIPLE_POINT_PRESSURE * 0.5;
    let at = parahydrogen_solid::TRIPLE_POINT_PRESSURE;

    // Below it the gas root exists and the liquid one is not asked for; at it the liquid root
    // is the one the residual uses, and both are finite.
    let gas = freezing_point::residual(t_tp, below, freezing_point::FluidRoot::Gas, &calibration)
        .expect("a gas root below the triple-point pressure");
    let liquid =
        freezing_point::residual(t_tp, at, freezing_point::FluidRoot::Liquid, &calibration)
            .expect("a liquid root at the triple-point pressure");
    assert!(gas.is_finite() && liquid.is_finite());
    assert!(
        (gas - liquid).abs() > 1e-3,
        "the two roots give the same residual, so the rule is not being applied"
    );
}
