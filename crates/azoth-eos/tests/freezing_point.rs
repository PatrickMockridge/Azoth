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
        let result = freezing_point(
            &names(case),
            case.vector("z").expect("the case declares z"),
            common::input_str(case, "solid"),
            pascals(common::input(case, "P")),
        )
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

/// A candidate that is not one of the fluid's own components is refused.
#[test]
fn a_candidate_the_fluid_does_not_have_is_refused() {
    let error = freezing_point(
        &["methane".to_string()],
        &[1.0],
        "para-hydrogen",
        pascals(1.0e5),
    )
    .expect_err("the fluid has no para-hydrogen in it");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

/// A `z` that does not match the component list is refused rather than zipped short.
#[test]
fn a_composition_of_the_wrong_length_is_refused() {
    let error = freezing_point(
        &["methane".to_string()],
        &[1.0, 0.0],
        "methane",
        pascals(1.0e5),
    )
    .expect_err("one component and two fractions");
    assert!(
        matches!(error, AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

/// **Methane does not freeze, and that is NeqSim's own guard rather than an absence.**
///
/// `ComponentSolid.fugcoef` returns `1e30` for methane before any arithmetic, so its solid is
/// infinitely volatile and the residual has no sign change. Measured: the search refuses at 1,
/// 10 and 50 bara rather than reporting methane's triple point of 90.69 K, which is what the
/// tabulated route would otherwise solve to.
#[test]
fn a_methane_candidate_is_refused() {
    for p_bar in [1.0, 10.0, 50.0] {
        let error = freezing_point(
            &["methane".to_string()],
            &[1.0],
            "methane",
            pascals(p_bar * 1.0e5),
        )
        .expect_err("methane has no solid in NeqSim");
        assert!(
            matches!(error, AzothError::SolverNotConverged { .. }),
            "{p_bar} bar: {error:?}"
        );
    }
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

/// The feed NeqSim's own freezing test uses, as mole fractions.
const LNG: [&str; 5] = ["CO2", "nitrogen", "methane", "ethane", "propane"];

const LNG_Z: [f64; 5] = [
    0.089_484_367_947_023_1,
    0.579_634_022_102_985,
    0.170_546_677_734_326,
    0.144_227_745_985_202,
    0.016_107_186_230_464_2,
];

/// **The tabulated route, against NeqSim's own answer rather than its own test.**
///
/// `FreezingPointTemperatureFlashTest.testLNGFreezingPointFlashAfterFluidOnlyTPFlash` asserts
/// only that the answer lies between 90 and 220 K; the capture prints it. Measured, this
/// reproduces `184.710072436643` K at 5 bara to `3.4e-9` relative and `194.118909895905` K at
/// 50 bara to `8.2e-8`.
#[test]
fn the_tabulated_route_reproduces_the_capture() {
    let names: Vec<String> = LNG.iter().map(|name| (*name).to_string()).collect();
    for (p_bar, want) in [(5.0, 184.710_072_436_643), (50.0, 194.118_909_895_905)] {
        let result = freezing_point(&names, &LNG_Z, "CO2", pascals(p_bar * 1.0e5))
            .unwrap_or_else(|e| panic!("{p_bar} bar: {e:?}"));
        assert_eq!(
            result.component, "CO2",
            "{p_bar} bar: the controlling component"
        );
        assert!(
            (result.temperature.value / want - 1.0).abs() < 1.0e-7,
            "{p_bar} bar: {} against NeqSim's {want}",
            result.temperature.value
        );
    }
}

/// **At 20 bara NeqSim refuses and this answers, and the answer is suspect.**
///
/// NeqSim's own failure is `the freezing-point bracket for CO2 collapsed without satisfying
/// Gibbs equilibrium; the fluid residual is discontinuous or the requested density root is
/// unavailable`. This solves to `194.8203` K there - which is **above** the 50 bara answer of
/// `194.1189`, and a freezing point should rise with pressure, so this is a second root of the
/// residual rather than the physical one. Recorded as a divergence rather than pinned as an
/// answer: a case that asserted it would be asserting a number this test distrusts.
#[test]
fn the_20_bara_state_neqsim_refuses_is_recorded_and_not_trusted() {
    let names: Vec<String> = LNG.iter().map(|name| (*name).to_string()).collect();
    let at_20 = freezing_point(&names, &LNG_Z, "CO2", pascals(20.0 * 1.0e5))
        .expect("this solves where NeqSim's bracket collapses");
    let at_50 =
        freezing_point(&names, &LNG_Z, "CO2", pascals(50.0 * 1.0e5)).expect("50 bara solves");
    assert!(
        (at_20.temperature.value / 194.820_296_026_635_28 - 1.0).abs() < 1.0e-7,
        "the 20 bara answer is now {}",
        at_20.temperature.value
    );
    assert!(
        at_20.temperature.value > at_50.temperature.value,
        "the 20 bara answer is {} and the 50 bara one {}, so the 20 bara root is no longer          above it and the divergence this test records has changed shape",
        at_20.temperature.value,
        at_50.temperature.value
    );
}
