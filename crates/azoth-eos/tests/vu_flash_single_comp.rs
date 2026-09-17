//! Spec-driven tests for the `eos.vu_flash_single_comp` model.

use azoth_core::units::{cubic_meters_per_mole, joules_per_mole, pascals};
use azoth_eos::databank;
use azoth_eos::vu_flash_single_comp::vu_flash_single_comp;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.vu_flash_single_comp";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::VuFlashSingleCompResult {
    let (mixture, ideal_gas) =
        databank::mixture_of(case.list("components").expect("components"), None)
            .unwrap_or_else(|e| panic!("case `{}` should resolve but failed: {e}", case.id));
    vu_flash_single_comp(
        &mixture,
        &ideal_gas,
        pascals(common::input(case, "P")),
        cubic_meters_per_mole(common::input(case, "V")),
        joules_per_mole(common::input(case, "U")),
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
        common::assert_close(
            result.beta,
            common::expected(case, "beta"),
            case.tolerance,
            &format!("{context} (beta)"),
        );
        common::assert_consistent(&result, context);
        assert!(
            result.warnings.is_empty(),
            "{context}: the cases are states whose volume is consistent, so none warns"
        );
    }
}

/// The split is the lever rule, at any internal energy between the two saturated ends.
///
/// The case where `beta` is one half could pass on a model that always answers one half,
/// so this walks the span. The two ends are not inputs anyone knows - they are what the
/// model computes - so they are recovered from two calls and then used to predict a
/// third: `beta` is linear in `U` and the volume is linear in `beta`, and a model whose
/// answer is a *rule* reproduces both to the last bit.
#[test]
fn the_split_is_the_fraction_of_the_span_the_energy_is() {
    let (mixture, ideal_gas) = databank::mixture_of(&["propane"], None).expect("propane");
    let p = pascals(1.0e6);
    let at = |u: f64| {
        vu_flash_single_comp(
            &mixture,
            &ideal_gas,
            p,
            cubic_meters_per_mole(0.001),
            joules_per_mole(u),
        )
        .expect("a state inside the two-phase region has a split")
    };

    // Two probes, well inside, from which the span follows.
    let (u_a, u_b) = (-10_000.0, -5_000.0);
    let (a, b) = (at(u_a), at(u_b));
    let slope = (u_b - u_a) / (b.beta - a.beta);
    let u_liq = u_a - a.beta * slope;
    let u_vap = u_liq + slope;

    // *Just* inside each end, which is as close as the model may be asked: exactly at
    // an end the span recovered from two probes overshoots by its own round-off, and
    // the model refuses a state outside it - correctly. `v_liq` and `v_vap` are then
    // what the linear rule extrapolates to, which is the claim the volumes below test.
    let inset = 1.0e-6 * (u_vap - u_liq);
    for (u, wanted) in [(u_liq + inset, 1.0e-6), (u_vap - inset, 1.0 - 1.0e-6)] {
        let result = at(u);
        assert!(
            (result.beta - wanted).abs() < 1e-5,
            "at U = {u} the split is {} rather than {wanted}",
            result.beta
        );
    }
    let (v_liq, v_vap) = (at(u_liq + inset).v.value, at(u_vap - inset).v.value);

    for fraction in [0.1, 0.25, 0.5, 0.75, 0.9] {
        let u = u_liq + fraction * (u_vap - u_liq);
        let result = at(u);
        assert!(
            (result.beta - fraction).abs() < 1e-9,
            "U is {fraction} of the span but the split is {}",
            result.beta
        );
        let wanted = (1.0 - fraction) * v_liq + fraction * v_vap;
        // The endpoints themselves were recovered from probes at `1e-6` of the span,
        // so the linear rule carries that much of the span's width - `2e-9 m**3/mol`
        // here - and the tolerance is that rather than the arithmetic's last bit.
        assert!(
            (result.v.value - wanted).abs() < 1e-8,
            "the volume at {fraction} is {} but the lever rule gives {wanted}",
            result.v.value
        );
    }
}

/// A pressure at or above the critical one, and an energy outside the span, are refused.
///
/// Both are states this model has no split for rather than ones it failed to find, and
/// the error types say which: a caller asking a pure component for a two-phase state
/// above its critical pressure has asked a question with no answer, not hit a solver
/// limit.
#[test]
fn a_state_the_saturation_line_does_not_reach_is_refused() {
    let (mixture, ideal_gas) = databank::mixture_of(&["propane"], None).expect("propane");
    let pc = mixture.components()[0].pc.value;

    let above_critical = vu_flash_single_comp(
        &mixture,
        &ideal_gas,
        pascals(pc * 1.1),
        cubic_meters_per_mole(0.001),
        joules_per_mole(-5000.0),
    )
    .expect_err("a pressure above the critical one has no saturation temperature");
    assert!(
        matches!(above_critical, azoth_core::AzothError::OutOfRange { .. }),
        "got {above_critical:?}"
    );

    let subcooled = vu_flash_single_comp(
        &mixture,
        &ideal_gas,
        pascals(1.0e6),
        cubic_meters_per_mole(0.001),
        joules_per_mole(-1.0e6),
    )
    .expect_err("an energy below the saturated liquid's has no split");
    assert!(
        matches!(subcooled, azoth_core::AzothError::OutOfRange { .. }),
        "got {subcooled:?}"
    );
}

/// A mixture is not this model's subject, and it says which model is.
#[test]
fn a_mixture_is_refused_by_name() {
    let (mixture, ideal_gas) =
        databank::mixture_of(&["methane", "n-butane"], None).expect("the pair");
    let error = vu_flash_single_comp(
        &mixture,
        &ideal_gas,
        pascals(1.0e6),
        cubic_meters_per_mole(0.001),
        joules_per_mole(-5000.0),
    )
    .expect_err("a mixture has no single saturation line");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "got {error:?}"
    );
}
