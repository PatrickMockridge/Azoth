//! Spec-driven tests for the `eos.tp_solid_flash` model.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::{model_gen, tp_solid_flash};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.tp_solid_flash";

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty());

    for case in spec.cases {
        let names = case.list("components").expect("components");
        let result = tp_solid_flash(
            names,
            common::input_str(case, "solid"),
            kelvins(common::input(case, "T")),
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
            common::input_str(case, "eos"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        for (field, actual, expected) in [
            (
                "solid_fraction",
                result.solid_fraction,
                common::expected(case, "solid_fraction"),
            ),
            (
                "phase_count",
                f64::from(result.phase_count),
                common::expected(case, "phase_count"),
            ),
        ] {
            common::assert_close(
                actual,
                expected,
                case.tolerance,
                &format!("{context} ({field})"),
            );
        }
        // **The split closes to the fraction solve's own tolerance and not to round-off.**
        // The solid's share is the material balance's *remainder* at fractions the Newton
        // stopped on - `eos.tp_solid_flash`'s `algorithm.tolerance`, `1e-8` - so the sum
        // carries that and nothing finer. NeqSim's own states carry the same: `2e-11` at
        // 273.15 K, `4e-9` at the carbon-dioxide state.
        let total: f64 = result.beta.iter().sum();
        common::assert_close(
            total,
            1.0,
            result.residual.max(1.0e-09),
            &format!("{context} (the split sums to one)"),
        );
        for (index, row) in result.x.iter().enumerate() {
            common::assert_close(
                row.iter().sum::<f64>(),
                1.0,
                1.0e-09,
                &format!("{context} (phase {index} is a composition)"),
            );
        }
    }
}

/// **The three states the probe captured, layer by layer**, against NeqSim's own numbers.
///
/// The case pins the solid's share; this pins the gas the solid was taken from, which is the
/// half a wrong `E_solid` would leave right. The tolerance is `1e-6` for the reason the case
/// states: the outer loop stops on a `1e-3` change and takes at least four steps, so neither
/// implementation's answer is a converged one.
#[test]
fn the_gas_agrees_with_the_two_route_capture() {
    let names = ["water", "methane"];
    for (temperature, solid, water) in [
        (273.15, 0.499_785_776_882_074, 0.000_428_262_703_006_822),
        (253.15, 0.499_967_499_063_729, 6.499_760_686_357_46e-05),
    ] {
        let result = tp_solid_flash(
            &names,
            "water",
            kelvins(temperature),
            pascals(1.0e6),
            &[0.5, 0.5],
            "srk",
        )
        .expect("computes");
        common::assert_close(
            result.solid_fraction,
            solid,
            1.0e-06,
            &format!("the solid at {temperature} K"),
        );
        common::assert_close(
            result.x[0][0],
            water,
            1.0e-06,
            &format!("the gas's water at {temperature} K"),
        );
        // The gas is one phase here and the solid is the other: the aqueous the fluid-only
        // flash seeded was removed by the solve, and `phase_count` is where that shows.
        assert_eq!(result.phase_count, 2, "at {temperature} K");
    }
}

/// A feed the named substance is not in is refused rather than answered.
#[test]
fn a_solid_the_feed_does_not_carry_is_refused() {
    let error = tp_solid_flash(
        &["water", "methane"],
        "benzene",
        kelvins(273.15),
        pascals(1.0e6),
        &[0.5, 0.5],
        "srk",
    )
    .expect_err("no benzene, no benzene ice");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

/// **A component whose melt data the table does not carry is refused, where NeqSim computes
/// with the zero.**
///
/// `acetic acid`'s `HEATOFFUSION` is `0.0`, and zero is how this table spells an absence
/// everywhere else. NeqSim reads it as a value - the fusion term vanishes and the solid's
/// coefficient becomes the reference liquid's - so a substance that does not melt is
/// indistinguishable there from one whose heat of fusion nobody entered. The databank's own
/// rule is to refuse, and this is where that shows.
#[test]
fn a_substance_the_table_states_no_melt_for_is_refused() {
    let error = tp_solid_flash(
        &["water", "acetic acid"],
        "acetic acid",
        kelvins(273.15),
        pascals(1.0e6),
        &[0.5, 0.5],
        "srk",
    )
    .expect_err("no heat of fusion, no melt");
    assert!(
        matches!(error, azoth_core::AzothError::InvalidInput { .. }),
        "{error:?}"
    );
}

/// **A feed above the named component's triple point is answered, not refused**: the screen
/// finds a negative amount and the fluid flash is the whole of the answer.
///
/// Methane's triple point is `90.69` K, so a methane/methane-free feed at 273.15 K is well
/// above it. NeqSim's `checkAndAddSolidPhase` skips such a component outright; the same state
/// is reached here through the screen, which is the one place the two can disagree.
#[test]
fn above_the_triple_point_the_fluid_flash_is_the_answer() {
    let result = tp_solid_flash(
        &["methane", "ethane"],
        "methane",
        kelvins(273.15),
        pascals(1.0e6),
        &[0.5, 0.5],
        "srk",
    )
    .expect("computes");
    assert_eq!(result.solid_fraction, 0.0);
    assert_eq!(result.phase_count, 2);
}

/// **The solid's share is monotone in temperature, and the fugacity is what drives it.**
///
/// Ice's fugacity coefficient falls with temperature, so the amount that freezes must rise;
/// this checks the shape rather than the number, which the case already pins.
#[test]
fn colder_freezes_more() {
    let names = ["water", "methane"];
    let mut previous = 0.0;
    let mut previous_coefficient = f64::INFINITY;
    for step in 0..5 {
        let temperature = 273.15 - f64::from(step);
        let result = tp_solid_flash(
            &names,
            "water",
            kelvins(temperature),
            pascals(1.0e6),
            &[0.5, 0.5],
            "srk",
        )
        .expect("computes");
        assert!(
            result.solid_fraction > previous,
            "at {temperature} K the solid is {} and at {} K it was {previous}",
            result.solid_fraction,
            temperature + 1.0
        );
        assert!(
            result.solid_fugacity_coefficient < previous_coefficient,
            "at {temperature} K the coefficient is {}",
            result.solid_fugacity_coefficient
        );
        previous = result.solid_fraction;
        previous_coefficient = result.solid_fugacity_coefficient;
    }
}
