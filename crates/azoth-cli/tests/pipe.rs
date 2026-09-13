//! Integration tests for the CLI's composition.
//!
//! The interesting logic in `azoth pipe` is not the argument parsing: it is
//! turning a flow and a bore into a velocity, choosing a friction factor, adding
//! the fitting loss by a *different* method from the straight-pipe one, and
//! deciding which warnings survive to the report. That is what these exercise.

use azoth_cli::pipe::{self, FlowUnit, FrictionMethod};
use azoth_core::WarningCode;

/// The command from the brief, as a test.
fn headline() -> pipe::PipeResult {
    pipe::compute(
        "water",
        20.0,
        10.0,
        FlowUnit::CubicMetresPerHour,
        0.1,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &["90_elbow".to_string(), "gate_valve_open".to_string()],
        FrictionMethod::Colebrook,
    )
    .expect("the headline example should compute")
}

#[test]
fn the_headline_example_computes() {
    let result = headline();

    // Velocity from 10 m3/h through a 100 mm bore: A = 7.854e-3 m2, so
    // v = 0.0027778 / 7.854e-3 = 0.35368 m/s.
    assert!(
        (result.velocity.value - 0.353678).abs() < 1e-5,
        "velocity was {}",
        result.velocity.value
    );
    // Re = rho v D / mu.
    assert!((result.re - 35_233.6).abs() < 1.0, "Re was {}", result.re);
    assert!(result.dp_total > 0.0);
    assert!(result.dp_straight > 0.0);
    assert!(result.dp_fittings > 0.0);
}

#[test]
fn total_is_the_sum_of_the_two_contributions() {
    // The composition is the thing this tool does that no single calc does, so
    // it is worth asserting the arithmetic rather than just that a number came
    // out.
    let result = headline();
    let expected = result.dp_straight + result.dp_fittings;
    assert!(
        (result.dp_total - expected).abs() < 1e-9,
        "total {} is not straight {} + fittings {}",
        result.dp_total,
        result.dp_straight,
        result.dp_fittings
    );
}

#[test]
fn no_fittings_means_no_fitting_loss() {
    let with = headline();
    let without = pipe::compute(
        "water",
        20.0,
        10.0,
        FlowUnit::CubicMetresPerHour,
        0.1,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &[],
        FrictionMethod::Colebrook,
    )
    .unwrap();

    assert_eq!(without.dp_fittings, 0.0);
    assert_eq!(without.k_total, 0.0);
    // And the straight-pipe term is unaffected by whether fittings were listed.
    assert_eq!(without.dp_straight, with.dp_straight);
    assert!(without.dp_total < with.dp_total);
}

#[test]
fn mass_and_volumetric_flow_describe_the_same_state() {
    // 10 m3/h of water at 998.2 kg/m3 is 2.7778 kg/s. The two must give the same
    // velocity, or the kg/s branch is converting wrongly.
    let volumetric = headline();
    let mass = pipe::compute(
        "water",
        20.0,
        10.0 * 998.2 / 3600.0,
        FlowUnit::KilogramsPerSecond,
        0.1,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &[],
        FrictionMethod::Colebrook,
    )
    .unwrap();

    assert!(
        (mass.velocity.value - volumetric.velocity.value).abs() < 1e-12,
        "mass flow gave {} m/s, volumetric gave {}",
        mass.velocity.value,
        volumetric.velocity.value
    );
}

#[test]
fn the_explicit_and_implicit_friction_factors_agree_within_the_published_claim() {
    // The CLI offers a choice, so the two must not silently disagree by more
    // than the approximation is documented to cost.
    let colebrook = headline();
    let swamee = pipe::compute(
        "water",
        20.0,
        10.0,
        FlowUnit::CubicMetresPerHour,
        0.1,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &[],
        FrictionMethod::SwameeJain,
    )
    .unwrap();

    let relative =
        (swamee.friction_factor - colebrook.friction_factor).abs() / colebrook.friction_factor;
    assert!(relative < 0.01, "friction factors differ by {relative:.4}");
    assert!(
        swamee.iterations.is_none(),
        "the explicit form must not report iterations"
    );
    assert!(colebrook.iterations.is_some());
}

#[test]
fn warnings_are_not_repeated() {
    // Each contributing calc checks its own inputs, so the same condition can be
    // raised more than once. A list that repeats itself is a list people skim.
    let result = headline();
    let mut seen = std::collections::HashSet::new();
    for warning in &result.warnings {
        assert!(
            seen.insert((warning.code, warning.field.clone())),
            "duplicate warning {:?} about {:?}",
            warning.code,
            warning.field
        );
    }
}

#[test]
fn supplying_viscosity_means_the_regime_is_actually_checked() {
    // The CLI knows the viscosity, so it must pass it on. Handing `None` to
    // darcy_weisbach would make that calc correctly report RANGE_CHECK_SKIPPED -
    // which would be misleading here, because this tool has just computed the
    // Reynolds number itself.
    let result = headline();
    assert!(
        !result
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::RangeCheckSkipped),
        "the CLI knows the viscosity; nothing should report as unchecked: {:?}",
        result.warnings
    );
}

#[test]
fn transitional_flow_is_flagged_rather_than_hidden() {
    // Air at 2 L/s in a 50 mm pipe lands around Re = 3400.
    let result = pipe::compute(
        "air",
        20.0,
        2.0,
        FlowUnit::LitresPerSecond,
        0.05,
        20.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &[],
        FrictionMethod::Colebrook,
    )
    .unwrap();

    assert!(
        result.re > 2000.0 && result.re < 4000.0,
        "expected transitional flow, got Re = {}",
        result.re
    );
    assert!(
        result
            .warnings
            .iter()
            .any(|w| w.code == WarningCode::TransitionalFlow),
        "transitional flow must be flagged: {:?}",
        result.warnings.iter().map(|w| w.code).collect::<Vec<_>>()
    );
}

#[test]
fn unknown_fittings_and_unknown_fluids_are_rejected() {
    let bad_fitting = pipe::compute(
        "water",
        20.0,
        10.0,
        FlowUnit::CubicMetresPerHour,
        0.1,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &["no_such_fitting".to_string()],
        FrictionMethod::Colebrook,
    );
    assert!(bad_fitting.is_err());

    let bad_fluid = pipe::compute(
        "unobtainium",
        20.0,
        10.0,
        FlowUnit::CubicMetresPerHour,
        0.1,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &[],
        FrictionMethod::Colebrook,
    );
    assert!(bad_fluid.is_err());
}

#[test]
fn a_temperature_outside_the_table_is_refused_not_extrapolated() {
    // Water's viscosity changes by a factor of six across 0-100 C, so a
    // straight-line extension past either end would be a confident wrong number.
    for temperature in [-10.0, 150.0] {
        let result = pipe::compute(
            "water",
            temperature,
            10.0,
            FlowUnit::CubicMetresPerHour,
            0.1,
            100.0,
            pipe::DEFAULT_ROUGHNESS_M,
            &[],
            FrictionMethod::Colebrook,
        );
        assert!(result.is_err(), "{temperature} C should have been refused");
    }
}

#[test]
fn zero_diameter_is_refused_before_dividing() {
    let result = pipe::compute(
        "water",
        20.0,
        10.0,
        FlowUnit::CubicMetresPerHour,
        0.0,
        100.0,
        pipe::DEFAULT_ROUGHNESS_M,
        &[],
        FrictionMethod::Colebrook,
    );
    assert!(result.is_err());
}

#[test]
fn flow_unit_and_method_names_parse_strictly() {
    // An unrecognised unit must not fall back to a default: reading m3/h as L/s
    // is a factor of 3600, and reading it as "whatever the default was" is
    // worse, because nothing looks wrong.
    assert!(FlowUnit::parse("m3/h").is_ok());
    assert!(FlowUnit::parse("L/s").is_ok());
    assert!(FlowUnit::parse("kg/s").is_ok());
    assert!(FlowUnit::parse("gpm").is_err());

    assert!(FrictionMethod::parse("colebrook").is_ok());
    assert!(FrictionMethod::parse("swamee-jain").is_ok());
    assert!(FrictionMethod::parse("blasius").is_err());
}
