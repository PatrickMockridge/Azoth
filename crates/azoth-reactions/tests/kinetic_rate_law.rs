//! The kinetic rate law, against the probe that reads it off a fluid's own reactions.
//!
//! The oracle is `validation/neqsim/KineticsProbe.java` and its capture
//! `captures/kinetics_probe.tsv`. Three reactions of a CO2/water chemical system are printed
//! there with their laws and their rate factors at two temperatures, and the same three are
//! switched to the Arrhenius branch afterwards - so both laws are pinned on the reactions a
//! fluid actually carries rather than on a construction of the test's own.
//!
//! **The legacy law ignores the reaction.** Both temperatures come back identical across all
//! three reactions, which is what a law written from literals does; and the stored
//! `rateFactor` of `0.003` is not read by it. The test asserts that sameness deliberately: it is
//! the property that makes the Arrhenius branch a migration rather than a variant.

use azoth_reactions::kinetic_rate_law::{
    KineticRateLaw, ReferenceKinetics, kinetic_rate_law, rate_factor,
};
use azoth_test_support as common;

/// From `captures/kinetics_probe.tsv`, read rather than transcribed: the capture is the oracle
/// and a number copied out of it is a number that can drift from it.
fn captured(key: &str, occurrence: usize) -> f64 {
    let capture = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../validation/neqsim/captures/kinetics_probe.tsv"),
    )
    .expect("the capture is committed");
    capture
        .lines()
        .filter_map(|line| line.trim().strip_prefix(&format!("{key}=")))
        .nth(occurrence)
        .unwrap_or_else(|| panic!("no occurrence {occurrence} of `{key}`"))
        .parse()
        .expect("a number")
}

/// **The legacy correlation is one number per temperature.** All three of the fluid's reactions
/// return the same value at `298.15 K` and the same at `373.15 K`, and the reaction's own
/// stored rate factor - `0.003` on every row of the capture - is not read.
#[test]
fn the_legacy_law_is_the_same_for_every_reaction() {
    // `legacy_rate_factor_at_298.15` appears three times - once per reaction - and
    // `legacy_rate_factor_at_373.15` three times after it.
    let cold = captured("legacy_rate_factor_at_298.15", 0);
    let hot = captured("legacy_rate_factor_at_373.15", 0);

    let stored = ReferenceKinetics {
        law: KineticRateLaw::LegacyTemperatureCorrelation,
        reference_rate: 0.003,
        activation_energy: 0.0,
        reference_temperature: 298.15,
    };
    assert!(
        (rate_factor(&stored, 298.15).expect("a rate") - cold).abs() / cold < 1.0e-12,
        "298.15 K: {} against the capture's {cold}",
        rate_factor(&stored, 298.15).expect("a rate")
    );
    assert!(
        (rate_factor(&stored, 373.15).expect("a rate") - hot).abs() / hot < 1.0e-12,
        "373.15 K: {} against the capture's {hot}",
        rate_factor(&stored, 373.15).expect("a rate")
    );

    // The parameters are not read: a different rate, energy and reference temperature give the
    // same two numbers, which is what "the same for every reaction" means.
    let other = ReferenceKinetics {
        reference_rate: 12.0,
        activation_energy: 90_000.0,
        reference_temperature: 400.0,
        ..stored
    };
    assert_eq!(
        rate_factor(&other, 298.15).expect("a rate"),
        rate_factor(&stored, 298.15).expect("a rate"),
        "the legacy branch ignores the reaction's parameters"
    );
}

/// **The Arrhenius branch, on the same reactions.** The probe sets a reference rate of `1e-3`,
/// an activation energy of `50000 J/mol` and a reference temperature of `298.15 K`, and reads
/// the factor at the same two temperatures.
#[test]
fn the_arrhenius_law_reads_its_three_parameters() {
    let cold = captured("arrhenius_rate_factor_at_298.15", 0);
    let hot = captured("arrhenius_rate_factor_at_373.15", 0);

    let kinetics = ReferenceKinetics {
        law: KineticRateLaw::ReferenceArrhenius,
        reference_rate: 1.0e-3,
        activation_energy: 50_000.0,
        reference_temperature: 298.15,
    };
    // At the reference temperature the exponential is one, so the rate factor is the rate.
    assert!(
        (rate_factor(&kinetics, 298.15).expect("a rate") - cold).abs() / cold < 1.0e-12,
        "at the reference temperature the factor is the rate: {cold}"
    );
    assert!(
        (rate_factor(&kinetics, 373.15).expect("a rate") - hot).abs() / hot < 1.0e-12,
        "373.15 K: {} against the capture's {hot}",
        rate_factor(&kinetics, 373.15).expect("a rate")
    );
}

/// **A zero reference rate short-circuits**, before the exponential, as the class does it.
#[test]
fn a_zero_rate_is_zero() {
    let kinetics = ReferenceKinetics {
        law: KineticRateLaw::ReferenceArrhenius,
        reference_rate: 0.0,
        activation_energy: 50_000.0,
        reference_temperature: 298.15,
    };
    assert_eq!(rate_factor(&kinetics, 400.0).expect("a rate"), 0.0);
}

/// The refusals, both of them: the class throws for a temperature that is not finite and
/// positive, and for reference parameters that are.
#[test]
fn the_two_refusals() {
    let kinetics = ReferenceKinetics {
        law: KineticRateLaw::ReferenceArrhenius,
        reference_rate: 1.0e-3,
        activation_energy: 50_000.0,
        reference_temperature: 298.15,
    };
    assert!(rate_factor(&kinetics, 0.0).is_err(), "zero kelvin");
    assert!(rate_factor(&kinetics, -1.0).is_err(), "below zero");
    let negative_rate = ReferenceKinetics {
        reference_rate: -1.0,
        ..kinetics
    };
    assert!(
        rate_factor(&negative_rate, 300.0).is_err(),
        "a negative rate"
    );
    // The legacy branch does not consult them, so it answers where the Arrhenius one refuses.
    let legacy = ReferenceKinetics {
        law: KineticRateLaw::LegacyTemperatureCorrelation,
        ..negative_rate
    };
    assert!(rate_factor(&legacy, 300.0).is_ok());
}

/// The selector's own fallback: an absent name is the legacy law, which is the state the
/// class's database-built reactions are in.
#[test]
fn the_selector_answers_the_legacy_law() {
    assert_eq!(
        KineticRateLaw::parse("legacy").expect("a law"),
        KineticRateLaw::LegacyTemperatureCorrelation
    );
    assert_eq!(
        KineticRateLaw::parse("reference_arrhenius").expect("a law"),
        KineticRateLaw::ReferenceArrhenius
    );
    assert!(KineticRateLaw::parse("something else").is_err());
}

// ---------------------------------------------------------------------------
// The registered model: `reactions.kinetic_rate_law`, against the cases its spec declares.
// ---------------------------------------------------------------------------

/// The model's id, so the table can be asked for it rather than the cases transcribed.
const MODEL_ID: &str = "reactions.kinetic_rate_law";

fn call_case(case: &azoth_core::spec::TestCase) -> azoth_reactions::KineticRateLawResult {
    kinetic_rate_law(
        case.string("law").expect("law"),
        common::input(case, "T"),
        common::input(case, "reference_rate"),
        common::input(case, "activation_energy"),
        common::input(case, "reference_temperature"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec =
        azoth_reactions::model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call_case(case);
        let context = &format!("{}::{}", spec.id, case.id);
        let expected = case
            .expected_value("rate_factor")
            .unwrap_or_else(|| panic!("the case `{}` declares rate_factor", case.id));
        common::assert_close(
            result.rate_factor,
            expected,
            case.tolerance,
            &format!("{context} (rate_factor)"),
        );
        common::assert_consistent(&result, context);
    }
}
