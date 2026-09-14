//! Every spec case of every unit operation, run against the **Rust** implementation alone.
//!
//! # Why this file exists
//!
//! `crates/azoth-process` was the only crate in the tree with no tests of its own. The
//! eight unit operations were exercised exclusively through Python calling the compiled
//! extension, which means the Rust side had never been run on its own - and a defect that
//! happened to be masked by how Python drives it would have been invisible.
//!
//! # One file, not eight
//!
//! The rest of this tree gives a model one test file per model. That is right when the
//! files differ. These would not: eight copies of the same loop over `spec.cases`, each
//! differing only in which function it calls and which fields it reads. The whole point
//! of the loop is that the *spec* decides what is checked, so the per-model part is a
//! dispatch and nothing more.
//!
//! # What is checked
//!
//! Every value the case declares as expected, and nothing else. A case that declares no
//! `beta` is not checked for one - the loop asks for each name in turn and skips what the
//! spec does not claim, so a spec that gains an expectation is checked the next time this
//! runs without an edit here.

mod common;

use azoth_core::spec::TestCase;
use azoth_core::units::{kelvins, pascals, watts};
use azoth_process::{
    compressor, expander, heater, mixer, model_gen, pump, separator, splitter, throttling_valve,
};
use azoth_test_support as support;

/// Compare a scalar against the case's expectation, if it declares one.
fn check(case: &TestCase, name: &str, actual: f64) {
    if let Some(expected) = case.expected_value(name) {
        support::assert_close(
            actual,
            expected,
            case.tolerance,
            &format!("{} ({name})", case.id),
        );
    }
}

/// Compare a vector against the case's expectation, entry by entry.
fn check_vector(case: &TestCase, name: &str, actual: &[f64]) {
    if let Some(expected) = case.expected_vector(name) {
        assert_eq!(
            actual.len(),
            expected.len(),
            "{}: {name} has a different length",
            case.id
        );
        for (index, (got, want)) in actual.iter().zip(expected).enumerate() {
            support::assert_close(
                *got,
                *want,
                case.tolerance,
                &format!("{} ({name}[{index}])", case.id),
            );
        }
    }
}

/// The active cases of one model, or a failure naming what is missing.
fn cases(id: &str) -> &'static [TestCase] {
    let spec = model_gen::model(id).unwrap_or_else(|| panic!("{id} is not in the model table"));
    assert!(!spec.cases.is_empty(), "{id} declares no cases");
    spec.cases
}

fn z_of(case: &TestCase) -> &'static [f64] {
    case.vector("z").expect("the case declares z")
}

#[test]
fn separator_cases() {
    for case in cases("process.separator") {
        if !case.is_active() {
            continue;
        }
        let mixture = common::mixture_from_case(case);
        let ideal_gas = common::ideal_gas_from_case(case);
        let result = separator(
            &mixture,
            &ideal_gas,
            kelvins(support::input(case, "T")),
            pascals(support::input(case, "P")),
            support::input(case, "n"),
            z_of(case),
            pascals(support::input(case, "pressure_drop")),
            watts(support::input(case, "heat_duty")),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        check(case, "T", result.temperature.value);
        check(case, "P", result.pressure.value);
        if let Some(beta) = result.beta {
            check(case, "beta", beta);
        }
        check(case, "gas_flow", result.gas_flow);
        check(case, "liquid_flow", result.liquid_flow);
        check_vector(case, "gas_z", &result.gas_z);
        check_vector(case, "liquid_z", &result.liquid_z);
    }
}

#[test]
fn mixer_cases() {
    for case in cases("process.mixer") {
        if !case.is_active() {
            continue;
        }
        let mixture = common::mixture_from_case(case);
        let ideal_gas = common::ideal_gas_from_case(case);
        let temperatures: Vec<_> = case
            .vector("T")
            .expect("the case declares T")
            .iter()
            .map(|value| kelvins(*value))
            .collect();
        let pressures: Vec<_> = case
            .vector("P")
            .expect("the case declares P")
            .iter()
            .map(|value| pascals(*value))
            .collect();
        let flows = case.vector("n").expect("the case declares n");
        let compositions = case.matrix("z").expect("the case declares z");

        let result = mixer(
            &mixture,
            &ideal_gas,
            &temperatures,
            &pressures,
            flows,
            compositions,
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        check(case, "T", result.temperature.value);
        check(case, "P", result.pressure.value);
        check(case, "flow", result.flow);
        if let Some(beta) = result.beta {
            check(case, "beta", beta);
        }
        check_vector(case, "z_out", &result.z_out);
    }
}

#[test]
fn splitter_cases() {
    for case in cases("process.splitter") {
        if !case.is_active() {
            continue;
        }
        let mixture = common::mixture_from_case(case);
        let result = splitter(
            &mixture,
            kelvins(support::input(case, "T")),
            pascals(support::input(case, "P")),
            support::input(case, "n"),
            z_of(case),
            case.vector("fractions")
                .expect("the case declares fractions"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        check(case, "T", result.temperature.value);
        check(case, "P", result.pressure.value);
        if let Some(beta) = result.beta {
            check(case, "beta", beta);
        }
        check_vector(case, "flows", &result.flows);
    }
}

#[test]
fn throttling_valve_cases() {
    for case in cases("process.throttling_valve") {
        if !case.is_active() {
            continue;
        }
        let mixture = common::mixture_from_case(case);
        let ideal_gas = common::ideal_gas_from_case(case);
        let result = throttling_valve(
            &mixture,
            &ideal_gas,
            kelvins(support::input(case, "T")),
            pascals(support::input(case, "P")),
            z_of(case),
            pascals(support::input(case, "pressure_drop")),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        check(case, "T", result.temperature.value);
        check(case, "P", result.pressure.value);
        if let Some(beta) = result.beta {
            check(case, "beta", beta);
        }
    }
}

#[test]
fn heater_cases() {
    for case in cases("process.heater") {
        if !case.is_active() {
            continue;
        }
        let mixture = common::mixture_from_case(case);
        let ideal_gas = common::ideal_gas_from_case(case);
        let result = heater(
            &mixture,
            &ideal_gas,
            kelvins(support::input(case, "T")),
            pascals(support::input(case, "P")),
            support::input(case, "n"),
            z_of(case),
            pascals(support::input(case, "pressure_drop")),
            watts(support::input(case, "heat_duty")),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        check(case, "T", result.temperature.value);
        check(case, "P", result.pressure.value);
        if let Some(beta) = result.beta {
            check(case, "beta", beta);
        }
    }
}

/// The three machines, which take the same arguments and share a procedure.
///
/// Dispatched rather than written three times for the reason the module doc gives, and
/// because their *sameness* is the point: `crates/azoth-process/src/isentropic.rs` is one
/// procedure and the difference between them is one direction.
macro_rules! machine_cases {
    ($name:ident, $id:literal, $call:path) => {
        #[test]
        fn $name() {
            for case in cases($id) {
                if !case.is_active() {
                    continue;
                }
                let mixture = common::mixture_from_case(case);
                let ideal_gas = common::ideal_gas_from_case(case);
                let result = $call(
                    &mixture,
                    &ideal_gas,
                    kelvins(support::input(case, "T")),
                    pascals(support::input(case, "P")),
                    support::input(case, "n"),
                    z_of(case),
                    pascals(support::input(case, "outlet_pressure")),
                    support::input(case, "efficiency"),
                )
                .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

                check(case, "T", result.temperature.value);
                check(case, "P", result.pressure.value);
                check(case, "power", result.power);
                check(
                    case,
                    "isentropic_temperature",
                    result.isentropic_temperature.value,
                );
                if let Some(beta) = result.beta {
                    check(case, "beta", beta);
                }
            }
        }
    };
}

machine_cases!(compressor_cases, "process.compressor", compressor);
machine_cases!(pump_cases, "process.pump", pump);
machine_cases!(expander_cases, "process.expander", expander);
