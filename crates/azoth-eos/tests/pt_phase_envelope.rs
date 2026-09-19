//! Spec-driven tests for the `eos.pt_phase_envelope` model.

use azoth_core::units::pascals;
use azoth_core::{AzothError, CalcResult};
use azoth_eos::mixture::Mixture;
use azoth_eos::{Cubic, databank, model_gen, pt_phase_envelope};
use azoth_test_support as common;

const MODEL_ID: &str = "eos.pt_phase_envelope";

fn methane_butane() -> Mixture {
    databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None)
        .expect("the pair resolves")
        .0
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    for case in spec.cases {
        let names = case
            .list("components")
            .expect("the case declares components");
        let mixture = databank::mixture_of(names, Cubic::Pr, None)
            .expect("the case's components resolve")
            .0;
        let result = pt_phase_envelope(
            &mixture,
            pascals(common::input(case, "P")),
            case.vector("z").expect("z"),
        )
        .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id));

        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.critical_temperature.value,
            common::expected(case, "critical_temperature"),
            case.tolerance,
            &format!("{context} (critical_temperature)"),
        );
        common::assert_close(
            result.critical_pressure.value,
            common::expected(case, "critical_pressure"),
            case.tolerance,
            &format!("{context} (critical_pressure)"),
        );
        for name in [
            "cricondenbar_temperature",
            "cricondenbar_pressure",
            "cricondentherm_temperature",
            "cricondentherm_pressure",
        ] {
            let actual = match name {
                "cricondenbar_temperature" => result.cricondenbar_temperature.value,
                "cricondenbar_pressure" => result.cricondenbar_pressure.value,
                "cricondentherm_temperature" => result.cricondentherm_temperature.value,
                _ => result.cricondentherm_pressure.value,
            };
            common::assert_close(
                actual,
                common::expected(case, name),
                case.tolerance,
                &format!("{context} ({name})"),
            );
        }
        assert_eq!(
            result.iterations as f64,
            common::expected(case, "iterations"),
            "{context}: iteration count"
        );
    }
}

/// Both branches are non-empty, ordered, and climb in pressure from the starting point.
#[test]
fn both_branches_trace_from_the_low_pressure() {
    let result = pt_phase_envelope(&methane_butane(), pascals(100000.0), &[0.5, 0.5]).unwrap();
    assert!(
        result.bubble_temperature.len() > 10,
        "bubble branch too short"
    );
    assert!(result.dew_temperature.len() > 10, "dew branch too short");
    assert!(
        (result.bubble_pressure[0] - 100000.0).abs() < 1000.0,
        "bubble branch should start near the starting pressure"
    );
    assert!(
        (result.dew_pressure[0] - 100000.0).abs() < 1000.0,
        "dew branch should start near the starting pressure"
    );
    // Each branch climbs overall: it starts near the low pressure and ends near the
    // critical point, which is above it in both temperature and pressure.
    for (t, p) in [
        (&result.bubble_temperature, &result.bubble_pressure),
        (&result.dew_temperature, &result.dew_pressure),
    ] {
        assert!(
            *t.last().unwrap() > t[0],
            "a branch should climb in temperature"
        );
        assert!(
            *p.last().unwrap() > p[0],
            "a branch should climb in pressure"
        );
    }
    assert!(
        result.is_clean(),
        "unexpected warnings: {:?}",
        result.warnings
    );
}

#[test]
fn a_single_component_is_refused() {
    let propane = databank::mixture_of(&["propane"], Cubic::Pr, None)
        .expect("propane resolves")
        .0;
    let err = pt_phase_envelope(&propane, pascals(1.0e5), &[1.0]).unwrap_err();
    assert!(matches!(err, AzothError::InvalidInput { .. }));
    assert!(err.to_string().contains("pure_saturation"));
}

#[test]
fn a_malformed_composition_is_refused() {
    let mixture = methane_butane();
    for z in [vec![0.5, 0.6], vec![0.6, -0.1], vec![0.5, 0.4, 0.1]] {
        assert!(
            pt_phase_envelope(&mixture, pascals(1.0e5), &z).is_err(),
            "z = {z:?}"
        );
    }
}

#[test]
fn a_non_positive_pressure_is_refused() {
    let mixture = methane_butane();
    for p in [0.0, -1.0] {
        assert!(
            pt_phase_envelope(&mixture, pascals(p), &[0.5, 0.5]).is_err(),
            "P={p}"
        );
    }
}
