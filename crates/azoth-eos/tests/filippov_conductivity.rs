//! Spec-driven tests for the `eos.filippov_conductivity` model.

use azoth_core::units::{kelvins, kilograms_per_mole, pascals};
use azoth_eos::filippov_conductivity::filippov_conductivity;
use azoth_test_support as common;

const MODEL_ID: &str = "eos.filippov_conductivity";

fn call(case: &azoth_core::spec::TestCase) -> azoth_eos::FilippovConductivityResult {
    let flat = case
        .matrix("liquid_conductivity")
        .expect("the case states the coefficient matrix");
    assert_eq!(flat.len() % 3, 0, "three coefficients per component");
    let rows: Vec<[f64; 3]> = flat.chunks(3).map(|row| [row[0], row[1], row[2]]).collect();
    let molar_mass: Vec<_> = case
        .vector("molar_mass")
        .expect("the case states the molar masses")
        .iter()
        .map(|value| kilograms_per_mole(*value))
        .collect();
    let tc: Vec<_> = case
        .vector("tc")
        .expect("the case states the critical temperatures")
        .iter()
        .map(|value| kelvins(*value))
        .collect();
    let pc: Vec<_> = case
        .vector("pc")
        .expect("the case states the critical pressures")
        .iter()
        .map(|value| pascals(*value))
        .collect();
    filippov_conductivity(
        &rows,
        &molar_mass,
        case.vector("z")
            .expect("the case states the mole fractions"),
        &tc,
        &pc,
        kelvins(case.input("T").expect("the case states T")),
        pascals(case.input("P").expect("the case states P")),
        case.input("filippov_coefficient")
            .expect("the case states the coefficient"),
        case.flag("use_pressure_correction")
            .expect("the case states the correction flag"),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = azoth_eos::model_gen::model(MODEL_ID).expect("the model should be in its table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = &format!("{}::{}", spec.id, case.id);
        common::assert_close(
            result.k.value,
            case.expected_value("k").expect("the case states k"),
            case.tolerance,
            &format!("{context} (k)"),
        );
        common::assert_close(
            result.pressure_correction,
            case.expected_value("pressure_correction")
                .expect("the case states pressure_correction"),
            case.tolerance,
            &format!("{context} (pressure_correction)"),
        );
        assert_eq!(
            result.pure_component_conductivity.len(),
            case.vector("z").expect("the case states z").len(),
            "{context}: one pure conductivity per component"
        );
        if let Some(expected) = case.expected_vector("pure_component_conductivity") {
            assert_eq!(
                result.pure_component_conductivity.len(),
                expected.len(),
                "{context}: the expected vector's length"
            );
            for (index, (value, want)) in result
                .pure_component_conductivity
                .iter()
                .zip(expected)
                .enumerate()
            {
                common::assert_close(
                    *value,
                    *want,
                    case.tolerance,
                    &format!("{context} (pure_component_conductivity[{index}])"),
                );
            }
        }
        common::assert_consistent(&result, context);
    }
}

/// **Methane's polynomial is negative at 300 K**, so the oil case takes the `1e-10` floor;
/// a port that skipped it would weight a negative conductivity into the mean.
#[test]
fn a_negative_polynomial_is_floored() {
    let result = filippov_conductivity(
        &[
            [0.290304, -4.720407e-4, -4.320339e-6],
            [0.24688, -4.719916e-4, -6.249512e-8],
        ],
        &[kilograms_per_mole(0.016043), kilograms_per_mole(0.058123)],
        &[0.09542082642782022, 0.9045791735721797],
        &[kelvins(190.56), kelvins(425.12)],
        &[pascals(4.599e6), pascals(3.796e6)],
        kelvins(300.0),
        pascals(2.0e6),
        0.72,
        true,
    )
    .expect("the oil state runs");
    assert!(
        (result.pure_component_conductivity[0] - 1.0e-10).abs() < 1e-20,
        "methane's polynomial is floored, not negative"
    );
    assert!(
        (result.k.value - 0.09486574751651838).abs() < 1e-15,
        "k = {} against the capture's 0.09486574751651838",
        result.k.value
    );
}

/// The Missenard correction is above one at 1000 bar and can be switched off; the capture
/// records both answers.
#[test]
fn the_pressure_correction_bites_and_can_be_switched_off() {
    let run = |correction: bool| {
        filippov_conductivity(
            &[[-0.384, 0.00525, -6.37e-6]],
            &[kilograms_per_mole(0.018015)],
            &[1.0],
            &[kelvins(647.3)],
            &[pascals(22.089e6)],
            kelvins(300.0),
            pascals(1.0e8),
            0.72,
            correction,
        )
        .expect("the water state runs")
    };
    let on = run(true);
    let off = run(false);
    assert!((on.k.value - 0.6203856060798368).abs() < 1e-15);
    assert!((on.pressure_correction - 1.004347751464848).abs() < 1e-15);
    assert!((off.k.value - 0.6177000000000002).abs() < 1e-15);
    assert_eq!(off.pressure_correction, 1.0);
}
