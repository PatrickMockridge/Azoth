//! Spec-driven tests for the `characterization.tbp_fraction_component` model.
//!
//! Every expected value is a field of
//! `validation/neqsim/captures/tbp_fraction_component_probe.tsv`, read from the component's own
//! getters after `addTBPfraction`.

use azoth_characterization::{TbpFractionComponentResult, model_gen, tbp_fraction_component};
use azoth_core::spec::TestCase;
use azoth_core::units::{kilograms_per_cubic_meter, kilograms_per_mole};
use azoth_test_support as common;

const MODEL_ID: &str = "characterization.tbp_fraction_component";

fn call(case: &TestCase) -> TbpFractionComponentResult {
    tbp_fraction_component(
        kilograms_per_mole(case.input("molar_mass").expect("every case states a mass")),
        kilograms_per_cubic_meter(case.input("density").expect("every case states a density")),
    )
    .unwrap_or_else(|e| panic!("case `{}` should compute but failed: {e}", case.id))
}

#[test]
fn every_case_in_the_spec() {
    let spec = model_gen::model(MODEL_ID).expect("the model should be in its own table");
    assert!(!spec.cases.is_empty(), "the model should have cases");

    for case in spec.cases {
        let result = call(case);
        let context = format!("{}::{}", spec.id, case.id);
        let pairs: [(&str, f64); 19] = [
            ("critical_temperature", result.critical_temperature.value),
            ("critical_pressure", result.critical_pressure.value),
            ("boiling_temperature", result.boiling_temperature.value),
            ("acentric_factor", result.acentric_factor),
            ("attraction_exponent", result.attraction_exponent),
            ("critical_volume", result.critical_volume.value),
            ("racket_z", result.racket_z),
            ("racket_z_cpa", result.racket_z_cpa),
            ("parachor", result.parachor),
            ("critical_viscosity", result.critical_viscosity.value),
            (
                "triple_point_temperature",
                result.triple_point_temperature.value,
            ),
            ("heat_of_fusion", result.heat_of_fusion.value),
            (
                "ideal_gas_enthalpy_of_formation",
                result.ideal_gas_enthalpy_of_formation.value,
            ),
            ("cp_a", result.cp_a),
            ("cp_b", result.cp_b),
            ("cp_c", result.cp_c),
            ("cp_d", result.cp_d),
            // The case records the declared unit, ångström; the kernel answers SI.
            (
                "lennard_jones_diameter",
                result.lennard_jones_diameter.value * 1.0e10,
            ),
            ("lennard_jones_energy", result.lennard_jones_energy.value),
        ];
        for (name, actual) in pairs {
            common::assert_close(
                actual,
                case.expected_value(name)
                    .unwrap_or_else(|| panic!("{context}: the case states no {name}")),
                case.tolerance,
                &format!("{context}::{name}"),
            );
        }
        common::assert_consistent(&result, &context);
    }
}

#[test]
fn a_cut_without_a_positive_mass_or_density_is_refused() {
    for (mass, density) in [(0.0, 730.0), (0.095, 0.0), (-0.1, 730.0)] {
        assert!(
            tbp_fraction_component(kilograms_per_mole(mass), kilograms_per_cubic_meter(density))
                .is_err(),
            "mass {mass} and density {density} should be refused"
        );
    }
}

#[test]
fn the_acentric_factor_is_the_fitted_coefficients_root() {
    // The model's Edmister value and the root `setm` writes are different numbers; the card
    // states the root, and the two disagree by more than any tolerance would hide.
    let result =
        tbp_fraction_component(kilograms_per_mole(0.095), kilograms_per_cubic_meter(730.0))
            .expect("computes");
    assert!((result.acentric_factor - 0.464_933_101_556_860_1).abs() < 1e-12);
    assert!((result.acentric_factor - 0.458_707_752_128_320_5).abs() > 1e-3);
}
