//! `thermal.boiling_heat_transfer`, against
//! `validation/neqsim/captures/boiling_heat_transfer_probe.tsv`.
//!
//! The capture drives `BoilingHeatTransfer` on the class's R134a-like test fluid, with the
//! composed state and one row per branch: no wall superheat, the quality clamps at both
//! ends, the horizontal stratified correction, the interval counts and a zero heat flux.

use azoth_test_support as common;
use azoth_thermal::boiling_heat_transfer::{BoilingRequest, rate};
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.boiling_heat_transfer";

fn request(case: &azoth_core::spec::TestCase) -> BoilingRequest {
    BoilingRequest {
        mass_flux: common::input(case, "mass_flux"),
        tube_inner_diameter: common::input(case, "tube_inner_diameter"),
        liquid_density: common::input(case, "liquid_density"),
        vapor_density: common::input(case, "vapor_density"),
        liquid_viscosity: common::input(case, "liquid_viscosity"),
        vapor_viscosity: common::input(case, "vapor_viscosity"),
        liquid_cp: common::input(case, "liquid_cp"),
        liquid_conductivity: common::input(case, "liquid_conductivity"),
        surface_tension: common::input(case, "surface_tension"),
        heat_of_vaporization: common::input(case, "heat_of_vaporization"),
        wall_superheat: common::input(case, "wall_superheat"),
        saturation_pressure_difference: common::input(case, "saturation_pressure_difference"),
        heat_flux: common::input(case, "heat_flux"),
        vapor_quality: common::input(case, "vapor_quality"),
        quality_in: common::input(case, "quality_in"),
        quality_out: common::input(case, "quality_out"),
        intervals: case.input("intervals"),
        horizontal: case.flag("horizontal"),
    }
}

#[test]
fn every_case_in_the_spec() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    common::assert_skips_are_explained(spec);

    let mut executed = 0;
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        if case.kind != "worked_example" && case.kind != "reference" {
            continue;
        }
        let out = rate(&request(case))
            .unwrap_or_else(|e| panic!("test `{}` should compute but failed: {e}", case.id));
        let context = &format!("{}::{}", spec.id, case.id);
        for (field, value, tolerance) in [
            ("chen_htc", out.chen_htc.value, case.tolerance),
            (
                "gungor_winterton_htc",
                out.gungor_winterton_htc.value,
                case.tolerance,
            ),
            (
                "gungor_winterton_corrected_htc",
                out.gungor_winterton_corrected_htc.value,
                case.tolerance,
            ),
            ("average_htc", out.average_htc.value, case.tolerance),
            (
                "martinelli_parameter",
                out.martinelli_parameter,
                case.tolerance,
            ),
            (
                "chen_enhancement_factor",
                out.chen_enhancement_factor,
                case.tolerance,
            ),
            (
                "chen_suppression_factor",
                out.chen_suppression_factor,
                case.tolerance,
            ),
        ] {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(value, expected, tolerance, &format!("{context}::{field}"));
            }
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_split_interval_count_and_zero_inputs_are_refused() {
    let base = BoilingRequest {
        mass_flux: 300.0,
        tube_inner_diameter: 0.01905,
        liquid_density: 1147.0,
        vapor_density: 50.1,
        liquid_viscosity: 1.65e-4,
        vapor_viscosity: 1.25e-5,
        liquid_cp: 1480.0,
        liquid_conductivity: 0.078,
        surface_tension: 0.006,
        heat_of_vaporization: 163000.0,
        wall_superheat: 5.0,
        saturation_pressure_difference: 5.0e4,
        heat_flux: 2.0e4,
        vapor_quality: 0.3,
        quality_in: 0.1,
        quality_out: 0.9,
        intervals: Some(7.5),
        horizontal: Some(false),
    };
    assert!(rate(&base).is_err());

    let zero_flux = BoilingRequest {
        mass_flux: 0.0,
        ..base.clone()
    };
    assert!(rate(&zero_flux).is_err());

    let negative_superheat = BoilingRequest {
        wall_superheat: -1.0,
        ..base
    };
    assert!(rate(&negative_superheat).is_err());
}
