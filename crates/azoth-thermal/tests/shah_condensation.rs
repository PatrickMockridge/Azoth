//! `thermal.shah_condensation`, against `validation/neqsim/captures/shah_condensation_probe.tsv`.
//!
//! The capture drives `ShahCondensation` on the class's R134a-like test fluid, with the
//! composed state and one row per branch: quality one, quality below zero, both
//! reduced-pressure clamps, an odd interval count, the gravity branch, a liquid lighter
//! than the vapour, and the three validity bands.

use azoth_test_support as common;
use azoth_thermal::shah_condensation::{ShahRequest, rate};
use azoth_thermal::spec_gen;

const CALC_ID: &str = "thermal.shah_condensation";

fn request(case: &azoth_core::spec::TestCase) -> ShahRequest {
    ShahRequest {
        mass_flux: common::input(case, "mass_flux"),
        tube_inner_diameter: common::input(case, "tube_inner_diameter"),
        liquid_density: common::input(case, "liquid_density"),
        liquid_viscosity: common::input(case, "liquid_viscosity"),
        liquid_cp: common::input(case, "liquid_cp"),
        liquid_conductivity: common::input(case, "liquid_conductivity"),
        vapor_density: common::input(case, "vapor_density"),
        reduced_pressure: common::input(case, "reduced_pressure"),
        vapor_quality: common::input(case, "vapor_quality"),
        quality_in: common::input(case, "quality_in"),
        quality_out: common::input(case, "quality_out"),
        intervals: case.input("intervals"),
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
            ("liquid_only_htc", out.liquid_only_htc.value, case.tolerance),
            (
                "liquid_only_reynolds",
                out.liquid_only_reynolds,
                case.tolerance,
            ),
            ("local_htc", out.local_htc.value, case.tolerance),
            ("average_htc", out.average_htc.value, case.tolerance),
            ("vertical_htc", out.vertical_htc.value, case.tolerance),
        ] {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(value, expected, tolerance, &format!("{context}::{field}"));
            }
        }
        if let Some(expected) = case.flag("in_valid_range") {
            assert_eq!(out.in_valid_range, expected, "{context}::in_valid_range");
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_split_interval_count_and_zero_inputs_are_refused() {
    let base = ShahRequest {
        mass_flux: 300.0,
        tube_inner_diameter: 0.01905,
        liquid_density: 1147.0,
        liquid_viscosity: 1.65e-4,
        liquid_cp: 1480.0,
        liquid_conductivity: 0.078,
        vapor_density: 50.1,
        reduced_pressure: 0.1,
        vapor_quality: 0.5,
        quality_in: 0.9,
        quality_out: 0.1,
        intervals: Some(7.5),
    };
    assert!(rate(&base).is_err());

    let zero_flux = ShahRequest {
        mass_flux: 0.0,
        ..base.clone()
    };
    assert!(rate(&zero_flux).is_err());

    let zero_pressure = ShahRequest {
        reduced_pressure: 0.0,
        ..base
    };
    assert!(rate(&zero_pressure).is_err());
}
