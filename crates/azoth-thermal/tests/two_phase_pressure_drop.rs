//! `thermal.two_phase_pressure_drop`, against
//! `validation/neqsim/captures/two_phase_pressure_drop_probe.tsv`.
//!
//! The capture drives `TwoPhasePressureDrop` on the class's R134a-like test fluid, with the
//! composed state and one row per branch: both quality endpoints, a zero surface tension and
//! a zero flux, both sweep directions, and the odd and below-minimum interval counts.

use azoth_test_support as common;
use azoth_thermal::spec_gen;
use azoth_thermal::two_phase_pressure_drop::{
    TwoPhasePressureDropRequest, fanning_friction, friedel_gradient, gravitational_gradient,
    muller_steinhagen_heck_gradient, rate,
};

const CALC_ID: &str = "thermal.two_phase_pressure_drop";

fn request(case: &azoth_core::spec::TestCase) -> TwoPhasePressureDropRequest {
    TwoPhasePressureDropRequest {
        mass_flux: common::input(case, "mass_flux"),
        vapor_quality: common::input(case, "vapor_quality"),
        tube_inner_diameter: common::input(case, "tube_inner_diameter"),
        tube_length: common::input(case, "tube_length"),
        liquid_density: common::input(case, "liquid_density"),
        vapor_density: common::input(case, "vapor_density"),
        liquid_viscosity: common::input(case, "liquid_viscosity"),
        vapor_viscosity: common::input(case, "vapor_viscosity"),
        surface_tension: common::input(case, "surface_tension"),
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
        for (field, value) in [
            ("friedel_gradient", out.friedel_gradient.value),
            ("friedel_pressure_drop", out.friedel_pressure_drop.value),
            (
                "friedel_average_pressure_drop",
                out.friedel_average_pressure_drop.value,
            ),
            (
                "muller_steinhagen_heck_gradient",
                out.muller_steinhagen_heck_gradient.value,
            ),
            ("gravitational_gradient", out.gravitational_gradient.value),
            (
                "acceleration_pressure_drop",
                out.acceleration_pressure_drop.value,
            ),
        ] {
            if let Some(expected) = case.expected_value(field) {
                common::assert_close(
                    value,
                    expected,
                    case.tolerance,
                    &format!("{context}::{field}"),
                );
            }
        }
        executed += 1;
    }
    assert!(executed > 0, "no case ran");
}

#[test]
fn the_two_friction_correlations_meet_at_the_quality_endpoints() {
    let (mass_flux, tube_id, liquid_density, vapor_density) = (300.0, 0.01905, 1147.0, 50.1);
    let (liquid_viscosity, vapor_viscosity) = (1.65e-4, 1.25e-5);
    let gradient = |surface_tension, x| {
        friedel_gradient(
            mass_flux,
            x,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            vapor_viscosity,
            surface_tension,
        )
    };
    let msh = |x| {
        muller_steinhagen_heck_gradient(
            mass_flux,
            x,
            tube_id,
            liquid_density,
            vapor_density,
            liquid_viscosity,
            vapor_viscosity,
        )
    };
    common::assert_close(
        gradient(0.006, 0.0),
        msh(0.0),
        1e-12,
        "liquid-only endpoint",
    );
    common::assert_close(
        gradient(0.006, 1.0),
        msh(1.0),
        1e-12,
        "vapour-only endpoint",
    );
    common::assert_close(msh(0.0), 47.7043657235868, 1e-9, "the all-liquid gradient");
    common::assert_close(msh(1.0), 572.9811504639314, 1e-9, "the all-vapour gradient");
}

#[test]
fn the_fanning_factor_switches_at_the_laminar_boundary() {
    assert_eq!(fanning_friction(0.0), 0.0);
    assert_eq!(fanning_friction(-5.0), 0.0);
    common::assert_close(fanning_friction(1000.0), 0.016, 1e-12, "laminar");
    common::assert_close(
        fanning_friction(10_000.0),
        0.079 * 10_000.0f64.powf(-0.25),
        1e-12,
        "turbulent",
    );
}

#[test]
fn the_gravitational_endpoints_are_the_pure_densities() {
    common::assert_close(
        gravitational_gradient(0.0, 1147.0, 50.1),
        1147.0 * 9.81,
        1e-12,
        "liquid",
    );
    common::assert_close(
        gravitational_gradient(1.0, 1147.0, 50.1),
        50.1 * 9.81,
        1e-12,
        "vapour",
    );
}

#[test]
fn a_non_integral_interval_count_and_a_non_positive_length_are_refused() {
    let base = TwoPhasePressureDropRequest {
        mass_flux: 300.0,
        vapor_quality: 0.5,
        tube_inner_diameter: 0.01905,
        tube_length: 5.0,
        liquid_density: 1147.0,
        vapor_density: 50.1,
        liquid_viscosity: 1.65e-4,
        vapor_viscosity: 1.25e-5,
        surface_tension: 0.006,
        quality_in: 0.1,
        quality_out: 0.9,
        intervals: Some(7.5),
    };
    assert!(rate(&base).is_err());

    let zero_length = TwoPhasePressureDropRequest {
        tube_length: 0.0,
        intervals: None,
        ..base.clone()
    };
    assert!(rate(&zero_length).is_err());

    let negative_surface_tension = TwoPhasePressureDropRequest {
        surface_tension: -0.001,
        intervals: None,
        ..base.clone()
    };
    assert!(rate(&negative_surface_tension).is_err());

    let zero_diameter = TwoPhasePressureDropRequest {
        tube_inner_diameter: 0.0,
        intervals: None,
        ..base
    };
    assert!(rate(&zero_diameter).is_err());
}

#[test]
fn a_zero_flux_is_an_answer_not_a_refusal() {
    let request = TwoPhasePressureDropRequest {
        mass_flux: 0.0,
        vapor_quality: 0.5,
        tube_inner_diameter: 0.01905,
        tube_length: 5.0,
        liquid_density: 1147.0,
        vapor_density: 50.1,
        liquid_viscosity: 1.65e-4,
        vapor_viscosity: 1.25e-5,
        surface_tension: 0.006,
        quality_in: 0.1,
        quality_out: 0.9,
        intervals: None,
    };
    let out = rate(&request).expect("a zero flux computes zero");
    assert_eq!(out.friedel_gradient.value, 0.0);
    assert_eq!(out.friedel_pressure_drop.value, 0.0);
    assert_eq!(out.friedel_average_pressure_drop.value, 0.0);
    assert_eq!(out.muller_steinhagen_heck_gradient.value, 0.0);
    assert_eq!(out.acceleration_pressure_drop.value, 0.0);
    common::assert_close(
        out.gravitational_gradient.value,
        941.8239194720564,
        1e-9,
        "gravity",
    );
}
