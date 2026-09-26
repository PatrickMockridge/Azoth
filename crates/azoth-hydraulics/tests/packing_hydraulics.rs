//! `hydraulics.packing_hydraulics`, against `validation/neqsim/captures/packing_probe.tsv`.
//!
//! The capture drives `PackingHydraulicsCalculator` directly on four states, printing the packing
//! it resolved, every input the correlations read and all sixteen outputs - so the cases hold the
//! numbers and this file holds the two verdicts, the table's own resolution and the two rows a
//! case cannot state.

use azoth_core::units::{
    kilograms_per_cubic_meter, kilograms_per_second, meters, newtons_per_meter, pascal_seconds,
    pascals,
};
use azoth_hydraulics::packing::{
    critical_surface_tension, display_name, normalize, packing, packing_or_default,
};
use azoth_hydraulics::packing_hydraulics::{PackingState, packing_hydraulics};
use azoth_hydraulics::spec_gen;
use azoth_test_support as common;

const CALC_ID: &str = "hydraulics.packing_hydraulics";

fn state(case: &azoth_core::spec::TestCase) -> PackingState {
    PackingState {
        column_diameter: meters(common::input(case, "column_diameter")),
        packed_height: common::input(case, "packed_height"),
        vapor_mass_flow: kilograms_per_second(common::input(case, "vapor_mass_flow")),
        liquid_mass_flow: kilograms_per_second(common::input(case, "liquid_mass_flow")),
        vapor_density: kilograms_per_cubic_meter(common::input(case, "vapor_density")),
        liquid_density: kilograms_per_cubic_meter(common::input(case, "liquid_density")),
        vapor_viscosity: pascal_seconds(common::input(case, "vapor_viscosity")),
        liquid_viscosity: pascal_seconds(common::input(case, "liquid_viscosity")),
        surface_tension: newtons_per_meter(common::input(case, "surface_tension")),
        vapor_diffusivity: common::input(case, "vapor_diffusivity"),
        liquid_diffusivity: common::input(case, "liquid_diffusivity"),
        hydraulic_capacity_factor: common::input(case, "hydraulic_capacity_factor"),
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
        match case.kind {
            "worked_example" | "reference" => {
                let out = packing_hydraulics(common::input_str(case, "packing"), state(case))
                    .unwrap_or_else(|e| {
                        panic!("test `{}` should compute but failed: {e}", case.id)
                    });
                let context = &format!("{}::{}", spec.id, case.id);
                for (field, value, tolerance) in [
                    ("flooding_velocity", out.flooding_velocity, case.tolerance),
                    ("vapor_velocity", out.vapor_velocity, case.tolerance),
                    ("liquid_velocity", out.liquid_velocity, case.tolerance),
                    ("f_factor", out.f_factor, case.tolerance),
                    ("percent_flood", out.percent_flood, case.tolerance),
                    (
                        "pressure_drop_per_meter",
                        out.pressure_drop_per_meter.value,
                        case.tolerance,
                    ),
                    (
                        "total_pressure_drop",
                        out.total_pressure_drop.value,
                        case.tolerance,
                    ),
                    ("wetted_area", out.wetted_area, case.tolerance),
                    ("k_ga", out.k_ga, case.tolerance),
                    ("k_la", out.k_la, case.tolerance),
                    ("htu_g", out.htu_g, case.tolerance),
                    ("htu_l", out.htu_l, case.tolerance),
                    ("htu_og", out.htu_og, case.tolerance),
                    ("hetp", out.hetp.value, case.tolerance),
                    ("theoretical_stages", out.theoretical_stages, case.tolerance),
                    (
                        "specific_surface_area",
                        out.specific_surface_area,
                        case.tolerance,
                    ),
                    ("void_fraction", out.void_fraction, case.tolerance),
                    ("packing_factor", out.packing_factor, case.tolerance),
                    ("wetting_rate", out.wetting_rate, case.tolerance),
                    (
                        "minimum_wetting_rate",
                        out.minimum_wetting_rate,
                        case.tolerance,
                    ),
                ] {
                    if let Some(expected) = case.expected_value(field) {
                        common::assert_close(
                            value,
                            expected,
                            tolerance,
                            &format!("{context} ({field})"),
                        );
                    }
                }
                common::assert_consistent(&out, context);
            }
            "property" => match case.property {
                Some("unit_round_trip") => {}
                other => panic!("{}::{}: unhandled property {other:?}", spec.id, case.id),
            },
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 2,
        "expected several active cases, ran {executed}"
    );
}

/// **The file's row replaces the built-in, and the resolved geometry is a number.**
///
/// `PackingSpecificationLibrary` registers its 22 built-ins and then loads the CSV, so a row
/// whose normalized name collides replaces the earlier registration: `Pall-Ring-50` resolves to
/// the file's *plastic* row - 111.1 m**2/m**3, void 0.919, factor 180 - rather than the built-in
/// metal one (120.0, 0.96, 66). A port that loaded the file first would answer the flood from the
/// other geometry, which is a difference the capture's own first row measures.
#[test]
fn the_file_row_replaces_the_built_in() {
    let pall = packing_or_default("Pall-Ring-50");
    assert_eq!(pall.name, "Pall-Ring-50");
    assert_eq!(pall.category, "random");
    assert!((pall.specific_surface_area - 111.1).abs() < 1e-12);
    assert!((pall.void_fraction - 0.919).abs() < 1e-12);
    assert!((pall.packing_factor - 180.0).abs() < 1e-12);
    // The Billet constants the file's own columns derive: `Cp` and `Ch/6`.
    assert!((pall.billet_liquid_constant - 0.698).abs() < 1e-12);
    assert!((pall.billet_gas_constant - 2.725 / 6.0).abs() < 1e-12);
    // The critical surface tension follows the material, and this row is plastic.
    assert!((pall.critical_surface_tension - 0.033).abs() < 1e-12);
    assert!((critical_surface_tension("ceramic") - 0.061).abs() < 1e-12);
    assert!((critical_surface_tension("metal") - 0.075).abs() < 1e-12);

    // A built-in the file does not touch keeps its own values, and the alias lookup is forgiving.
    let mellapak = packing_or_default("Mellapak-250Y");
    assert_eq!(mellapak.category, "structured");
    assert!((mellapak.specific_surface_area - 250.0).abs() < 1e-12);
    assert_eq!(normalize("PALL RING 50"), normalize("Pall-Ring-50"));
    assert!(packing("pallring 50").is_some(), "the forgiving key");
    // And a name nothing carries answers the default, which is the class's own fallback.
    assert_eq!(packing_or_default("nothing-like-this").name, "Pall-Ring-50");
    assert_eq!(display_name("pallring", "random", 50.0), "Pall-Ring-50");
    assert_eq!(display_name("Flexipac", "structured", 220.0), "Flexipac");
}

/// **The two verdicts, which are the class's own thresholds and not a case's field.**
///
/// The first captured state is below the minimum wetting rate for a random packing (`5e-5`
/// m**3/(m**2 s)) and at 9.67 per cent of flood, so both verdicts are false; the third, at 12
/// kg/s of liquid, wets and is still outside the 40-to-80-per-cent design window. Both are the
/// measurement - the class's own warnings say the same - and a port that reported them true would
/// be certifying a bed the class refuses.
#[test]
fn the_verdicts_are_the_classs_thresholds() {
    let run = |liquid_kg_per_s: f64| {
        packing_hydraulics(
            "Pall-Ring-50",
            PackingState {
                column_diameter: meters(1.0),
                packed_height: 5.0,
                vapor_mass_flow: kilograms_per_second(0.35),
                liquid_mass_flow: kilograms_per_second(liquid_kg_per_s),
                vapor_density: kilograms_per_cubic_meter(45.0),
                liquid_density: kilograms_per_cubic_meter(990.0),
                vapor_viscosity: pascal_seconds(1.8e-5),
                liquid_viscosity: pascal_seconds(6.5e-4),
                surface_tension: newtons_per_meter(0.072),
                vapor_diffusivity: 3.3e-7,
                liquid_diffusivity: 1.9e-9,
                hydraulic_capacity_factor: 1.0,
            },
        )
        .expect("the state computes")
    };
    let dry = run(3.5);
    assert!(!dry.wetting_ok, "3.5 kg/s is under the minimum");
    assert!(
        !dry.design_ok,
        "and 9.67 per cent of flood is under the window"
    );
    assert!((dry.minimum_wetting_rate - 5.0e-5).abs() < 1e-15);

    let wet = run(12.0);
    assert!(wet.wetting_ok, "12 kg/s wets the bed");
    assert!(
        !wet.design_ok,
        "and 17.5 per cent of flood is still oversized"
    );
    assert_eq!(
        pascals(wet.total_pressure_drop.value).value,
        wet.total_pressure_drop.value
    );
}
