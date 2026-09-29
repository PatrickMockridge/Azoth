//! `hydraulics.tray_hydraulics`, against
//! `validation/neqsim/captures/tray_hydraulics_probe.tsv`.
//!
//! The capture drives `TrayHydraulicsCalculator` directly on eight states, printing the sixteen
//! inputs it was given and all twenty-four answers - so the cases hold the numbers and this file
//! holds the four verdicts a case cannot carry.
//!
//! **And a twenty-fifth line the spec declares no output for**: `sizeColumnDiameter` mutates the
//! object it is read from, so the probe prints it last and the sizing is a *companion* of the
//! calculation rather than one of its results. This file holds its eight values.

use azoth_core::units::{
    kilograms_per_cubic_meter, kilograms_per_second, meters, millimeters, newtons_per_meter,
    pascal_seconds,
};
use azoth_hydraulics::spec_gen;
use azoth_hydraulics::tray_hydraulics::{
    TrayHydraulicsState, size_column_diameter, tray_hydraulics,
};
use azoth_test_support as common;

const CALC_ID: &str = "hydraulics.tray_hydraulics";

/// **The four verdicts per captured row**, which `TestCase` cannot carry.
///
/// `gen_registry`'s `collect_numbers` excludes flags deliberately - "so that a reader of the
/// generated table cannot tell a flag from a quantity" - and `case.flag` reads an *input*'s
/// flag, not an expectation's. So a boolean expectation has no field to arrive in, and the eight
/// rows are held here instead. The order is the capture's own, and the values are its.
const VERDICTS: &[(&str, bool, bool, bool, bool)] = &[
    // label, weeping_ok, entrainment_ok, downcommer_backup_ok, design_ok
    ("sieve_default", true, true, true, false),
    ("valve_tray", true, true, true, false),
    ("bubble_cap_tray", true, true, true, true),
    ("sieve_high_vapor", true, true, true, false),
    ("sieve_low_vapor", false, true, true, false),
    ("valve_low_vapor", false, true, true, false),
    ("sieve_stated_weir", true, true, true, false),
    ("sieve_wide", true, true, true, false),
    ("sieve_default_worked_example", true, true, true, false),
];

fn state(case: &azoth_core::spec::TestCase) -> TrayHydraulicsState {
    TrayHydraulicsState {
        tray_type: common::input_str(case, "tray_type").to_string(),
        column_diameter: meters(common::input(case, "column_diameter")),
        tray_spacing: meters(common::input(case, "tray_spacing")),
        weir_height: meters(common::input(case, "weir_height")),
        weir_length: meters(common::input(case, "weir_length")),
        downcommer_area_fraction: common::input(case, "downcommer_area_fraction"),
        // The spec declares the hole diameter in millimetres, and so does the class; the state
        // takes a length, as every other length here does.
        hole_diameter: millimeters(common::input(case, "hole_diameter")),
        hole_area_fraction: common::input(case, "hole_area_fraction"),
        design_flood_fraction: common::input(case, "design_flood_fraction"),
        vapor_mass_flow: kilograms_per_second(common::input(case, "vapor_mass_flow")),
        liquid_mass_flow: kilograms_per_second(common::input(case, "liquid_mass_flow")),
        vapor_density: kilograms_per_cubic_meter(common::input(case, "vapor_density")),
        liquid_density: kilograms_per_cubic_meter(common::input(case, "liquid_density")),
        liquid_viscosity: pascal_seconds(common::input(case, "liquid_viscosity")),
        surface_tension: newtons_per_meter(common::input(case, "surface_tension")),
        relative_volatility: common::input(case, "relative_volatility"),
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
                let out = tray_hydraulics(state(case)).unwrap_or_else(|e| {
                    panic!("test `{}` should compute but failed: {e}", case.id)
                });
                let context = &format!("{}::{}", spec.id, case.id);
                for (field, value, tolerance) in [
                    ("flooding_velocity", out.flooding_velocity, case.tolerance),
                    (
                        "actual_vapor_velocity",
                        out.actual_vapor_velocity,
                        case.tolerance,
                    ),
                    ("percent_flood", out.percent_flood, case.tolerance),
                    (
                        "minimum_vapor_velocity",
                        out.minimum_vapor_velocity,
                        case.tolerance,
                    ),
                    ("fs_factor", out.fs_factor, case.tolerance),
                    ("entrainment", out.entrainment, case.tolerance),
                    (
                        "downcommer_backup",
                        out.downcommer_backup.value,
                        case.tolerance,
                    ),
                    (
                        "downcommer_backup_fraction",
                        out.downcommer_backup_fraction,
                        case.tolerance,
                    ),
                    (
                        "total_tray_pressure_drop",
                        out.total_tray_pressure_drop.value,
                        case.tolerance,
                    ),
                    (
                        "total_tray_pressure_drop_mbar",
                        out.total_tray_pressure_drop_mbar,
                        case.tolerance,
                    ),
                    (
                        "dry_tray_pressure_drop",
                        out.dry_tray_pressure_drop.value,
                        case.tolerance,
                    ),
                    (
                        "liquid_head_pressure_drop",
                        out.liquid_head_pressure_drop.value,
                        case.tolerance,
                    ),
                    (
                        "residual_head_pressure_drop",
                        out.residual_head_pressure_drop.value,
                        case.tolerance,
                    ),
                    ("tray_efficiency", out.tray_efficiency, case.tolerance),
                    ("turndown_ratio", out.turndown_ratio, case.tolerance),
                    (
                        "calculated_weir_length",
                        out.calculated_weir_length.value,
                        case.tolerance,
                    ),
                    ("active_area", out.active_area.value, case.tolerance),
                    ("total_area", out.total_area.value, case.tolerance),
                    ("hole_area", out.hole_area.value, case.tolerance),
                    ("downcommer_area", out.downcommer_area.value, case.tolerance),
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

                // **The four verdicts, from the table above rather than from the case.**
                let (_, weeping, entrainment_ok, backup_ok, design) = VERDICTS
                    .iter()
                    .find(|(label, ..)| *label == case.id)
                    .unwrap_or_else(|| panic!("no captured verdicts for `{}`", case.id));
                assert_eq!(out.weeping_ok, *weeping, "{context}: weeping_ok");
                assert_eq!(
                    out.entrainment_ok, *entrainment_ok,
                    "{context}: entrainment_ok"
                );
                assert_eq!(
                    out.downcommer_backup_ok, *backup_ok,
                    "{context}: downcommer_backup_ok"
                );
                assert_eq!(out.design_ok, *design, "{context}: design_ok");
            }
            other => panic!("{}::{}: unknown test kind {other:?}", spec.id, case.id),
        }
        executed += 1;
    }
    assert!(
        executed >= 8,
        "expected the eight captured rows, ran {executed}"
    );
}

/// **The tray type's four branches, on one state**, which is what makes the type's fallback a
/// measurement rather than a reading: `sieve` and `valve` differ in the hole area, the flooding
/// factor and the orifice coefficient at once, and `bubble-cap` takes the `else` of the first two
/// while being the only type the weeping check returns for.
/// **`sizeColumnDiameter` on the capture's eight states**, in its own order.
///
/// The sizing is a *companion* of the calculation rather than one of its outputs: it mutates the
/// object it is read from - writing the trial `1.0` m and leaving the sized value behind - so the
/// probe prints it last, after every other line, and the spec declares no output for it.
///
/// **Every value is a standard table entry and no two rows agree by accident.** The eight cover
/// the table's `0.5`, `0.8`, `0.9`, `1.1` and `1.4`; the two `0.5` rows are the two low-vapour
/// ones, whose flooding velocity is `0.089` times the default's, and the two `0.8` rows include
/// the stated-weir row - because a stated weir length moves nothing, which the capture says in
/// the same breath.
const SIZED_DIAMETERS_M: &[(&str, f64)] = &[
    ("sieve_default", 0.8),
    ("valve_tray", 0.8),
    ("bubble_cap_tray", 0.9),
    ("sieve_high_vapor", 1.1),
    ("sieve_low_vapor", 0.5),
    ("valve_low_vapor", 0.5),
    ("sieve_stated_weir", 0.8),
    ("sieve_wide", 1.4),
    // The worked example is the first captured row, so it sizes with it.
    ("sieve_default_worked_example", 0.8),
];

#[test]
fn the_sized_diameter_is_the_standard_table_over_the_captures_rows() {
    let spec = common::spec(spec_gen::specs(), CALC_ID);
    for case in spec.all_tests() {
        if !case.is_active() {
            continue;
        }
        let Some((_, expected)) = SIZED_DIAMETERS_M
            .iter()
            .find(|(label, _)| *label == case.id)
        else {
            panic!("no captured sizing for `{}`", case.id);
        };
        let sized = size_column_diameter(&state(case))
            .unwrap_or_else(|e| panic!("`{}` should size but failed: {e}", case.id));
        assert_eq!(
            sized, *expected,
            "{}::{}: the sized diameter is {sized} against the capture's {expected}",
            spec.id, case.id
        );
    }
}

#[test]
fn the_tray_type_moves_four_things_at_one_state() {
    let run = |tray_type: &str| {
        tray_hydraulics(TrayHydraulicsState {
            tray_type: tray_type.to_string(),
            column_diameter: meters(1.0),
            tray_spacing: meters(0.6),
            weir_height: meters(0.05),
            weir_length: meters(-1.0),
            downcommer_area_fraction: 0.1,
            hole_diameter: millimeters(12.7),
            hole_area_fraction: 0.1,
            design_flood_fraction: 0.8,
            vapor_mass_flow: kilograms_per_second(1.0),
            liquid_mass_flow: kilograms_per_second(3.0),
            vapor_density: kilograms_per_cubic_meter(2.0),
            liquid_density: kilograms_per_cubic_meter(800.0),
            liquid_viscosity: pascal_seconds(1.0e-3),
            surface_tension: newtons_per_meter(0.02),
            relative_volatility: 2.0,
        })
        .expect("the state computes")
    };
    let sieve = run("sieve");
    let valve = run("valve");
    let bubble_cap = run("bubble-cap");

    // The hole area: the stated fraction, `0.14` and `0.12` - so `valve` is the largest and
    // `sieve` the smallest, its stated `0.1` being the least of the three.
    assert!(valve.hole_area.value > bubble_cap.hole_area.value);
    assert!(bubble_cap.hole_area.value > sieve.hole_area.value);
    // The flooding velocity's factor: `1.0`, `1.05` and `0.85`.
    assert!(valve.flooding_velocity > sieve.flooding_velocity);
    assert!(bubble_cap.flooding_velocity < sieve.flooding_velocity);
    // And the weeping check, where only the bubble-cap returns early.
    assert_eq!(bubble_cap.minimum_vapor_velocity, 0.0);
    assert_eq!(bubble_cap.turndown_ratio, 0.0);
    assert!(bubble_cap.weeping_ok);
}
