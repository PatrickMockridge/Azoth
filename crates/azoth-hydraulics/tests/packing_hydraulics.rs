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
use azoth_hydraulics::packing::{display_name, normalize, packing};
use azoth_hydraulics::packing_hydraulics::{
    PackingState, packing_hydraulics, packing_hydraulics_with,
};
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

/// **The built-ins are the reference, and the file is material-qualified.**
///
/// `PackingSpecificationLibrary` registers its 22 built-ins under their name and a `name-material`
/// alias, and a CSV row's display name carries its material, so a qualified row that already
/// exists is skipped: `Pall-Ring-50` is the built-in metal estimate (120.0 m**2/m**3, 0.96, 66)
/// and the file's metal 50 mm row (112.6, 0.951, 85) does not displace it. The file's plastic
/// 50 mm ring is reachable as `Pall-Ring-50-plastic` - the identity upstream PR #4405 (#4380,
/// #4392) preserved, where the file's order used to let the plastic row overwrite the built-in
/// geometry.
#[test]
fn the_built_in_is_the_reference_and_the_file_is_material_qualified() {
    let pall = packing("Pall-Ring-50").expect("the built-in metal estimate");
    assert_eq!(pall.name, "Pall-Ring-50");
    assert_eq!(pall.category, "random");
    assert_eq!(pall.material, "metal");
    assert!((pall.specific_surface_area - 120.0).abs() < 1e-12);
    assert!((pall.void_fraction - 0.96).abs() < 1e-12);
    assert!((pall.packing_factor - 66.0).abs() < 1e-12);
    assert!((pall.critical_surface_tension - 0.075).abs() < 1e-12);
    // The built-in's qualified alias is the same estimate.
    assert_eq!(packing("Pall-Ring-50-metal").unwrap(), pall);

    // The file's plastic 50 mm ring is the distinct entry the fix preserved.
    let plastic = packing("Pall-Ring-50-plastic").expect("the file's plastic row");
    assert_eq!(plastic.material, "plastic");
    assert!((plastic.specific_surface_area - 111.1).abs() < 1e-12);
    assert!((plastic.void_fraction - 0.919).abs() < 1e-12);
    assert!((plastic.packing_factor - 180.0).abs() < 1e-12);
    // The Billet constants the file's own columns derive: `Cp` and `Ch/6`.
    assert!((plastic.billet_liquid_constant - 0.698).abs() < 1e-12);
    assert!((plastic.billet_gas_constant - 2.725 / 6.0).abs() < 1e-12);
    assert!((plastic.critical_surface_tension - 0.033).abs() < 1e-12);

    // A file row with no built-in of its own is reachable unqualified.
    let ring35 = packing("Pall-Ring-35").expect("the file's 35 mm ring");
    assert!((ring35.specific_surface_area - 139.4).abs() < 1e-12);

    // The four placeholder rows the file carried are not packings.
    assert!(packing("Mellapak250X").is_none());
    assert!(packing("PallRingPacking").is_none());

    // A built-in the file does not touch keeps its own values, and the alias lookup is forgiving.
    let mellapak = packing("Mellapak-250Y").expect("the built-in");
    assert_eq!(mellapak.category, "structured");
    assert!((mellapak.specific_surface_area - 250.0).abs() < 1e-12);
    assert_eq!(normalize("PALL RING 50"), normalize("Pall-Ring-50"));
    assert!(packing("pallring 50").is_some(), "the forgiving key");
    // A name nothing carries is refused: upstream #4380 replaced the `Pall-Ring-50` fallback
    // with `require`, so a misspelt packing is an error rather than a silent geometry change.
    assert!(packing("nothing-like-this").is_none());
    assert_eq!(display_name("pallring", "random", 50.0), "Pall-Ring-50");
    assert_eq!(display_name("Flexipac", "structured", 220.0), "Flexipac");
}

/// **The two verdicts, which are the class's own thresholds and not a case's field.**
///
/// The first captured state is below the minimum wetting rate for a random packing (`5e-5`
/// m**3/(m**2 s)) and at 5.86 per cent of flood, so both verdicts are false; the third, at 12
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
        "and 5.86 per cent of flood is under the window"
    );
    assert!((dry.minimum_wetting_rate - 5.0e-5).abs() < 1e-15);

    let wet = run(12.0);
    assert!(wet.wetting_ok, "12 kg/s wets the bed");
    assert!(
        !wet.design_ok,
        "and 10.59 per cent of flood is still oversized"
    );
    assert_eq!(
        pascals(wet.total_pressure_drop.value).value,
        wet.total_pressure_drop.value
    );
}

/// **The category override reaches every place the category is read, and `None` is the name's own
/// row.**
///
/// `PackedColumn.setStructuredPacking` states the category on the *column* rather than selecting a
/// differently-named table row by it, so a port that could only state it by name would need one row
/// per packing per category. The override is held to the capture's own structured row
/// (`mellapak_250y_absorber`): `Mellapak-250Y` at `Some(true)` reproduces it, `Some(false)` does
/// not, and the geometry is the same row either way - what moves is the category's four reads.
///
/// **`None` being bit-identical is the load-bearing assertion.** It is what says the two functions
/// cannot state different numbers, so the four existing oracles cannot tell which one ran and the
/// new path costs them nothing.
#[test]
fn the_category_override_reaches_every_read_and_none_is_the_names_own_row() {
    let state = || PackingState {
        column_diameter: meters(1.2),
        packed_height: 8.0,
        vapor_mass_flow: kilograms_per_second(0.28),
        liquid_mass_flow: kilograms_per_second(4.0),
        vapor_density: kilograms_per_cubic_meter(30.0),
        liquid_density: kilograms_per_cubic_meter(1000.0),
        vapor_viscosity: pascal_seconds(1.7e-5),
        liquid_viscosity: pascal_seconds(5.5e-4),
        surface_tension: newtons_per_meter(0.065),
        vapor_diffusivity: 4.0e-7,
        liquid_diffusivity: 2.2e-9,
        hydraulic_capacity_factor: 1.0,
    };

    // `None` is the plain call, to the bit.
    let plain = packing_hydraulics("Mellapak-250Y", state()).expect("the state computes");
    let none = packing_hydraulics_with("Mellapak-250Y", state(), None).expect("the state computes");
    assert_eq!(plain, none, "`None` must be the name's own row exactly");

    // The name's own row is structured already, so `Some(true)` is the capture's row.
    let structured =
        packing_hydraulics_with("Mellapak-250Y", state(), Some(true)).expect("the state computes");
    assert_eq!(structured.packing_category, "structured");
    assert!((structured.wetted_area - 137.257_721_106_356_9).abs() < 1e-9);
    assert!((structured.hetp.value - 1.0).abs() < 1e-12);

    // `Some(false)` is the same geometry under the random category's coefficients.
    let random =
        packing_hydraulics_with("Mellapak-250Y", state(), Some(false)).expect("the state computes");
    assert_eq!(random.packing_category, "random");
    assert!(
        (random.specific_surface_area - structured.specific_surface_area).abs() < 1e-12,
        "the override moves the category, not the row's geometry"
    );
    // The three reads that branch, each named: the wetting minimum, the Leva exponent through the
    // pressure drop, and the HETP's rule of thumb.
    assert!(
        (random.minimum_wetting_rate - structured.minimum_wetting_rate).abs() > 1e-12,
        "the wetting minimum branches on the category"
    );
    assert!(
        (random.total_pressure_drop.value - structured.total_pressure_drop.value).abs() > 1e-12,
        "the Leva exponent branches on the category"
    );
    assert!(
        (random.hetp.value - structured.hetp.value).abs() > 1e-12,
        "the HETP's rule of thumb branches on the category"
    );
}
