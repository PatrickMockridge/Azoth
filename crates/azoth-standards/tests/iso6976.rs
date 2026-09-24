//! `standards.iso6976`'s own arithmetic, against the numbers its spec works out.
//!
//! The cross-implementation check is the case machinery's
//! (`python/tests/models/test_model_cases.py` runs every case through both backends); what
//! is here is the worked example the spec carries, plus the two refusals a gas can earn.

use azoth_core::units::kelvins;
use azoth_standards::iso6976;

fn names(list: &[&str]) -> Vec<String> {
    list.iter().map(|name| (*name).to_string()).collect()
}

#[test]
fn pure_methane_matches_the_worked_example() {
    let r = iso6976(
        &names(&["methane"]),
        &[1.0],
        kelvins(288.15),
        kelvins(298.15),
    )
    .expect("methane is in the standard's table");

    println!(
        "M={} Z={} d={} rho_i={} rho_r={} Hsup={} Hinf={}",
        r.molar_mass.value,
        r.compression_factor.value,
        r.relative_density.value,
        r.density_ideal.value,
        r.density_real.value,
        r.superior_calorific_value.value,
        r.inferior_calorific_value.value
    );

    // The spec's worked example, to the last bit it states.
    assert!((r.molar_mass.value - 0.016043).abs() < 1e-15);
    assert!((r.compression_factor.value - 0.99800191).abs() < 1e-15);
    assert!((r.superior_calorific_value.value - 890630.0).abs() < 1e-9);
    assert!((r.inferior_calorific_value.value - 802600.0).abs() < 1e-9);
    // The densities follow from the molar mass and the reference state.
    assert!((r.density_ideal.value - 0.6784954068722059).abs() < 1e-12);
    assert!((r.density_real.value - 0.6798538159833841).abs() < 1e-12);
    assert!((r.relative_density.value - 0.5547423719626039).abs() < 1e-12);
}

/// **A mixture is the sum of its rows, and the reference temperatures choose which ones.**
///
/// The two calorific values are stated at different reference temperatures, so the same gas
/// has more than one pair - which is why the temperature is an input rather than a constant
/// inside.
#[test]
fn a_mixture_sums_its_rows_and_the_reference_temperatures_bite() {
    let components = names(&["methane", "n-butane"]);
    let z = [0.9, 0.1];
    let at_25 = iso6976(&components, &z, kelvins(288.15), kelvins(298.15)).expect("a mixture");
    let at_0 = iso6976(&components, &z, kelvins(273.15), kelvins(273.15)).expect("a mixture");
    let at_60f = iso6976(&components, &z, kelvins(288.7), kelvins(288.7)).expect("a mixture");

    println!(
        "25 C: Hsup={} Hinf={} Z={} d={}",
        at_25.superior_calorific_value.value,
        at_25.inferior_calorific_value.value,
        at_25.compression_factor.value,
        at_25.relative_density.value
    );
    println!(
        "0 C:  Hsup={} Hinf={} Z={}",
        at_0.superior_calorific_value.value,
        at_0.inferior_calorific_value.value,
        at_0.compression_factor.value
    );
    println!(
        "60 F: Hsup={} Hinf={}",
        at_60f.superior_calorific_value.value, at_60f.inferior_calorific_value.value
    );

    // The sums are the weighted rows: the superior value is dominated by methane's, and the
    // heavier butane's is larger per mole - so the mixture's is between them.
    assert!(at_25.superior_calorific_value.value > 890630.0);
    assert!(at_25.superior_calorific_value.value < 2_877_000.0);
    // A reference temperature moves the numbers: the standard tabulates them, so they differ.
    assert!(at_25.superior_calorific_value.value != at_0.superior_calorific_value.value);
    assert!(at_0.compression_factor.value != at_25.compression_factor.value);
    assert!(at_60f.superior_calorific_value.value != at_25.superior_calorific_value.value);
}

/// **The two refusals: a gas the standard has no row for, and a reference temperature it
/// does not tabulate.**
#[test]
fn a_gas_outside_the_standard_and_a_temperature_outside_its_set_are_refused() {
    // Propylene is ISO 6976's and not the component databank's, so it is not in the compiled
    // table - the manifest records all sixteen such rows.
    let refused = iso6976(
        &names(&["propylene"]),
        &[1.0],
        kelvins(288.15),
        kelvins(298.15),
    );
    assert!(
        refused.is_err(),
        "propylene has no row in the compiled table"
    );

    // 30 C is not one of the standard's combustion temperatures, and the class would
    // silently use 25 instead.
    let refused = iso6976(
        &names(&["methane"]),
        &[1.0],
        kelvins(288.15),
        kelvins(303.15),
    );
    assert!(
        refused.is_err(),
        "30 C is not a reference temperature the table carries"
    );

    // 17 C is inside the volumetric range but is not a column of the table.
    let refused = iso6976(
        &names(&["methane"]),
        &[1.0],
        kelvins(290.15),
        kelvins(298.15),
    );
    assert!(
        refused.is_err(),
        "17 C is not a reference temperature the table carries"
    );
}
