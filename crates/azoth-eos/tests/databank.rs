//! The component databank: names in, a mixture out.

use azoth_eos::databank;

/// The numbers NeqSim's `COMP.csv` carries for methane, checked against the file
/// rather than against this crate.
///
/// Retyped here deliberately: a test that read the same table it is testing would
/// agree with any table at all. These are the values a reader can find in the vendored
/// source, and the point is that a calculation no longer has to be given them.
const METHANE_TC: f64 = 190.56;
const METHANE_PC: f64 = 4_599_000.0;
const METHANE_OMEGA: f64 = 0.0115;

#[test]
fn a_name_resolves_to_the_constants_neqsim_ships() {
    let methane = databank::entry("methane").expect("methane is in the databank");
    assert!((methane.tc - METHANE_TC).abs() < 1e-6, "Tc: {}", methane.tc);
    assert!((methane.pc - METHANE_PC).abs() < 1e-3, "Pc: {}", methane.pc);
    assert!(
        (methane.omega - METHANE_OMEGA).abs() < 1e-9,
        "omega: {}",
        methane.omega
    );

    // The five `Cp` coefficients, which are what a spec used to retype as four
    // dimensionless numbers against a scale of its own.
    assert!(
        (methane.cp[0] - 37.978_352).abs() < 1e-6,
        "cp_a: {}",
        methane.cp[0]
    );
    assert!(
        (methane.cp[1] - -0.074_618_15).abs() < 1e-9,
        "cp_b: {}",
        methane.cp[1]
    );
    assert!(methane.cp[4] > 0.0, "cp_e: {}", methane.cp[4]);
}

#[test]
fn a_name_is_matched_without_regard_to_case_or_surrounding_space() {
    let lower = databank::entry("methane").expect("methane");
    let padded = databank::entry("  METHANE  ").expect("METHANE");
    assert_eq!(lower.name, padded.name);
}

#[test]
fn a_name_the_databank_does_not_have_is_refused() {
    let error = databank::entry("unobtainium").unwrap_err();
    let message = error.to_string();
    assert!(
        message.contains("unobtainium"),
        "the failure should name the substance it could not find: {message}"
    );
}

#[test]
fn a_mixture_comes_back_with_its_ideal_gas_model() {
    let (mixture, ideal_gas) = databank::mixture_of(&["methane", "n-butane"]).expect("a mixture");
    assert_eq!(mixture.len(), 2);
    // Both vectors are one entry per component, which is the property the parallel
    // vectors in the specs could only assert.
    assert_eq!(ideal_gas.cp_a.len(), 2);
    assert_eq!(ideal_gas.cp_e.len(), 2);
    // The pair's interaction parameter is in the databank, not the ideal-mixture zero.
    assert!(
        mixture.kij(0, 1) != 0.0,
        "methane/n-butane should carry a fitted interaction parameter"
    );
}

#[test]
fn a_pair_the_databank_does_not_carry_falls_back_to_the_ideal_mixture() {
    // NeqSim's own reader substitutes zero for an absent pair, and so does this: an
    // absent parameter is the ideal-mixture default rather than a failure.
    assert_eq!(databank::kij("methane", "krypton"), 0.0);
}

#[test]
fn an_empty_mixture_is_refused() {
    assert!(databank::mixture_of(&[]).is_err());
}

#[test]
fn the_databank_carries_the_components_the_calculations_use() {
    let names = databank::names();
    assert!(
        names.len() > 100,
        "expected a full table, found {}",
        names.len()
    );
    for name in [
        "methane", "ethane", "propane", "n-butane", "co2", "water", "nitrogen",
    ] {
        assert!(
            names.contains(&name.to_string()),
            "the databank has no `{name}`"
        );
    }
}
