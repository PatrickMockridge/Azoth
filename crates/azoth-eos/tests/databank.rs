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
    let methane = databank::entry("methane", None).expect("methane is in the databank");
    assert!((methane.tc - METHANE_TC).abs() < 1e-6, "Tc: {}", methane.tc);
    assert!((methane.pc - METHANE_PC).abs() < 1e-3, "Pc: {}", methane.pc);
    assert!(
        (methane.omega - METHANE_OMEGA).abs() < 1e-9,
        "omega: {}",
        methane.omega
    );

    // The five `Cp` coefficients, which are what a spec used to retype as four
    // dimensionless numbers against a scale of its own. `Some` for everything the
    // table carries: it holds the polynomial for every row it has.
    let cp = methane
        .cp
        .expect("the table carries a heat-capacity polynomial for every row");
    assert!((cp[0] - 37.978_352).abs() < 1e-6, "cp_a: {}", cp[0]);
    assert!((cp[1] - -0.074_618_15).abs() < 1e-9, "cp_b: {}", cp[1]);
    assert!(cp[4] > 0.0, "cp_e: {}", cp[4]);
}

#[test]
fn the_transport_fields_carry_what_neqsim_ships() {
    // Molar mass and critical volume feed the density a pipe kernel assembles; the
    // dipole feeds the gas viscosity. Water is the non-zero-dipole case, so a default
    // of zero cannot pass.
    let methane = databank::entry("methane", None).expect("methane");
    let water = databank::entry("water", None).expect("water");
    let molar_mass = methane.molar_mass.expect("the table carries molar mass");
    assert!(
        (molar_mass - 0.016_043).abs() < 1e-6,
        "molar mass: {molar_mass}"
    );
    assert!(
        (methane.critical_volume.expect("critical volume") - 9.9e-5).abs() < 1e-9,
        "critical volume: {:?}",
        methane.critical_volume
    );
    assert!(
        (water.dipole.expect("dipole") - 1.8).abs() < 1e-9,
        "dipole: {:?}",
        water.dipole
    );
}

#[test]
fn the_liquid_transport_fields_carry_what_neqsim_ships() {
    // The liquid viscosity and conductivity correlations read the model selector and
    // its coefficients. n-butane is model 2, which a selector defaulted to the common
    // model 3 would silently get wrong.
    let methane = databank::entry("methane", None).expect("methane");
    let butane = databank::entry("n-butane", None).expect("n-butane");
    assert_eq!(methane.liqviscmodel, Some(3));
    assert_eq!(butane.liqviscmodel, Some(2));
    let v = methane.liqvisc.expect("the table carries liquid viscosity");
    assert!((v[0] - -26.87).abs() < 1e-9, "liqvisc1: {}", v[0]);
    assert!((v[1] - 1150.0).abs() < 1e-9, "liqvisc2: {}", v[1]);
    let c = methane
        .liquid_conductivity
        .expect("the table carries liquid conductivity");
    assert!(
        (c[0] - 0.290_304).abs() < 1e-9,
        "liquidconductivity1: {}",
        c[0]
    );
}

#[test]
fn a_name_is_matched_without_regard_to_case_or_surrounding_space() {
    let lower = databank::entry("methane", None).expect("methane");
    let padded = databank::entry("  METHANE  ", None).expect("METHANE");
    assert_eq!(lower.name, padded.name);
}

#[test]
fn a_name_the_databank_does_not_have_is_refused() {
    let error = databank::entry("unobtainium", None).unwrap_err();
    let message = error.to_string();
    assert!(
        message.contains("unobtainium"),
        "the failure should name the substance it could not find: {message}"
    );
}

#[test]
fn a_mixture_comes_back_with_its_ideal_gas_model() {
    let (mixture, ideal_gas) =
        databank::mixture_of(&["methane", "n-butane"], None).expect("a mixture");
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
    assert_eq!(databank::kij("methane", "krypton", None), 0.0);
}

#[test]
fn the_nrtl_matrices_resolve_from_names() {
    // NRTLALPHA/NRTLGIJ/NRTLGJI from INTER.csv for methanol/water, flattened row-major.
    let dij = databank::nrtl_dij(&["methanol", "water"]);
    let alpha = databank::nrtl_alpha(&["methanol", "water"]);
    assert_eq!(dij, vec![0.0, -48.68, 610.6, 0.0]);
    assert_eq!(alpha, vec![0.0, 0.303, 0.303, 0.0]);
}

#[test]
fn the_nrtl_energy_is_directional() {
    // `g_ij != g_ji`: reversing the name order swaps the off-diagonal energy, while the
    // symmetric `alpha` is unchanged.
    assert_eq!(
        databank::nrtl_dij(&["water", "methanol"]),
        vec![0.0, 610.6, -48.68, 0.0]
    );
    assert_eq!(
        databank::nrtl_alpha(&["water", "methanol"]),
        vec![0.0, 0.303, 0.303, 0.0]
    );
}

#[test]
fn the_unifac_parameters_resolve_from_names() {
    // The resolved inputs for methanol/water: two subgroups, one per component.
    let p = databank::unifac_parameters(&["methanol", "water"]).expect("methanol and water");
    assert_eq!(p.groups, vec![1.0, 0.0, 0.0, 1.0]);
    assert_eq!(p.group_r, vec![1.4311, 0.92]);
    assert_eq!(p.group_q, vec![1.432, 1.4]);
    assert_eq!(p.aij, vec![0.0, -181.0, 289.6, 0.0]);
}

#[test]
fn a_name_without_unifac_groups_is_refused() {
    // The databank ships UNIFAC groups for a subset of its substances; a name without
    // one is a caller error, not a silent ideal-mixture default.
    assert!(databank::unifac_parameters(&["methane", "unobtainium"]).is_err());
}

#[test]
fn an_empty_mixture_is_refused() {
    assert!(databank::mixture_of(&[], None).is_err());
}

#[test]
fn the_databank_carries_the_components_the_calculations_use() {
    let names = databank::names(None);
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

// ---------------------------------------------------------------------------
// The overlay: a keycard, as a value
// ---------------------------------------------------------------------------
//
// One test per rule, because the rules are what a keycard is. Every one of them is
// sabotage-verified: the comment on each names the break that makes it fail, and each
// break is a plausible way to write the merge wrong rather than a contrived one.

/// An overlay naming only `omega`, which is the per-parameter rule's test case.
fn omega_only() -> databank::Overlay {
    let mut overlay = databank::Overlay::new();
    overlay.set_component(
        "methane",
        databank::ComponentOverride {
            omega: Some(0.5),
            ..Default::default()
        },
    );
    overlay
}

#[test]
fn a_card_overrides_a_shipped_name() {
    // Sabotage: ignore the overlay in `entry`'s `(Some, Some)` arm.
    let base = databank::entry("methane", None).unwrap();
    let carded = databank::entry("methane", Some(&omega_only())).unwrap();
    assert!(
        (carded.omega - 0.5).abs() < 1e-12,
        "omega: {}",
        carded.omega
    );
    assert!(
        (carded.omega - base.omega).abs() > 1e-3,
        "the shipped omega and the card's must differ, or this proves nothing"
    );
}

#[test]
fn an_override_keeps_the_parameters_it_does_not_name() {
    // Parameter by parameter. Sabotage: `over.tc.unwrap_or(0.0)` instead of
    // `unwrap_or(base.tc)` - which is exactly what whole-record replacement does, and
    // it would make a user correcting one value restate the others.
    let base = databank::entry("methane", None).unwrap();
    let carded = databank::entry("methane", Some(&omega_only())).unwrap();
    assert_eq!(
        carded.tc, base.tc,
        "Tc was not named, so it is the shipped one"
    );
    assert_eq!(
        carded.pc, base.pc,
        "Pc was not named, so it is the shipped one"
    );
    assert_eq!(
        carded.cp, base.cp,
        "the polynomial is not a card's to state"
    );
}

#[test]
fn a_zero_kij_from_a_card_is_not_absent() {
    // Overriding a fitted pair back to ideal mixing is a caller stating something, and
    // a lookup that reads a zero as "no opinion" undoes it with no symptom.
    // Sabotage: `.filter(|v| *v != 0.0)` anywhere in the lookup.
    let fitted = databank::kij("methane", "n-butane", None);
    assert!(
        fitted.abs() > 1e-6,
        "the fixture needs a fitted non-zero pair, got {fitted}"
    );

    let mut overlay = databank::Overlay::new();
    overlay.set_kij("methane", "n-butane", 0.0).unwrap();
    assert_eq!(
        databank::kij("methane", "n-butane", Some(&overlay)),
        0.0,
        "a card's explicit zero did not beat the fitted value"
    );
}

#[test]
fn a_kij_override_is_paired_whichever_order_it_is_read_in() {
    // Sabotage: store only the ordering the caller wrote.
    let mut overlay = databank::Overlay::new();
    overlay.set_kij("n-butane", "methane", 0.5).unwrap();
    assert_eq!(databank::kij("methane", "n-butane", Some(&overlay)), 0.5);
    assert_eq!(databank::kij("n-butane", "methane", Some(&overlay)), 0.5);
}

#[test]
fn a_component_paired_with_itself_is_refused() {
    // Sabotage: accept it. `Mixture::new` catches a non-zero diagonal only for a pair
    // whose *both* names are in the mixture being built, so a self-pair on a name
    // nothing builds would be stored and read by nothing.
    let mut overlay = databank::Overlay::new();
    let refused = overlay.set_kij("methane", "methane", 0.1);
    assert!(refused.is_err(), "a self-pair was accepted");
}

#[test]
fn a_card_added_substance_has_no_polynomial() {
    // Sabotage: default `cp` to `[0.0; 5]`, which is a zero heat capacity wearing the
    // shape of a polynomial.
    let mut overlay = databank::Overlay::new();
    overlay.set_component(
        "unobtainium",
        databank::ComponentOverride {
            tc: Some(500.0),
            pc: Some(2.0e6),
            omega: Some(0.3),
            ..Default::default()
        },
    );
    let added = databank::entry("unobtainium", Some(&overlay)).expect("the card adds it");
    assert!(added.cp.is_none(), "an added substance has no polynomial");
    assert_eq!(added.name, "unobtainium");
}

#[test]
fn a_partial_new_component_is_refused() {
    // Sabotage: complete it from the nearest shipped substance, which invents data.
    let mut overlay = databank::Overlay::new();
    overlay.set_component(
        "unobtainium",
        databank::ComponentOverride {
            tc: Some(500.0),
            ..Default::default()
        },
    );
    let refused = databank::entry("unobtainium", Some(&overlay));
    assert!(refused.is_err(), "a partial substance was completed");
}

#[test]
fn a_name_in_neither_source_is_refused_by_class() {
    // The class matters: Python raises `PropertyUnavailableError` for this, and a
    // resolution failure inside Rust has to be the same kind of thing or the two
    // languages report one lookup two ways. Sabotage: leave it as `InvalidInput`.
    let missing = databank::entry("unobtainium", Some(&databank::Overlay::new()));
    assert!(
        matches!(
            missing,
            Err(azoth_core::AzothError::PropertyUnavailable { .. })
        ),
        "expected PropertyUnavailable, got {missing:?}"
    );
}

#[test]
fn the_overlay_s_components_are_available() {
    // Sabotage: drop the union in `names`.
    let mut overlay = databank::Overlay::new();
    overlay.set_component("unobtainium", databank::ComponentOverride::default());
    let names = databank::names(Some(&overlay));
    assert!(names.contains(&"unobtainium".to_string()));
    assert!(
        names.contains(&"methane".to_string()),
        "the table is still there"
    );
}

#[test]
fn the_table_is_the_file() {
    // `all_entries` and `all_kij` are the file, and a card's rows are not in the file:
    // `python/tests/test_data_agreement.py` compares their output against
    // `data/components/*.csv` byte for byte, and an overlay-aware table would break the
    // one test that proves both languages read one databank.
    // Sabotage: make either overlay-aware.
    let mut overlay = databank::Overlay::new();
    overlay.set_component("unobtainium", databank::ComponentOverride::default());
    overlay.set_kij("methane", "ethane", 0.9).unwrap();

    let names: Vec<&str> = databank::all_entries()
        .iter()
        .map(|e| e.name.as_str())
        .collect();
    assert!(
        !names.contains(&"unobtainium"),
        "the table gained a card's row"
    );

    let pairs = databank::all_kij();
    assert!(
        !pairs
            .iter()
            .any(|(a, b, _)| a == "methane" && b == "ethane"),
        "the table gained a card's pair"
    );
}

#[test]
fn a_mixture_needs_a_polynomial() {
    // Sabotage: skip the `cp.is_none()` check, which yields an `IdealGasModel` of zeros
    // and an enthalpy that is silently wrong.
    let mut overlay = databank::Overlay::new();
    overlay.set_component(
        "unobtainium",
        databank::ComponentOverride {
            tc: Some(500.0),
            pc: Some(2.0e6),
            omega: Some(0.3),
            ..Default::default()
        },
    );
    let refused = databank::mixture_of(&["unobtainium"], Some(&overlay));
    assert!(refused.is_err(), "a mixture was built without a polynomial");
}

#[test]
fn a_card_component_with_a_polynomial_has_an_enthalpy() {
    let mut overlay = databank::Overlay::new();
    overlay.set_component(
        "unobtainium",
        databank::ComponentOverride {
            tc: Some(500.0),
            pc: Some(2.0e6),
            omega: Some(0.3),
            cp: Some([20.0, 0.1, 0.0, 0.0, 0.0]),
        },
    );
    let (_, ideal_gas) = databank::mixture_of(&["unobtainium"], Some(&overlay))
        .expect("a card component with a polynomial has an enthalpy");
    assert_eq!(ideal_gas.cp_a, vec![20.0]);
    assert_eq!(ideal_gas.cp_b, vec![0.1]);
}

#[test]
fn an_overlay_changes_the_mixture_it_builds() {
    // The end of the path: a card reaches a `Mixture`, not only a lookup. Both
    // directions, because a test that checked only the first would pass with the
    // overlay ignored.
    let (base, _) = databank::mixture_of(&["methane", "n-butane"], None).unwrap();
    let (carded, _) = databank::mixture_of(&["methane", "n-butane"], Some(&omega_only())).unwrap();

    assert!(
        (carded.components()[0].omega - 0.5).abs() < 1e-12,
        "the card did not reach the mixture"
    );
    assert_ne!(
        carded.components()[0].omega,
        base.components()[0].omega,
        "the shipped and carded mixtures must differ"
    );
    assert_eq!(
        carded.components()[0].tc,
        base.components()[0].tc,
        "Tc was not named, so it is the shipped one"
    );
}

#[test]
fn an_empty_overlay_is_the_shipped_data() {
    // The case that keeps the comparison honest: if this one fails, every other
    // difference above is unreadable.
    let mut empty = databank::Overlay::new();
    assert!(empty.is_empty());
    empty.set_kij("methane", "n-butane", 0.0).unwrap();
    assert!(!empty.is_empty(), "an overlay holding a pair is not empty");
    assert!(databank::Overlay::default().is_empty());
}
