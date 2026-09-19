//! The component databank: names in, a mixture out.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::mixture::RootSide;

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
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("a mixture");
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
    assert_eq!(databank::kij("methane", "krypton", Cubic::Pr, None), 0.0);
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
    assert!(databank::mixture_of(&[], Cubic::Pr, None).is_err());
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
    let fitted = databank::kij("methane", "n-butane", Cubic::Pr, None);
    assert!(
        fitted.abs() > 1e-6,
        "the fixture needs a fitted non-zero pair, got {fitted}"
    );

    let mut overlay = databank::Overlay::new();
    overlay.set_kij("methane", "n-butane", 0.0).unwrap();
    assert_eq!(
        databank::kij("methane", "n-butane", Cubic::Pr, Some(&overlay)),
        0.0,
        "a card's explicit zero did not beat the fitted value"
    );
}

#[test]
fn a_kij_override_is_paired_whichever_order_it_is_read_in() {
    // Sabotage: store only the ordering the caller wrote.
    let mut overlay = databank::Overlay::new();
    overlay.set_kij("n-butane", "methane", 0.5).unwrap();
    assert_eq!(
        databank::kij("methane", "n-butane", Cubic::Pr, Some(&overlay)),
        0.5
    );
    assert_eq!(
        databank::kij("n-butane", "methane", Cubic::Pr, Some(&overlay)),
        0.5
    );
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
            .any(|(a, b, _, _)| a == "methane" && b == "ethane"),
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
    let refused = databank::mixture_of(&["unobtainium"], Cubic::Pr, Some(&overlay));
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
            association: None,
        },
    );
    let (_, ideal_gas) = databank::mixture_of(&["unobtainium"], Cubic::Pr, Some(&overlay))
        .expect("a card component with a polynomial has an enthalpy");
    assert_eq!(ideal_gas.cp_a, vec![20.0]);
    assert_eq!(ideal_gas.cp_b, vec![0.1]);
}

#[test]
fn an_overlay_changes_the_mixture_it_builds() {
    // The end of the path: a card reaches a `Mixture`, not only a lookup. Both
    // directions, because a test that checked only the first would pass with the
    // overlay ignored.
    let (base, _) = databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).unwrap();
    let (carded, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, Some(&omega_only())).unwrap();

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

/// The association parameters, read from the table rather than restated from it.
///
/// Water's are CPA's published set: `eps/R` of 2003.1 K against the literature's
/// 2003.4, and `kappa_AB` of 0.0692 exactly, which is how the energy's unit - J/mol,
/// not kelvin - is fixed.
#[test]
fn an_associating_component_carries_its_scheme_and_parameters() {
    let water = databank::entry("water", None).expect("water is in the databank");
    let a = water.association.expect("water names a scheme");
    assert_eq!(a.scheme, azoth_eos::association::SiteScheme::FourC);
    assert_eq!(a.sites, 4);
    assert!((a.energy - 16655.0).abs() < 1e-9, "eps: {}", a.energy);
    assert!(
        (a.volume_srk - 0.0692).abs() < 1e-12,
        "beta: {}",
        a.volume_srk
    );

    // The fitted covolume is NeqSim's internal scale, and it is *not* the cubic's own
    // `0.08664 R Tc/Pc` - which is 2.11 in the same scale. A CPA mixture's equation of
    // state is not determined by `Tc` and `Pc`, so a loader that quietly used the cubic's
    // value would be wrong by 45% on water.
    assert!(
        (a.b_srk - 1.4515).abs() < 1e-9,
        "the fitted covolume: {}",
        a.b_srk
    );
    assert!(
        (a.b_srk / 2.1127 - 1.0).abs() > 0.3,
        "and it must not be the cubic's own"
    );
}

#[test]
fn a_component_the_table_gives_no_scheme_has_none() {
    for name in ["methane", "ethane", "n-butane", "nitrogen"] {
        let entry = databank::entry(name, None).expect("in the databank");
        assert!(
            entry.association.is_none(),
            "{name} carries `0` in the scheme column, which is the table's marker for \
             a component with no scheme"
        );
    }
}

/// The four schemes the table names, each mapped to the one the kernel knows.
#[test]
fn every_scheme_the_table_names_is_carried() {
    let cases = [
        ("water", azoth_eos::association::SiteScheme::FourC),
        ("methanol", azoth_eos::association::SiteScheme::TwoB),
        ("co2", azoth_eos::association::SiteScheme::TwoA),
        ("benzene", azoth_eos::association::SiteScheme::OneA),
    ];
    for (name, scheme) in cases {
        let entry = databank::entry(name, None).expect("in the databank");
        assert_eq!(
            entry.association.map(|a| a.scheme),
            Some(scheme),
            "{name}'s scheme"
        );
    }
}

/// Two rows carry a blank `associationboundingvolume_pr`, and both are components whose
/// every sibling association cell is already zero.
///
/// This pins the loader's blank-means-zero rule to the data it was written for. If
/// upstream fills those cells - or blanks a different one - this fails, rather than the
/// rule silently becoming wrong.
#[test]
fn the_only_blank_association_cells_are_pr_volumes_on_zero_rows() {
    for name in ["h2so4", "hno3"] {
        let entry = databank::entry(name, None).expect("in the databank");
        let a = entry.association.expect("names a scheme");
        assert_eq!(a.volume_pr, 0.0);
        assert_eq!(a.energy, 0.0);
        assert_eq!(a.a_pr, 0.0);
        assert_eq!(a.b_pr, 0.0);
        assert!(!a.scheme.self_bonds(), "{name} is a 1A component");
    }
}

/// The association parameters reach the thing a model is handed, not just the table.
///
/// The path is table -> `Entry` -> `Component` -> `Mixture`, and a break anywhere in it
/// leaves a CPA model silently running on critical constants - which for water is 45%
/// wrong on the covolume and is not a difference any test downstream would attribute to
/// the loader.
#[test]
fn a_mixture_carries_its_components_association_parameters() {
    use azoth_eos::association::SiteScheme;

    let (mixture, _) =
        databank::mixture_of(&["water", "methane"], Cubic::Pr, None).expect("a mixture");
    let water = mixture.components()[0]
        .association
        .as_ref()
        .expect("water names a scheme");
    assert_eq!(water.scheme, SiteScheme::FourC);
    assert!((water.energy - 16655.0).abs() < 1e-9);
    assert!((water.b_srk - 1.4515).abs() < 1e-9);

    // And a component the table gives no scheme to stays without one, rather than
    // inheriting a neighbour's.
    assert!(
        mixture.components()[1].association.is_none(),
        "methane is not associating"
    );
}

/// A mixture with no associating component says so, rather than failing.
#[test]
fn a_mixture_of_non_associating_components_has_none() {
    let (mixture, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("a mixture");
    assert!(
        mixture.components().iter().all(|c| c.association.is_none()),
        "neither methane nor n-butane names a scheme"
    );
}

/// The databank resolves the association parameters NeqSim's own kernel reads.
///
/// The check is at NeqSim's state rather than this library's: `CpaProbe` prints the site
/// fractions and `dFCPAdN` for water/methanol 0.6/0.4 at 300 K and 100 bar, with a total
/// volume of `2.62032674828829e-5 m3/mol`. Feeding those through the mixture's own
/// resolved parameters is what proves the *selection* is right - the SRK family and not
/// the PR one, SI and not NeqSim's internal scale, and the fitted covolumes rather than
/// the cubic's.
#[test]
fn an_associating_mixture_resolves_neqsims_parameters() {
    use azoth_eos::Cubic;
    // `SystemSrkCPA`, so the cubic is Soave's - and the family matters: the association
    // parameters are per cubic, and water's `kappa_AB` is 0.0692 for SRK against
    // 0.046473789 for PR. `Cubic::default()` is PR, so a mixture that did not say which
    // it is would read the other family's fluid.
    let (mixture, _) =
        databank::mixture_of(&["water", "methanol"], Cubic::Srk, None).expect("a mixture");
    let mixture = mixture
        .with_cubic(Cubic::Srk)
        .with_association()
        .expect("water and methanol bond");
    let association = mixture.association().expect("an associating mixture");

    // `bcpa_srk` in SI, which is the column times 1e-5.
    let covolumes = [1.4515e-5, 3.0978e-5];
    let state = association
        .solve(&covolumes, &[0.6, 0.4], 2.620_326_748_288_29e-5, 300.0)
        .expect("a solvable state");
    for (site, want) in [(0, 0.101_330_316_299_160), (4, 0.031_557_041_194_167_5)] {
        assert!(
            (state.fractions[site] / want - 1.0).abs() < 1.0e-10,
            "xsite[{site}]: {} vs NeqSim {want}",
            state.fractions[site]
        );
    }
    for (i, want) in [(0, -9.807_081_619_647_33), (1, -8.298_303_821_393_14)] {
        assert!(
            (state.ln_phi[i] / want - 1.0).abs() < 1.0e-10,
            "dFCPAdN[{i}]: {} vs NeqSim {want}",
            state.ln_phi[i]
        );
    }
}

/// The associating interaction column is its own column, and a card may state it.
///
/// Water/methanol is `-0.153` in `cpakij_SRK` against `KIJPR`'s `-0.0789` - a factor of
/// two - so an override of one is not an override of the other, and a mixture that read
/// the wrong one would be a different fluid.
#[test]
fn the_associating_interaction_column_is_its_own() {
    use azoth_eos::association::AssociationCubic;

    let names = ["water", "methanol"];
    let table = databank::cpa_kij(&names, AssociationCubic::Srk, None);
    assert!(
        (table[1] - -0.153).abs() < 1.0e-12,
        "the CPA column: {}",
        table[1]
    );

    let mut overlay = databank::Overlay::new();
    overlay
        .set_cpa_kij("water", "methanol", AssociationCubic::Srk, -0.08)
        .expect("a pair");
    let carded = databank::cpa_kij(&names, AssociationCubic::Srk, Some(&overlay));
    assert!((carded[1] - -0.08).abs() < 1.0e-12, "the card's value wins");
    // The family is part of the key, so an SRK override is not a PR one.
    let pr_table = databank::cpa_kij(&names, AssociationCubic::Pr, None);
    let pr_carded = databank::cpa_kij(&names, AssociationCubic::Pr, Some(&overlay));
    assert!(
        (pr_carded[1] - pr_table[1]).abs() < 1.0e-12,
        "the PR column is its own: {} against the table's {}",
        pr_carded[1],
        pr_table[1]
    );

    // And the classical column is untouched by a CPA override, in both directions.
    let (classical, _) =
        databank::mixture_of(&names, Cubic::Pr, Some(&overlay)).expect("a mixture");
    assert!(
        (classical.kij(0, 1) - -0.0789).abs() < 1.0e-9,
        "the classical pair: {}",
        classical.kij(0, 1)
    );
}

/// The one call resolves what the four-step recipe did, and the cubic it is given is the
/// one all three of its uses read.
#[test]
fn an_associating_mixture_is_one_call() {
    use azoth_eos::{Cubic, RootSide};

    let names = ["water", "methanol"];
    let (srk, _) = databank::associating_mixture_of(&names, Cubic::Srk, None).expect("SRK-CPA");
    let (pr, _) = databank::associating_mixture_of(&names, Cubic::Pr, None).expect("PR-CPA");

    // The family is read from the cubic, so the two differ - water's `kappa_AB` is
    // 0.0692 for SRK against 0.046473789 for PR, which is a different fluid.
    let state = |mixture: &azoth_eos::Mixture| {
        let reduced = mixture
            .reduced_parameters(
                azoth_core::units::kelvins(300.0),
                azoth_core::units::pascals(1.0e7),
            )
            .expect("reduced parameters");
        mixture
            .phase_state(&reduced, &[0.6, 0.4], RootSide::Liquid)
            .expect("a liquid root")
    };
    assert!(
        (state(&srk).z / state(&pr).z - 1.0).abs() > 1.0e-4,
        "the cubic selects the association family, so the two roots differ"
    );
}

/// The CPA root, against NeqSim's, and the three pieces that get it there.
///
/// This test began as the record of a 47% divergence. It is now an equality: the liquid
/// root agrees with `CpaProbe`'s `0.105050962879418` for water/methanol 0.6/0.4 at 300 K
/// and 100 bar to better than one part in a million, and the reduced parameters it is
/// solved from agree to `1e-6`.
///
/// Three things were wrong, and they were found in this order:
///
/// **The association carries a pressure, so the root is not the cubic's.** `PhaseSrkCPA`
/// `.molarVolume` solves `BonV - (B/n) dFdV() - P B/(n R T) = 0`, and
/// `PhaseSrkCPA.dFdV()` is `super.dFdV() + dFCPAdV()` - the cubic *plus* the association.
/// NeqSim's root is where the *total* pressure equals the specified one, and the cubic's
/// own root at NeqSim's parameters is `0.15229232`, 47% away. `Mixture::associating_root`
/// solves the residual `f(Z) = Z - Z/(Z-B) + A Z/((Z+d1 B)(Z+d2 B)) - Pa Z/P`, whose first
/// three terms are the cubic's own.
///
/// **The interaction matrix is `cpakij_SRK`, not `KIJPR`.** `SystemSrkCPA` runs
/// `setMixingRule(10)`, a CPA rule, and it reads the `cpa` columns. Water/methanol is
/// `-0.153` there against `KIJPR`'s `-0.0789`, which is a 3.4% difference in the mixture's
/// `A` - and the wrong `A` was *cancelling* a second error, which is why the two had to be
/// fixed together.
///
/// **The CPA volume translation is not applied here, and applying it was the second
/// error.** `ComponentSrk.getVolumeCorrection` is
/// `0.40768 (0.29441 - Z_RA) R Tc/Pc`, non-zero for water (`racketZCPA` is `0.296941807`)
/// and zero for methanol. Adding it moved the root *away* by 0.58%; without it the root is
/// NeqSim's to seven figures. **Why NeqSim returns a non-zero correction whose effect is
/// not in its reported root is not resolved** - either `PhaseSrkEos.molarVolume` does not
/// apply it on this path, or `getZ` is computed from the uncorrected volume. That is a
/// read of the method's tail, not a guess, and it is the next question if a fluid ever
/// needs the translation.
///
/// **Two causes were guessed from the size of the gap before any of this and both were
/// wrong.** What settled it was extending `CpaProbe` to print the phase's own `A`, `B` and
/// `Z` and then reading `molarVolume`. When a port disagrees, print the source's
/// intermediates - the reduced parameters below are asserted for exactly that reason.
#[test]
fn the_cpa_root_matches_neqsim() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide};

    let names = ["water", "methanol"];
    let (mixture, _) = databank::associating_mixture_of(&names, Cubic::Srk, None)
        .expect("water and methanol bond");
    let reduced = mixture
        .reduced_parameters(kelvins(300.0), pascals(1.0e7))
        .expect("reduced parameters");
    let state = mixture
        .phase_state(&reduced, &[0.6, 0.4], RootSide::Liquid)
        .expect("a liquid root");

    // The fluid is the same one before the root is compared at all: the probe reports a
    // dimensional `B` of `2.11002000000000`, which reduces to `0.0845923635`.
    let b_mix: f64 = reduced.b[0] * 0.6 + reduced.b[1] * 0.4;
    assert!(
        (b_mix / 0.084_592_363_5 - 1.0).abs() < 1.0e-9,
        "the reduced B: azoth {b_mix} vs NeqSim 0.0845923635"
    );

    assert!(
        (state.z / 0.105_050_962_879_418 - 1.0).abs() < 1.0e-3,
        "the liquid root: azoth {} vs NeqSim 0.105050962879418 - a gap here means the \
         volume root, the CPA interaction matrix or the association has regressed",
        state.z
    );
}

/// The association's derivative surface, against finite differences of the phase state.
///
/// `phase_derivatives` is what the second-order flash, both saturation operators and the
/// envelope solve on, and for an associating mixture every one of its three families
/// crosses a surface the cubic's does not: the root comes from `associating_root`, so the
/// volume responds to `n_j`, `T` and `P` through the association's pressure as well as the
/// cubic's.
///
/// **The oracle is the phase state, not the derivation.** `phase_state` re-solves the
/// associating root at each perturbed state, so this is the one check that sees a term
/// missing from the root's sensitivities - which is where the difference lives. A check
/// that shared the derivation could not: the same identity assembled twice agrees with
/// itself, which is what cost this tranche a session on the pure `dV^2` companion.
#[test]
fn an_associating_mixtures_derivative_surface_matches_the_phase_state() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide};

    let names = ["water", "methanol"];
    let (mixture, _) = databank::associating_mixture_of(&names, Cubic::Srk, None)
        .expect("water and methanol bond");
    let (t, p) = (320.0, 2.0e6);
    let x = [0.674_266_566_871_552, 0.325_733_433_128_448];
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("reduced parameters");
    let state = mixture
        .phase_state(&reduced, &x, RootSide::Liquid)
        .expect("a liquid root");
    let d = mixture
        .phase_derivatives(&reduced, &x, state.z)
        .expect("the derivative surface");

    // The surface's own `ln phi` is the state's, which is the guard that makes the three
    // comparisons below comparisons of derivatives rather than of different quantities.
    for i in 0..2 {
        assert!(
            (d.ln_phi[i] - state.ln_phi[i]).abs() < 1e-12,
            "ln phi_{i}: surface {} vs state {}",
            d.ln_phi[i],
            state.ln_phi[i]
        );
    }

    let ln_phi_at = |temperature: f64, pressure: f64, composition: &[f64]| -> Vec<f64> {
        let reduced = mixture
            .reduced_parameters(kelvins(temperature), pascals(pressure))
            .expect("reduced parameters");
        mixture
            .phase_state(&reduced, composition, RootSide::Liquid)
            .expect("a liquid root")
            .ln_phi
    };

    // Composition, at constant T and P. The perturbation is on a mole number at unit
    // total, so the mole fractions move with it: `x_i = n_i/(sum n)` is what makes this
    // the normalised partial NeqSim reports, and the one `d_ln_phi_dn` is written in.
    let h = 1.0e-7;
    for j in 0..2 {
        let mut up = x;
        up[j] += h;
        let mut down = x;
        down[j] -= h;
        let up = up.map(|value| value / (1.0 + h));
        let down = down.map(|value| value / (1.0 - h));
        let (up, down) = (ln_phi_at(t, p, &up), ln_phi_at(t, p, &down));
        for i in 0..2 {
            let numerical = (up[i] - down[i]) / (2.0 * h);
            assert!(
                (d.d_ln_phi_dn[i][j] / numerical - 1.0).abs() < 1e-5,
                "d ln phi_{i}/dn_{j}: analytic {} vs numerical {numerical}",
                d.d_ln_phi_dn[i][j]
            );
        }
    }

    // Temperature at constant P and composition, and pressure at constant T.
    let step = 1.0e-3;
    let (up, down) = (ln_phi_at(t + step, p, &x), ln_phi_at(t - step, p, &x));
    for i in 0..2 {
        let numerical = (up[i] - down[i]) / (2.0 * step);
        assert!(
            (d.d_ln_phi_dt[i] / numerical - 1.0).abs() < 1e-5,
            "d ln phi_{i}/dT: analytic {} vs numerical {numerical}",
            d.d_ln_phi_dt[i]
        );
    }
    let step = 100.0;
    let (up, down) = (ln_phi_at(t, p + step, &x), ln_phi_at(t, p - step, &x));
    for i in 0..2 {
        let numerical = (up[i] - down[i]) / (2.0 * step);
        assert!(
            (d.d_ln_phi_dp[i] / numerical - 1.0).abs() < 1e-5,
            "d ln phi_{i}/dP: analytic {} vs numerical {numerical}",
            d.d_ln_phi_dp[i]
        );
    }
}

/// The shipped keycard is the databank, at the level a model reads it.
///
/// `databank/keycard.toml` is generated from the same table the library ships, so
/// resolving an associating mixture *through* it must give the same fluid as resolving it
/// with no card at all. That is the property the generator exists for, and the association
/// is what makes it non-trivial: the card states the fitted attraction and covolume in SI
/// and the table carries them in NeqSim's internal scale, so this fails by a factor of a
/// hundred thousand if either crossing is missing.
#[test]
fn the_shipped_card_is_the_databank() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide};

    let path =
        std::path::PathBuf::from(env!("CARGO_MANIFEST_DIR")).join("../../databank/keycard.toml");
    let card = azoth_eos::card::Card::from_path(&path).expect("the shipped card reads");

    let names = ["water", "methanol"];
    let (shipped, _) = databank::associating_mixture_of(&names, Cubic::Srk, None).expect("SRK-CPA");
    let (carded, _) = databank::associating_mixture_of(&names, Cubic::Srk, Some(card.overlay()))
        .expect("SRK-CPA through the card");

    let root = |mixture: &azoth_eos::Mixture| {
        let reduced = mixture
            .reduced_parameters(kelvins(300.0), pascals(1.0e7))
            .expect("reduced parameters");
        mixture
            .phase_state(&reduced, &[0.6, 0.4], RootSide::Liquid)
            .expect("a liquid root")
            .z
    };
    assert!(
        (root(&carded) / root(&shipped) - 1.0).abs() < 1.0e-12,
        "the card's fluid: {} against the table's {}",
        root(&carded),
        root(&shipped)
    );
}

/// The associating root at low pressure, against NeqSim's own flash.
///
/// **This is the regime the root finder used to get wrong.** At 1 bar the substituted
/// cubic's attraction is small enough that it has one real root, so `z_min` and `z_max`
/// are the same number - the vapour root - and Newton from it converges there and stops.
/// The liquid root exists only because the association's pressure creates it, so it has
/// to be found rather than seeded, and before the geometric walk it was not: every
/// low-pressure CPA flash saw the same root twice and reported a single phase.
///
/// The oracle is NeqSim's `TPflash`, which splits this feed. Its gas phase is at
/// `Z = 0.949619265923351` and its aqueous phase at `0.000876327166634292`, and the two
/// tests below take one each.
#[test]
fn the_cpa_vapour_root_matches_neqsim_at_low_pressure() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide};

    let (mixture, _) = databank::associating_mixture_of(&["water", "methanol"], Cubic::Srk, None)
        .expect("water and methanol bond");
    let (t, p) = (356.0, 1.0e5);
    let x = [0.6, 0.4];
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("reduced parameters");
    let state = mixture
        .phase_state(&reduced, &x, RootSide::Vapour)
        .expect("a vapour root");
    assert!(
        (state.z / 0.949_619_265_923_351 - 1.0).abs() < 1.0e-9,
        "the vapour root: azoth {} vs NeqSim 0.949619265923351",
        state.z
    );

    // The kernel at that state too, which `CpaProbe` prints: `FCPA` is
    // `A_assoc/(R T)` and `dFCPAdN[i]` is the association's `ln phi_i`.
    let association = mixture.association().expect("an associating mixture");
    let covolumes: Vec<f64> = reduced.b.iter().map(|b| b * 8.3144621 * t / p).collect();
    let kernel = association
        .solve(&covolumes, &x, state.z * 8.3144621 * t / p, t)
        .expect("a solvable state");
    assert!(
        (kernel.helmholtz_rt / -0.049_437_345_696_229_4 - 1.0).abs() < 1.0e-9,
        "FCPA: {} vs NeqSim -0.0494373456962294",
        kernel.helmholtz_rt
    );
    for (i, want) in [(0, -0.085_911_899_466_216_2), (1, -0.113_824_565_989_815)] {
        assert!(
            (kernel.ln_phi[i] / want - 1.0).abs() < 1.0e-9,
            "dFCPAdN[{i}]: {} vs NeqSim {want}",
            kernel.ln_phi[i]
        );
    }
}

/// The liquid root at the same state, at NeqSim's own aqueous composition.
///
/// `CpaProbe` prints the flashed aqueous phase's `Z` and mole fractions, so this compares
/// the root azoth finds against the root NeqSim's flash settled on rather than against a
/// hand-derived value.
#[test]
fn the_cpa_liquid_root_is_found_at_low_pressure() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide};

    let (mixture, _) = databank::associating_mixture_of(&["water", "methanol"], Cubic::Srk, None)
        .expect("water and methanol bond");
    let (t, p) = (356.0, 1.0e5);
    let x = [0.674_266_566_871_552, 0.325_733_433_128_448];
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("reduced parameters");
    let liquid = mixture
        .phase_state(&reduced, &x, RootSide::Liquid)
        .expect("a liquid root");
    let vapour = mixture
        .phase_state(&reduced, &x, RootSide::Vapour)
        .expect("a vapour root");

    assert!(
        (liquid.z / 0.000_876_327_166_634_292 - 1.0).abs() < 1.0e-3,
        "the liquid root: azoth {} vs NeqSim 0.000876327166634292",
        liquid.z
    );
    // And it is the *lower* root, which is what makes this a test of the branch rather
    // than of a number: a finder that returned the vapour root twice would pass the
    // comparison above only by coincidence, and cannot pass this one.
    assert!(
        liquid.z < vapour.z / 100.0,
        "the two sides found the same root: {} against {}",
        liquid.z,
        vapour.z
    );
}

/// **The cubic `ln phi` is the derivative of the cubic Helmholtz energy at the volume it is
/// given, and the closed form is that derivative only at the cubic's own root.**
///
/// The closed form - `(b_i/B)(Z-1) - ln(Z-B) - (A/(B(d1-d2)))(2*abar_i/a - b_i/B)*ln(...)` -
/// carries a `Z - 1`. That `Z - 1` is the *equation of state's* `Z - 1`, written as
/// `B/(Z-B) - A Z/((Z+d1 B)(Z+d2 B))`, and the two agree only where the volume solves the
/// SRK equation. For a pure component at 350 K and 50 bar: at the cubic's own root the closed
/// form is `-10.4890220129` against a finite difference of that Helmholtz energy of
/// `-10.48902202`, and at `Z = 0.5` it is `-1.94579104285` against `-3.47184011449`.
///
/// An associating mixture's volume is not the cubic's root - the association carries a
/// pressure - so the two are not interchangeable here, and they differ by `0.0328` in
/// component 0 and `0.0699` in component 1 at 356 K and 1 bara. `Cubic::eos_z_minus_one` is
/// the equation of state's, and the model now uses it. **That gap was recorded as NeqSim's
/// for a session**, on the strength of a probe that subtracted `FV()` where the chain rule
/// needs `dFdV()` - `FV() + dFCPAdV()`, and the association's `dFCPAdV` is seventeen times
/// the cubic's `FV()` here. `validation/neqsim/CpaFdProbe.java` now prints both, and
/// `CpaSweep`'s `euler_*` keys say the same without a finite difference in it.
///
/// The derivative is computed three ways below: NeqSim's `dFdN` agrees with this model's
/// cubic part, a finite difference of *this model's own* Helmholtz energy at a fixed volume
/// agrees with it, and the closed form does not.
#[test]
fn the_cubic_ln_phi_is_the_derivative_at_the_association_shifted_root() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide};

    let (mixture, _) = databank::associating_mixture_of(&["water", "methanol"], Cubic::Srk, None)
        .expect("water and methanol bond");
    let (t, p) = (356.0, 1.0e5);
    let x = [0.6, 0.4];
    let reduced = mixture
        .reduced_parameters(kelvins(t), pascals(p))
        .expect("reduced parameters");
    let state = mixture
        .phase_state(&reduced, &x, RootSide::Vapour)
        .expect("a vapour root");
    let association = mixture.association().expect("an associating mixture");
    let covolumes: Vec<f64> = reduced.b.iter().map(|b| b * 8.3144621 * t / p).collect();
    let kernel = association
        .solve(&covolumes, &x, state.z * 8.3144621 * t / p, t)
        .expect("a solvable state");

    let (z, b_red) = (state.z, state.b_mix);
    let (d1, d2) = (Cubic::Srk.delta1(), Cubic::Srk.delta2());
    let i_term = ((z + d1 * b_red) / (z + d2 * b_red)).ln();
    let coefficient = state.a_mix / ((d1 - d2) * b_red);

    // NeqSim's `lnPhi[i]` and `dFCPAdN[i]` at this state, from `CpaSweep 356 1 0.6`. Their
    // difference is NeqSim's cubic part, which `CpaFdProbe` measures to be `dF/dn_i` there.
    let neqsim_cubic = [
        -0.038_499_945_545_132_2 - -0.085_911_899_466_216_2,
        -0.069_416_678_113_187_5 - -0.113_824_565_989_815,
    ];

    // `d(F/RT)/dn_i` at constant `V`, from this model's own Helmholtz energy. Its
    // compressibility argument is the state's `Z`, and holding it is what holds the volume:
    // the energy is homogeneous of degree one in `(V, n)`, so `n_i` moving with `Z` fixed
    // moves the composition and the total together and the volume not at all.
    let h = 1.0e-7;
    let helmholtz_at = |n: [f64; 2]| {
        mixture
            .helmholtz_energy(&reduced, &n, z)
            .expect("a Helmholtz energy")
    };

    for (i, &neqsim) in neqsim_cubic.iter().enumerate() {
        let abar: f64 = (0..2)
            .map(|j| x[j] * (1.0 - mixture.kij(i, j)) * (reduced.a[i] * reduced.a[j]).sqrt())
            .sum();
        let b_ratio = reduced.b[i] / b_red;
        let factor = 2.0 * abar / state.a_mix - b_ratio;
        let closed_form = b_ratio * (z - 1.0) - (z - b_red).ln() - coefficient * factor * i_term;

        let mut up = x;
        up[i] += h;
        let mut down = x;
        down[i] -= h;
        // The energy now carries the association, so its gradient is the *whole*
        // fugacity coefficient - `ln phi` itself, once `ln z` is taken off - and not the
        // cubic's share of it. That is the identity `helmholtz_energy`'s own test pins,
        // arriving here as the check that the cubic's share is what the two routes agree
        // about once the association is added back.
        let finite_difference =
            (helmholtz_at(up) - helmholtz_at(down)) / (2.0 * h) - z.ln() - kernel.ln_phi[i];

        let azoth_cubic = state.ln_phi[i] - kernel.ln_phi[i];
        assert!(
            (azoth_cubic / neqsim - 1.0).abs() < 1.0e-8,
            "component {i}: this model's cubic ln phi is {azoth_cubic} and NeqSim's \
             `dFdN` less its association term is {neqsim}"
        );
        assert!(
            (azoth_cubic / finite_difference - 1.0).abs() < 1.0e-6,
            "component {i}: this model's cubic ln phi is {azoth_cubic} but a finite \
             difference of its own Helmholtz energy at a fixed volume, less the \
             association, is {finite_difference}"
        );
        assert!(
            (closed_form / azoth_cubic - 1.0).abs() > 0.01,
            "component {i}: the closed form {closed_form} has met the derivative \
             {azoth_cubic} - which would mean this volume IS the cubic's root, so the \
             association has stopped carrying a pressure"
        );
    }
}

/// The PC-SAFT columns are read, and zero is how the table spells absence.
///
/// The three are vendored and unparsed until a model reads them, and 47 of the 286 rows
/// carry `0` in all three rather than a blank - so "the table has no PC-SAFT set for this
/// substance" is a *value*, and a model that computed with `m = 0` would be solving for a
/// fluid with no segments rather than refusing one it has no parameters for.
///
/// Water's numbers are the check against the vendored file rather than against this crate:
/// `m = 1.0656`, `sigma = 3.0007 Å`, `epsilon/k = 366.51 K`, which is what NeqSim's
/// `COMP.csv` carries and what `Component.java:547-549` reads, with the ångström division
/// already done by the generator.
#[test]
fn the_saft_vr_mie_columns_are_read_and_zero_means_absent() {
    // Methane's published set: `lambda_r = 12.65`, `lambda_a = 6`, `m = 1`, `sigma =
    // 3.7412` angstrom, `epsilon/k = 153.36`.
    let methane = databank::entry("methane", None).expect("methane");
    assert!(
        (methane.lambda_r_mie - 12.65).abs() < 1e-9,
        "lambda_r: {}",
        methane.lambda_r_mie
    );
    assert!(
        (methane.lambda_a_mie - 6.0).abs() < 1e-9,
        "lambda_a: {}",
        methane.lambda_a_mie
    );
    assert!((methane.m_mie - 1.0).abs() < 1e-12, "m: {}", methane.m_mie);
    assert!(
        (methane.sigma_mie - 3.7412e-10).abs() < 1e-14,
        "sigma: {} - the table is SI, so this is metres and not angstrom",
        methane.sigma_mie
    );
    assert!(
        (methane.epsik_mie - 153.36).abs() < 1e-9,
        "epsilon/k: {}",
        methane.epsik_mie
    );

    // A substance with no SAFT-VR-Mie set is absent as a zero in all five, which is the
    // shape a model has to refuse rather than compute with. **Twelve of the table's 286
    // rows carry one** - the light alkanes, `co2`, nitrogen and water - and the count is
    // asserted rather than left to a sample, so a generator that stopped emitting the
    // columns fails here rather than at a model. Counting it off the CSV with a
    // comma-splitting tool gives 58, because the file quotes fields that contain commas.
    let carried = databank::names(None)
        .into_iter()
        .filter(|name| databank::entry(name, None).is_ok_and(|e| e.m_mie > 0.0))
        .count();
    assert_eq!(carried, 12, "rows carrying a SAFT-VR-Mie set");

    // **`m` is the absence marker and `lambda_r` is not.** The table carries the standard
    // `12`/`6` on every row, so methanol - which has no set - reports a repulsive exponent
    // of 12 beside a segment number of zero. A model that keyed absence on `lambda_r`
    // would solve for a fluid with no segments on 274 of the 286 rows.
    let absent = databank::entry("methanol", None).expect("methanol");
    assert_eq!(absent.m_mie, 0.0, "methanol has no SAFT-VR-Mie set");
    assert_eq!(
        absent.lambda_r_mie, 12.0,
        "but its lambda_r is the table's default"
    );
    assert_eq!(absent.sigma_mie, 0.0, "and its sigma is absent, like its m");
}

/// already done by the generator.
#[test]
fn the_pcsaft_columns_are_read_and_zero_means_absent() {
    let water = databank::entry("water", None).expect("water");
    assert!((water.m_saft - 1.0656).abs() < 1e-9, "m: {}", water.m_saft);
    assert!(
        (water.sigma_saft - 3.0007e-10).abs() < 1e-14,
        "sigma: {} - the table is SI, so this is metres and not ångström",
        water.sigma_saft
    );
    assert!(
        (water.epsik_saft - 366.51).abs() < 1e-9,
        "epsilon/k: {}",
        water.epsik_saft
    );

    let methane = databank::entry("methane", None).expect("methane");
    assert!(
        (methane.m_saft - 1.0).abs() < 1e-12,
        "m: {}",
        methane.m_saft
    );

    // A substance the table carries no PC-SAFT set for: absent as a zero, in all three,
    // which is the shape a model has to refuse rather than compute with.
    let absent = databank::names(None)
        .into_iter()
        .filter(|name| {
            databank::entry(name, None).is_ok_and(|e| e.m_saft == 0.0 && e.sigma_saft == 0.0)
        })
        .count();
    assert_eq!(
        absent, 47,
        "rows with no PC-SAFT set; the vendored table has 47"
    );
}

/// The PC-SAFT interaction column is a third fit, and a sparse one.
///
/// Methane/n-butane is `0.022` here against the `0.01289789` its SRK and PR columns share,
/// which is the check that it is read from `KIJPCSAFT` rather than from a column that
/// happened to be nearby - the two are close enough that a slip would look plausible.
///
/// **NeqSim reads it in the only configuration it can run.** Its default mixing rule leaves
/// the interaction matrix null, and a mixture throws a NullPointerException before reaching
/// a `k_ij`; the classic rule is what a caller must set, and its branch is keyed on the
/// PC-SAFT phase class and reads this column. So there is no usable NeqSim PC-SAFT that
/// ignores it.
#[test]
fn the_pcsaft_interaction_column_is_read() {
    let pair = databank::pcsaft_kij(&["methane", "n-butane"]);
    assert!(
        (pair[1] - 0.022).abs() < 1e-12,
        "methane/n-butane: {} - the table's `KIJPCSAFT`, not its `kijsrk`",
        pair[1]
    );
    assert!(
        (pair[1] - databank::kij("methane", "n-butane", Cubic::Pr, None)).abs() > 1e-6,
        "the two columns must differ here, or this proves nothing"
    );
    // Symmetric, and zero for a pair the table does not carry - NeqSim's ideal-mixture
    // default rather than a failure.
    assert_eq!(pair[1], pair[2]);
    let absent = databank::pcsaft_kij(&["water", "methanol"]);
    assert_eq!(absent[1], 0.0);
}

/// **The cubic's interaction column is the cubic's, and the table says so 76 times.**
///
/// NeqSim selects on the phase class: `phase.getClass().getName().equals(
/// "neqsim.thermo.phase.PhasePrEos")` reads `KIJPR` and every other phase - SRK, RK, and
/// even `PhaseUMRCPA`, which extends `PhasePrEos` without being it - reads `KIJSRK`.
///
/// **azoth read `KIJPR` for every cubic until this.** No test could see it, because every
/// pair the validation cases use has the *same* value in both columns - methane/n-butane
/// is `0.01289789` in each, methane/propane `0.00747722`, water/methanol `-0.0789` - so
/// the whole SRK validation ran where reading the wrong column changes nothing. The
/// sweep below is the check that covers the other 76.
#[test]
fn the_interaction_column_follows_the_cubic() {
    let pairs = databank::all_kij();
    let differing: Vec<_> = pairs.iter().filter(|(_, _, pr, srk)| pr != srk).collect();
    assert_eq!(
        differing.len(),
        76,
        "the number of pairs whose `KIJSRK` and `KIJPR` differ, measured over the 516 \
         in-scope pairs"
    );

    for (a, b, pr, srk) in differing {
        assert_eq!(
            databank::kij(a, b, Cubic::Pr, None),
            *pr,
            "{a}/{b}: the PR cubic should read `KIJPR`"
        );
        assert_eq!(
            databank::kij(a, b, Cubic::Srk, None),
            *srk,
            "{a}/{b}: the SRK cubic should read `KIJSRK`"
        );
    }

    // And the resolution, not only the lookup: `mixture_of` must carry the column of the
    // cubic it was asked for, because that is what every model actually calls.
    for (cubic, expected) in [(Cubic::Pr, 0.135_000_01), (Cubic::Srk, 0.1018)] {
        let (mixture, _) =
            databank::mixture_of(&["propane", "co2"], cubic, None).expect("the pair resolves");
        assert_eq!(mixture.kij(0, 1), expected, "{cubic:?}");
    }

    // The oracle: NeqSim 3.20.0's `SystemSrkEos` at 350 K and 30 bar with `z = 0.5/0.5`,
    // from `validation/neqsim/SrkKijProbe.java`. azoth reproduces it to fifteen digits
    // with `KIJSRK` and is 0.57% out with `KIJPR`.
    let (mixture, _) = databank::mixture_of(&["propane", "co2"], Cubic::Srk, None).unwrap();
    let reduced = mixture
        .reduced_parameters(kelvins(350.0), pascals(3.0e6))
        .expect("a state");
    let state = mixture
        .phase_state(&reduced, &[0.5, 0.5], RootSide::Vapour)
        .expect("a root");
    assert!(
        (state.z / 0.823_765_416_605_303 - 1.0).abs() < 1.0e-12,
        "NeqSim's SRK gives 0.823765416605303 and this gives {}",
        state.z
    );
}
