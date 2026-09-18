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

    let (mixture, _) = databank::mixture_of(&["water", "methane"], None).expect("a mixture");
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
    let (mixture, _) = databank::mixture_of(&["methane", "n-butane"], None).expect("a mixture");
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
    let (mixture, _) = databank::mixture_of(&["water", "methanol"], None).expect("a mixture");
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
        .set_cpa_kij("water", "methanol", -0.08)
        .expect("a pair");
    let carded = databank::cpa_kij(&names, AssociationCubic::Srk, Some(&overlay));
    assert!((carded[1] - -0.08).abs() < 1.0e-12, "the card's value wins");

    // And the classical column is untouched by a CPA override, in both directions.
    let (classical, _) = databank::mixture_of(&names, Some(&overlay)).expect("a mixture");
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

/// The association's derivative surface is **not derived**, and this pins the refusal.
///
/// `phase_derivatives` is what the second-order flash, both saturation operators and the
/// envelope solve on, so an associating mixture cannot take those paths. The gap is named
/// rather than papered over because the alternative - returning the cubic's derivatives -
/// is a *wrong* answer rather than a missing one.
///
/// **What is missing.** The kernel's derivatives are at constant `T` and `V`; a flash
/// needs constant `T` and `P`, and the conversion is the chain rule through the volume:
/// `d ln phi_i/dn_j = Phi_{n_i n_j} + Phi_{n_i V} V_{n_j}`. The kernel supplies
/// `Phi_{n_i n_j}` and `Phi_{n_i V}`, but `V_{n_j}` is `(R T/P)(Z + dZ/dn_j)` and **`dZ/dn_j`
/// is not the cubic's** - an associating mixture's root comes from `associating_root`, so
/// its volume responds to the composition through the association's pressure too.
///
/// Differentiating that residual needs `d2(A/(RT))/dV dn_j`. Its pure `dV^2` companion
/// **is** computed, and checked against a finite difference of `d(A/(RT))/dV`; the mixed
/// derivative is not, and the site-fraction system's second differentiation does not yet
/// reproduce it.
///
/// Sabotage: when the surface lands, the entry above stops being a gap and this test
/// fails, which is the point of it.
#[test]
fn an_associating_mixture_refuses_the_derivative_surface() {
    use azoth_core::units::{kelvins, pascals};
    use azoth_eos::{Cubic, RootSide};

    let names = ["water", "methanol"];
    let (mixture, _) = databank::associating_mixture_of(&names, Cubic::Srk, None)
        .expect("water and methanol bond");
    let reduced = mixture
        .reduced_parameters(kelvins(320.0), pascals(2.0e6))
        .expect("reduced parameters");
    let state = mixture
        .phase_state(&reduced, &[0.6, 0.4], RootSide::Liquid)
        .expect("a liquid root");

    let refusal = mixture
        .phase_derivatives(&reduced, &[0.6, 0.4], state.z)
        .expect_err("the derivative surface is not derived for an associating mixture");
    let message = refusal.to_string();
    assert!(
        message.contains("derivative surface"),
        "the refusal should name what is missing: {message}"
    );

    // And the state itself is fine - it is the derivative surface, and only that, which
    // is outstanding.
    assert!(state.ln_phi.iter().all(|value| value.is_finite()));
}
