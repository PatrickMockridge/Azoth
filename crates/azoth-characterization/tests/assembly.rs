//! Tests for `mixture_from_characterization`, the row-to-fluid step.
//!
//! The assembler has no oracle of its own: NeqSim's route is a mutation sequence on a system,
//! and the probe records the rows rather than the fluid. So each side is held to a reference it
//! already has - the databank for a base row and for a pseudo row carrying the databank's own
//! numbers - which is what makes the derivation (the Peneloux shifts, the Soave coefficient) a
//! claim rather than a restatement.

use azoth_characterization::{FluidComponent, mixture_from_characterization};
use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank;
use azoth_eos::{Cubic, Mixture, srk_kappa};

fn base_row(fluid: f64, name: &str, moles: f64) -> FluidComponent {
    FluidComponent {
        fluid,
        name: name.to_string(),
        pseudo_component: false,
        moles,
        molar_mass: None,
        density: None,
        normal_boiling_point: None,
        critical_temperature: None,
        critical_pressure: None,
        acentric_factor: None,
        critical_volume: None,
        racket_z: None,
        racket_z_cpa: None,
        parachor: None,
        critical_viscosity: None,
        triple_point_temperature: None,
        heat_of_fusion: None,
        ideal_gas_enthalpy_of_formation: None,
        cp_a: None,
        cp_b: None,
        cp_c: None,
        cp_d: None,
        attractive_m: None,
    }
}

/// A pseudo row carrying a databank substance's own numbers, as a cut characterised to match it.
fn pseudo_row(name: &str, moles: f64) -> FluidComponent {
    let entry = databank::entry(name, None).expect("the reference substance resolves");
    let cp = entry.cp.expect("the databank carries the polynomial");
    FluidComponent {
        fluid: 0.0,
        name: format!("{name}_PC"),
        pseudo_component: true,
        moles,
        molar_mass: entry.molar_mass,
        density: None,
        normal_boiling_point: None,
        critical_temperature: Some(entry.tc),
        critical_pressure: Some(entry.pc),
        acentric_factor: Some(entry.omega),
        critical_volume: entry.critical_volume,
        racket_z: Some(entry.rackett_z),
        racket_z_cpa: None,
        parachor: Some(entry.parachor),
        critical_viscosity: None,
        triple_point_temperature: None,
        heat_of_fusion: None,
        ideal_gas_enthalpy_of_formation: None,
        cp_a: Some(cp[0]),
        cp_b: Some(cp[1]),
        cp_c: Some(cp[2]),
        cp_d: Some(cp[3]),
        attractive_m: None,
    }
}

fn assert_reduced_parameters_match(ours: &Mixture, theirs: &Mixture, at: &str) {
    let (left, right) = (
        ours.reduced_parameters(kelvins(320.0), pascals(4.0e6))
            .expect("a state"),
        theirs
            .reduced_parameters(kelvins(320.0), pascals(4.0e6))
            .expect("a state"),
    );
    for (label, left, right) in [
        ("a", left.a.as_slice(), right.a.as_slice()),
        ("b", left.b.as_slice(), right.b.as_slice()),
        ("psi", left.psi.as_slice(), right.psi.as_slice()),
        ("psi_t", left.psi_t.as_slice(), right.psi_t.as_slice()),
    ] {
        for (index, (left, right)) in left.iter().zip(right).enumerate() {
            assert!(
                (left - right).abs() <= 1e-12 * right.abs().max(f64::MIN_POSITIVE),
                "{at}::{label}[{index}]: {left} against {right}"
            );
        }
    }
}

#[test]
fn a_base_only_table_is_the_databank_mixture() {
    let rows = [
        base_row(0.0, "methane", 2.0),
        base_row(0.0, "ethane", 1.0),
        base_row(1.0, "propane", 0.5),
    ];
    let (ours, ideal_gas) =
        mixture_from_characterization(&rows, Cubic::Srk).expect("the names resolve");
    let (theirs, their_polynomial) =
        databank::mixture_of(&["methane", "ethane", "propane"], Cubic::Srk, None)
            .expect("the names resolve");

    assert_eq!(ours.len(), theirs.len());
    for (index, (left, right)) in ours
        .components()
        .iter()
        .zip(theirs.components())
        .enumerate()
    {
        assert_eq!(left.tc.value, right.tc.value, "component {index}: Tc");
        assert_eq!(left.pc.value, right.pc.value, "component {index}: Pc");
        assert_eq!(left.omega, right.omega, "component {index}: omega");
        assert_eq!(
            left.molar_mass, right.molar_mass,
            "component {index}: molar mass"
        );
        assert_eq!(
            left.volume_shift, right.volume_shift,
            "component {index}: the PR shift"
        );
        assert_eq!(
            left.volume_shift_srk, right.volume_shift_srk,
            "component {index}: the SRK shift"
        );
        assert_eq!(
            left.parachor, right.parachor,
            "component {index}: the parachor"
        );
    }
    assert_reduced_parameters_match(&ours, &theirs, "base_only");

    assert_eq!(ideal_gas.cp_a, their_polynomial.cp_a);
    assert_eq!(ideal_gas.cp_b, their_polynomial.cp_b);
    assert_eq!(ideal_gas.cp_c, their_polynomial.cp_c);
    assert_eq!(ideal_gas.cp_d, their_polynomial.cp_d);
    assert_eq!(ideal_gas.cp_e, their_polynomial.cp_e);
}

#[test]
fn a_pseudo_row_carries_the_databank_derived_shifts_and_alpha() {
    let entry = databank::entry("n-heptane", None).expect("the reference substance resolves");
    let rows = [pseudo_row("n-heptane", 1.0)];
    let (ours, _) = mixture_from_characterization(&rows, Cubic::Srk).expect("the row is complete");
    let component = &ours.components()[0];

    let pr = azoth_eos::pr_peneloux_shift(
        entry.omega,
        kelvins(entry.tc),
        pascals(entry.pc),
        Some(entry.rackett_z),
    )
    .expect("a shift");
    let srk = azoth_eos::srk_peneloux_shift(
        entry.omega,
        kelvins(entry.tc),
        pascals(entry.pc),
        Some(entry.rackett_z),
    )
    .expect("a shift");
    assert_eq!(component.volume_shift, pr.c.value);
    assert_eq!(component.volume_shift_srk, srk.c.value);
    // The row states no fitted coefficient, so Soave's own from the acentric factor.
    assert_eq!(
        component.alpha_params,
        vec![srk_kappa(entry.omega).expect("a coefficient").kappa]
    );

    let (theirs, _) = databank::mixture_of(&["n-heptane"], Cubic::Srk, None)
        .expect("the reference substance resolves");
    assert_reduced_parameters_match(&ours, &theirs, "pseudo_row");
}

#[test]
fn a_pseudo_row_without_its_heat_capacity_is_refused() {
    let mut row = pseudo_row("n-heptane", 1.0);
    row.cp_c = None;
    let error = mixture_from_characterization(&[row], Cubic::Srk)
        .expect_err("an enthalpy needs the polynomial");
    assert!(
        error.to_string().contains("components[0].cp_c"),
        "the refusal should name the row's field: {error}"
    );
}

#[test]
fn a_base_name_the_databank_does_not_carry_is_refused() {
    let rows = [base_row(0.0, "unobtainium", 1.0)];
    let error = mixture_from_characterization(&rows, Cubic::Srk).expect_err("the name has no row");
    assert!(
        error.to_string().contains("unobtainium"),
        "the refusal should name the substance: {error}"
    );
}
