//! `Alpha::SoreideWhitson` - the one alpha that is not a function of the component alone.
//!
//! The model that consumes it is not registered yet, so these are the tests that keep the
//! variant honest in the meantime: it is refused where it cannot be resolved, and it moves
//! water's attraction where it can.

use azoth_core::AzothError;
use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::alpha_term::Alpha;
use azoth_eos::mixing_rule::{MixingRule, SoreideWhitsonRole};
use azoth_eos::mixture::RootSide;

const NAMES: [&str; 3] = ["CO2", "methane", "water"];

fn roles() -> Vec<SoreideWhitsonRole> {
    vec![
        SoreideWhitsonRole::CarbonDioxide,
        SoreideWhitsonRole::Hydrocarbon,
        SoreideWhitsonRole::Water,
    ]
}

/// The rule on its own, so the matrix can be reached without a mixture - the field it
/// lives in is private and there is no accessor for it.
fn rule(salinity: f64) -> MixingRule {
    MixingRule::SoreideWhitson {
        kij: vec![0.0; 9],
        roles: roles(),
        salinity,
    }
}

fn mixture(salinity: f64) -> azoth_eos::mixture::Mixture {
    let (base, _) = azoth_eos::databank::mixture_of(&NAMES, Cubic::Pr, None).unwrap();
    base.with_alpha(Alpha::SoreideWhitson)
        .with_mixing_rule(rule(salinity))
}

/// **Refused beside any other rule.** The variant reads which component is water from the
/// rule's roles and the salinity from the same place, so a classic `kij` beside it leaves
/// it nothing to read - and resolving against a zero would be a water alpha quietly
/// computed at fresh water.
#[test]
fn it_is_refused_beside_another_rule() {
    let (base, _) = azoth_eos::databank::mixture_of(&NAMES, Cubic::Pr, None).unwrap();
    let broken = base.with_alpha(Alpha::SoreideWhitson);
    let err = broken
        .reduced_parameters(kelvins(313.15), pascals(5_000_000.0))
        .expect_err("the classic rule carries no roles and no salinity");
    assert!(matches!(err, AzothError::InvalidInput { .. }), "{err:?}");
    assert_eq!(err.field(), Some("mixing_rule"));
}

/// **The salinity moves water's attraction and nothing else's.**
///
/// `A_i = omega_a alpha_i Pr_i / Tr_i^2`, so the entry is proportional to the alpha: a port
/// that had left the brine out would have a water `A` that does not move - which is exactly
/// what NeqSim's own does, at every concentration up to `75.96 mol/kg`.
#[test]
fn the_salinity_moves_water_and_only_water() {
    let fresh = mixture(0.0)
        .reduced_parameters(kelvins(313.15), pascals(5_000_000.0))
        .unwrap();
    let salty = mixture(4.0)
        .reduced_parameters(kelvins(313.15), pascals(5_000_000.0))
        .unwrap();

    let water = 2;
    assert!(
        (salty.a[water] - fresh.a[water]).abs() > 1.0e-3,
        "water's A is {} fresh and {} at 4 mol/kg",
        fresh.a[water],
        salty.a[water]
    );
    for other in [0, 1] {
        assert_eq!(
            fresh.a[other], salty.a[other],
            "component {other} is PR78 and takes no salinity"
        );
    }

    // And the aqueous half of the matrix moves with it, which is the *other* place the
    // salinity reaches - the two are separate paths and this is the one the rule owns.
    //
    // **`x_water > 0.8` or the correlation does not apply at all**, which is NeqSim's own
    // gate and the first composition this test used (water at 0.5) hit it and returned the
    // base matrix unchanged.
    let x = [0.02, 0.08, 0.90];
    let (tr, omega) = (&[0.9, 1.1, 0.5][..], &[0.1, 0.01, 0.34][..]);
    let fresh_kij = rule(0.0).phase_kij(&[0.0; 9], tr, omega, &x);
    let salty_kij = rule(4.0).phase_kij(&[0.0; 9], tr, omega, &x);
    assert_ne!(fresh_kij, salty_kij);
    assert_eq!(
        fresh_kij[1], salty_kij[1],
        "a pair that is not water's keeps the base value"
    );
}

/// The rule and the alpha together produce a phase, which is the whole point of the pair.
#[test]
fn the_mixture_computes_a_phase() {
    let reduced = mixture(0.0)
        .reduced_parameters(kelvins(313.15), pascals(5_000_000.0))
        .unwrap();
    let state = mixture(0.0)
        .phase_state(&reduced, &[0.2, 0.3, 0.5], RootSide::Vapour)
        .unwrap();
    assert!(state.z > 0.0 && state.z < 1.0, "z = {}", state.z);
    assert_eq!(state.ln_phi.len(), 3);
    assert!(state.ln_phi.iter().all(|value| value.is_finite()));
}
