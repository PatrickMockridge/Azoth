//! The mixing-rule seam: the classic rule is the default, and the temperature-
//! dependent rule resolves its `kij` once per state.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::Cubic;
use azoth_eos::databank;
use azoth_eos::mixing_rule::{MixingRule, SoreideWhitsonRole};

/// A symmetric 2 x 2 interaction matrix with `k01 = k10 = k`.
fn pair_matrix(k: f64) -> Vec<f64> {
    vec![0.0, k, k, 0.0]
}

#[test]
fn classic_t_resolves_kij_inversely_or_linearly() {
    let kij_t = pair_matrix(10.0);
    let inverse = MixingRule::ClassicT {
        kij: pair_matrix(0.05),
        kij_t: kij_t.clone(),
        inverse_temperature: true,
    };
    assert_eq!(inverse.effective_kij(300.0)[1], 0.05 + 10.0 / 300.0);

    let linear = MixingRule::ClassicT {
        kij: pair_matrix(0.05),
        kij_t,
        inverse_temperature: false,
    };
    assert_eq!(
        linear.effective_kij(300.0)[1],
        0.05 + 10.0 * (300.0 / 273.15 - 1.0)
    );
}

#[test]
fn classic_t2_resolves_each_pair_against_its_own_form() {
    let rule = MixingRule::ClassicT2 {
        kij: pair_matrix(0.05),
        kij_t: pair_matrix(10.0),
        inverse_temperature: vec![false, false, true, true],
    };
    let effective = rule.effective_kij(300.0);
    // The (0,1) interaction is flagged linear: `k0 + kt * T`, plain `T`.
    assert_eq!(effective[1], 0.05 + 10.0 * 300.0);
    // The (1,0) interaction is flagged inverse: `k0 + kt / T`.
    assert_eq!(effective[2], 0.05 + 10.0 / 300.0);
    // The diagonal is zero under either form.
    assert_eq!(effective[0], 0.0);
    assert_eq!(effective[3], 0.0);
}

/// The Soreide-Whitson aqueous correlations, against NeqSim 3.20.0's
/// `getkijWhitsonSoreideAqueous` (LEGACY) at S = 2 mol/kg, T = 300 K.
///
/// The four components are water, methane, nitrogen and CO2 in that order, with the
/// reduced temperatures and acentric factors NeqSim resolves for them. The expected
/// values are NeqSim's output, recorded rather than recomputed.
#[test]
fn soreide_whitson_reproduces_neqsims_aqueous_kij() {
    use SoreideWhitsonRole::{CarbonDioxide, Hydrocarbon, Nitrogen, Water};

    let base = vec![
        0.0, 0.485, 0.4778, 0.1896, 0.485, 0.0, 0.0311, 0.107, 0.4778, 0.0311, 0.0, 0.0, 0.1896,
        0.107, 0.0, 0.0,
    ];
    let rule = MixingRule::SoreideWhitson {
        kij: base.clone(),
        roles: vec![Water, Hydrocarbon, Nitrogen, CarbonDioxide],
        salinity: 2.0,
    };
    let tr = [
        0.463_463_618_105_978_7,
        1.574_307_304_785_894_4,
        2.379_064_234_734_338_3,
        0.986_225_714_191_788_1,
    ];
    let omega = [0.344, 0.0115, 0.0403, 0.2276];

    // A non-aqueous phase (water.x < 0.8) keeps the base matrix untouched.
    let gas = rule.phase_kij(&base, &tr, &omega, &[0.2, 0.3, 0.2, 0.3]);
    assert_eq!(gas, base);

    // The water-rich phase applies the correlation only where the *second* component
    // is water - NeqSim's asymmetry - so (gas, water) entries differ from (water, gas).
    let kij = rule.phase_kij(&base, &tr, &omega, &[0.9, 0.05, 0.03, 0.02]);
    assert!(
        (kij[4] - (-0.175_225_767_447_587_53)).abs() < 1e-12,
        "methane-water HC"
    );
    assert!(
        (kij[8] - (-0.574_890_598_658_567_8)).abs() < 1e-12,
        "nitrogen-water N2"
    );
    assert!(
        (kij[12] - (-0.081_272_539_490_276_38)).abs() < 1e-12,
        "CO2-water CO2"
    );
    // The water-first entries keep the base kij.
    assert_eq!(kij[1], 0.485, "water-methane");
    assert_eq!(kij[2], 0.4778, "water-nitrogen");
    assert_eq!(kij[3], 0.1896, "water-CO2");
    // The water-water and the gas-gas pairs are zero and the base respectively.
    assert_eq!(kij[0], 0.0);
    assert_eq!(kij[5], 0.0);
    assert_eq!(kij[6], 0.0311, "methane-nitrogen");
}

#[test]
fn classic_rule_is_unchanged_by_temperature() {
    let classic = MixingRule::Classic {
        kij: pair_matrix(0.05),
    };
    assert_eq!(classic.effective_kij(300.0)[1], 0.05);
    assert_eq!(classic.effective_kij(200.0)[1], 0.05);
}

#[test]
fn a_mixture_resolves_the_rule_at_reduced_parameters() {
    let (base, _) =
        databank::mixture_of(&["methane", "n-butane"], Cubic::Pr, None).expect("the pair resolves");
    let kij0 = base.kij(0, 1);
    let n = base.len();
    let mut kij = vec![0.0; n * n];
    kij[1] = kij0;
    kij[n] = kij0;
    let mut kij_t = vec![0.0; n * n];
    kij_t[1] = 1.0;
    kij_t[n] = 1.0;

    let mixture = base.with_mixing_rule(MixingRule::ClassicT {
        kij,
        kij_t,
        inverse_temperature: true,
    });
    let reduced = mixture
        .reduced_parameters(kelvins(300.0), pascals(1_000_000.0))
        .expect("the state reduces");
    assert!((reduced.kij[1] - (kij0 + 1.0 / 300.0)).abs() < 1e-15);
}

/// The UMR rule's two inputs, against NeqSim 3.20.0's `SystemUMRCPAEoS`.
///
/// Methane/water 0.98/0.02 at 298.15 K, from `validation/neqsim/UmrCpaProbe.java`: the
/// state NeqSim's own `TPflashUMRCPADehydrationLifecycleTest` pins its water-in-gas
/// envelope at. The expected values are NeqSim's output, recorded rather than recomputed.
///
/// **The UNIFAC set is the `_umrmc` one**, because the UMR-CPA model pairs the rule
/// with a Mathias-Copeman attraction term, and NeqSim's
/// `PhaseGEUnifacUMRPRU.useMcInteractionParameters` selects that table for attractive
/// terms 13, 19 and 22. The `_umr` tables give a different `ln gamma` here, so this key
/// is what says which was read.
#[test]
fn the_unifac_umrpru_set_reproduces_neqsims_activity_coefficients() {
    use azoth_eos::databank::UmrpruSet;

    let params = databank::unifac_umrpru_parameters(&["methane", "water"], UmrpruSet::Umrmc)
        .expect("the pair resolves");
    let result = azoth_eos::unifac_umrpru_activity_coefficients(&params, 298.15, &[0.98, 0.02])
        .expect("the coefficients evaluate");

    // The probe's `lnGamma[0]`/`lnGamma[1]`.
    for (i, expected) in [0.000_245_756_047_888_7_f64, 0.887_875_427_408_146]
        .iter()
        .enumerate()
    {
        assert!(
            (result.ln_gamma[i] - expected).abs() < 1e-10,
            "ln_gamma[{i}] = {}, NeqSim's {expected}",
            result.ln_gamma[i]
        );
    }
}

/// The UMR rule's `alpha_mix`, against the probe's own `alpha_mix`.
///
/// The inputs are the probe's **printed** `qPure_aT_over_bRT` and `lnGamma`, so this
/// checks the rule's arithmetic and its constant against a measured NeqSim number
/// without depending on a databank column. `hwfc` is `-1/0.53` there because the GE
/// model is `UNIFAC_UMRPRU`; the Huron-Vidal rule beside it takes the cubic's own
/// `hv_constant()`, and nothing about a single component distinguishes the two.
#[test]
fn the_umr_rule_reproduces_neqsims_alpha_mix() {
    use azoth_eos::mixing_rule::umr_ader;

    let qpure = [3.071_869_686_431_51, 6.766_461_302_779_70];
    let ln_gamma = [0.000_245_756_047_888_7, 0.887_875_427_408_146];
    let ader = umr_ader(&qpure, &ln_gamma);
    let alpha_mix = 0.98 * ader[0] + 0.02 * ader[1];

    assert!(
        (alpha_mix - 3.111_802_368_805_46).abs() < 1e-13,
        "alpha_mix = {alpha_mix}, NeqSim's 3.11180236880546"
    );
}
