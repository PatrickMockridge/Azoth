//! The mixing-rule seam: the classic rule is the default, and the temperature-
//! dependent rule resolves its `kij` once per state.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::databank;
use azoth_eos::mixing_rule::MixingRule;

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
        databank::mixture_of(&["methane", "n-butane"], None).expect("the pair resolves");
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
