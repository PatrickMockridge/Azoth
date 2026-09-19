//! The SAFT-VR-Mie flash, against NeqSim's `TPflashSAFT`.
//!
//! The oracle is `validation/neqsim/SaftVrMieFlashProbe.java`, which runs the class three
//! ways - dispatched, directly, and directly on a system that was initialised first - and
//! **only the third splits anything**. The dispatched route reads the feed composition before
//! its own `init(0)`, so it solves a zero feed: the numbers below are the third column's.
//! The whole finding is written up at `~/Desktop/neqsim-tpflashsaft-reads-a-zero-feed.md`.

use azoth_core::units::{kelvins, pascals};
use azoth_eos::tp_flash_saft::tp_flash_saft;

fn names() -> Vec<String> {
    vec!["methane".to_string(), "n-butane".to_string()]
}

fn flash(t: f64, p_bar: f64) -> azoth_eos::results::SaftFlashResult {
    tp_flash_saft(&names(), kelvins(t), pascals(p_bar * 1.0e5), &[0.6, 0.4]).expect("a flash")
}

/// **Three states, and every one of them splits** - where NeqSim's own dispatch reports a
/// single gas phase at all three.
///
/// The tolerance is `1e-5` and the reason is the kernel's: this model's `eta` derivatives are
/// central differences, so the K-values carry that difference's arithmetic and the loop's own
/// residual plateaus at about `2e-6`.
#[test]
fn the_flash_splits_where_neqsim_does() {
    for (t, p, beta) in [
        (350.0, 30.0, 0.970_118_166_758_440),
        (250.0, 30.0, 0.535_808_206_580_083),
        (200.0, 30.0, 0.406_167_928_288_498),
    ] {
        let result = flash(t, p);
        assert_eq!(result.phase, azoth_eos::Phase::TwoPhase, "at {t} K");
        let solved = result.beta.expect("a split has a vapour fraction");
        assert!(
            (solved / beta - 1.0).abs() < 1.0e-5,
            "at {t} K the vapour fraction is {solved} against NeqSim's {beta}"
        );
    }
}

/// The compositions and K-values at 250 K, which the probe prints for its `init(0)`-first run:
/// `x = [0.163196975571159, 0.836803024428841]`, `y = [0.978419594354518, 0.0215804056454816]`,
/// `K = [5.99533083741967, 0.0257891174857343]`.
///
/// **`1e-4`, because the loop stops on a relative change of `1e-5`**: the iterate it returns
/// carries that change by construction, so asking the answer for more than the stopping rule
/// asked the iterate for would be measuring the tolerance rather than the flash.
#[test]
fn the_compositions_are_neqsims() {
    let result = flash(250.0, 30.0);
    for (i, (x, y, k)) in [
        (
            0.163_196_975_571_159,
            0.978_419_594_354_518,
            5.995_330_837_419_67,
        ),
        (
            0.836_803_024_428_841,
            0.021_580_405_645_481_6,
            0.025_789_117_485_734_3,
        ),
    ]
    .into_iter()
    .enumerate()
    {
        assert!(
            (result.x[i] / x - 1.0).abs() < 1.0e-4,
            "liquid {i}: {} against {x}",
            result.x[i]
        );
        assert!(
            (result.y[i] / y - 1.0).abs() < 1.0e-4,
            "vapour {i}: {} against {y}",
            result.y[i]
        );
        assert!(
            (result.k[i] / k - 1.0).abs() < 1.0e-4,
            "K {i}: {} against {k}",
            result.k[i]
        );
    }
    assert!(
        result.iterations > 0 && result.iterations < 50,
        "the loop should converge rather than exhaust its steps: {}",
        result.iterations
    );
}

/// The seed is Wilson's, and the databank's constants give NeqSim's own.
///
/// The probe prints `K = [5.58119075962650, 0.0138173256591864]` for methane/n-butane at
/// 250 K and 30 bara. **The ratio of pressures is the only place a unit could hide**, because
/// it is dimensionless - so the databank's pascals and NeqSim's bar give the same seed, and
/// this checks that rather than assuming it.
#[test]
fn the_seed_is_wilsons() {
    let methane = azoth_eos::databank::entry("methane", None).expect("methane");
    let butane = azoth_eos::databank::entry("n-butane", None).expect("n-butane");
    let k = azoth_eos::tp_flash_saft::wilson_k(
        &[methane.pc, butane.pc],
        &[methane.omega, butane.omega],
        &[methane.tc, butane.tc],
        250.0,
        30.0e5,
    );
    assert!(
        (k[0] / 5.581_190_759_626_50 - 1.0).abs() < 1.0e-12,
        "methane's Wilson K is {} against NeqSim's 5.58119075962650",
        k[0]
    );
    assert!(
        (k[1] / 0.013_817_325_659_186_4 - 1.0).abs() < 1.0e-12,
        "n-butane's is {} against 0.0138173256591864",
        k[1]
    );
}
