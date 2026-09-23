//! The DIIS accelerator, against the probe that feeds it directly.
//!
//! The oracle is `validation/neqsim/ReactiveFlashProbe.java`'s `fluid=diis-accelerator` block,
//! whose capture is `captures/reactive_flash_probe.tsv`. The class is reached only from inside
//! the RAND solve, so its extrapolation never appears in a flash's numbers - and it changes
//! them: the forced-one-phase 300 K water-gas shift accepts 19 of its extrapolated steps, which
//! the same capture records against that fluid. Feeding the accelerator a scripted sequence is
//! what makes the Pulay system, the rolling buffer's wrap and the `null` on a singular matrix
//! all readable.
//!
//! The sequence is three entries wide and seven long against a history of four, so five of the
//! seven extrapolations run on a buffer that has wrapped.

use azoth_reactions::diis::DiisAccelerator;

/// `extrapolated[step]` from the capture, for the two steps where a combination exists.
const CAPTURED: [(usize, [f64; 3]); 6] = [
    (
        1,
        [
            1.592_150_170_648_464_9,
            3.184_300_341_296_929_7,
            4.776_450_511_945_395,
        ],
    ),
    (
        2,
        [
            0.499_999_999_999_399_6,
            0.999_999_999_998_799_2,
            1.499_999_999_998_198_8,
        ],
    ),
    (
        3,
        [
            0.500_000_000_000_008,
            1.000_000_000_000_016,
            1.500_000_000_000_021_3,
        ],
    ),
    (
        4,
        [
            0.499_999_999_999_989_8,
            0.999_999_999_999_979_6,
            1.499_999_999_999_977,
        ],
    ),
    (
        5,
        [
            0.500_000_000_000_852_7,
            1.000_000_000_001_705_3,
            1.500_000_000_002_572_2,
        ],
    ),
    (
        6,
        [
            0.500_000_000_001_357_1,
            1.000_000_000_002_714_3,
            1.500_000_000_004_064_3,
        ],
    ),
];

/// The scripted sequence the probe feeds: `iterate_i = 0.5 (step + 1) (i + 1)` and
/// `residual_i = (i + 1) / (step + 1) + 0.1 step`, three components wide.
fn entry(step: usize) -> ([f64; 3], [f64; 3]) {
    let mut iterate = [0.0; 3];
    let mut residual = [0.0; 3];
    for i in 0..3 {
        let component = i as f64 + 1.0;
        iterate[i] = 0.5 * (step as f64 + 1.0) * component;
        residual[i] = component / (step as f64 + 1.0) + 0.1 * step as f64;
    }
    (iterate, residual)
}

/// **The Pulay combination, on a buffer that wraps.** Two pairs is the class's minimum, and the
/// history holds four, so the later steps extrapolate over a rolling window - which is what
/// makes `bufferIndex`'s wrap a measurement rather than an assumption. The first entry returns
/// `None` and the combination appears on the second.
#[test]
fn the_pulay_combination_is_the_classes() {
    let mut diis = DiisAccelerator::new(3, 4);
    assert_eq!(diis.count(), 0, "nothing is stored at construction");
    assert!(!diis.can_extrapolate(), "and nothing can be combined");

    for step in 0..7 {
        let (iterate, residual) = entry(step);
        diis.add_entry(&iterate, &residual)
            .expect("the shapes agree");
        assert_eq!(diis.count(), (step + 1).min(4), "count at step {step}");
        assert_eq!(
            diis.can_extrapolate(),
            step > 0,
            "can_extrapolate at step {step}"
        );

        let extrapolated = diis.extrapolate();
        if step == 0 {
            assert!(extrapolated.is_none(), "one pair is not a combination");
            continue;
        }
        let extrapolated = extrapolated.expect("two pairs are");
        let (_, captured) = CAPTURED
            .iter()
            .find(|(index, _)| *index == step)
            .expect("every step past the first is in the capture");
        for (component, (got, want)) in extrapolated.iter().zip(captured).enumerate() {
            assert!(
                (got - want).abs() / want.abs() < 1.0e-9,
                "step {step}, component {component}: {got} against the capture's {want}"
            );
        }
    }
}

/// **A singular Pulay system is `null`, not a combination.** Two pairs whose residuals are
/// equal make every entry of the overlap matrix the same number, so the elimination leaves a
/// zero pivot: the class reports that it *can* extrapolate - it has the two entries it asks
/// for - and then returns `null`.
#[test]
fn identical_residuals_are_refused() {
    let mut singular = DiisAccelerator::new(3, 4);
    singular
        .add_entry(&[1.0, 1.0, 1.0], &[1.0, 2.0, 3.0])
        .expect("the shapes agree");
    singular
        .add_entry(&[2.0, 2.0, 2.0], &[1.0, 2.0, 3.0])
        .expect("the shapes agree");
    assert!(singular.can_extrapolate());
    assert!(
        singular.extrapolate().is_none(),
        "the capture has this one at null"
    );
}

/// **`reset` empties the history and the accelerator refuses again**, which is what the
/// capture's three lines after the loop show.
#[test]
fn reset_discards_the_history() {
    let mut diis = DiisAccelerator::new(3, 4);
    let (iterate, residual) = entry(0);
    diis.add_entry(&iterate, &residual)
        .expect("the shapes agree");
    let (iterate, residual) = entry(1);
    diis.add_entry(&iterate, &residual)
        .expect("the shapes agree");
    assert!(diis.extrapolate().is_some());

    diis.reset();
    assert_eq!(diis.count(), 0);
    assert!(!diis.can_extrapolate());
    assert!(diis.extrapolate().is_none());
}

/// The class throws where this refuses: a vector that is not the buffer's width would leave a
/// slot half-written, and a wrapped slot keeps the values it had, so the entry is rejected
/// instead of silently truncated.
#[test]
fn a_short_vector_is_refused() {
    let mut diis = DiisAccelerator::new(3, 4);
    assert!(
        diis.add_entry(&[1.0, 2.0], &[1.0, 2.0, 3.0]).is_err(),
        "two entries against a vector length of three"
    );
    assert_eq!(diis.count(), 0, "and nothing was stored");
}
